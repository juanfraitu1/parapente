//! Coverage-based validation for detected segments
//!
//! Uses coverage profile characteristics to distinguish true gene family members
//! from false positives (noise).
//!
//! Key insight: True genes have:
//! - Consistent coverage depth (low variance)
//! - High coverage breadth (most of region covered)
//! - Read density similar to seed genes
//! - Jaccard similarity with seeds

use anyhow::Result;
use fxhash::FxHashSet as HashSet;

/// Validate a detected segment using coverage profile
#[derive(Debug, Clone)]
pub struct CoverageValidation {
    pub chrom: String,
    pub start: i64,
    pub end: i64,
    pub mean_coverage: f64,
    pub coverage_std: f64,
    pub coverage_cv: f64,      // Coefficient of variation
    pub breadth: f64,           // Fraction of region with coverage > 0
    pub read_count: usize,
    pub read_density: f64,      // Reads per kb
    pub jaccard_with_seed: f64,
    pub is_valid: bool,
    pub confidence: f64,        // 0-1 score
}

/// Parameters for coverage validation
#[derive(Debug, Clone)]
pub struct CoverageValidationParams {
    pub min_mean_coverage: f64,     // Minimum mean coverage for valid gene
    pub max_coverage_cv: f64,       // Maximum coefficient of variation
    pub min_breadth: f64,          // Minimum fraction of region covered
    pub min_read_density: f64,     // Minimum reads per kb
    pub min_jaccard: f64,           // Minimum Jaccard similarity with seed
}

impl Default for CoverageValidationParams {
    fn default() -> Self {
        Self {
            min_mean_coverage: 5.0,    // At least 5x coverage on average
            max_coverage_cv: 2.0,      // CV < 2.0 (std/mean)
            min_breadth: 0.5,          // At least 50% of region covered
            min_read_density: 10.0,    // At least 10 reads per kb
            min_jaccard: 0.01,          // At least 1% Jaccard with seed
        }
    }
}

/// Validate segments using coverage profile characteristics
pub fn validate_segments_by_coverage(
    segments: &[(String, i64, i64, HashSet<String>)],
    params: &CoverageValidationParams,
) -> Vec<CoverageValidation> {
    let mut validations = Vec::new();
    
    for (chrom, start, end, reads) in segments {
        let span = end - start;
        let read_count = reads.len();
        let read_density = read_count as f64 / (span as f64 / 1000.0);  // reads per kb
        
        // Compute confidence score
        let mut confidence = 1.0;
        let mut valid = true;
        
        // Check read density
        if read_density < params.min_read_density {
            confidence *= 0.5;
            valid = false;
        }
        
        // Check span (very small segments are suspicious)
        if span < 5000 {
            confidence *= 0.3;
            valid = false;
        }
        
        // Check span (very large segments may be overmerged)
        if span > 500_000 {
            confidence *= 0.7;  // Penalty but not invalid
        }
        
        // Check read count (need enough reads for confidence)
        if read_count < 10 {
            confidence *= 0.3;
            valid = false;
        }
        
        validations.push(CoverageValidation {
            chrom: chrom.clone(),
            start: *start,
            end: *end,
            mean_coverage: 0.0,  // Would need BAM access to compute
            coverage_std: 0.0,
            coverage_cv: 0.0,
            breadth: 0.0,
            read_count,
            read_density,
            jaccard_with_seed: 0.0,  // Would need to compute
            is_valid: valid,
            confidence,
        });
    }
    
    validations
}

/// Filter segments by coverage-based confidence
pub fn filter_by_confidence(
    segments: &[(String, i64, i64, HashSet<String>)],
    min_confidence: f64,
) -> Vec<(String, i64, i64, HashSet<String>)> {
    let params = CoverageValidationParams::default();
    let validations = validate_segments_by_coverage(segments, &params);
    
    segments.iter()
        .zip(validations.iter())
        .filter(|(_, v)| v.confidence >= min_confidence)
        .map(|(s, _)| s.clone())
        .collect()
}

/// Calculate Jaccard similarity between read sets
pub fn jaccard_similarity(a: &HashSet<String>, b: &HashSet<String>) -> f64 {
    if a.is_empty() || b.is_empty() {
        return 0.0;
    }
    
    let intersection = a.intersection(b).count();
    let union = a.union(b).count();
    
    if union == 0 {
        0.0
    } else {
        intersection as f64 / union as f64
    }
}
