//! Enhanced coverage-based validation for gene family detection
//!
//! This module provides multi-layer validation using coverage features:
//! 1. Peak analysis (height, width, symmetry)
//! 2. Slope analysis (boundary sharpness)
//! 3. Pattern matching (gene-like coverage profile)

use anyhow::Result;

/// Quality metrics derived from coverage profile
#[derive(Debug, Clone)]
pub struct CoverageQuality {
    pub peak_height: f64,           // Maximum coverage in region
    pub peak_width: i64,            // Width at half-maximum (gene size estimate)
    pub symmetry_score: f64,        // 0=perfectly symmetric, 1=highly asymmetric
    pub left_boundary_score: f64,   // 5' boundary steepness (0-1)
    pub right_boundary_score: f64,  // 3' boundary steepness (0-1)
    pub overall_score: f64,         // Combined quality score (0-1)
}

impl CoverageQuality {
    /// Returns true if coverage pattern looks like a real gene
    pub fn is_gene_like(&self, min_score: f64) -> bool {
        self.overall_score >= min_score
            && self.peak_height >= 5.0  // At least 5 reads
            && self.symmetry_score < 0.5  // Not too asymmetric
    }
    
    /// Returns confidence level based on quality metrics
    pub fn confidence(&self) -> ValidationConfidence {
        if self.overall_score >= 0.8 
            && self.symmetry_score < 0.3 
            && self.left_boundary_score > 0.5 
            && self.right_boundary_score > 0.5 {
            ValidationConfidence::High
        } else if self.overall_score >= 0.5 {
            ValidationConfidence::Medium
        } else {
            ValidationConfidence::Low
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ValidationConfidence {
    High,
    Medium,
    Low,
}

/// Analyze coverage profile quality
/// 
/// Computes peak width at half maximum, symmetry, and boundary slopes
pub fn analyze_coverage_quality(bins: &[f64], bin_size: i64) -> CoverageQuality {
    if bins.is_empty() {
        return CoverageQuality {
            peak_height: 0.0,
            peak_width: 0,
            symmetry_score: 1.0,
            left_boundary_score: 0.0,
            right_boundary_score: 0.0,
            overall_score: 0.0,
        };
    }
    
    let peak_height = bins.iter().copied().fold(0.0f64, f64::max);
    let half_max = peak_height / 2.0;
    
    // Find peak center (index of max)
    let peak_idx = bins.iter()
        .enumerate()
        .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
        .map(|(i, _)| i)
        .unwrap_or(0);
    
    // Calculate peak width at half maximum (FWHM approximation)
    let mut left_half = peak_idx;
    let mut right_half = peak_idx;
    
    for i in (0..peak_idx).rev() {
        if bins[i] < half_max {
            left_half = i;
            break;
        }
    }
    
    for i in peak_idx..bins.len() {
        if bins[i] < half_max {
            right_half = i;
            break;
        }
    }
    
    let peak_width = (right_half - left_half) as i64 * bin_size;
    
    // Calculate symmetry: compare area left vs right of peak
    let left_sum: f64 = bins[..peak_idx].iter().sum();
    let right_sum: f64 = bins[peak_idx+1..].iter().sum();
    let total_sum = left_sum + right_sum;
    
    let symmetry_score = if total_sum > 0.0 {
        (left_sum - right_sum).abs() / total_sum
    } else {
        1.0
    };
    
    // Calculate boundary slopes (steepness at edges)
    let slope_window = 3.min(bins.len() / 4);
    
    let left_slope = if slope_window > 0 && peak_idx > slope_window {
        let start_cov = bins[..slope_window].iter().sum::<f64>() / slope_window as f64;
        let peak_cov = bins[peak_idx];
        (peak_cov - start_cov) / peak_idx as f64
    } else {
        0.0
    };
    
    let right_slope = if slope_window > 0 && peak_idx + slope_window < bins.len() {
        let end_cov = bins[bins.len()-slope_window..].iter().sum::<f64>() / slope_window as f64;
        let peak_cov = bins[peak_idx];
        (peak_cov - end_cov) / (bins.len() - peak_idx) as f64
    } else {
        0.0
    };
    
    // Normalize slopes to 0-1 range (assuming max reasonable slope is 5.0)
    let left_boundary_score = (left_slope / 5.0).clamp(0.0, 1.0);
    let right_boundary_score = (right_slope / 5.0).clamp(0.0, 1.0);
    
    // Overall quality score
    let overall_score = if peak_height > 0.0 {
        let height_score = (peak_height / 50.0).min(1.0);  // Normalize to 50 reads max
        let boundary_avg = (left_boundary_score + right_boundary_score) / 2.0;
        let symmetry_penalty = 1.0 - symmetry_score;  // Higher is better
        
        (height_score * 0.3) + (boundary_avg * 0.4) + (symmetry_penalty * 0.3)
    } else {
        0.0
    };
    
    CoverageQuality {
        peak_height,
        peak_width,
        symmetry_score,
        left_boundary_score,
        right_boundary_score,
        overall_score,
    }
}

/// Multi-layer validation result
#[derive(Debug, Clone)]
pub struct ValidationResult {
    pub coverage_confidence: ValidationConfidence,
    pub jaccard_score: f64,
    pub isoseq_score: Option<f64>,
    pub multi_seed_confirmed: bool,
    pub overall_confidence: ValidationConfidence,
    pub is_valid: bool,
}

/// Comprehensive validation using all available signals
pub fn validate_gene_core(
    coverage_quality: &CoverageQuality,
    jaccard_score: f64,
    isoseq_score: Option<f64>,
    found_by_n_seeds: usize,
    params: &ValidationParams,
) -> ValidationResult {
    // Layer 1: Coverage quality
    let coverage_confidence = coverage_quality.confidence();
    let coverage_pass = coverage_quality.overall_score >= params.min_coverage_score;
    
    // Layer 2: Multi-mapping (Jaccard)
    let jaccard_pass = jaccard_score >= params.min_jaccard;
    
    // Layer 3: Iso-Seq (if available)
    let isoseq_pass = match isoseq_score {
        Some(score) => score >= params.min_isoseq_sim,
        None => true, // Skip if not available
    };
    
    // Layer 4: Multi-seed confirmation
    let multi_seed_confirmed = found_by_n_seeds >= params.min_seeds_for_confirmation;
    
    // Calculate overall confidence
    let mut confidence_score = 0.0;
    let mut weights = 0.0;
    
    // Coverage quality contribution
    let coverage_weight = match coverage_confidence {
        ValidationConfidence::High => (1.0, 0.3),
        ValidationConfidence::Medium => (0.6, 0.3),
        ValidationConfidence::Low => (0.2, 0.3),
    };
    confidence_score += coverage_weight.0 * coverage_weight.1;
    weights += coverage_weight.1;
    
    // Jaccard contribution
    let jaccard_weight = 0.2;
    confidence_score += jaccard_score.min(1.0) * jaccard_weight;
    weights += jaccard_weight;
    
    // Iso-Seq contribution (if available)
    if let Some(score) = isoseq_score {
        let isoseq_weight = 0.3;
        confidence_score += score * isoseq_weight;
        weights += isoseq_weight;
    }
    
    // Multi-seed bonus
    if multi_seed_confirmed {
        confidence_score += 0.2;
    }
    
    // Normalize
    let normalized_score = if weights > 0.0 {
        confidence_score / (weights + if multi_seed_confirmed { 0.0 } else { 0.2 })
    } else {
        0.0
    };
    
    let overall_confidence = if normalized_score >= 0.8 {
        ValidationConfidence::High
    } else if normalized_score >= 0.5 {
        ValidationConfidence::Medium
    } else {
        ValidationConfidence::Low
    };
    
    // Final decision
    let is_valid = coverage_pass && jaccard_pass && isoseq_pass;
    
    ValidationResult {
        coverage_confidence,
        jaccard_score,
        isoseq_score,
        multi_seed_confirmed,
        overall_confidence,
        is_valid,
    }
}

/// Parameters for validation
#[derive(Debug, Clone)]
pub struct ValidationParams {
    pub min_coverage_score: f64,       // Minimum coverage quality (0-1)
    pub min_jaccard: f64,              // Minimum Jaccard similarity
    pub min_isoseq_sim: f64,           // Minimum Iso-Seq similarity
    pub min_seeds_for_confirmation: usize, // Seeds that must find this locus
}

impl Default for ValidationParams {
    fn default() -> Self {
        Self {
            min_coverage_score: 0.4,
            min_jaccard: 0.01,
            min_isoseq_sim: 0.6,
            min_seeds_for_confirmation: 2,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_perfect_gene_profile() {
        // Symmetric peak: 0,1,2,3,4,5,6,5,4,3,2,1,0
        let bins: Vec<f64> = vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0];
        let quality = analyze_coverage_quality(&bins, 100);
        
        assert!(quality.peak_height > 5.0);
        assert!(quality.symmetry_score < 0.2);  // Very symmetric
        assert!(quality.left_boundary_score > 0.3);
        assert!(quality.right_boundary_score > 0.3);
        assert!(quality.overall_score > 0.5);
        assert_eq!(quality.confidence(), ValidationConfidence::High);
    }
    
    #[test]
    fn test_asymmetric_profile() {
        // Asymmetric: 0,1,2,3,4,5,6,2,1,0,0,0,0
        let bins: Vec<f64> = vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0];
        let quality = analyze_coverage_quality(&bins, 100);
        
        assert!(quality.symmetry_score > 0.3);  // Asymmetric
        assert!(quality.confidence() != ValidationConfidence::High);
    }
}
