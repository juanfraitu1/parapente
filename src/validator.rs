//! Validation module for distinguishing true genes from false positives

use anyhow::Result;
use noodles::bam;
use noodles::core::Position;
use noodles::sam::alignment::Record as SamRecord;
use std::collections::HashSet;

/// Validation evidence for a detected locus
#[derive(Debug, Clone)]
pub struct ValidationEvidence {
    /// Coverage-based standard deviation (higher = more gene-like)
    pub coverage_std: f64,
    /// Number of spliced reads
    pub spliced_read_count: usize,
    /// Number of unspliced reads
    pub unspliced_read_count: usize,
    /// Transcript length variability (coefficient of variation)
    pub length_cv: f64,
    /// Mean transcript length
    pub mean_length: f64,
    /// Consistency with seed's splice pattern (0-1)
    pub splice_consistency: f64,
    /// Whether this locus is found by multiple seeds
    pub multi_seed_support: bool,
}

impl ValidationEvidence {
    /// Calculate a composite confidence score (0-1)
    pub fn confidence_score(&self) -> f64 {
        // Coverage variability is key - genes have high std due to exons/introns
        let coverage_score = if self.coverage_std > 50.0 {
            1.0
        } else if self.coverage_std > 20.0 {
            0.7
        } else if self.coverage_std > 10.0 {
            0.4
        } else {
            0.1
        };

        // Spliced reads indicate mature mRNA (gene-like)
        let total_reads = self.spliced_read_count + self.unspliced_read_count;
        let splicing_score = if total_reads > 0 {
            let ratio = self.spliced_read_count as f64 / total_reads as f64;
            // Both very high and very low can be suspicious
            if ratio > 0.3 && ratio < 0.95 {
                0.9 // Good balance
            } else if ratio >= 0.95 {
                0.6 // Might be pre-mRNA contamination
            } else {
                0.4 // Low splicing
            }
        } else {
            0.0
        };

        // Length consistency
        let length_score = if self.length_cv < 0.3 {
            0.9 // Very consistent
        } else if self.length_cv < 0.5 {
            0.7
        } else if self.length_cv < 0.7 {
            0.4
        } else {
            0.2 // Highly variable
        };

        // Multi-seed support is strong evidence
        let seed_score = if self.multi_seed_support { 1.0 } else { 0.5 };

        // Weighted combination
        let score = coverage_score * 0.35
            + splicing_score * 0.25
            + length_score * 0.20
            + seed_score * 0.20;

        score
    }

    /// Classify as true gene, uncertain, or false positive
    pub fn classify(&self) -> LocusClassification {
        let score = self.confidence_score();
        
        if score > 0.7 {
            LocusClassification::TrueGene
        } else if score > 0.4 {
            LocusClassification::Uncertain
        } else {
            LocusClassification::FalsePositive
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum LocusClassification {
    TrueGene,
    Uncertain,
    FalsePositive,
}

impl std::fmt::Display for LocusClassification {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            LocusClassification::TrueGene => write!(f, "true_gene"),
            LocusClassification::Uncertain => write!(f, "uncertain"),
            LocusClassification::FalsePositive => write!(f, "false_positive"),
        }
    }
}

/// Collect coverage statistics for a region
pub fn analyze_coverage(
    bam_path: &str,
    chrom: &str,
    start: usize,
    end: usize,
    bin_size: usize,
) -> Result<(f64, Vec<usize>)> {
    let mut reader = bam::io::indexed_reader::Builder::default()
        .build_from_path(bam_path)?;

    let header = reader.read_header()?;

    let start_pos = Position::try_from(start)?;
    let end_pos = Position::try_from(end)?;
    let region = noodles::core::Region::new(chrom, start_pos..=end_pos);

    let num_bins = (end - start) / bin_size + 1;
    let mut coverage = vec![0usize; num_bins];

    for result in reader.query(&header, &region)? {
        let record = result?;

        if record.flags().is_unmapped() {
            continue;
        }

        let read_start = record.alignment_start()
            .transpose()?
            .map(|p| usize::from(p))
            .unwrap_or(0);
        let read_end = record.alignment_end()
            .transpose()?
            .map(|p| usize::from(p))
            .unwrap_or(read_start + 1);

        // Clip to region
        let rs = read_start.max(start);
        let re = read_end.min(end);

        if rs >= re {
            continue;
        }

        let start_bin = (rs - start) / bin_size;
        let end_bin = (re - start) / bin_size;

        for b in start_bin..=end_bin.min(num_bins - 1) {
            coverage[b] += 1;
        }
    }

    // Calculate standard deviation
    let mean = coverage.iter().sum::<usize>() as f64 / coverage.len() as f64;
    let variance = coverage
        .iter()
        .map(|&x| (x as f64 - mean).powi(2))
        .sum::<f64>()
        / coverage.len() as f64;
    let std = variance.sqrt();

    Ok((std, coverage))
}

/// Analyze transcript characteristics
pub fn analyze_transcripts(
    bam_path: &str,
    chrom: &str,
    start: usize,
    end: usize,
) -> Result<(usize, usize, Vec<usize>)> {
    let mut reader = bam::io::indexed_reader::Builder::default()
        .build_from_path(bam_path)?;

    let header = reader.read_header()?;

    let start_pos = Position::try_from(start)?;
    let end_pos = Position::try_from(end)?;
    let region = noodles::core::Region::new(chrom, start_pos..=end_pos);

    let mut spliced = 0;
    let mut unspliced = 0;
    let mut lengths = Vec::new();

    for result in reader.query(&header, &region)? {
        let record = result?;

        if record.flags().is_unmapped() {
            continue;
        }

        // Check for splice junctions (CIGAR N operation)
        let has_splice = record.cigar().iter().any(|op| {
            use noodles::sam::alignment::record::cigar::Op;
            matches!(op, Op::Skip(_))
        });

        if has_splice {
            spliced += 1;
        } else {
            unspliced += 1;
        }

        // Get transcript length
        let read_start = record.alignment_start()
            .transpose()?
            .map(|p| usize::from(p))
            .unwrap_or(0);
        let read_end = record.alignment_end()
            .transpose()?
            .map(|p| usize::from(p))
            .unwrap_or(read_start + 1);
        
        lengths.push(read_end - read_start);
    }

    Ok((spliced, unspliced, lengths))
}

/// Validate a locus using multiple criteria
pub fn validate_locus(
    bam_path: &str,
    chrom: &str,
    start: i64,
    end: i64,
    multi_seed_support: bool,
) -> Result<ValidationEvidence> {
    let start_usize = start.max(1) as usize;
    let end_usize = end.max(1) as usize;

    // Coverage analysis
    let (coverage_std, _) = analyze_coverage(bam_path, chrom, start_usize, end_usize, 1000)?;

    // Transcript analysis
    let (spliced, unspliced, lengths) = analyze_transcripts(bam_path, chrom, start_usize, end_usize)?;

    // Calculate length statistics
    let (mean_length, length_cv) = if !lengths.is_empty() {
        let mean = lengths.iter().sum::<usize>() as f64 / lengths.len() as f64;
        let variance = lengths
            .iter()
            .map(|&x| (x as f64 - mean).powi(2))
            .sum::<f64>()
            / lengths.len() as f64;
        let std = variance.sqrt();
        let cv = if mean > 0.0 { std / mean } else { 0.0 };
        (mean, cv)
    } else {
        (0.0, 0.0)
    };

    Ok(ValidationEvidence {
        coverage_std,
        spliced_read_count: spliced,
        unspliced_read_count: unspliced,
        length_cv,
        mean_length,
        splice_consistency: 0.0, // Would need seed comparison
        multi_seed_support,
    })
}

/// Batch validate multiple loci
pub fn validate_loci_batch(
    bam_path: &str,
    loci: &[(String, i64, i64)],
    multi_seed_hits: Option<&HashSet<String>>, // chrom:start-end strings
) -> Result<Vec<ValidationEvidence>> {
    loci.iter()
        .map(|(chrom, start, end)| {
            let key = format!("{}:{}-{}", chrom, start, end);
            let multi_seed = multi_seed_hits
                .map(|set| set.contains(&key))
                .unwrap_or(false);
            validate_locus(bam_path, chrom, *start, *end, multi_seed)
        })
        .collect()
}
