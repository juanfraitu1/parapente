//! Advanced Iso-Seq validation using read quality and structural features

use anyhow::Result;
use noodles::bam;
use noodles::core::Position;
use noodles::sam::alignment::Record;
use std::collections::HashMap;

/// Advanced transcript features for validation
/// NOTE: MAPQ is intentionally NOT used for filtering since multi-mapping reads (MAPQ=0)
/// are essential for gene family detection
#[derive(Debug, Clone)]
pub struct AdvancedTranscriptFeatures {
    pub total_transcripts: usize,
    pub mean_length: f64,
    pub length_cv: f64,
    pub strand_plus: usize,
    pub strand_minus: usize,
    pub coverage_std: f64,
    // Iso-Seq specific features (MAPQ excluded - multi-mapping reads are vital)
    pub mean_divergence: f64,
    pub mean_alignment_score: f64,
    pub mean_identity: f64,           // from cm tag (matches / total CIGAR ops)
    pub n_full_length: usize,
    pub n_partial: usize,
    pub intron_count_distribution: HashMap<usize, usize>,
    pub mean_introns: f64,
    pub strand_consistency: f64,      // 0-1, higher = more consistent
    pub mean_read_quality: f64,       // from rq tag (PacBio read quality)
    pub mean_num_passes: f64,         // from np tag (PacBio CCS passes)
    pub chimeric_read_ratio: f64,     // reads with supplementary alignments / total
    pub secondary_alignment_ratio: f64, // reads with secondary alignments / total
}

impl AdvancedTranscriptFeatures {
    /// Calculate similarity score to another feature set (0-1)
    /// Uses Iso-Seq specific features for validation
    pub fn similarity(&self, other: &AdvancedTranscriptFeatures) -> f64 {
        // Strand consistency (25% weight) - genes have consistent strand
        let strand_score = if (self.strand_plus > self.strand_minus) == 
                               (other.strand_plus > other.strand_minus) { 1.0 } else { 0.0 };
        
        // Length similarity (20% weight) - paralogs have similar transcript lengths
        let len_ratio = other.mean_length / self.mean_length;
        let length_score = if len_ratio >= 0.7 && len_ratio <= 1.3 { 1.0 } 
                          else if len_ratio >= 0.5 && len_ratio <= 2.0 { 0.5 } else { 0.0 };
        
        // Length consistency - low CV = gene-like (15% weight)
        let cv_score = if other.length_cv < 0.5 { 1.0 } else { 0.5 };
        
        // Coverage variability - high std = exon/intron structure (15% weight)
        let cov_score = if other.coverage_std > 50.0 { 1.0 } else { 0.5 };
        
        // Intron structure similarity (10% weight)
        let intron_diff = (self.mean_introns - other.mean_introns).abs();
        let intron_score = if intron_diff <= 2.0 { 1.0 } 
                          else if intron_diff <= 5.0 { 0.5 } else { 0.0 };
        
        // Divergence similarity - paralogs have similar divergence patterns (10% weight)
        let div_diff = (self.mean_divergence - other.mean_divergence).abs();
        let divergence_score = if div_diff < 0.02 { 1.0 }
                               else if div_diff < 0.05 { 0.7 }
                               else if div_diff < 0.10 { 0.4 } else { 0.0 };
        
        // Identity similarity (5% weight)
        let id_diff = (self.mean_identity - other.mean_identity).abs();
        let identity_score = if id_diff < 0.02 { 1.0 }
                             else if id_diff < 0.05 { 0.7 }
                             else if id_diff < 0.10 { 0.4 } else { 0.0 };
        
        strand_score * 0.25 + 
        length_score * 0.20 + 
        cv_score * 0.15 + 
        cov_score * 0.15 +
        intron_score * 0.10 +
        divergence_score * 0.10 +
        identity_score * 0.05
    }
}

/// Count introns from CIGAR string
fn count_introns(record: &bam::Record) -> usize {
    let cigar = record.cigar();
    let mut intron_count = 0;
    
    for result in cigar.iter() {
        if let Ok(op) = result {
            if op.kind() == noodles::sam::alignment::record::cigar::op::Kind::Skip {
                intron_count += 1;
            }
        }
    }
    
    intron_count
}

/// Extract read length from auxiliary data (rl tag - PacBio)
fn get_read_length(record: &bam::Record) -> Option<i64> {
    for result in record.data().iter() {
        if let Ok((tag, value)) = result {
            if tag.as_ref() == b"rl" {
                match value {
                    noodles::sam::alignment::record::data::field::Value::Int32(v) => return Some(v as i64),
                    noodles::sam::alignment::record::data::field::Value::Int16(v) => return Some(v as i64),
                    _ => {}
                }
            }
        }
    }
    None
}

/// Extract read quality from rq tag (PacBio)
fn get_read_quality(record: &bam::Record) -> Option<f64> {
    for result in record.data().iter() {
        if let Ok((tag, value)) = result {
            if tag.as_ref() == b"rq" {
                match value {
                    noodles::sam::alignment::record::data::field::Value::Float(v) => return Some(v as f64),
                    _ => {}
                }
            }
        }
    }
    None
}

/// Extract number of passes from np tag (PacBio CCS)
fn get_num_passes(record: &bam::Record) -> Option<i32> {
    for result in record.data().iter() {
        if let Ok((tag, value)) = result {
            if tag.as_ref() == b"np" {
                match value {
                    noodles::sam::alignment::record::data::field::Value::Int32(v) => return Some(v),
                    noodles::sam::alignment::record::data::field::Value::Int16(v) => return Some(v as i32),
                    noodles::sam::alignment::record::data::field::Value::Int8(v) => return Some(v as i32),
                    _ => {}
                }
            }
        }
    }
    None
}

/// Extract number of CIGAR matches from cm tag (used to calculate identity)
fn get_cigar_matches(record: &bam::Record) -> Option<i32> {
    for result in record.data().iter() {
        if let Ok((tag, value)) = result {
            if tag.as_ref() == b"cm" {
                match value {
                    noodles::sam::alignment::record::data::field::Value::Int32(v) => return Some(v),
                    noodles::sam::alignment::record::data::field::Value::Int16(v) => return Some(v as i32),
                    noodles::sam::alignment::record::data::field::Value::Int8(v) => return Some(v as i32),
                    _ => {}
                }
            }
        }
    }
    None
}

/// Calculate identity from cm tag and CIGAR operations
fn calculate_identity(record: &bam::Record) -> Option<f64> {
    let cm = get_cigar_matches(record)?;
    
    // Count total reference-consuming CIGAR operations
    let cigar = record.cigar();
    let mut total_ref_ops = 0i32;
    
    for result in cigar.iter() {
        if let Ok(op) = result {
            use noodles::sam::alignment::record::cigar::op::Kind;
            match op.kind() {
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch | 
                Kind::Deletion | Kind::Skip => {
                    total_ref_ops += op.len() as i32;
                }
                _ => {}
            }
        }
    }
    
    if total_ref_ops > 0 {
        Some(cm as f64 / total_ref_ops as f64)
    } else {
        None
    }
}

/// Extract divergence (de tag) from record
fn get_divergence(record: &bam::Record) -> Option<f64> {
    for result in record.data().iter() {
        if let Ok((tag, value)) = result {
            if tag.as_ref() == b"de" {
                if let noodles::sam::alignment::record::data::field::Value::Float(v) = value {
                    return Some(v as f64);
                }
            }
        }
    }
    None
}

/// Extract alignment score (AS tag) from record
fn get_alignment_score(record: &bam::Record) -> Option<i64> {
    for result in record.data().iter() {
        if let Ok((tag, value)) = result {
            if tag.as_ref() == b"AS" {
                match value {
                    noodles::sam::alignment::record::data::field::Value::Int32(v) => return Some(v as i64),
                    noodles::sam::alignment::record::data::field::Value::Int16(v) => return Some(v as i64),
                    noodles::sam::alignment::record::data::field::Value::Int8(v) => return Some(v as i64),
                    _ => {}
                }
            }
        }
    }
    None
}

/// Analyze transcripts with advanced Iso-Seq features
/// NOTE: MAPQ filtering is intentionally disabled - multi-mapping reads (MAPQ=0) 
/// are essential for gene family detection
pub fn analyze_advanced_features(
    bam_path: &str,
    chrom: &str,
    start: i64,
    end: i64,
    _min_mapq: u8,  // Kept for API compatibility but ignored
    max_divergence: f64,
) -> Result<AdvancedTranscriptFeatures> {
    let mut reader = bam::io::indexed_reader::Builder::default()
        .build_from_path(bam_path)?;

    let header = reader.read_header()?;
    
    let start_pos = Position::try_from(start.max(1) as usize)?;
    let end_pos = Position::try_from(end.max(1) as usize)?;
    let region = noodles::core::Region::new(chrom, start_pos..=end_pos);

    let mut total = 0usize;
    let mut plus_strand = 0usize;
    let mut minus_strand = 0usize;
    let mut lengths = Vec::new();
    let mut divergences = Vec::new();
    let mut alignment_scores = Vec::new();
    let mut identities = Vec::new();
    let mut read_qualities = Vec::new();
    let mut num_passes = Vec::new();
    let mut full_length = 0usize;
    let mut partial = 0usize;
    let mut intron_counts: HashMap<usize, usize> = HashMap::new();
    let mut chimeric_count = 0usize;
    let mut secondary_count = 0usize;

    // For coverage calculation
    let bin_size = 1000usize;
    let num_bins = ((end - start) as usize / bin_size) + 1;
    let mut coverage = vec![0usize; num_bins];

    for result in reader.query(&header, &region)? {
        let record = result?;

        if record.flags().is_unmapped() {
            continue;
        }

        // Track secondary and supplementary alignments (these are VALID for gene family detection)
        if record.flags().is_secondary() {
            secondary_count += 1;
        }
        if record.flags().is_supplementary() {
            chimeric_count += 1;
            // For chimeric reads, we still want to analyze them but note the chimeric status
        }

        // Filter by divergence only (NOT by MAPQ - multi-mapping reads are essential)
        let div = get_divergence(&record).unwrap_or(1.0);
        if div > max_divergence {
            continue;
        }

        total += 1;

        // Get positions
        let read_start = match record.alignment_start() {
            Some(Ok(p)) => p.get() as i64,
            _ => continue,
        };
        
        let read_end = match record.alignment_end() {
            Some(Ok(p)) => p.get() as i64,
            _ => read_start + 1000, // fallback
        };

        let length = (read_end - read_start) as f64;
        lengths.push(length);

        // Strand
        let flags_val = record.flags().bits();
        if flags_val & 0x10 != 0 {
            minus_strand += 1;
        } else {
            plus_strand += 1;
        }

        // Iso-Seq quality metrics (NOT including MAPQ)
        divergences.push(div);
        
        if let Some(ascore) = get_alignment_score(&record) {
            alignment_scores.push(ascore as f64);
        }
        
        // Identity from cm tag
        if let Some(ident) = calculate_identity(&record) {
            identities.push(ident);
        }
        
        // Read quality (rq tag - PacBio)
        if let Some(rq) = get_read_quality(&record) {
            read_qualities.push(rq);
        }
        
        // Number of passes (np tag - PacBio CCS)
        if let Some(np) = get_num_passes(&record) {
            num_passes.push(np as f64);
        }

        // Full-length vs partial (heuristic: reads with 3+ introns and good quality = full-length)
        let n_introns = count_introns(&record);
        *intron_counts.entry(n_introns).or_default() += 1;
        
        // Full-length: has introns (spliced) AND good divergence AND reasonable identity
        let has_good_identity = calculate_identity(&record)
            .map(|i| i > 0.95)
            .unwrap_or(true);
        
        if n_introns >= 2 && div < 0.02 && has_good_identity {
            full_length += 1;
        } else {
            partial += 1;
        }

        // Coverage
        let start_bin = ((read_start - start).max(0) as usize / bin_size).min(num_bins - 1);
        let end_bin = ((read_end - start).max(0) as usize / bin_size).min(num_bins - 1);
        for bin in start_bin..=end_bin {
            coverage[bin] += 1;
        }
    }

    // Calculate statistics
    let mean_length = if !lengths.is_empty() {
        lengths.iter().sum::<f64>() / lengths.len() as f64
    } else { 0.0 };

    let length_std = if !lengths.is_empty() {
        let mean = mean_length;
        let variance = lengths.iter().map(|l| (l - mean).powi(2)).sum::<f64>() / lengths.len() as f64;
        variance.sqrt()
    } else { 0.0 };

    let length_cv = if mean_length > 0.0 { length_std / mean_length } else { 0.0 };

    let mean_divergence = if !divergences.is_empty() {
        divergences.iter().sum::<f64>() / divergences.len() as f64
    } else { 1.0 };

    let mean_alignment_score = if !alignment_scores.is_empty() {
        alignment_scores.iter().sum::<f64>() / alignment_scores.len() as f64
    } else { 0.0 };

    let mean_identity = if !identities.is_empty() {
        identities.iter().sum::<f64>() / identities.len() as f64
    } else { 0.0 };

    let mean_read_quality = if !read_qualities.is_empty() {
        read_qualities.iter().sum::<f64>() / read_qualities.len() as f64
    } else { 0.0 };

    let mean_num_passes = if !num_passes.is_empty() {
        num_passes.iter().sum::<f64>() / num_passes.len() as f64
    } else { 0.0 };

    let coverage_std = if !coverage.is_empty() {
        let mean_cov = coverage.iter().sum::<usize>() as f64 / coverage.len() as f64;
        let variance = coverage.iter()
            .map(|&c| (c as f64 - mean_cov).powi(2))
            .sum::<f64>() / coverage.len() as f64;
        variance.sqrt()
    } else { 0.0 };

    // Calculate mean introns
    let total_introns: usize = intron_counts.iter()
        .map(|(n, count)| n * count)
        .sum();
    let mean_introns = if total > 0 { 
        total_introns as f64 / total as f64 
    } else { 0.0 };

    // Strand consistency (fraction of dominant strand)
    let strand_consistency = if total > 0 {
        plus_strand.max(minus_strand) as f64 / total as f64
    } else { 0.0 };

    // Chimeric and secondary ratios (these are NORMAL for gene families)
    let total_observed = total + secondary_count + chimeric_count;
    let chimeric_ratio = if total_observed > 0 {
        chimeric_count as f64 / total_observed as f64
    } else { 0.0 };
    
    let secondary_ratio = if total_observed > 0 {
        secondary_count as f64 / total_observed as f64
    } else { 0.0 };

    Ok(AdvancedTranscriptFeatures {
        total_transcripts: total,
        mean_length,
        length_cv,
        strand_plus: plus_strand,
        strand_minus: minus_strand,
        coverage_std,
        mean_divergence,
        mean_alignment_score,
        mean_identity,
        n_full_length: full_length,
        n_partial: partial,
        intron_count_distribution: intron_counts,
        mean_introns,
        strand_consistency,
        mean_read_quality,
        mean_num_passes,
        chimeric_read_ratio: chimeric_ratio,
        secondary_alignment_ratio: secondary_ratio,
    })
}

/// Validate a locus against seed using advanced features
pub fn validate_advanced(
    bam_path: &str,
    seed_chrom: &str,
    seed_start: usize,
    seed_end: usize,
    target_chrom: &str,
    target_start: i64,
    target_end: i64,
    min_mapq: u8,
    max_divergence: f64,
) -> Result<(f64, AdvancedTranscriptFeatures, AdvancedTranscriptFeatures)> {
    let seed_features = analyze_advanced_features(
        bam_path, seed_chrom, seed_start as i64, seed_end as i64,
        min_mapq, max_divergence
    )?;
    
    let target_features = analyze_advanced_features(
        bam_path, target_chrom, target_start, target_end,
        min_mapq, max_divergence
    )?;
    
    let similarity = seed_features.similarity(&target_features);
    
    Ok((similarity, seed_features, target_features))
}

/// Filter loci by advanced Iso-Seq features
/// NOTE: MAPQ is intentionally NOT used - multi-mapping reads are essential for gene family detection
pub fn filter_by_advanced_features(
    features: &AdvancedTranscriptFeatures,
    min_transcripts: usize,
    max_divergence: f64,
    min_strand_consistency: f64,
) -> bool {
    if features.total_transcripts < min_transcripts {
        return false;
    }
    
    // Filter by divergence (but NOT by MAPQ!)
    if features.mean_divergence > max_divergence {
        return false;
    }
    
    // Genes should have consistent strand
    if features.strand_consistency < min_strand_consistency {
        return false;
    }
    
    // Additional Iso-Seq specific filters
    // High identity to reference is expected for true genes
    if features.mean_identity < 0.90 && features.mean_identity > 0.0 {
        return false;
    }
    
    true
}
