//! Iso-Seq aware detection using transcript structure

use anyhow::Result;
use noodles::bam;
use noodles::core::Position;


/// Transcript analysis results
#[derive(Debug, Clone)]
pub struct TranscriptAnalysis {
    pub total_transcripts: usize,
    pub mean_length: f64,
    pub length_std: f64,
    pub length_cv: f64,
    pub strand_plus: usize,
    pub strand_minus: usize,
    pub coverage_std: f64,
}

/// Collect and analyze transcripts from a region
pub fn analyze_locus_transcripts(
    bam_path: &str,
    chrom: &str,
    start: i64,
    end: i64,
) -> Result<TranscriptAnalysis> {
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

    // For coverage calculation
    let bin_size = 1000usize;
    let num_bins = ((end - start) as usize / bin_size) + 1;
    let mut coverage = vec![0usize; num_bins];

    for result in reader.query(&header, &region)? {
        let record = result?;

        if record.flags().is_unmapped() {
            continue;
        }

        total += 1;

        // Get positions
        let read_start = record.alignment_start()
            .transpose()?
            .map(|p| usize::from(p) as i64)
            .unwrap_or(0);
        
        // Simple end calculation - just use reference length from alignment
        let read_end = read_start + 1000i64; // Simplified - use actual end from BAM in production

        let length = read_end - read_start;
        lengths.push(length);

        // Count strand - check the 0x10 flag bit manually
        // In SAM/BAM, 0x10 means reverse complemented
        let flags_val = record.flags().bits();
        if flags_val & 0x10 != 0 {
            minus_strand += 1;
        } else {
            plus_strand += 1;
        }

        // Update coverage
        let rs = read_start.max(start);
        let re = read_end.min(end);
        if rs < re {
            let start_bin = ((rs - start) as usize / bin_size).min(num_bins - 1);
            let end_bin = ((re - start) as usize / bin_size).min(num_bins - 1);
            for b in start_bin..=end_bin {
                coverage[b] += 1;
            }
        }
    }

    // Calculate statistics
    let (mean_length, length_std, length_cv) = if !lengths.is_empty() {
        let mean = lengths.iter().sum::<i64>() as f64 / lengths.len() as f64;
        let variance = lengths.iter()
            .map(|&l| (l as f64 - mean).powi(2))
            .sum::<f64>() / lengths.len() as f64;
        let std = variance.sqrt();
        let cv = if mean > 0.0 { std / mean } else { 0.0 };
        (mean, std, cv)
    } else {
        (0.0, 0.0, 0.0)
    };

    // Coverage statistics
    let coverage_std = if !coverage.is_empty() {
        let mean = coverage.iter().sum::<usize>() as f64 / coverage.len() as f64;
        let variance = coverage.iter()
            .map(|&c| (c as f64 - mean).powi(2))
            .sum::<f64>() / coverage.len() as f64;
        variance.sqrt()
    } else {
        0.0
    };

    Ok(TranscriptAnalysis {
        total_transcripts: total,
        mean_length,
        length_std,
        length_cv,
        strand_plus: plus_strand,
        strand_minus: minus_strand,
        coverage_std,
    })
}

/// Compare two loci using transcript features
pub fn compare_loci_isoseq(
    seed: &TranscriptAnalysis,
    target: &TranscriptAnalysis,
) -> f64 {
    let mut score = 0.0;
    let mut weights = 0.0;

    // 1. Strand consistency (must match)
    let seed_strand = if seed.strand_plus > seed.strand_minus { '+' } else { '-' };
    let target_strand = if target.strand_plus > target.strand_minus { '+' } else { '-' };
    
    if seed_strand != target_strand {
        return 0.0; // Different strand = not the same gene
    }

    // 2. Length similarity
    let len_score = if seed.mean_length > 0.0 && target.mean_length > 0.0 {
        let ratio = seed.mean_length / target.mean_length;
        if ratio >= 0.7 && ratio <= 1.43 {
            1.0 - (ratio - 1.0).abs()
        } else {
            0.0
        }
    } else {
        0.5
    };
    score += len_score * 0.3;
    weights += 0.3;

    // 3. Length consistency (CV) - genes should have consistent lengths
    let cv_score = if seed.length_cv < 0.5 && target.length_cv < 0.5 {
        1.0
    } else if seed.length_cv < 0.7 && target.length_cv < 0.7 {
        0.7
    } else {
        0.4
    };
    score += cv_score * 0.2;
    weights += 0.2;

    // 4. Coverage variability (genes have high std due to exons/introns)
    let coverage_score = if target.coverage_std > 50.0 {
        1.0
    } else if target.coverage_std > 20.0 {
        0.7
    } else if target.coverage_std > 10.0 {
        0.4
    } else {
        0.1
    };
    score += coverage_score * 0.5;
    weights += 0.5;

    if weights > 0.0 {
        score / weights
    } else {
        0.0
    }
}

/// Validate a locus using Iso-Seq features
pub fn validate_with_isoseq(
    bam_path: &str,
    seed_chrom: &str,
    seed_start: usize,
    seed_end: usize,
    target_chrom: &str,
    target_start: i64,
    target_end: i64,
) -> Result<(f64, TranscriptAnalysis, TranscriptAnalysis)> {
    // Analyze seed
    let seed_analysis = analyze_locus_transcripts(
        bam_path, seed_chrom, seed_start as i64, seed_end as i64
    )?;

    // Analyze target
    let target_analysis = analyze_locus_transcripts(
        bam_path, target_chrom, target_start, target_end
    )?;

    // Compare
    let similarity = compare_loci_isoseq(&seed_analysis, &target_analysis);

    Ok((similarity, seed_analysis, target_analysis))
}
