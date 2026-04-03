//! Boundary refinement for gene family detection
//!
//! Post-processing module to trim over-extended detected regions to
//! more accurate gene boundaries using:
//! 1. Iso-Seq transcript boundaries
//! 2. Jaccard-based flank trimming

use anyhow::Result;
use fxhash::FxHashSet as HashSet;
use noodles::bam;
use noodles::core::Position;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::Record as SamRecord;

/// How `--refine-boundaries` derives coordinates from BAM alignments.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, clap::ValueEnum)]
pub enum IsoBoundaryMode {
    /// Percentiles on raw alignment start/end (can mix opposite-strand noise).
    Legacy,
    /// Majority strand, robust percentiles on 5'/3' termini, fused with CIGAR exon outer bounds.
    #[default]
    StrandExon,
}

// Strand-specific coverage valley detection implemented below

/// Refined locus with trimmed boundaries
#[derive(Debug, Clone)]
pub struct RefinedLocus {
    pub chrom: String,
    pub start: i64,
    pub end: i64,
    /// Original detected region
    pub original_start: i64,
    pub original_end: i64,
    /// Read names supporting this locus
    pub reads: HashSet<String>,
    /// Jaccard similarity after trimming
    pub jaccard_with_seed: f64,
    /// Number of Iso-Seq transcripts
    pub n_transcripts: usize,
    /// Transcript density (transcripts per kb)
    pub transcript_density: f64,
    /// Trim amount upstream
    pub trim_5p: i64,
    /// Trim amount downstream
    pub trim_3p: i64,
    /// Refinement method used
    pub refinement_method: RefinementMethod,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum RefinementMethod {
    IsoSeqBoundaries,
    JaccardTrim,
    CoverageDrop,
    None,
}

impl std::fmt::Display for RefinementMethod {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RefinementMethod::IsoSeqBoundaries => write!(f, "IsoSeq"),
            RefinementMethod::JaccardTrim => write!(f, "Jaccard"),
            RefinementMethod::CoverageDrop => write!(f, "Coverage"),
            RefinementMethod::None => write!(f, "None"),
        }
    }
}

/// Parameters for boundary refinement
#[derive(Debug, Clone)]
pub struct RefinementParams {
    /// Window size for Jaccard calculation (bp)
    pub window_size: i64,
    /// Step size for sliding window (bp)
    pub step_size: i64,
    /// Minimum Jaccard to keep a window
    pub min_jaccard_window: f64,
    /// Jaccard drop threshold for boundary detection
    pub jaccard_drop_threshold: f64,
    /// Minimum transcript density for Iso-Seq boundaries (transcripts per kb)
    pub min_transcript_density: f64,
    /// Maximum trim fraction (don't trim more than this)
    pub max_trim_fraction: f64,
    /// Minimum region size after trimming
    pub min_region_size: i64,
    /// Keep all significant segments instead of just the best one
    pub keep_all_segments: bool,
    /// Minimum Jaccard for a segment to be considered significant
    pub min_segment_jaccard: f64,
    /// Iso-Seq coordinate model for transcript envelope
    pub iso_boundary_mode: IsoBoundaryMode,
}

impl Default for RefinementParams {
    fn default() -> Self {
        Self {
            window_size: 1000,
            step_size: 500,
            min_jaccard_window: 0.01,
            jaccard_drop_threshold: 0.5,
            min_transcript_density: 0.5,
            max_trim_fraction: 0.8,
            min_region_size: 5000,
            keep_all_segments: true,  // Default to true for multi-seed compatibility
            min_segment_jaccard: 0.05,
            iso_boundary_mode: IsoBoundaryMode::default(),
        }
    }
}

/// Per-alignment geometry for boundary inference
#[derive(Debug, Clone)]
struct TranscriptGeom {
    /// Genomic 5' end (TSS side in reference coordinates)
    t5: i64,
    /// Genomic 3' end (TES side in reference coordinates)
    t3: i64,
    strand: char,
    /// First spliced block start on reference (equals alignment start if unspliced)
    exon_lo: i64,
    /// Last spliced block end on reference (exclusive end coordinate, same convention as aln_end)
    exon_hi: i64,
}

/// Exon blocks from CIGAR (reference intervals); returns outer shell [first block start, last block end).
fn exon_outer_from_cigar(record: &bam::Record) -> Option<(i64, i64)> {
    if record.flags().is_unmapped() {
        return None;
    }
    let aln_start = record
        .alignment_start()
        .transpose()
        .ok()??
        .get() as i64;
    let mut ref_pos = aln_start;
    let mut block_start: Option<i64> = None;
    let mut blocks: Vec<(i64, i64)> = Vec::new();

    for result in record.cigar().iter() {
        let op = result.ok()?;
        let len = op.len() as i64;
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch | Kind::Deletion => {
                if block_start.is_none() {
                    block_start = Some(ref_pos);
                }
                ref_pos += len;
            }
            Kind::Skip => {
                if let Some(bs) = block_start.take() {
                    if ref_pos > bs {
                        blocks.push((bs, ref_pos));
                    }
                }
                ref_pos += len;
            }
            _ => {}
        }
    }
    if let Some(bs) = block_start {
        if ref_pos > bs {
            blocks.push((bs, ref_pos));
        }
    }
    if blocks.is_empty() {
        return None;
    }
    let ex_lo = blocks.first().map(|b| b.0)?;
    let ex_hi = blocks.last().map(|b| b.1)?;
    Some((ex_lo, ex_hi))
}

/// Collect transcript geometry (strand, 5'/3' termini, exon shell) from a region
fn collect_transcript_geometry(
    bam_path: &str,
    chrom: &str,
    start: i64,
    end: i64,
) -> Result<Vec<TranscriptGeom>> {
    let mut reader = bam::io::indexed_reader::Builder::default()
        .build_from_path(bam_path)?;
    let header = reader.read_header()?;

    let region_start = Position::try_from(start.max(1) as usize)?;
    let region_end = Position::try_from(end.max(1) as usize)?;
    let region = noodles::core::Region::new(chrom, region_start..=region_end);

    let mut out = Vec::new();

    for result in reader.query(&header, &region)? {
        let record = result?;

        if record.flags().is_unmapped() {
            continue;
        }

        let aln_start = record
            .alignment_start()
            .transpose()?
            .map(|p: Position| p.get() as i64)
            .unwrap_or(0);
        let aln_end = record
            .alignment_end()
            .transpose()?
            .map(|p: Position| p.get() as i64)
            .unwrap_or(aln_start + 1);

        let strand = if record.flags().is_reverse_complemented() {
            '-'
        } else {
            '+'
        };

        let (t5, t3) = if strand == '+' {
            (aln_start, aln_end)
        } else {
            (aln_end, aln_start)
        };

        let (exon_lo, exon_hi) = exon_outer_from_cigar(&record).unwrap_or((aln_start, aln_end));

        out.push(TranscriptGeom {
            t5,
            t3,
            strand,
            exon_lo,
            exon_hi,
        });
    }

    Ok(out)
}

fn percentile_i64(values: &mut [i64], q: f64) -> Option<i64> {
    if values.is_empty() {
        return None;
    }
    let q = q.clamp(0.0, 1.0);
    values.sort_unstable();
    let n = values.len();
    if n == 1 {
        return Some(values[0]);
    }
    let pos = (n - 1) as f64 * q;
    let lo = pos.floor() as usize;
    let hi = (pos.ceil() as usize).min(n - 1);
    if lo >= hi {
        return Some(values[lo]);
    }
    let w = pos - lo as f64;
    Some(((values[lo] as f64) * (1.0 - w) + (values[hi] as f64) * w).round() as i64)
}

/// Strand-aware bounds: robust envelope on dominant strand, fused with CIGAR exon outer shell.
fn find_strand_exon_iso_bounds(geoms: &[TranscriptGeom], region_start: i64, region_end: i64) -> Option<(i64, i64)> {
    if geoms.is_empty() {
        return None;
    }
    let plus_n = geoms.iter().filter(|g| g.strand == '+').count();
    let minus_n = geoms.len() - plus_n;
    let dom = if plus_n >= minus_n { '+' } else { '-' };

    let mut los: Vec<i64> = Vec::new();
    let mut his: Vec<i64> = Vec::new();
    let mut ex_los: Vec<i64> = Vec::new();
    let mut ex_his: Vec<i64> = Vec::new();

    for g in geoms.iter().filter(|g| g.strand == dom) {
        let lo = g.t5.min(g.t3);
        let hi = g.t5.max(g.t3);
        los.push(lo);
        his.push(hi);
        ex_los.push(g.exon_lo.min(g.exon_hi));
        ex_his.push(g.exon_lo.max(g.exon_hi));
    }

    if los.is_empty() {
        for g in geoms {
            let lo = g.t5.min(g.t3);
            let hi = g.t5.max(g.t3);
            los.push(lo);
            his.push(hi);
            ex_los.push(g.exon_lo.min(g.exon_hi));
            ex_his.push(g.exon_lo.max(g.exon_hi));
        }
    }

    let p25_lo = percentile_i64(&mut los, 0.25)?;
    let p75_hi = percentile_i64(&mut his, 0.75)?;
    let p25_ex_lo = percentile_i64(&mut ex_los, 0.25)?;
    let p75_ex_hi = percentile_i64(&mut ex_his, 0.75)?;

    let tx_start = p25_lo.max(p25_ex_lo).clamp(region_start, region_end);
    let tx_end = p75_hi.min(p75_ex_hi).clamp(region_start, region_end);

    if tx_end > tx_start {
        Some((tx_start, tx_end))
    } else {
        None
    }
}

/// Legacy envelope: p10 on alignment starts, p90 on alignment ends (all reads, strand-mixed).
fn find_legacy_percentile_bounds(
    geoms: &[TranscriptGeom],
    region_start: i64,
    region_end: i64,
) -> (i64, i64) {
    if geoms.is_empty() {
        return (region_start, region_end);
    }
    let mut starts: Vec<i64> = geoms.iter().map(|g| g.t5.min(g.t3)).collect();
    let mut ends: Vec<i64> = geoms.iter().map(|g| g.t5.max(g.t3)).collect();
    starts.sort_unstable();
    ends.sort_unstable();
    let n = starts.len();
    let p10_idx = (n / 10).min(n.saturating_sub(1));
    let p90_idx = ((n * 9) / 10).min(n.saturating_sub(1));
    (starts[p10_idx], ends[p90_idx])
}

/// Calculate Jaccard similarity for reads in a window
fn calculate_window_jaccard(
    window_reads: &HashSet<String>,
    seed_reads: &HashSet<String>,
) -> f64 {
    if window_reads.is_empty() || seed_reads.is_empty() {
        return 0.0;
    }

    let intersection: usize = window_reads.iter().filter(|x| seed_reads.contains(*x)).count();
    let union = window_reads.len() + seed_reads.len() - intersection;

    if union > 0 {
        intersection as f64 / union as f64
    } else {
        0.0
    }
}

/// Find reads in a specific window
fn collect_reads_in_window(
    bam_path: &str,
    chrom: &str,
    window_start: i64,
    window_end: i64,
) -> Result<HashSet<String>> {
    let mut reader = bam::io::indexed_reader::Builder::default()
        .build_from_path(bam_path)?;
    let header = reader.read_header()?;

    let region_start = Position::try_from(window_start.max(1) as usize)?;
    let region_end = Position::try_from(window_end.max(1) as usize)?;
    let region = noodles::core::Region::new(chrom, region_start..=region_end);

    let mut reads: HashSet<String> = HashSet::default();

    for result in reader.query(&header, &region)? {
        let record = result?;

        if record.flags().is_unmapped() {
            continue;
        }

        let name = match record.name() {
            Some(n) => match std::str::from_utf8(n) {
                Ok(s) => s,
                Err(_) => continue,
            },
            None => continue,
        };

        reads.insert(name.to_string());
    }

    Ok(reads)
}

/// Jaccard-based flank trimming
/// 
/// Slide windows from edges inward and find where Jaccard drops significantly
fn jaccard_based_trimming(
    bam_path: &str,
    chrom: &str,
    region_start: i64,
    region_end: i64,
    seed_reads: &HashSet<String>,
    params: &RefinementParams,
) -> Result<(i64, i64)> {
    let span = region_end - region_start;
    
    // Don't trim if region is already small
    if span < params.min_region_size * 2 {
        return Ok((region_start, region_end));
    }

    // Scan from 5' end (upstream)
    let mut new_start = region_start;
    let mut prev_jaccard = 1.0;

    for pos in (region_start..region_end).step_by(params.step_size as usize) {
        let window_end = (pos + params.window_size).min(region_end);
        if window_end - pos < params.window_size / 2 {
            break;
        }

        let window_reads = collect_reads_in_window(bam_path, chrom, pos, window_end)?;
        let jaccard = calculate_window_jaccard(&window_reads, seed_reads);

        // Check for significant drop
        if jaccard < params.min_jaccard_window && 
           prev_jaccard > jaccard * (1.0 + params.jaccard_drop_threshold) {
            new_start = pos;
            break;
        }

        prev_jaccard = jaccard.max(params.min_jaccard_window);
    }

    // Scan from 3' end (downstream)
    let mut new_end = region_end;
    prev_jaccard = 1.0;

    for pos in ((region_start..region_end).rev()).step_by(params.step_size as usize) {
        let window_start = (pos - params.window_size).max(region_start);
        if pos - window_start < params.window_size / 2 {
            break;
        }

        let window_reads = collect_reads_in_window(bam_path, chrom, window_start, pos)?;
        let jaccard = calculate_window_jaccard(&window_reads, seed_reads);

        if jaccard < params.min_jaccard_window && 
           prev_jaccard > jaccard * (1.0 + params.jaccard_drop_threshold) {
            new_end = pos;
            break;
        }

        prev_jaccard = jaccard.max(params.min_jaccard_window);
    }

    // Ensure we don't trim too much
    let max_trim = (span as f64 * params.max_trim_fraction) as i64;
    let trim_5p = new_start - region_start;
    let trim_3p = region_end - new_end;

    if trim_5p > max_trim {
        new_start = region_start + max_trim;
    }
    if trim_3p > max_trim {
        new_end = region_end - max_trim;
    }

    // Ensure minimum size
    if new_end - new_start < params.min_region_size {
        let center = (region_start + region_end) / 2;
        new_start = center - params.min_region_size / 2;
        new_end = center + params.min_region_size / 2;
    }

    Ok((new_start, new_end))
}

/// Build strand-specific coverage profile
fn build_strand_coverage_profile(
    bam_path: &str,
    chrom: &str,
    start: i64,
    end: i64,
    bin_size: i64,
) -> Result<(Vec<f64>, Vec<f64>)> {
    let mut reader = bam::io::indexed_reader::Builder::default()
        .build_from_path(bam_path)?;
    let header = reader.read_header()?;

    let region_start = Position::try_from(start.max(1) as usize)?;
    let region_end = Position::try_from(end.max(1) as usize)?;
    let region = noodles::core::Region::new(chrom, region_start..=region_end);

    let n_bins = ((end - start) / bin_size) as usize + 1;
    let mut plus_coverage = vec![0.0f64; n_bins];
    let mut minus_coverage = vec![0.0f64; n_bins];

    for result in reader.query(&header, &region)? {
        let record = result?;
        
        if record.flags().is_unmapped() {
            continue;
        }

        let aln_start = record.alignment_start()
            .transpose()?
            .map(|p: Position| p.get() as i64)
            .unwrap_or(start);
        let aln_end = record.alignment_end()
            .transpose()?
            .map(|p: Position| p.get() as i64)
            .unwrap_or(aln_start + 1);

        // Clip to region
        let clip_start = aln_start.max(start);
        let clip_end = aln_end.min(end);

        if clip_end <= clip_start {
            continue;
        }

        // Determine strand
        let is_minus = record.flags().is_reverse_complemented();

        // Add coverage to bins
        let start_bin = ((clip_start - start) / bin_size) as usize;
        let end_bin = ((clip_end - start) / bin_size) as usize;

        for bin_idx in start_bin..=end_bin.min(n_bins - 1) {
            if is_minus {
                minus_coverage[bin_idx] += 1.0;
            } else {
                plus_coverage[bin_idx] += 1.0;
            }
        }
    }

    Ok((plus_coverage, minus_coverage))
}

/// Find valleys in a coverage profile
fn find_valleys_in_profile(
    coverage: &[f64],
    min_prominence_frac: f64,
) -> Vec<usize> {
    let max_cov = coverage.iter().copied().fold(0.0f64, f64::max);
    if max_cov <= 0.0 {
        return vec![];
    }

    let mut valleys = Vec::new();

    for i in 1..coverage.len() - 1 {
        let prev = coverage[i - 1];
        let curr = coverage[i];
        let next = coverage[i + 1];

        // Local minimum
        if curr < prev && curr < next {
            // Check prominence
            let left_peak = (0..i).rev()
                .map(|j| coverage[j])
                .find(|&c| c > curr)
                .unwrap_or(curr);
            let right_peak = (i + 1..coverage.len())
                .map(|j| coverage[j])
                .find(|&c| c > curr)
                .unwrap_or(curr);
            let peak_height = left_peak.min(right_peak);
            
            if peak_height > 0.0 {
                let prominence = (peak_height - curr) / peak_height;
                if prominence >= min_prominence_frac {
                    valleys.push(i);
                }
            }
        }
    }

    valleys
}

/// Calculate strand ratio (plus/minus) for a region
fn get_strand_ratio(
    bam_path: &str,
    chrom: &str,
    start: i64,
    end: i64,
) -> Result<f64> {
    let mut reader = bam::io::indexed_reader::Builder::default()
        .build_from_path(bam_path)?;
    let header = reader.read_header()?;

    let region_start = Position::try_from(start.max(1) as usize)?;
    let region_end = Position::try_from(end.max(1) as usize)?;
    let region = noodles::core::Region::new(chrom, region_start..=region_end);

    let mut plus = 0usize;
    let mut minus = 0usize;

    for result in reader.query(&header, &region)? {
        let record = result?;
        
        if record.flags().is_unmapped() {
            continue;
        }

        if record.flags().is_reverse_complemented() {
            minus += 1;
        } else {
            plus += 1;
        }
    }

    if minus == 0 {
        Ok(99.0) // All plus
    } else {
        Ok(plus as f64 / minus as f64)
    }
}

/// Split a locus using hierarchical valley detection
/// 
/// Tier 1: Major valleys (high prominence >= 0.5) - gene boundaries
/// Tier 2: Subvalleys (lower prominence >= 0.15) - only if no major valleys
fn split_by_coverage_valleys(
    bam_path: &str,
    chrom: &str,
    start: i64,
    end: i64,
    seed_reads: &HashSet<String>,
) -> Result<Vec<(i64, i64)>> {
    let bin_size = 500i64;
    let min_segment_bp = 10000i64;
    
    // Tier 1: Major valleys (gene boundaries) - high prominence
    let major_prominence = 0.50;
    // Tier 2: Subvalleys (introns/internal) - lower prominence
    let sub_prominence = 0.15;

    // Build strand-specific coverage
    let (plus_cov, minus_cov) = build_strand_coverage_profile(
        bam_path, chrom, start, end, bin_size
    )?;

    // Try MAJOR valleys first (gene boundaries)
    let plus_major = find_valleys_in_profile(&plus_cov, major_prominence);
    let minus_major = find_valleys_in_profile(&minus_cov, major_prominence);
    
    let mut all_valley_bins: Vec<usize> = plus_major.clone();
    all_valley_bins.extend(minus_major.clone());
    all_valley_bins.sort_unstable();
    all_valley_bins.dedup();

    let valley_type = if all_valley_bins.is_empty() {
        // No major valleys - try subvalleys
        let plus_sub = find_valleys_in_profile(&plus_cov, sub_prominence);
        let minus_sub = find_valleys_in_profile(&minus_cov, sub_prominence);
        
        // Filter to valleys NOT in major list
        let plus_sub_filtered: Vec<usize> = plus_sub.into_iter()
            .filter(|v| !plus_major.contains(v))
            .collect();
        let minus_sub_filtered: Vec<usize> = minus_sub.into_iter()
            .filter(|v| !minus_major.contains(v))
            .collect();
        
        all_valley_bins = plus_sub_filtered;
        all_valley_bins.extend(minus_sub_filtered);
        all_valley_bins.sort_unstable();
        all_valley_bins.dedup();
        
        "subvalleys"
    } else {
        "MAJOR valleys"
    };

    println!("      Using {}: {} total (major: {}+, {}-)", 
             valley_type,
             all_valley_bins.len(),
             plus_major.len(),
             minus_major.len());

    if all_valley_bins.is_empty() {
        return Ok(vec![(start, end)]);
    }

    // Split into segments
    let mut segments = Vec::new();
    let mut seg_start = start;

    for &bin_idx in &all_valley_bins {
        let valley_pos = start + (bin_idx as i64) * bin_size;
        
        // Check segment size
        if valley_pos - seg_start >= min_segment_bp {
            segments.push((seg_start, valley_pos));
            seg_start = valley_pos;
        }
    }

    // Add final segment
    if end - seg_start >= min_segment_bp {
        segments.push((seg_start, end));
    } else if !segments.is_empty() {
        // Extend last segment
        let last_idx = segments.len() - 1;
        segments[last_idx] = (segments[last_idx].0, end);
    } else {
        segments.push((start, end));
    }

    println!("      Coverage valleys: {} → {} segments", 
             all_valley_bins.len(), segments.len());

    // Merge small segments to prevent oversplitting
    let merged = merge_small_segments(&segments, 15000);
    if merged.len() != segments.len() {
        println!("      Merged small segments: {} → {} segments", 
                 segments.len(), merged.len());
    }

    Ok(merged)
}

/// Merge segments that are too small (likely coverage artifacts)
/// Simple greedy merge: consecutive small segments get merged together
fn merge_small_segments(
    segments: &[(i64, i64)],
    min_size: i64,
) -> Vec<(i64, i64)> {
    if segments.len() < 2 {
        return segments.to_vec();
    }

    let mut merged: Vec<(i64, i64)> = Vec::new();
    let mut current_start = segments[0].0;
    let mut current_end = segments[0].1;
    let mut current_size = current_end - current_start;

    for i in 1..segments.len() {
        let (next_start, next_end) = segments[i];
        let next_size = next_end - next_start;

        if current_size < min_size || next_size < min_size {
            // Merge: extend current to include next
            current_end = next_end;
            current_size = current_end - current_start;
        } else {
            // Both large enough, finalize current and start new
            merged.push((current_start, current_end));
            current_start = next_start;
            current_end = next_end;
            current_size = next_size;
        }
    }

    // Add final segment
    merged.push((current_start, current_end));

    merged
}

/// Refine a single locus using multiple methods
/// Returns a vector of refined segments (usually 1, but can be multiple for tandem arrays)
pub fn refine_locus(
    bam_path: &str,
    chrom: &str,
    start: i64,
    end: i64,
    seed_reads: &HashSet<String>,
    params: &RefinementParams,
) -> Result<Vec<RefinedLocus>> {
    let original_start = start;
    let original_end = end;
    let span = end - start;
    let mut refined_loci = Vec::new();

    // Method 1: For large regions, try coverage valley splitting
    if span > 100000 {
        println!("      Large region detected ({} kb), trying coverage valley splitting...", span / 1000);
        match split_by_coverage_valleys(bam_path, chrom, start, end, seed_reads) {
            Ok(segments) if segments.len() > 1 => {
                // Calculate Jaccard for each segment
                let mut segment_scores: Vec<((i64, i64), f64, HashSet<String>)> = Vec::new();
                
                for (seg_start, seg_end) in &segments {
                    let seg_reads = collect_reads_in_window(bam_path, chrom, *seg_start, *seg_end)?;
                    let jaccard = calculate_window_jaccard(&seg_reads, seed_reads);
                    segment_scores.push(((*seg_start, *seg_end), jaccard, seg_reads));
                }

                if params.keep_all_segments {
                    // Keep all segments that meet the minimum Jaccard threshold
                    println!("      Keeping all segments with Jaccard >= {:.3}", params.min_segment_jaccard);
                    
                    for ((seg_start, seg_end), seg_jaccard, seg_reads) in &segment_scores {
                        if *seg_jaccard >= params.min_segment_jaccard {
                            // Apply further refinement to each segment
                            let refined = refine_single_segment(
                                bam_path, chrom, *seg_start, *seg_end, 
                                original_start, original_end,
                                seg_reads, seed_reads, params
                            )?;
                            refined_loci.push(refined);
                            println!("      Kept segment: {}-{} (Jaccard: {:.3})", 
                                     seg_start, seg_end, seg_jaccard);
                        } else {
                            println!("      Dropped segment: {}-{} (Jaccard: {:.3} < {:.3})", 
                                     seg_start, seg_end, seg_jaccard, params.min_segment_jaccard);
                        }
                    }
                    
                    if !refined_loci.is_empty() {
                        return Ok(refined_loci);
                    }
                    // If no segments passed threshold, fall through to use best segment
                }
                
                // Default: use only the best segment
                let mut best_segment = segments[0];
                let mut best_jaccard = 0.0;
                let mut best_reads = HashSet::default();

                for ((seg_start, seg_end), jaccard, seg_reads) in &segment_scores {
                    if *jaccard > best_jaccard {
                        best_jaccard = *jaccard;
                        best_segment = (*seg_start, *seg_end);
                        best_reads = seg_reads.clone();
                    }
                }

                // Now apply other refinement methods to this segment
                let seg_start = best_segment.0;
                let seg_end = best_segment.1;
                println!("      Selected segment: {}-{} (Jaccard: {:.3})", seg_start, seg_end, best_jaccard);
                
                let strand_ratio = get_strand_ratio(bam_path, chrom, seg_start, seg_end)?;
                println!("      Strand ratio: {:.1}", strand_ratio);
                println!("      Final segment: {}-{} ({} kb)", seg_start, seg_end, (seg_end-seg_start)/1000);
                
                let refined = refine_single_segment(
                    bam_path, chrom, seg_start, seg_end,
                    original_start, original_end,
                    &best_reads, seed_reads, params
                )?;
                return Ok(vec![refined]);
            }
            _ => {}
        }
    }

    // No splitting needed or splitting failed - refine the whole locus
    let locus_reads = collect_reads_in_window(bam_path, chrom, start, end)?;
    let refined = refine_single_segment(
        bam_path, chrom, start, end,
        original_start, original_end,
        &locus_reads, seed_reads, params
    )?;
    Ok(vec![refined])
}

/// Refine a single segment (helper function)
fn refine_single_segment(
    bam_path: &str,
    chrom: &str,
    start: i64,
    end: i64,
    original_start: i64,
    original_end: i64,
    locus_reads: &HashSet<String>,
    seed_reads: &HashSet<String>,
    params: &RefinementParams,
) -> Result<RefinedLocus> {
    let span = original_end - original_start;

    // Method 1: Iso-Seq transcript boundaries (legacy vs strand + exon shell)
    let transcripts = collect_transcript_geometry(bam_path, chrom, start, end)?;
    let (tx_start, tx_end) = match params.iso_boundary_mode {
        IsoBoundaryMode::Legacy => find_legacy_percentile_bounds(&transcripts, start, end),
        IsoBoundaryMode::StrandExon => find_strand_exon_iso_bounds(&transcripts, start, end)
            .unwrap_or_else(|| find_legacy_percentile_bounds(&transcripts, start, end)),
    };

    let tx_span = tx_end - tx_start;
    let tx_density = transcripts.len() as f64 / (span as f64 / 1000.0);

    let (refined_start, refined_end, method) = if tx_density >= params.min_transcript_density &&
                                                  tx_span >= params.min_region_size &&
                                                  tx_span <= (span as f64 * 0.8) as i64 {
        // Use transcript boundaries if density is good and size is reasonable
        (tx_start, tx_end, RefinementMethod::IsoSeqBoundaries)
    } else {
        // Method 2: Jaccard-based trimming
        let (jac_start, jac_end) = jaccard_based_trimming(
            bam_path, chrom, start, end, seed_reads, params
        )?;
        
        if (jac_start > start + params.step_size) || 
           (jac_end < end - params.step_size) {
            (jac_start, jac_end, RefinementMethod::JaccardTrim)
        } else {
            (start, end, RefinementMethod::None)
        }
    };

    // Collect reads for refined region
    let refined_reads = collect_reads_in_window(bam_path, chrom, refined_start, refined_end)?;
    let jaccard = calculate_window_jaccard(&refined_reads, seed_reads);

    Ok(RefinedLocus {
        chrom: chrom.to_string(),
        start: refined_start,
        end: refined_end,
        original_start,
        original_end,
        reads: refined_reads,
        jaccard_with_seed: jaccard,
        n_transcripts: transcripts.len(),
        transcript_density: tx_density,
        trim_5p: refined_start - original_start,
        trim_3p: original_end - refined_end,
        refinement_method: method,
    })
}

/// Refine multiple loci
pub fn refine_loci(
    bam_path: &str,
    loci: &[(String, i64, i64, HashSet<String>)],
    seed_reads: &HashSet<String>,
    params: &RefinementParams,
) -> Result<Vec<RefinedLocus>> {
    let mut refined = Vec::new();
    let mut total_orig_span = 0i64;
    let mut total_new_span = 0i64;

    for (i, (chrom, start, end, _locus_reads)) in loci.iter().enumerate() {
        println!("  Refining locus {}: {}:{}-{}...", i + 1, chrom, start, end);
        let locus_orig_span = end - start;
        total_orig_span += locus_orig_span;

        match refine_locus(bam_path, chrom, *start, *end, seed_reads, params) {
            Ok(rl_vec) => {
                let locus_new_span: i64 = rl_vec.iter().map(|rl| rl.end - rl.start).sum();
                total_new_span += locus_new_span;
                
                let reduction = (1.0 - locus_new_span as f64 / locus_orig_span as f64) * 100.0;
                
                if rl_vec.len() == 1 {
                    let rl = &rl_vec[0];
                    println!("    {} → {} bp ({:.1}% reduction, method: {})",
                             locus_orig_span, locus_new_span, reduction, rl.refinement_method);
                    if rl.trim_5p > 0 || rl.trim_3p > 0 {
                        println!("    Trimmed: {}bp 5', {}bp 3'", rl.trim_5p, rl.trim_3p);
                    }
                } else {
                    println!("    Split into {} segments ({} → {} bp, {:.1}% reduction)",
                             rl_vec.len(), locus_orig_span, locus_new_span, reduction);
                    for (j, rl) in rl_vec.iter().enumerate() {
                        println!("      Segment {}: {}-{} ({} bp, Jaccard: {:.3})",
                                 j + 1, rl.start, rl.end, rl.end - rl.start, rl.jaccard_with_seed);
                    }
                }

                refined.extend(rl_vec);
            }
            Err(e) => {
                eprintln!("    Warning: Refinement failed: {}", e);
            }
        }
    }
    
    // Print overall summary
    if total_orig_span > 0 {
        let avg_reduction = (1.0 - total_new_span as f64 / total_orig_span as f64) * 100.0;
        println!("\nTotal span reduction: {:.1}% ({} bp → {} bp)", 
                 avg_reduction, total_orig_span, total_new_span);
    }

    Ok(refined)
}

/// Write refined loci to BED file
pub fn write_refined_bed(
    output_path: &str,
    loci: &[RefinedLocus],
    seed_reads: &HashSet<String>,
) -> Result<()> {
    use std::fs::File;
    use std::io::Write;

    let mut file = File::create(output_path)?;

    writeln!(file, "# Refined Gene Family Detection Results")?;
    writeln!(file, "# Columns: chrom, start, end, name, jaccard, transcripts, tx_density, trim_5p, trim_3p, method, original_span")?;
    writeln!(file, "#")?;

    for (i, locus) in loci.iter().enumerate() {
        let name = format!("refined_locus_{}", i + 1);
        let orig_span = locus.original_end - locus.original_start;
        let new_span = locus.end - locus.start;

        writeln!(
            file,
            "{}\t{}\t{}\t{}\t{:.4}\t{}\t{:.2}\t{}\t{}\t{}\t{}",
            locus.chrom,
            locus.start,
            locus.end,
            name,
            locus.jaccard_with_seed,
            locus.n_transcripts,
            locus.transcript_density,
            locus.trim_5p,
            locus.trim_3p,
            locus.refinement_method,
            orig_span
        )?;
    }

    println!("\nWrote {} refined loci to {}", loci.len(), output_path);
    
    // Print summary
    let total_orig: i64 = loci.iter().map(|l| l.original_end - l.original_start).sum();
    let total_new: i64 = loci.iter().map(|l| l.end - l.start).sum();
    let avg_reduction = (1.0 - total_new as f64 / total_orig as f64) * 100.0;
    
    println!("Total span reduction: {:.1}% ({} bp → {} bp)", 
             avg_reduction, total_orig, total_new);

    Ok(())
}
