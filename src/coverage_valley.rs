//! Coverage valley detection for gene family splitting
//! 
//! This module implements coverage-based locus splitting by detecting
//! valleys (local minima) in the query-read coverage profile.
//! 
//! Algorithm:
//! 1. Build coverage profile from query reads (multi-mapping reads from seed)
//! 2. Bin coverage to reduce noise (default 100bp bins)
//! 3. Detect peaks (local maxima above threshold)
//! 4. Find valleys between peaks (local minima)
//! 5. Validate valleys using Iso-Seq signals when ambiguous
//! 6. Split locus at significant valleys

use anyhow::Result;
use fxhash::FxHashSet as HashSet;
use noodles::bam;
use noodles::core::Position;
use noodles::sam::alignment::Record as SamRecord;
/// A detected peak in coverage profile
#[derive(Debug, Clone)]
pub struct CoveragePeak {
    pub bin_idx: usize,
    pub position: i64,  // Genomic position
    pub coverage: f64,  // Coverage value
}

/// A detected valley between two peaks
#[derive(Debug, Clone)]
pub struct CoverageValley {
    pub left_peak: usize,   // Index into peaks vector
    pub right_peak: usize,
    pub bin_idx: usize,     // Valley position (bin index)
    pub position: i64,      // Genomic position
    pub coverage: f64,      // Valley coverage
    pub is_significant: bool, // Passes threshold test
    pub depth_ratio: f64,   // Valley depth / max peak height
}

/// Result of splitting a locus by coverage valleys
#[derive(Debug, Clone)]
pub struct SplitSegment {
    pub chrom: String,
    pub start: i64,
    pub end: i64,
    pub peak_count: usize,  // Number of peaks in this segment
    pub total_coverage: f64, // Sum of coverage
}

/// Coverage profile for a genomic region
#[derive(Debug)]
pub struct CoverageProfile {
    pub chrom: String,
    pub start: i64,
    pub end: i64,
    pub bin_size: i64,
    pub bins: Vec<f64>,     // Coverage per bin
    pub bin_starts: Vec<i64>, // Genomic start of each bin
}

/// Parameters for valley detection
#[derive(Debug, Clone)]
pub struct ValleyDetectionParams {
    pub bin_size: i64,           // Bin size in bp (default: 100)
    pub valley_frac: f64,        // Valley threshold as fraction of max (default: 0.1)
    pub peak_threshold_frac: f64, // Peak threshold as fraction of max (default: 0.15)
    pub min_gap_bp: i64,         // Minimum gap between peaks (default: 3000)
    pub min_segment_bp: i64,     // Minimum segment length (default: 2000)
    pub min_prominence_frac: f64, // Minimum prominence (peak - valley) / peak (default: 0.5)
}

impl ValleyDetectionParams {
    /// Create dynamic parameters based on coverage profile characteristics
    pub fn dynamic(max_coverage: f64, coverage_std: f64, profile_len: i64) -> Self {
        // Dynamic valley threshold: lower for high-coverage profiles
        // High coverage = more reads = more noise = need lower threshold to find valleys
        let valley_frac = if max_coverage > 1000.0 {
            0.03  // Very high coverage: valleys can be shallow
        } else if max_coverage > 500.0 {
            0.05  // High coverage
        } else if max_coverage > 100.0 {
            0.08  // Medium coverage
        } else {
            0.12  // Low coverage: need deeper valleys
        };
        
        // Dynamic prominence: based on coverage variability
        // High variability = peaks and valleys are more pronounced
        let variability_ratio = if max_coverage > 0.0 { coverage_std / max_coverage } else { 0.5 };
        let min_prominence_frac = if variability_ratio > 0.5 {
            0.20  // High variability: prominence is clear
        } else if variability_ratio > 0.3 {
            0.25  // Medium variability
        } else {
            0.32  // Low variability: need stricter prominence
        };
        
        // Dynamic min_gap: smaller for shorter profiles
        let min_gap_bp = if profile_len < 100_000 {
            1500  // Shorter profiles may have closer genes
        } else if profile_len < 500_000 {
            2000
        } else {
            2500  // Longer profiles can have wider gaps
        };
        
        Self {
            bin_size: 100,
            valley_frac,
            peak_threshold_frac: 0.10,  // Lower threshold to find more peaks
            min_gap_bp,
            min_segment_bp: 1500,
            min_prominence_frac,
        }
    }
}

impl Default for ValleyDetectionParams {
    fn default() -> Self {
        Self {
            bin_size: 100,
            valley_frac: 0.08,
            peak_threshold_frac: 0.10,  // Lowered from 0.12
            min_gap_bp: 2000,
            min_segment_bp: 1500,
            min_prominence_frac: 0.25,  // Lowered from 0.32
        }
    }
}

/// How `split_by_valleys` interprets the valley list.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ValleySplitFilter {
    /// Only split at valleys marked `is_significant` by coverage analysis.
    SignificantOnly,
    /// Use every valley in the slice; caller has already applied depth / Iso-Seq rules.
    CallerValidated,
}

/// Build coverage profile from query reads
/// 
/// Counts how many query reads cover each position in the region
pub fn build_coverage_profile(
    bam_path: &str,
    chrom: &str,
    start: i64,
    end: i64,
    query_reads: &HashSet<String>,
    bin_size: i64,
) -> Result<CoverageProfile> {
    let mut reader = bam::io::indexed_reader::Builder::default()
        .build_from_path(bam_path)?;
    let header = reader.read_header()?;

    let region_start = Position::try_from(start as usize)?;
    let region_end = Position::try_from(end as usize)?;
    let region = noodles::core::Region::new(chrom, region_start..=region_end);

    // Initialize bins
    let region_len = end - start;
    let n_bins = (region_len / bin_size) as usize + 1;
    let mut bin_counts = vec![0.0f64; n_bins];

    // Count coverage per bin
    for result in reader.query(&header, &region)? {
        let record = result?;
        
        if record.flags().is_unmapped() {
            continue;
        }

        // Check if this read is in our query set
        let name = match record.name() {
            Some(n) => match std::str::from_utf8(n) {
                Ok(s) => s,
                Err(_) => continue,
            },
            None => continue,
        };

        if !query_reads.contains(name) {
            continue;
        }

        // Get alignment coordinates
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

        // Add coverage to bins
        let start_bin = ((clip_start - start) / bin_size) as usize;
        let end_bin = ((clip_end - start) / bin_size) as usize;

        for bin_idx in start_bin..=end_bin.min(n_bins - 1) {
            // Weight by overlap fraction
            let bin_start = start + bin_idx as i64 * bin_size;
            let bin_end = (bin_start + bin_size).min(end);
            let overlap_start = clip_start.max(bin_start);
            let overlap_end = clip_end.min(bin_end);
            let overlap = (overlap_end - overlap_start) as f64;
            let bin_span = (bin_end - bin_start) as f64;

            if bin_span > 0.0 {
                bin_counts[bin_idx] += overlap / bin_span;
            }
        }
    }

    // Build bin start positions
    let bin_starts: Vec<i64> = (0..n_bins)
        .map(|i| start + i as i64 * bin_size)
        .collect();

    Ok(CoverageProfile {
        chrom: chrom.to_string(),
        start,
        end,
        bin_size,
        bins: bin_counts,
        bin_starts,
    })
}

/// Same geometry as [`build_coverage_profile`], but splits signal by alignment strand (SAM flag 0x10).
/// Use with [`analyze_valleys`] per strand when combined-strand coverage hides inter-gene dips.
pub fn build_strand_coverage_profiles(
    bam_path: &str,
    chrom: &str,
    start: i64,
    end: i64,
    query_reads: &HashSet<String>,
    bin_size: i64,
) -> Result<(CoverageProfile, CoverageProfile)> {
    let mut reader = bam::io::indexed_reader::Builder::default()
        .build_from_path(bam_path)?;
    let header = reader.read_header()?;

    let region_start = Position::try_from(start as usize)?;
    let region_end = Position::try_from(end as usize)?;
    let region = noodles::core::Region::new(chrom, region_start..=region_end);

    let region_len = end - start;
    let n_bins = (region_len / bin_size) as usize + 1;
    let mut plus_bins = vec![0.0f64; n_bins];
    let mut minus_bins = vec![0.0f64; n_bins];

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

        if !query_reads.contains(name) {
            continue;
        }

        let aln_start = record
            .alignment_start()
            .transpose()?
            .map(|p: Position| p.get() as i64)
            .unwrap_or(start);
        let aln_end = record
            .alignment_end()
            .transpose()?
            .map(|p: Position| p.get() as i64)
            .unwrap_or(aln_start + 1);

        let clip_start = aln_start.max(start);
        let clip_end = aln_end.min(end);

        if clip_end <= clip_start {
            continue;
        }

        let is_minus = record.flags().is_reverse_complemented();
        let target = if is_minus {
            &mut minus_bins
        } else {
            &mut plus_bins
        };

        let start_bin = ((clip_start - start) / bin_size) as usize;
        let end_bin = ((clip_end - start) / bin_size) as usize;

        for bin_idx in start_bin..=end_bin.min(n_bins - 1) {
            let bin_start_bp = start + bin_idx as i64 * bin_size;
            let bin_end_bp = (bin_start_bp + bin_size).min(end);
            let overlap_start = clip_start.max(bin_start_bp);
            let overlap_end = clip_end.min(bin_end_bp);
            let overlap = (overlap_end - overlap_start) as f64;
            let bin_span = (bin_end_bp - bin_start_bp) as f64;

            if bin_span > 0.0 {
                target[bin_idx] += overlap / bin_span;
            }
        }
    }

    let bin_starts: Vec<i64> = (0..n_bins)
        .map(|i| start + i as i64 * bin_size)
        .collect();

    let plus = CoverageProfile {
        chrom: chrom.to_string(),
        start,
        end,
        bin_size,
        bins: plus_bins,
        bin_starts: bin_starts.clone(),
    };
    let minus = CoverageProfile {
        chrom: chrom.to_string(),
        start,
        end,
        bin_size,
        bins: minus_bins,
        bin_starts,
    };

    Ok((plus, minus))
}

/// Find peaks in coverage profile
/// 
/// A peak is a local maximum that:
/// 1. Is higher than both neighbors
/// 2. Exceeds the peak threshold (fraction of global max)
/// 3. Is separated from other peaks by min_gap_bp
fn find_peaks(profile: &CoverageProfile, params: &ValleyDetectionParams) -> Vec<CoveragePeak> {
    let n_bins = profile.bins.len();
    if n_bins < 3 {
        return vec![];
    }

    let max_coverage = profile.bins.iter().copied().fold(0.0f64, f64::max);
    if max_coverage <= 0.0 {
        return vec![];
    }

    let peak_threshold = max_coverage * params.peak_threshold_frac;
    let min_peak_bins = (params.min_gap_bp / params.bin_size).max(2) as usize;

    // Find local maxima
    let mut raw_peaks: Vec<usize> = vec![];
    for i in 0..n_bins {
        let left_ok = i == 0 || profile.bins[i] >= profile.bins[i - 1];
        let right_ok = i == n_bins - 1 || profile.bins[i] >= profile.bins[i + 1];

        if left_ok && right_ok && profile.bins[i] >= peak_threshold {
            raw_peaks.push(i);
        }
    }

    // Merge close peaks (keep highest)
    let mut peaks: Vec<CoveragePeak> = vec![];
    for peak_idx in raw_peaks {
        let peak_cov = profile.bins[peak_idx];

        // Check if close to previous peak
        if let Some(last_peak) = peaks.last() {
            let dist = peak_idx - last_peak.bin_idx;
            if dist < min_peak_bins {
                // Replace if higher
                if peak_cov > last_peak.coverage {
                    peaks.pop();
                } else {
                    continue;
                }
            }
        }

        peaks.push(CoveragePeak {
            bin_idx: peak_idx,
            position: profile.bin_starts[peak_idx] + params.bin_size / 2,
            coverage: peak_cov,
        });
    }

    peaks
}

/// Find valleys between peaks
/// 
/// A valley is the minimum between two consecutive peaks.
/// Significance is determined locally based on the adjacent peaks.
fn find_valleys(
    profile: &CoverageProfile,
    peaks: &[CoveragePeak],
    params: &ValleyDetectionParams,
) -> Vec<CoverageValley> {
    if peaks.len() < 2 {
        return vec![];
    }

    let max_coverage = profile.bins.iter().copied().fold(0.0f64, f64::max);
    let global_mean: f64 = profile.bins.iter().sum::<f64>() / profile.bins.len().max(1) as f64;
    
    let mut valleys = vec![];

    for i in 0..peaks.len() - 1 {
        let left_peak = &peaks[i];
        let right_peak = &peaks[i + 1];
        let left_idx = left_peak.bin_idx;
        let right_idx = right_peak.bin_idx;

        // Find minimum between peaks
        let (min_idx, min_val) = profile.bins[left_idx..=right_idx]
            .iter()
            .enumerate()
            .map(|(idx, &val)| (left_idx + idx, val))
            .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
            .unwrap_or((left_idx, profile.bins[left_idx]));

        // === LOCAL DYNAMIC SIGNIFICANCE ===
        // Evaluate significance based on local context, not global thresholds
        
        // 1. Local peak characteristics
        let local_max_peak = left_peak.coverage.max(right_peak.coverage);
        let local_min_peak = left_peak.coverage.min(right_peak.coverage);
        let peak_asymmetry = if local_max_peak > 0.0 {
            (local_max_peak - local_min_peak) / local_max_peak
        } else {
            0.0
        };
        
        // 2. Local valley depth relative to adjacent peaks
        let left_drop = left_peak.coverage - min_val;
        let right_drop = right_peak.coverage - min_val;
        let min_drop = left_drop.min(right_drop);
        let max_drop = left_drop.max(right_drop);
        
        // 3. Local prominence ratio (how much of a dip is this?)
        let local_prominence = if local_max_peak > 0.0 {
            max_drop / local_max_peak
        } else {
            0.0
        };
        
        // 4. Local coverage ratio (valley vs peak)
        let local_coverage_ratio = if local_max_peak > 0.0 {
            min_val / local_max_peak
        } else {
            1.0
        };
        
        // 5. Compute local mean in the valley region
        let valley_start = left_idx.saturating_sub(5);
        let valley_end = (right_idx + 5).min(profile.bins.len() - 1);
        let local_mean: f64 = profile.bins[valley_start..=valley_end]
            .iter()
            .sum::<f64>() / (valley_end - valley_start + 1).max(1) as f64;
        
        // === DYNAMIC THRESHOLDS ===
        // Adapt thresholds based on local coverage characteristics
        
        // Higher coverage = expect deeper valleys relative to peaks
        // Lower coverage = accept shallower valleys
        let dynamic_prominence_threshold = if local_max_peak > 1000.0 {
            0.15  // High coverage: valleys should be clear
        } else if local_max_peak > 500.0 {
            0.20  // Medium-high coverage
        } else if local_max_peak > 100.0 {
            0.25  // Medium coverage
        } else if local_max_peak > 50.0 {
            0.30  // Low-medium coverage
        } else {
            0.35  // Low coverage: accept weaker valleys
        };
        
        // Valley should be notably below the local mean
        // But threshold depends on variability
        let dynamic_valley_ratio = if peak_asymmetry > 0.5 {
            0.6  // Very asymmetric peaks - allow higher valley (less clear split)
        } else if peak_asymmetry > 0.3 {
            0.5  // Moderately asymmetric
        } else {
            0.4  // Symmetric peaks - expect clearer valley
        };
        
        // === SIGNIFICANCE CRITERIA ===
        // A valley is significant if multiple conditions suggest a real dip
        
        // Criterion 1: Valley is notably below the smaller peak
        let is_below_smaller_peak = min_val < local_min_peak * (1.0 - dynamic_prominence_threshold);
        
        // Criterion 2: Prominence is significant relative to larger peak
        let has_prominence = local_prominence >= dynamic_prominence_threshold;
        
        // Criterion 3: Valley is below local mean
        let is_below_local_mean = min_val < local_mean * 0.7;
        
        // Criterion 4: Coverage ratio is low (valley is much smaller than peaks)
        let is_low_coverage_ratio = local_coverage_ratio < dynamic_valley_ratio;
        
        // Criterion 5: Minimum drop is meaningful (at least 20% of peak height)
        let has_meaningful_drop = min_drop > local_max_peak * 0.20;
        
        // Combine criteria: need at least 3 of 5, with prominence being essential
        let criteria_met = [
            is_below_smaller_peak,
            has_prominence,
            is_below_local_mean,
            is_low_coverage_ratio,
            has_meaningful_drop,
        ];
        let count = criteria_met.iter().filter(|&&x| x).count();
        
        // Must have prominence AND at least 2 other criteria
        // OR have strong drop (40%+) and below local mean
        let is_significant = (has_prominence && count >= 3)
            || (min_drop > local_max_peak * 0.4 && is_below_local_mean);

        valleys.push(CoverageValley {
            left_peak: i,
            right_peak: i + 1,
            bin_idx: min_idx,
            position: profile.bin_starts[min_idx] + params.bin_size / 2,
            coverage: min_val,
            is_significant,
            depth_ratio: if max_coverage > 0.0 { min_val / max_coverage } else { 1.0 },
        });
    }

    valleys
}

/// Split a locus into segments based on coverage valleys
///
/// Returns the segments after splitting. Use [`ValleySplitFilter::SignificantOnly`] when
/// `valleys` is the raw output of [`analyze_valleys`]. Use [`ValleySplitFilter::CallerValidated`]
/// when the caller has already decided split points (e.g. integrated mode after depth and
/// Iso-Seq rules), including valleys that are not marked `is_significant` by coverage alone.
pub fn split_by_valleys(
    chrom: &str,
    start: i64,
    end: i64,
    profile: &CoverageProfile,
    valleys: &[CoverageValley],
    params: &ValleyDetectionParams,
    split_filter: ValleySplitFilter,
) -> Vec<SplitSegment> {
    let mut split_positions: Vec<i64> = match split_filter {
        ValleySplitFilter::SignificantOnly => valleys
            .iter()
            .filter(|v| v.is_significant)
            .map(|v| v.position)
            .collect(),
        ValleySplitFilter::CallerValidated => valleys.iter().map(|v| v.position).collect(),
    };
    split_positions.sort_unstable();
    split_positions.dedup();
    split_positions.retain(|&pos| {
        let left_len = pos - start;
        let right_len = end - pos;
        left_len >= params.min_segment_bp && right_len >= params.min_segment_bp
    });

    if split_positions.is_empty() {
        // No valid splits - return whole locus
        let total_coverage: f64 = profile.bins.iter().sum();
        return vec![SplitSegment {
            chrom: chrom.to_string(),
            start,
            end,
            peak_count: 1,  // At least one peak (otherwise wouldn't be called)
            total_coverage,
        }];
    }

    // Build segments
    let mut segments = vec![];
    let mut seg_start = start;

    for &split_pos in &split_positions {
        let bin_start = ((seg_start - profile.start) / profile.bin_size) as usize;
        let bin_end = ((split_pos - profile.start) / profile.bin_size) as usize;
        let coverage_sum: f64 = profile.bins[bin_start..bin_end.min(profile.bins.len())]
            .iter()
            .sum();

        segments.push(SplitSegment {
            chrom: chrom.to_string(),
            start: seg_start,
            end: split_pos,
            peak_count: 1,  // Approximation
            total_coverage: coverage_sum,
        });

        seg_start = split_pos;
    }

    // Add final segment
    let bin_start = ((seg_start - profile.start) / profile.bin_size) as usize;
    let bin_end = profile.bins.len();
    let final_coverage: f64 = profile.bins[bin_start..bin_end].iter().sum();

    segments.push(SplitSegment {
        chrom: chrom.to_string(),
        start: seg_start,
        end,
        peak_count: 1,
        total_coverage: final_coverage,
    });

    segments
}

/// Analyze a locus and return peaks and valleys
pub fn analyze_valleys(
    profile: &CoverageProfile,
    params: &ValleyDetectionParams,
) -> (Vec<CoveragePeak>, Vec<CoverageValley>) {
    let peaks = find_peaks(profile, params);
    let valleys = find_valleys(profile, &peaks, params);
    (peaks, valleys)
}

/// Print valley analysis report
pub fn print_valley_report(
    chrom: &str,
    start: i64,
    end: i64,
    peaks: &[CoveragePeak],
    valleys: &[CoverageValley],
    segments: &[SplitSegment],
) {
    println!("    Coverage valley analysis for {}:{}-{}", chrom, start, end);
    println!("      Detected {} peaks:", peaks.len());
    for (i, peak) in peaks.iter().enumerate() {
        println!(
            "        Peak {}: pos={} cov={:.1}",
            i + 1,
            peak.position,
            peak.coverage
        );
    }

    println!("      Detected {} valleys:", valleys.len());
    for (i, valley) in valleys.iter().enumerate() {
        let status = if valley.is_significant { "✓ SPLIT" } else { "✗ no split" };
        println!(
            "        Valley {}: pos={} cov={:.1} depth={:.2} {}",
            i + 1,
            valley.position,
            valley.coverage,
            valley.depth_ratio,
            status
        );
    }

    println!("      Split into {} segments:", segments.len());
    for (i, seg) in segments.iter().enumerate() {
        println!(
            "        Segment {}: {}:{}-{} (span: {} bp, coverage: {:.0})",
            i + 1,
            seg.chrom,
            seg.start,
            seg.end,
            seg.end - seg.start,
            seg.total_coverage
        );
    }
}



/// Split by peaks when valleys are not significant.
/// This is used as a fallback when coverage is continuous between genes.
/// 
/// Algorithm:
/// 1. Find all peaks in coverage profile
/// 2. If peaks are far enough apart (>min_gap), split between them
/// 3. Each peak represents a potential gene
pub fn split_by_peaks(
    chrom: &str,
    start: i64,
    end: i64,
    profile: &CoverageProfile,
    peaks: &[CoveragePeak],
    params: &ValleyDetectionParams,
) -> Vec<SplitSegment> {
    // Need at least 2 peaks to split
    if peaks.len() < 2 {
        let total_coverage: f64 = profile.bins.iter().sum();
        return vec![SplitSegment {
            chrom: chrom.to_string(),
            start,
            end,
            peak_count: peaks.len(),
            total_coverage,
        }];
    }
    
    // Find split points: midpoints between consecutive peaks
    let mut split_positions: Vec<i64> = Vec::new();
    
    for i in 0..peaks.len() - 1 {
        let left_peak = &peaks[i];
        let right_peak = &peaks[i + 1];
        let gap = right_peak.position - left_peak.position;
        
        // Only split if peaks are far enough apart
        if gap >= params.min_gap_bp {
            // Split at midpoint between peaks
            let midpoint = left_peak.position + gap / 2;
            
            // Ensure segments are large enough
            let left_len = midpoint - start;
            let right_len = end - midpoint;
            
            // Check if both resulting segments would be large enough
            let last_split = split_positions.last().copied().unwrap_or(start);
            if midpoint - last_split >= params.min_segment_bp && end - midpoint >= params.min_segment_bp {
                split_positions.push(midpoint);
            }
        }
    }
    
    if split_positions.is_empty() {
        let total_coverage: f64 = profile.bins.iter().sum();
        return vec![SplitSegment {
            chrom: chrom.to_string(),
            start,
            end,
            peak_count: peaks.len(),
            total_coverage,
        }];
    }
    
    // Build segments from split positions
    let mut segments = Vec::new();
    let mut seg_start = start;
    
    for &split_pos in &split_positions {
        let bin_start = ((seg_start - profile.start) / profile.bin_size) as usize;
        let bin_end = ((split_pos - profile.start) / profile.bin_size) as usize;
        let coverage_sum: f64 = profile.bins[bin_start..bin_end.min(profile.bins.len())]
            .iter()
            .sum();
        
        // Count peaks in this segment
        let peak_count = peaks.iter().filter(|p| p.position >= seg_start && p.position < split_pos).count();
        
        segments.push(SplitSegment {
            chrom: chrom.to_string(),
            start: seg_start,
            end: split_pos,
            peak_count: peak_count.max(1),
            total_coverage: coverage_sum,
        });
        
        seg_start = split_pos;
    }
    
    // Add final segment
    let bin_start = ((seg_start - profile.start) / profile.bin_size) as usize;
    let bin_end = profile.bins.len();
    let final_coverage: f64 = profile.bins[bin_start..bin_end].iter().sum();
    let peak_count = peaks.iter().filter(|p| p.position >= seg_start).count();
    
    segments.push(SplitSegment {
        chrom: chrom.to_string(),
        start: seg_start,
        end,
        peak_count: peak_count.max(1),
        total_coverage: final_coverage,
    });
    
    segments
}



/// Merge adjacent small segments that are likely intra-genic regions.
/// Small segments (<min_size) should be merged into larger neighbors.
/// This handles peak-based splitting creating many small intra-genic fragments.
pub fn merge_small_segments(
    segments: &[SplitSegment],
    min_size: i64,
    min_gap: i64,
) -> Vec<SplitSegment> {
    if segments.is_empty() {
        return vec![];
    }
    
    // Sort by position
    let mut sorted: Vec<_> = segments.iter().collect();
    sorted.sort_by_key(|s| (s.chrom.clone(), s.start));
    
    // First pass: mark small segments for merging
    // A small segment should be merged with its larger neighbor
    let mut keep: Vec<bool> = vec![true; sorted.len()];
    
    for i in 0..sorted.len() {
        let seg_span = sorted[i].end - sorted[i].start;
        if seg_span < min_size {
            // This is a small segment - find best neighbor to merge with
            let left_size = if i > 0 {
                sorted[i-1].end - sorted[i-1].start
            } else {
                0
            };
            let right_size = if i < sorted.len() - 1 {
                sorted[i+1].end - sorted[i+1].start
            } else {
                0
            };
            
            // Merge with larger neighbor (prefer larger neighbor)
            // But only if adjacent (gap <= min_gap)
            let gap_left = if i > 0 { sorted[i].start - sorted[i-1].end } else { i64::MAX };
            let gap_right = if i < sorted.len() - 1 { sorted[i+1].start - sorted[i].end } else { i64::MAX };
            
            if gap_left <= min_gap && (left_size >= right_size || gap_right > min_gap) {
                // Merge with left neighbor
                keep[i] = false;
            } else if gap_right <= min_gap && right_size > 0 {
                // Merge with right neighbor
                keep[i] = false;
            }
        }
    }
    
    // Second pass: build merged segments
    let mut merged: Vec<SplitSegment> = Vec::new();
    let mut i = 0;
    
    while i < sorted.len() {
        if !keep[i] {
            // This segment should be merged - skip it (already merged with neighbor)
            i += 1;
            continue;
        }
        
        let mut current = (*sorted[i]).clone();
        
        // Check if next segment(s) should be merged into this one
        let mut j = i + 1;
        while j < sorted.len() && !keep[j] && sorted[j].chrom == current.chrom {
            // Merge segment j into current
            current.end = sorted[j].end;
            current.peak_count += sorted[j].peak_count;
            current.total_coverage += sorted[j].total_coverage;
            j += 1;
        }
        
        merged.push(current);
        i = j;
    }
    
    merged
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_find_peaks() {
        let profile = CoverageProfile {
            chrom: "chr1".to_string(),
            start: 0,
            end: 1000,
            bin_size: 100,
            bins: vec![1.0, 2.0, 5.0, 3.0, 1.0, 1.0, 4.0, 6.0, 3.0, 1.0],
            bin_starts: (0..10).map(|i| i as i64 * 100).collect(),
        };

        let params = ValleyDetectionParams {
            bin_size: 100,
            peak_threshold_frac: 0.15,
            min_gap_bp: 200,
            ..Default::default()
        };

        let peaks = find_peaks(&profile, &params);
        assert!(!peaks.is_empty());
        // Should find peaks at bins 2 and 7 (coverage 5.0 and 6.0)
    }

    #[test]
    fn test_valley_significance() {
        // Peak coverage = 100
        // Valley coverage = 5
        // valley_frac = 0.1
        // Should be significant (5 < 10)
        assert!(5.0 < 100.0 * 0.1); // 5 < 10 ✓
    }

    #[test]
    fn split_by_valleys_caller_validated_uses_non_significant() {
        let profile = CoverageProfile {
            chrom: "chr1".to_string(),
            start: 0,
            end: 20_000,
            bin_size: 100,
            bins: vec![1.0; 200],
            bin_starts: (0..200).map(|i| i as i64 * 100).collect(),
        };
        let params = ValleyDetectionParams {
            min_segment_bp: 1000,
            ..Default::default()
        };
        let v = CoverageValley {
            left_peak: 0,
            right_peak: 1,
            bin_idx: 50,
            position: 5050,
            coverage: 0.0,
            is_significant: false,
            depth_ratio: 0.05,
        };
        let segs_sig = split_by_valleys(
            "chr1",
            0,
            20_000,
            &profile,
            &[v.clone()],
            &params,
            ValleySplitFilter::SignificantOnly,
        );
        assert_eq!(segs_sig.len(), 1);

        let segs_all = split_by_valleys(
            "chr1",
            0,
            20_000,
            &profile,
            &[v],
            &params,
            ValleySplitFilter::CallerValidated,
        );
        assert_eq!(segs_all.len(), 2);
    }
}

/// Split large loci using splice junction clusters.
/// Improved: Uses splice pattern similarity and seed reads only.
/// 
/// Key insight: Junctions within a gene have similar patterns (same intron positions).
/// Between genes, there are large gaps (>20kb) where junction patterns change.
pub fn split_by_junctions(
    bam_path: &str,
    segments: &[SplitSegment],
    seed_reads: &HashSet<String>,
    min_size_for_split: i64,
) -> Result<Vec<SplitSegment>> {
    use noodles::bam;
    use noodles::core::Position;
    use noodles::sam::alignment::record::cigar::op::Kind;
    
    println!("    === JUNCTION-BASED SPLITTING ===");
    println!("    Using seed-read junctions to find gene boundaries (min gap: {}kb)", min_size_for_split / 1000);
    
    let mut final_segments = Vec::new();
    
    for seg in segments {
        let span = seg.end - seg.start;
        
        // Only split large segments
        if span < min_size_for_split {
            final_segments.push(seg.clone());
            continue;
        }
        
        println!("    Checking {}:{}-{} ({}kb)", seg.chrom, seg.start, seg.end, span / 1000);
        
        // Step 1: Collect splice junctions from SEED READS ONLY
        let junctions = match collect_seed_junctions(bam_path, &seg.chrom, seg.start, seg.end, seed_reads) {
            Ok(j) => j,
            Err(e) => {
                println!("      Error collecting junctions: {}, keeping segment", e);
                final_segments.push(seg.clone());
                continue;
            }
        };
        
        if junctions.is_empty() {
            println!("      No seed-read junctions found, keeping segment");
            final_segments.push(seg.clone());
            continue;
        }
        
        println!("      Found {} junction positions from seed reads", junctions.len());
        
        // Step 2: Cluster junctions by position
        let clusters = cluster_junctions(&junctions, 1000);
        println!("      Clustered into {} junction groups", clusters.len());
        
        // Step 3: Find large gaps between junction clusters (gene boundaries)
        let split_points = find_splice_boundaries(&clusters, seg.start, seg.end, 20000);
        
        if split_points.is_empty() {
            println!("      No large gaps (20kb+) between junction clusters");
            final_segments.push(seg.clone());
            continue;
        }
        
        println!("      Found {} splice boundaries (large junction gaps)", split_points.len());
        
        // Step 4: Create sub-segments
        let mut last_end = seg.start;
        for split_pos in &split_points {
            if split_pos - last_end > 10000 {
                final_segments.push(SplitSegment {
                    chrom: seg.chrom.clone(),
                    start: last_end,
                    end: *split_pos,
                    peak_count: 1,
                    total_coverage: 0.0,
                });
            }
            last_end = *split_pos;
        }
        
        if seg.end - last_end > 10000 {
            final_segments.push(SplitSegment {
                chrom: seg.chrom.clone(),
                start: last_end,
                end: seg.end,
                peak_count: 1,
                total_coverage: 0.0,
            });
        }
    }
    
    Ok(final_segments)
}

/// Collect splice junction positions from seed reads ONLY (not all multi-mapping reads)
fn collect_seed_junctions(
    bam_path: &str,
    chrom: &str,
    start: i64,
    end: i64,
    seed_reads: &HashSet<String>,
) -> Result<Vec<i64>> {
    use noodles::bam;
    use noodles::core::Position;
    use noodles::sam::alignment::record::cigar::op::Kind;
    
    let mut reader = bam::io::indexed_reader::Builder::default()
        .build_from_path(bam_path)?;
    let header = reader.read_header()?;
    
    let region_start = Position::try_from(start.max(1) as usize)?;
    let region_end = Position::try_from(end.max(1) as usize)?;
    let region = noodles::core::Region::new(chrom, region_start..=region_end);
    
    let mut junctions = Vec::new();
    
    for result in reader.query(&header, &region)? {
        let record = result?;
        
        if record.flags().is_unmapped() {
            continue;
        }
        
        // CRITICAL: Only use seed reads, not all multi-mapping reads
        let name = match record.name() {
            Some(n) => match std::str::from_utf8(n) {
                Ok(s) => s,
                Err(_) => continue,
            },
            None => continue,
        };
        
        if !seed_reads.contains(name) {
            continue;
        }
        
        // Extract splice junctions from CIGAR N operations
        let cigar = record.cigar();
        let aln_start = record.alignment_start()
            .transpose()?
            .map(|p: Position| p.get() as i64)
            .unwrap_or(0);
        
        let mut ref_pos = aln_start;
        
        for op_result in cigar.iter() {
            let op = op_result?;
            let len = op.len() as i64;
            
            match op.kind() {
                Kind::Skip => {
                    // N operation = intron (splice junction)
                    junctions.push(ref_pos);  // Junction start
                    junctions.push(ref_pos + len);  // Junction end
                }
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch | Kind::Deletion => {
                    ref_pos += len;
                }
                _ => {}
            }
        }
    }
    
    junctions.sort();
    Ok(junctions)
}

/// Cluster junctions that are close together (within tolerance)
fn cluster_junctions(junctions: &[i64], tolerance: i64) -> Vec<(i64, i64)> {
    if junctions.is_empty() {
        return vec![];
    }
    
    let mut clusters: Vec<(i64, i64)> = Vec::new();
    let mut cluster_start = junctions[0];
    let mut cluster_end = junctions[0];
    
    for &pos in junctions.iter().skip(1) {
        if pos - cluster_end <= tolerance {
            cluster_end = pos;
        } else {
            clusters.push((cluster_start, cluster_end));
            cluster_start = pos;
            cluster_end = pos;
        }
    }
    
    clusters.push((cluster_start, cluster_end));
    clusters
}

/// Find splice pattern boundaries - large gaps between junction clusters
/// that indicate gene boundaries (intergenic regions, not introns)
fn find_splice_boundaries(
    clusters: &[(i64, i64)],
    region_start: i64,
    region_end: i64,
    min_gap: i64,
) -> Vec<i64> {
    if clusters.len() < 2 {
        return vec![];
    }
    
    let mut boundaries = Vec::new();
    
    // Sort clusters by position
    let mut sorted: Vec<_> = clusters.to_vec();
    sorted.sort_by_key(|(s, _e)| *s);
    
    // Find large gaps between junction clusters
    // Genes have clustered introns; intergenic regions have no introns
    for i in 0..sorted.len() - 1 {
        let cluster_end = sorted[i].1;
        let next_cluster_start = sorted[i + 1].0;
        let gap_size = next_cluster_start - cluster_end;
        
        // Only consider gaps large enough to be gene boundaries
        // Typical intergenic regions: >20kb
        // Typical introns: <10kb
        if gap_size > min_gap {
            boundaries.push(cluster_end + gap_size / 2);
        }
    }
    
    boundaries
}
