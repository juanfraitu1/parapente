//! Split merged read clouds using Iso-Seq properties
//!
//! Strategy:
//! 1. Find valleys in read coverage (gaps between genes)
//! 2. Use transcript boundary clusters (5' and 3' ends)
//! 3. Check strand consistency changes
//! 4. Look for intron structure changes

use anyhow::Result;
use fxhash::FxHashSet as HashSet;
use noodles::bam;
use noodles::core::Position;
use noodles::sam::alignment::Record as SamRecord;
use std::collections::HashMap;

/// A split point suggestion with confidence score
#[derive(Debug, Clone)]
struct SplitPoint {
    position: i64,
    confidence: f64,
    evidence_type: String,
}

/// Coverage bin for valley detection
#[derive(Debug, Clone)]
struct CoverageBin {
    start: i64,
    end: i64,
    count: usize,
}

/// Calculate coverage profile with given bin size
fn calculate_coverage_profile(
    bam_path: &str,
    chrom: &str,
    region_start: i64,
    region_end: i64,
    bin_size: i64,
    target_reads: Option<&HashSet<String>>,
) -> Result<Vec<CoverageBin>> {
    let mut reader = bam::io::indexed_reader::Builder::default()
        .build_from_path(bam_path)?;
    let header = reader.read_header()?;
    
    let start_pos = Position::try_from(region_start.max(1) as usize)?;
    let end_pos = Position::try_from(region_end.max(1) as usize)?;
    let region = noodles::core::Region::new(chrom, start_pos..=end_pos);
    
    let num_bins = ((region_end - region_start) as usize / bin_size as usize) + 1;
    let mut bin_counts = vec![0usize; num_bins];
    
    for result in reader.query(&header, &region)? {
        let record = result?;
        
        if record.flags().is_unmapped() {
            continue;
        }
        
        // Filter by read names if specified
        if let Some(targets) = target_reads {
            if let Some(name) = record.name() {
                if let Ok(name_str) = std::str::from_utf8(name) {
                    if !targets.contains(name_str) {
                        continue;
                    }
                }
            }
        }
        
        let start = record.alignment_start()
            .transpose()?
            .map(|p| p.get() as i64)
            .unwrap_or(0);
        let end = record.alignment_end()
            .transpose()?
            .map(|p| p.get() as i64)
            .unwrap_or(start + 1);
        
        // Update bins
        let start_bin = ((start - region_start).max(0) as usize / bin_size as usize)
            .min(num_bins - 1);
        let end_bin = ((end - region_start).max(0) as usize / bin_size as usize)
            .min(num_bins - 1);
        
        for bin_idx in start_bin..=end_bin {
            bin_counts[bin_idx] += 1;
        }
    }
    
    // Convert to bins
    let bins: Vec<CoverageBin> = bin_counts
        .iter()
        .enumerate()
        .map(|(i, &count)| CoverageBin {
            start: region_start + (i as i64 * bin_size),
            end: region_start + ((i + 1) as i64 * bin_size),
            count,
        })
        .collect();
    
    Ok(bins)
}

/// Find valleys in coverage (potential split points)
fn find_coverage_valleys(
    bins: &[CoverageBin],
    window_size: usize,
    min_valley_depth: f64,
) -> Vec<SplitPoint> {
    if bins.len() < window_size * 2 + 1 {
        return vec![];
    }
    
    let mut valleys = Vec::new();
    
    for i in window_size..bins.len() - window_size {
        let center_count = bins[i].count as f64;
        
        // Calculate mean coverage in surrounding windows
        let left_sum: usize = bins[i - window_size..i].iter().map(|b| b.count).sum();
        let right_sum: usize = bins[i + 1..i + 1 + window_size].iter().map(|b| b.count).sum();
        
        let left_mean = left_sum as f64 / window_size as f64;
        let right_mean = right_sum as f64 / window_size as f64;
        let surrounding_mean = (left_mean + right_mean) / 2.0;
        
        // Valley if center is much lower than surroundings
        if surrounding_mean > 0.0 {
            let depth = 1.0 - (center_count / surrounding_mean);
            if depth >= min_valley_depth && center_count < 5.0 {
                valleys.push(SplitPoint {
                    position: (bins[i].start + bins[i].end) / 2,
                    confidence: depth,
                    evidence_type: "coverage_valley".to_string(),
                });
            }
        }
    }
    
    valleys
}

/// Transcript boundary information
#[derive(Debug, Clone)]
struct TranscriptBoundary {
    position: i64,
    is_start: bool, // true = 5' end (TSS), false = 3' end (TES)
    strand: char,
}

/// Collect transcript boundaries from region
fn collect_transcript_boundaries(
    bam_path: &str,
    chrom: &str,
    region_start: i64,
    region_end: i64,
    target_reads: Option<&HashSet<String>>,
) -> Result<Vec<TranscriptBoundary>> {
    let mut reader = bam::io::indexed_reader::Builder::default()
        .build_from_path(bam_path)?;
    let header = reader.read_header()?;
    
    let start_pos = Position::try_from(region_start.max(1) as usize)?;
    let end_pos = Position::try_from(region_end.max(1) as usize)?;
    let region = noodles::core::Region::new(chrom, start_pos..=end_pos);
    
    let mut boundaries = Vec::new();
    
    for result in reader.query(&header, &region)? {
        let record = result?;
        
        if record.flags().is_unmapped() {
            continue;
        }
        
        // Filter by read names if specified
        if let Some(targets) = target_reads {
            if let Some(name) = record.name() {
                if let Ok(name_str) = std::str::from_utf8(name) {
                    if !targets.contains(name_str) {
                        continue;
                    }
                }
            }
        }
        
        let start = record.alignment_start()
            .transpose()?
            .map(|p| p.get() as i64)
            .unwrap_or(0);
        let end = record.alignment_end()
            .transpose()?
            .map(|p| p.get() as i64)
            .unwrap_or(start + 1);
        
        // Determine strand
        let flags_val = record.flags().bits();
        let strand = if flags_val & 0x10 != 0 { '-' } else { '+' };
        
        // 5' end (TSS) - start for + strand, end for - strand
        // 3' end (TES) - end for + strand, start for - strand
        if strand == '+' {
            boundaries.push(TranscriptBoundary {
                position: start,
                is_start: true,
                strand,
            });
            boundaries.push(TranscriptBoundary {
                position: end,
                is_start: false,
                strand,
            });
        } else {
            boundaries.push(TranscriptBoundary {
                position: end,
                is_start: true,
                strand,
            });
            boundaries.push(TranscriptBoundary {
                position: start,
                is_start: false,
                strand,
            });
        }
    }
    
    Ok(boundaries)
}

/// Cluster positions within one strand and end type (avoids mixing TSS with TES in one cluster).
fn cluster_positions_of_subset(
    boundaries: &[TranscriptBoundary],
    strand: char,
    want_start: bool,
    max_distance: i64,
    min_boundaries: usize,
    evidence_label: &str,
) -> Vec<SplitPoint> {
    let mut subset: Vec<&TranscriptBoundary> = boundaries
        .iter()
        .filter(|b| b.strand == strand && b.is_start == want_start)
        .collect();
    if subset.is_empty() {
        return vec![];
    }
    subset.sort_by_key(|b| b.position);

    let mut clusters: Vec<Vec<&TranscriptBoundary>> = Vec::new();
    for boundary in subset {
        let mut added = false;
        for cluster in &mut clusters {
            if let Some(last) = cluster.last() {
                if (boundary.position - last.position).abs() <= max_distance {
                    cluster.push(boundary);
                    added = true;
                    break;
                }
            }
        }
        if !added {
            clusters.push(vec![boundary]);
        }
    }

    clusters
        .into_iter()
        .filter(|c| c.len() >= min_boundaries)
        .map(|c| {
            let avg_pos = c.iter().map(|b| b.position).sum::<i64>() / c.len() as i64;
            SplitPoint {
                position: avg_pos,
                confidence: (c.len() as f64 / 10.0).min(1.0),
                evidence_type: evidence_label.to_string(),
            }
        })
        .collect()
}

/// Find clusters of transcript boundaries (strand-resolved TSS vs TES)
fn find_boundary_clusters(
    boundaries: &[TranscriptBoundary],
    max_distance: i64,
    min_boundaries: usize,
) -> Vec<SplitPoint> {
    if boundaries.is_empty() {
        return vec![];
    }

    let mut out = Vec::new();
    out.extend(cluster_positions_of_subset(
        boundaries,
        '+',
        true,
        max_distance,
        min_boundaries,
        "TSS_plus",
    ));
    out.extend(cluster_positions_of_subset(
        boundaries,
        '+',
        false,
        max_distance,
        min_boundaries,
        "TES_plus",
    ));
    out.extend(cluster_positions_of_subset(
        boundaries,
        '-',
        true,
        max_distance,
        min_boundaries,
        "TSS_minus",
    ));
    out.extend(cluster_positions_of_subset(
        boundaries,
        '-',
        false,
        max_distance,
        min_boundaries,
        "TES_minus",
    ));

    out.sort_by_key(|s| s.position);
    out
}

/// A sub-locus identified by splitting
#[derive(Debug, Clone)]
pub struct SubLocus {
    pub chrom: String,
    pub start: i64,
    pub end: i64,
    pub split_evidence: Vec<String>,
}

/// Split a merged cloud using multiple Iso-Seq signals
/// Split a merged cloud using multiple Iso-Seq signals (seed-aware)
pub fn split_cloud_by_isoseq_signals(
    bam_path: &str,
    chrom: &str,
    cloud_start: i64,
    cloud_end: i64,
    target_reads: Option<&HashSet<String>>,
    expected_genes: Option<usize>,
    seed_region: Option<(String, i64, i64)>, // (chrom, start, end) - seed location
) -> Result<Vec<SubLocus>> {
    println!("    Analyzing cloud {}:{}-{} for split points...", chrom, cloud_start, cloud_end);
    
    // Check if seed is within this cloud
    let seed_in_cloud = seed_region.as_ref().map(|(s_chrom, s_start, s_end)| {
        chrom == s_chrom && cloud_start < *s_end && cloud_end > *s_start
    }).unwrap_or(false);
    
    let (seed_start, seed_end) = if seed_in_cloud {
        let (_, s, e) = seed_region.unwrap();
        println!("      Seed region {}-{} is inside this cloud", s, e);
        (Some(s), Some(e))
    } else {
        (None, None)
    };
    
    let span = cloud_end - cloud_start;
    
    // 1. Calculate coverage profile
    let bin_size = (span / 100).max(1000); // ~100 bins or at least 1kb
    let coverage = calculate_coverage_profile(
        bam_path, chrom, cloud_start, cloud_end, bin_size, target_reads
    )?;
    
    // 2. Find coverage valleys
    let valleys = find_coverage_valleys(&coverage, 3, 0.5);
    println!("      Found {} coverage valleys", valleys.len());
    
    // 3. Collect transcript boundaries
    let boundaries = collect_transcript_boundaries(
        bam_path, chrom, cloud_start, cloud_end, target_reads
    )?;
    println!("      Found {} transcript boundaries", boundaries.len());
    
    // 4. Find boundary clusters
    let boundary_clusters = find_boundary_clusters(&boundaries, 5000, 20);
    println!("      Found {} boundary clusters", boundary_clusters.len());
    
    // 5. Combine all split points
    let mut all_splits: Vec<SplitPoint> = Vec::new();
    all_splits.extend(valleys);
    all_splits.extend(boundary_clusters);
    
    // Sort by position
    all_splits.sort_by_key(|s| s.position);
    
    // Merge nearby split points
    let mut merged_splits: Vec<SplitPoint> = Vec::new();
    for split in all_splits {
        if let Some(last) = merged_splits.last_mut() {
            if (split.position - last.position).abs() < 10000 {
                // Merge - take weighted average
                let total_conf = last.confidence + split.confidence;
                last.position = ((last.position as f64 * last.confidence 
                    + split.position as f64 * split.confidence) / total_conf) as i64;
                last.confidence = total_conf.min(1.0);
                last.evidence_type = format!("{}+{}", last.evidence_type, split.evidence_type);
            } else {
                merged_splits.push(split);
            }
        } else {
            merged_splits.push(split);
        }
    }
    
    println!("      Total unique split points: {}", merged_splits.len());
    
    // If we have expected gene count, filter to best N-1 split points
    if let Some(n) = expected_genes {
        if n > 1 && merged_splits.len() >= n - 1 {
            // Seed-aware selection: ensure split points are placed to isolate the seed
            if seed_in_cloud {
                // We need at least one split upstream and one downstream of seed
                let upstream_splits: Vec<_> = merged_splits.iter()
                    .filter(|s| seed_start.map(|ss| s.position < ss).unwrap_or(false))
                    .cloned()
                    .collect();
                let downstream_splits: Vec<_> = merged_splits.iter()
                    .filter(|s| seed_end.map(|se| s.position > se).unwrap_or(false))
                    .cloned()
                    .collect();
                
                println!("      Seed-aware selection: {} upstream, {} downstream candidates", 
                         upstream_splits.len(), downstream_splits.len());
                
                // Select best from each side
                let mut selected = Vec::new();
                if !upstream_splits.is_empty() {
                    let best = upstream_splits.iter()
                        .max_by(|a, b| a.confidence.partial_cmp(&b.confidence).unwrap())
                        .unwrap();
                    selected.push(best.clone());
                }
                if !downstream_splits.is_empty() {
                    let best = downstream_splits.iter()
                        .max_by(|a, b| a.confidence.partial_cmp(&b.confidence).unwrap())
                        .unwrap();
                    selected.push(best.clone());
                }
                
                // If we need more splits, add from either side by confidence
                let need_more = (n - 1).saturating_sub(selected.len());
                if need_more > 0 {
                    let mut remaining: Vec<_> = merged_splits.iter()
                        .filter(|s| !selected.iter().any(|sel| sel.position == s.position))
                        .cloned()
                        .collect();
                    remaining.sort_by(|a, b| b.confidence.partial_cmp(&a.confidence).unwrap());
                    selected.extend(remaining.into_iter().take(need_more));
                }
                
                merged_splits = selected;
                merged_splits.sort_by_key(|s| s.position);
                println!("      Selected {} seed-aware split points for {} genes", 
                         merged_splits.len(), n);
            } else {
                // No seed in cloud - use standard selection
                merged_splits.sort_by(|a, b| b.confidence.partial_cmp(&a.confidence).unwrap());
                merged_splits.truncate(n - 1);
                merged_splits.sort_by_key(|s| s.position);
                println!("      Selected top {} split points for {} genes", n - 1, n);
            }
        }
    }
    
    // Create sub-loci
    if merged_splits.is_empty() {
        // No splits - return whole cloud as one locus
        return Ok(vec![SubLocus {
            chrom: chrom.to_string(),
            start: cloud_start,
            end: cloud_end,
            split_evidence: vec!["no_clear_splits".to_string()],
        }]);
    }
    
    let mut subloci = Vec::new();
    let mut current_start = cloud_start;
    
    for (i, split) in merged_splits.iter().enumerate() {
        subloci.push(SubLocus {
            chrom: chrom.to_string(),
            start: current_start,
            end: split.position,
            split_evidence: vec![format!("split_{}: {}", i + 1, split.evidence_type)],
        });
        current_start = split.position;
    }
    
    // Add final locus
    subloci.push(SubLocus {
        chrom: chrom.to_string(),
        start: current_start,
        end: cloud_end,
        split_evidence: vec!["final_segment".to_string()],
    });
    
    println!("      Split into {} sub-loci:", subloci.len());
    for (i, sub) in subloci.iter().enumerate() {
        println!("        Locus {}: {}:{}-{} ({} bp)", 
                 i + 1, sub.chrom, sub.start, sub.end, sub.end - sub.start);
    }
    
    Ok(subloci)
}
