//! Splice-aware sub-cloud clustering
//!
//! Strategy:
//! 1. Find the initial broad read cloud (all paralogs)
//! 2. Extract splice junction patterns from each read
//! 3. Cluster reads by similar splice patterns (same gene family member)
//! 4. Create sub-clouds from splice-based clusters
//!
//! Key insight: Paralogs may share sequence similarity (multi-map) but have
//! different intron/exon structures due to different functional constraints.

use anyhow::Result;
use fxhash::FxHashSet as HashSet;
use noodles::bam;
use noodles::core::Position;
use noodles::sam::alignment::Record as SamRecord;
use std::collections::HashMap;

/// A splice junction: (start, end) positions
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct SpliceJunction {
    start: i64,
    end: i64,
}

/// Splice signature for a read - the set of all splice junctions
#[derive(Debug, Clone)]
struct SpliceSignature {
    read_name: String,
    junctions: Vec<SpliceJunction>,
    chrom: String,
    start: i64,
    end: i64,
}

impl SpliceSignature {
    /// Calculate similarity to another signature (0-1)
    fn similarity(&self, other: &SpliceSignature) -> f64 {
        if self.junctions.is_empty() || other.junctions.is_empty() {
            return 0.0;
        }
        
        // Count shared junctions (with some position tolerance)
        let tolerance = 10; // bp
        let mut shared = 0;
        
        for j1 in &self.junctions {
            for j2 in &other.junctions {
                if (j1.start - j2.start).abs() <= tolerance && 
                   (j1.end - j2.end).abs() <= tolerance {
                    shared += 1;
                    break;
                }
            }
        }
        
        // Jaccard-like similarity
        let union = self.junctions.len() + other.junctions.len() - shared;
        if union > 0 {
            shared as f64 / union as f64
        } else {
            0.0
        }
    }
    
    fn n_junctions(&self) -> usize {
        self.junctions.len()
    }
}

/// Extract splice junctions from a BAM record
fn extract_splice_junctions(record: &bam::Record, chrom: &str) -> Option<SpliceSignature> {
    if record.flags().is_unmapped() {
        return None;
    }
    
    let name = record.name()
        .and_then(|n| std::str::from_utf8(n).ok())
        .unwrap_or("unknown")
        .to_string();
    
    let start = record.alignment_start()
        .transpose()
        .ok()??
        .get() as i64;
    
    let end = record.alignment_end()
        .transpose()
        .ok()??
        .get() as i64;
    
    let cigar = record.cigar();
    let mut junctions = Vec::new();
    let mut ref_pos = start;
    
    for result in cigar.iter() {
        if let Ok(op) = result {
            use noodles::sam::alignment::record::cigar::op::Kind;
            match op.kind() {
                Kind::Skip => {
                    // N operation = intron
                    let intron_start = ref_pos;
                    let intron_end = ref_pos + op.len() as i64;
                    junctions.push(SpliceJunction {
                        start: intron_start,
                        end: intron_end,
                    });
                    ref_pos = intron_end;
                }
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch | 
                Kind::Deletion => {
                    ref_pos += op.len() as i64;
                }
                _ => {}
            }
        }
    }
    
    Some(SpliceSignature {
        read_name: name,
        junctions,
        chrom: chrom.to_string(),
        start,
        end,
    })
}

/// Collect splice signatures for all reads in a region
fn collect_splice_signatures(
    bam_path: &str,
    chrom: &str,
    region_start: i64,
    region_end: i64,
    target_reads: Option<&HashSet<String>>,
) -> Result<Vec<SpliceSignature>> {
    let mut reader = bam::io::indexed_reader::Builder::default()
        .build_from_path(bam_path)?;
    let header = reader.read_header()?;
    
    let start_pos = Position::try_from(region_start.max(1) as usize)?;
    let end_pos = Position::try_from(region_end.max(1) as usize)?;
    let region = noodles::core::Region::new(chrom, start_pos..=end_pos);
    
    let mut signatures = Vec::new();
    let mut seen_reads: HashSet<String> = HashSet::default();
    
    for result in reader.query(&header, &region)? {
        let record = result?;
        
        if record.flags().is_unmapped() {
            continue;
        }
        
        // Get read name
        let name = match record.name() {
            Some(n) => match std::str::from_utf8(n) {
                Ok(s) => s.to_string(),
                Err(_) => continue,
            },
            None => continue,
        };
        
        // Filter by target reads if specified
        if let Some(targets) = target_reads {
            if !targets.contains(&name) {
                continue;
            }
        }
        
        // Only process one alignment per read (prefer primary)
        if seen_reads.contains(&name) {
            continue;
        }
        seen_reads.insert(name.clone());
        
        if let Some(sig) = extract_splice_junctions(&record, chrom) {
            if !sig.junctions.is_empty() {
                signatures.push(sig);
            }
        }
    }
    
    Ok(signatures)
}

/// A sub-cloud based on splice pattern similarity
#[derive(Debug, Clone)]
pub struct SubCloud {
    pub chrom: String,
    pub start: i64,
    pub end: i64,
    pub read_names: HashSet<String>,
    pub n_spliced_reads: usize,
    pub representative_junctions: Vec<SpliceJunction>,
}

/// Cluster reads by splice pattern similarity
fn cluster_by_splice_patterns(
    signatures: &[SpliceSignature],
    min_similarity: f64,
    min_cluster_size: usize,
) -> Vec<Vec<&SpliceSignature>> {
    if signatures.is_empty() {
        return vec![];
    }
    
    let mut clusters: Vec<Vec<&SpliceSignature>> = Vec::new();
    
    for sig in signatures {
        let mut best_cluster: Option<usize> = None;
        let mut best_similarity = 0.0;
        
        // Find best matching cluster
        for (i, cluster) in clusters.iter().enumerate() {
            // Compare to cluster representative (first member)
            if let Some(rep) = cluster.first() {
                let sim = sig.similarity(rep);
                if sim >= min_similarity && sim > best_similarity {
                    best_similarity = sim;
                    best_cluster = Some(i);
                }
            }
        }
        
        if let Some(idx) = best_cluster {
            clusters[idx].push(sig);
        } else {
            clusters.push(vec![sig]);
        }
    }
    
    // Filter by minimum size
    clusters.into_iter()
        .filter(|c| c.len() >= min_cluster_size)
        .collect()
}

/// Create sub-clouds from splice-based clusters
fn create_sub_clouds(clusters: Vec<Vec<&SpliceSignature>>) -> Vec<SubCloud> {
    clusters.into_iter().map(|cluster| {
        let chrom = cluster[0].chrom.clone();
        let start = cluster.iter().map(|s| s.start).min().unwrap_or(0);
        let end = cluster.iter().map(|s| s.end).max().unwrap_or(0);
        
        let read_names: HashSet<String> = cluster
            .iter()
            .map(|s| s.read_name.clone())
            .collect();
        
        // Find representative junctions (present in >50% of reads)
        let mut junction_counts: HashMap<SpliceJunction, usize> = HashMap::default();
        for sig in &cluster {
            for junction in &sig.junctions {
                *junction_counts.entry(junction.clone()).or_default() += 1;
            }
        }
        
        let threshold = cluster.len() / 2;
        let representative: Vec<SpliceJunction> = junction_counts
            .into_iter()
            .filter(|(_, count)| *count >= threshold)
            .map(|(j, _)| j)
            .collect();
        
        SubCloud {
            chrom,
            start,
            end,
            n_spliced_reads: cluster.len(),
            read_names,
            representative_junctions: representative,
        }
    }).collect()
}

/// Main entry point: create splice-aware sub-clouds
pub fn create_splice_subclouds(
    bam_path: &str,
    chrom: &str,
    region_start: i64,
    region_end: i64,
    target_reads: Option<&HashSet<String>>,
    min_similarity: f64,
    min_cluster_size: usize,
) -> Result<Vec<SubCloud>> {
    println!("    Collecting splice signatures...");
    
    let signatures = collect_splice_signatures(
        bam_path, chrom, region_start, region_end, target_reads
    )?;
    
    println!("      Found {} reads with splice junctions", signatures.len());
    
    if signatures.is_empty() {
        return Ok(vec![]);
    }
    
    // Show distribution of junction counts
    let mut junction_dist: HashMap<usize, usize> = HashMap::default();
    for sig in &signatures {
        *junction_dist.entry(sig.n_junctions()).or_default() += 1;
    }
    
    let mut sorted_dist: Vec<_> = junction_dist.iter().collect();
    sorted_dist.sort_by_key(|(k, _)| *k);
    println!("      Junction count distribution:");
    for (n_junctions, count) in sorted_dist.iter().take(10) {
        println!("        {} junctions: {} reads", n_junctions, count);
    }
    
    println!("    Clustering by splice pattern similarity (min={})...", min_similarity);
    
    let clusters = cluster_by_splice_patterns(&signatures, min_similarity, min_cluster_size);
    
    println!("      Created {} splice-based clusters", clusters.len());
    
    let subclouds = create_sub_clouds(clusters);
    
    println!("      Sub-clouds:");
    for (i, cloud) in subclouds.iter().enumerate() {
        println!(
            "        {}: {}:{}-{} ({} reads, {} representative junctions)",
            i + 1,
            cloud.chrom,
            cloud.start,
            cloud.end,
            cloud.n_spliced_reads,
            cloud.representative_junctions.len()
        );
    }
    
    Ok(subclouds)
}

/// Merge overlapping sub-clouds and add unspliced reads
pub fn refine_sub_clouds(
    subclouds: Vec<SubCloud>,
    all_reads: &HashSet<String>,
    min_reads: usize,
) -> Vec<SubCloud> {
    // Sort by position
    let mut sorted = subclouds;
    sorted.sort_by_key(|c| c.start);
    
    // Merge overlapping
    let mut merged: Vec<SubCloud> = Vec::new();
    
    for cloud in sorted {
        if let Some(last) = merged.last_mut() {
            if cloud.start < last.end && cloud.chrom == last.chrom {
                // Overlap - merge
                last.end = last.end.max(cloud.end);
                last.read_names.extend(cloud.read_names);
                last.n_spliced_reads += cloud.n_spliced_reads;
            } else {
                merged.push(cloud);
            }
        } else {
            merged.push(cloud);
        }
    }
    
    // Filter by size
    merged.into_iter()
        .filter(|c| c.read_names.len() >= min_reads)
        .collect()
}
