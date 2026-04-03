//! Read cloud analysis: find dense clusters of overlapping multi-mapping reads
//! 
//! Core idea: Seed reads that are truly from a gene family will form dense "clouds"
//! at paralog locations, where multiple reads overlap each other.
//! Isolated mappings (single reads) are likely noise.

use anyhow::Result;
use bstr::BString;
use fxhash::FxHashSet as HashSet;
use noodles::bam;
use noodles::core::Position;
use noodles::sam::alignment::Record as SamRecord;
use std::collections::HashMap;

/// A read alignment with its full extent
#[derive(Debug, Clone)]
struct ReadAlignment {
    name: String,
    chrom: String,
    start: i64,
    end: i64,
    is_primary: bool,
    is_secondary: bool,
    is_supplementary: bool,
    strand: char,  // '+' or '-'
    soft_clip_5p: i64,  // Length of soft clip at 5' end
    soft_clip_3p: i64,  // Length of soft clip at 3' end
    num_passes: i32,    // PacBio CCS passes (from np tag)
    nm_tag: i32,        // Number of mismatches (NM tag) - higher = more diverged = more weight
}

/// A dense cluster of overlapping reads
#[derive(Debug, Clone)]
pub struct ReadCloud {
    pub chrom: String,
    pub start: i64,
    pub end: i64,
    pub read_names: HashSet<String>,  // FxHashSet for consistency
    pub n_primary: usize,
    pub n_secondary: usize,
    pub n_supplementary: usize,
    pub density_score: f64, // reads per kb
    pub strand: Option<char>, // '+' or '-' if strand-specific, None if mixed
}

impl ReadCloud {
    pub fn span(&self) -> i64 {
        self.end - self.start
    }
    
    pub fn n_unique_reads(&self) -> usize {
        self.read_names.len()
    }
}

/// Collect ALL alignments (primary + secondary + supplementary) for a set of reads
fn collect_all_alignments(
    bam_path: &str,
    target_reads: &HashSet<String>,
    chroms: &[String],
) -> Result<Vec<ReadAlignment>> {
    let mut all_alignments = Vec::new();
    
    for chrom in chroms {
        let mut reader = bam::io::indexed_reader::Builder::default()
            .build_from_path(bam_path)?;
        let header = reader.read_header()?;
        
        let chrom_name: BString = chrom.as_bytes().into();
        let ref_seqs = header.reference_sequences();
        let tid = ref_seqs.get_index_of(&chrom_name)
            .ok_or_else(|| anyhow::anyhow!("Chromosome {} not found", chrom))?;
        
        let ref_len = ref_seqs[tid].length().get();
        let start = Position::MIN;
        let end = Position::new(ref_len).unwrap_or(Position::MAX);
        let region = noodles::core::Region::new(chrom.clone(), start..=end);
        
        for result in reader.query(&header, &region)? {
            let record = result?;
            
            if record.flags().is_unmapped() {
                continue;
            }
            
            let name = match record.name() {
                Some(n) => match std::str::from_utf8(n) {
                    Ok(s) => s.to_string(),
                    Err(_) => continue,
                },
                None => continue,
            };
            
            if !target_reads.contains(&name) {
                continue;
            }
            
            let start = record.alignment_start()
                .transpose()?
                .map(|p| p.get() as i64)
                .unwrap_or(0);
            
            let end = record.alignment_end()
                .transpose()?
                .map(|p| p.get() as i64)
                .unwrap_or(start + 1);
            
            let flags = record.flags();
            
            // Determine strand
            let strand = if flags.is_reverse_complemented() { '-' } else { '+' };
            
            // Extract soft clip lengths from CIGAR
            let (soft_clip_5p, soft_clip_3p) = extract_soft_clips(&record);
            
            // Extract number of passes (PacBio np tag)
            let num_passes = extract_num_passes(&record);
            
            // Extract NM tag (number of mismatches)
            let nm_tag = extract_nm_tag(&record);
            
            all_alignments.push(ReadAlignment {
                name,
                chrom: chrom.clone(),
                start,
                end,
                is_primary: !flags.is_secondary() && !flags.is_supplementary(),
                is_secondary: flags.is_secondary(),
                is_supplementary: flags.is_supplementary(),
                strand,
                soft_clip_5p,
                soft_clip_3p,
                num_passes,
                nm_tag,
            });
        }
    }
    
    Ok(all_alignments)
}

/// Extract soft clip lengths from CIGAR string
/// Returns (5' soft clip length, 3' soft clip length)
fn extract_soft_clips(record: &bam::Record) -> (i64, i64) {
    crate::soft_clip::extract_soft_clip_5p_3p_lens(record)
}

/// Extract number of passes from PacBio np tag
fn extract_num_passes(record: &bam::Record) -> i32 {
    for result in record.data().iter() {
        if let Ok((tag, value)) = result {
            if tag.as_ref() == b"np" {
                match value {
                    noodles::sam::alignment::record::data::field::Value::Int32(v) => return v,
                    noodles::sam::alignment::record::data::field::Value::Int16(v) => return v as i32,
                    noodles::sam::alignment::record::data::field::Value::Int8(v) => return v as i32,
                    _ => {}
                }
            }
        }
    }
    0 // Default: unknown
}

/// Extract NM tag (number of mismatches) from BAM record
fn extract_nm_tag(record: &bam::Record) -> i32 {
    for result in record.data().iter() {
        if let Ok((tag, value)) = result {
            if tag.as_ref() == b"NM" {
                match value {
                    noodles::sam::alignment::record::data::field::Value::Int32(v) => return v,
                    noodles::sam::alignment::record::data::field::Value::Int16(v) => return v as i32,
                    noodles::sam::alignment::record::data::field::Value::Int8(v) => return v as i32,
                    noodles::sam::alignment::record::data::field::Value::UInt8(v) => return v as i32,
                    _ => {}
                }
            }
        }
    }
    0 // Default: no mismatches
}

/// Build a coverage profile and find dense peaks
/// With strand separation and full-length weighting
fn find_dense_clouds(
    alignments: &[ReadAlignment],
    min_overlap_reads: usize,
    min_cloud_span: i64,
    max_cloud_span: i64,
    strand_aware: bool,
) -> Vec<ReadCloud> {
    if alignments.is_empty() {
        return vec![];
    }
    
    // Group by chromosome
    let mut by_chrom: HashMap<String, Vec<&ReadAlignment>> = HashMap::new();
    for aln in alignments {
        by_chrom.entry(aln.chrom.clone()).or_default().push(aln);
    }
    
    println!("    Processing {} chromosomes...", by_chrom.len());
    for (chrom, aligns) in &by_chrom {
        println!("      {}: {} alignments", chrom, aligns.len());
    }
    
    let mut all_clouds = Vec::new();
    
    for (chrom, chrom_aligns) in by_chrom {
        if strand_aware {
            // Separate by strand and find clouds independently
            let plus_aligns: Vec<_> = chrom_aligns.iter()
                .filter(|a| a.strand == '+')
                .copied()
                .collect();
            let minus_aligns: Vec<_> = chrom_aligns.iter()
                .filter(|a| a.strand == '-')
                .copied()
                .collect();
            
            println!("      {}: {} (+) / {} (-) alignments by strand", 
                     chrom, plus_aligns.len(), minus_aligns.len());
            
            let plus_clouds = find_clouds_on_strand(&plus_aligns, min_overlap_reads, 
                                                     min_cloud_span, max_cloud_span, '+');
            let minus_clouds = find_clouds_on_strand(&minus_aligns, min_overlap_reads,
                                                     min_cloud_span, max_cloud_span, '-');
            
            all_clouds.extend(plus_clouds);
            all_clouds.extend(minus_clouds);
        } else {
            // Original behavior - combine strands
            let clouds = find_clouds_on_strand(&chrom_aligns, min_overlap_reads,
                                               min_cloud_span, max_cloud_span, '+');
            all_clouds.extend(clouds);
        }
    }
    
    println!("    DEBUG: Total clouds created before merge: {}", all_clouds.len());
    
    // Merge overlapping clouds
    let merged = merge_overlapping_clouds(all_clouds);
    
    // Extend clouds by soft-clipped bases
    let extended = extend_clouds_by_softclips(merged, alignments);
    
    extended
}

/// Find clouds on a single strand (or combined)
fn find_clouds_on_strand(
    alignments: &[&ReadAlignment],
    min_overlap_reads: usize,
    min_cloud_span: i64,
    max_cloud_span: i64,
    strand: char,
) -> Vec<ReadCloud> {
    if alignments.is_empty() {
        return vec![];
    }
    
    let mut sorted: Vec<_> = alignments.iter().copied().collect();
    sorted.sort_by_key(|a| a.start);
    
    let mut clouds = Vec::new();
    let mut active_reads: Vec<&&ReadAlignment> = Vec::new();
    let mut max_active = 0;
    
    for aln in &sorted {
        // Remove reads that no longer overlap
        active_reads.retain(|a| a.end > aln.start);
        active_reads.push(aln);
        
        // Calculate weighted unique reads
        // Weight = full-length bonus × divergence bonus
        // - Full-length reads (≥10 passes): 2x weight
        // - Diverged reads (NM ≥ 5): up to 2x additional weight (reads mapping to diverged paralogs)
        let mut weighted_count = 0.0;
        let mut unique_names: HashSet<String> = HashSet::default();
        
        for a in &active_reads {
            if unique_names.insert(a.name.clone()) {
                // Base weight from number of passes
                let full_length_weight = if a.num_passes >= 10 { 2.0 } else { 1.0 };
                
                // Divergence weight: reads with more mismatches are more likely to be
                // from diverged paralogs (helps distinguish highly similar genes)
                // NM 0-2: 1.0x, NM 3-5: 1.3x, NM 6-10: 1.6x, NM >10: 2.0x
                let divergence_weight = match a.nm_tag {
                    0..=2 => 1.0,
                    3..=5 => 1.3,
                    6..=10 => 1.6,
                    _ => 2.0,
                };
                
                weighted_count += full_length_weight * divergence_weight;
            }
        }
        
        max_active = max_active.max(unique_names.len());
        
        // Use weighted count for threshold
        if weighted_count >= min_overlap_reads as f64 {
            let start = active_reads.iter().map(|a| a.start).min().unwrap_or(0);
            let end = active_reads.iter().map(|a| a.end).max().unwrap_or(0);
            let span = end - start;
            
            if span >= min_cloud_span && span <= max_cloud_span {
                let n_primary = active_reads.iter().filter(|a| a.is_primary).count();
                let n_secondary = active_reads.iter().filter(|a| a.is_secondary).count();
                let n_supp = active_reads.iter().filter(|a| a.is_supplementary).count();
                
                let density = if span > 0 {
                    unique_names.len() as f64 / (span as f64 / 1000.0)
                } else {
                    0.0
                };
                
                clouds.push(ReadCloud {
                    chrom: active_reads[0].chrom.clone(),
                    start,
                    end,
                    read_names: unique_names,
                    n_primary,
                    n_secondary,
                    n_supplementary: n_supp,
                    density_score: density,
                    strand: Some(strand),
                });
            }
        }
    }
    
    println!("      DEBUG: Max overlapping reads at any position: {}", max_active);
    
    clouds
}

/// Extend clouds by soft-clipped bases at edges
fn extend_clouds_by_softclips(
    mut clouds: Vec<ReadCloud>,
    all_alignments: &[ReadAlignment],
) -> Vec<ReadCloud> {
    for cloud in &mut clouds {
        let mut max_extension_5p = 0i64;
        let mut max_extension_3p = 0i64;
        
        // Find alignments at cloud edges
        for aln in all_alignments {
            if aln.chrom != cloud.chrom {
                continue;
            }
            
            // Check if alignment is at 5' edge of cloud (within 1kb)
            let at_5p_edge = (aln.start - cloud.start).abs() < 1000;
            let at_3p_edge = (aln.end - cloud.end).abs() < 1000;
            
            if at_5p_edge && aln.soft_clip_5p > max_extension_5p {
                max_extension_5p = aln.soft_clip_5p;
            }
            
            if at_3p_edge && aln.soft_clip_3p > max_extension_3p {
                max_extension_3p = aln.soft_clip_3p;
            }
        }
        
        // Apply extensions (cap at reasonable limits)
        let max_extension = 5000i64; // Don't extend more than 5kb
        cloud.start -= max_extension_5p.min(max_extension);
        cloud.end += max_extension_3p.min(max_extension);
        
        if max_extension_5p > 0 || max_extension_3p > 0 {
            println!("      Extended cloud by {}bp 5' and {}bp 3' (soft clips)",
                     max_extension_5p, max_extension_3p);
        }
    }
    
    clouds
}

/// Merge clouds that significantly overlap
fn merge_overlapping_clouds(clouds: Vec<ReadCloud>) -> Vec<ReadCloud> {
    if clouds.is_empty() {
        return vec![];
    }
    
    let mut sorted = clouds;
    sorted.sort_by(|a, b| a.chrom.cmp(&b.chrom).then(a.start.cmp(&b.start)));
    
    let mut merged: Vec<ReadCloud> = Vec::new();
    let mut current = sorted[0].clone();
    
    for cloud in sorted.into_iter().skip(1) {
        if cloud.chrom == current.chrom && cloud.start < current.end {
            // Overlapping - merge
            current.end = current.end.max(cloud.end);
            current.read_names.extend(cloud.read_names);
            current.n_primary += cloud.n_primary;
            current.n_secondary += cloud.n_secondary;
            current.n_supplementary += cloud.n_supplementary;
            // Recalculate density
            let span = current.span();
            if span > 0 {
                current.density_score = current.read_names.len() as f64 / (span as f64 / 1000.0);
            }
        } else {
            // No overlap - save current and start new
            merged.push(current);
            current = cloud;
        }
    }
    merged.push(current);
    
    merged
}

/// Main entry point: find read clouds for a set of seed reads
pub fn find_read_clouds(
    bam_path: &str,
    seed_reads: &HashSet<String>,
    chroms: &[String],
    min_overlap_reads: usize,
    min_cloud_span: i64,
    max_cloud_span: i64,
    strand_aware: bool,
) -> Result<Vec<ReadCloud>> {
    println!("  Collecting all alignments for {} seed reads...", seed_reads.len());
    if strand_aware {
        println!("  Using strand-aware cloud detection");
    }
    
    let all_alignments = collect_all_alignments(bam_path, seed_reads, chroms)?;
    
    let primary_count = all_alignments.iter().filter(|a| a.is_primary).count();
    let secondary_count = all_alignments.iter().filter(|a| a.is_secondary).count();
    let supp_count = all_alignments.iter().filter(|a| a.is_supplementary).count();
    
    println!("  Found {} total alignments:", all_alignments.len());
    println!("    Primary: {}", primary_count);
    println!("    Secondary: {}", secondary_count);
    println!("    Supplementary: {}", supp_count);
    
    println!("  Finding dense read clouds (min {} overlapping reads)...", min_overlap_reads);
    
    let clouds = find_dense_clouds(
        &all_alignments,
        min_overlap_reads,
        min_cloud_span,
        max_cloud_span,
        strand_aware,
    );
    
    println!("  Found {} read clouds", clouds.len());
    
    for (i, cloud) in clouds.iter().enumerate() {
        println!(
            "    Cloud {}: {}:{}-{} (span: {} bp, {} unique reads, density: {:.2} reads/kb)",
            i + 1,
            cloud.chrom,
            cloud.start,
            cloud.end,
            cloud.span(),
            cloud.n_unique_reads(),
            cloud.density_score
        );
    }
    
    Ok(clouds)
}

/// Filter clouds by quality metrics
pub fn filter_clouds(
    clouds: Vec<ReadCloud>,
    min_unique_reads: usize,
    min_density: f64,
    require_secondary: bool,
) -> Vec<ReadCloud> {
    clouds.into_iter()
        .filter(|c| {
            if c.n_unique_reads() < min_unique_reads {
                return false;
            }
            if c.density_score < min_density {
                return false;
            }
            if require_secondary && c.n_secondary == 0 {
                return false;
            }
            true
        })
        .collect()
}
