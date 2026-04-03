//! Transcript boundary-based splitting for tandem arrays

use anyhow::Result;
use noodles::bam;
use noodles::core::Position;
use noodles::sam::alignment::Record;
use std::collections::HashMap;

/// A transcript with boundary information
#[derive(Debug, Clone)]
pub struct TranscriptBoundary {
    pub name: String,
    pub chrom: String,
    pub start: i64,
    pub end: i64,
    pub strand: char,
    pub length: i64,
}

/// Gene locus derived from transcript boundary clustering
#[derive(Debug, Clone)]
pub struct GeneLocus {
    pub chrom: String,
    pub start: i64,
    pub end: i64,
    pub strand: char,
    pub n_transcripts: usize,
    pub transcript_names: Vec<String>,
}

/// Extract transcript boundaries from a BAM record
fn extract_transcript_boundary(record: &bam::Record, chrom: &str) -> Option<TranscriptBoundary> {
    if record.flags().is_unmapped() {
        return None;
    }
    
    let name = record.name()
        .and_then(|n| std::str::from_utf8(n).ok())
        .unwrap_or("unknown")
        .to_string();
    
    let start = match record.alignment_start() {
        Some(Ok(p)) => p.get() as i64,
        _ => return None,
    };
    
    // Get end from alignment_end if available, otherwise parse CIGAR
    let end = match record.alignment_end() {
        Some(Ok(p)) => p.get() as i64,
        _ => {
            // Calculate end from CIGAR
            let mut ref_pos = start;
            let cigar = record.cigar();
            
            for result in cigar.iter() {
                let op = match result {
                    Ok(op) => op,
                    Err(_) => continue,
                };
                
                let (op_kind, len) = (op.kind(), op.len());
                
                match op_kind {
                    noodles::sam::alignment::record::cigar::op::Kind::Match |
                    noodles::sam::alignment::record::cigar::op::Kind::SequenceMismatch |
                    noodles::sam::alignment::record::cigar::op::Kind::SequenceMatch |
                    noodles::sam::alignment::record::cigar::op::Kind::Deletion |
                    noodles::sam::alignment::record::cigar::op::Kind::Skip => {
                        ref_pos += len as i64;
                    }
                    _ => {} // Others don't consume reference
                }
            }
            ref_pos
        }
    };
    
    // Determine strand from flags
    let flags_val = record.flags().bits();
    let strand = if flags_val & 0x10 != 0 { '-' } else { '+' };
    
    Some(TranscriptBoundary {
        name,
        chrom: chrom.to_string(),
        start,
        end,
        strand,
        length: end - start,
    })
}

/// Collect all transcript boundaries from a region
pub fn collect_transcript_boundaries(
    bam_path: &str,
    chrom: &str,
    start: i64,
    end: i64,
) -> Result<Vec<TranscriptBoundary>> {
    let mut reader = bam::io::indexed_reader::Builder::default()
        .build_from_path(bam_path)?;

    let header = reader.read_header()?;

    let start_pos = Position::try_from(start.max(1) as usize)?;
    let end_pos = Position::try_from(end.max(1) as usize)?;
    let region = noodles::core::Region::new(chrom, start_pos..=end_pos);

    let mut transcripts = Vec::new();

    for result in reader.query(&header, &region)? {
        let record = result?;

        if record.flags().is_unmapped() {
            continue;
        }
        
        // Skip secondary and supplementary alignments
        if record.flags().is_secondary() || record.flags().is_supplementary() {
            continue;
        }

        if let Some(tx) = extract_transcript_boundary(&record, chrom) {
            // Filter out transcripts with suspiciously large spans (>200kb)
            // These are likely pre-mRNA or chimeric reads
            if tx.length <= 200000 {
                transcripts.push(tx);
            }
        }
    }

    Ok(transcripts)
}

/// Cluster transcripts by 5' end position (tss)
fn cluster_by_tss(
    transcripts: Vec<TranscriptBoundary>,
    max_tss_distance: i64,
) -> Vec<Vec<TranscriptBoundary>> {
    if transcripts.is_empty() {
        return vec![];
    }
    
    let mut strand_groups: HashMap<char, Vec<TranscriptBoundary>> = HashMap::new();
    for tx in transcripts {
        strand_groups.entry(tx.strand).or_default().push(tx);
    }
    
    let mut all_clusters = Vec::new();
    
    for (strand, mut txs) in strand_groups {
        // Sort by 5' end position
        txs.sort_by_key(|t| if strand == '+' { t.start } else { -t.end });
        
        let mut clusters: Vec<Vec<TranscriptBoundary>> = Vec::new();
        
        for tx in txs {
            let tss = if strand == '+' { tx.start } else { tx.end };
            
            let mut found_cluster = false;
            for cluster in &mut clusters {
                let rep_tss = if strand == '+' { 
                    cluster[0].start 
                } else { 
                    cluster[0].end 
                };
                
                if (tss - rep_tss).abs() < max_tss_distance {
                    cluster.push(tx.clone());
                    found_cluster = true;
                    break;
                }
            }
            
            if !found_cluster {
                clusters.push(vec![tx]);
            }
        }
        
        all_clusters.extend(clusters);
    }
    
    all_clusters
}

/// Cluster transcripts by 3' end position (tes)
fn cluster_by_tes(
    transcripts: Vec<TranscriptBoundary>,
    max_tes_distance: i64,
) -> Vec<Vec<TranscriptBoundary>> {
    if transcripts.is_empty() {
        return vec![];
    }
    
    let mut strand_groups: HashMap<char, Vec<TranscriptBoundary>> = HashMap::new();
    for tx in transcripts {
        strand_groups.entry(tx.strand).or_default().push(tx);
    }
    
    let mut all_clusters = Vec::new();
    
    for (strand, mut txs) in strand_groups {
        // Sort by 3' end position
        txs.sort_by_key(|t| if strand == '+' { t.end } else { -t.start });
        
        let mut clusters: Vec<Vec<TranscriptBoundary>> = Vec::new();
        
        for tx in txs {
            let tes = if strand == '+' { tx.end } else { tx.start };
            
            let mut found_cluster = false;
            for cluster in &mut clusters {
                let rep_tes = if strand == '+' { 
                    cluster[0].end 
                } else { 
                    cluster[0].start 
                };
                
                if (tes - rep_tes).abs() < max_tes_distance {
                    cluster.push(tx.clone());
                    found_cluster = true;
                    break;
                }
            }
            
            if !found_cluster {
                clusters.push(vec![tx]);
            }
        }
        
        all_clusters.extend(clusters);
    }
    
    all_clusters
}

/// Split a merged region into individual genes using transcript boundaries
pub fn split_region_by_transcript_boundaries(
    bam_path: &str,
    chrom: &str,
    region_start: i64,
    region_end: i64,
    max_boundary_distance: i64,
    min_transcripts: usize,
    expected_genes: Option<usize>,
) -> Result<Vec<GeneLocus>> {
    println!("  Extracting transcript boundaries from region...");
    let transcripts = collect_transcript_boundaries(bam_path, chrom, region_start, region_end)?;
    
    println!("  Found {} primary alignments", transcripts.len());
    
    if transcripts.len() < min_transcripts {
        // Return the whole region as one locus
        return Ok(vec![GeneLocus {
            chrom: chrom.to_string(),
            start: region_start,
            end: region_end,
            strand: '+',
            n_transcripts: transcripts.len(),
            transcript_names: transcripts.iter().map(|t| t.name.clone()).collect(),
        }]);
    }
    
    // Try clustering by TSS first
    println!("  Clustering by transcription start sites (TSS)...");
    let tss_clusters = cluster_by_tss(transcripts.clone(), max_boundary_distance);
    println!("  Found {} TSS clusters", tss_clusters.len());
    
    // Filter clusters by min transcripts
    let mut valid_clusters: Vec<Vec<TranscriptBoundary>> = tss_clusters.into_iter()
        .filter(|c| c.len() >= min_transcripts)
        .collect();
    
    // If we don't have expected number of clusters, try TES clustering
    if let Some(n) = expected_genes {
        if valid_clusters.len() != n {
            println!("  Trying 3' end clustering...");
            let tes_clusters = cluster_by_tes(transcripts, max_boundary_distance);
            println!("  Found {} TES clusters", tes_clusters.len());
            
            let valid_tes: Vec<_> = tes_clusters.into_iter()
                .filter(|c| c.len() >= min_transcripts)
                .collect();
            
            if valid_tes.len().abs_diff(n) < valid_clusters.len().abs_diff(n) {
                valid_clusters = valid_tes;
            }
        }
    }
    
    // Convert clusters to gene loci
    let mut loci = Vec::new();
    
    for cluster in valid_clusters {
        let chrom = cluster[0].chrom.clone();
        let start = cluster.iter().map(|t| t.start).min().unwrap_or(region_start);
        let end = cluster.iter().map(|t| t.end).max().unwrap_or(region_end);
        let strand = cluster[0].strand;
        
        loci.push(GeneLocus {
            chrom,
            start,
            end,
            strand,
            n_transcripts: cluster.len(),
            transcript_names: cluster.iter().map(|t| t.name.clone()).collect(),
        });
    }
    
    // Sort by position
    loci.sort_by_key(|l| l.start);
    
    println!("  Split into {} gene loci", loci.len());
    for (i, locus) in loci.iter().enumerate() {
        println!("    Gene {}: {}:{}-{} ({} transcripts)",
                 i + 1, locus.chrom, locus.start, locus.end, locus.n_transcripts);
    }
    
    Ok(loci)
}

/// Analyze transcript boundary distribution for a locus
pub fn analyze_boundary_distribution(
    bam_path: &str,
    chrom: &str,
    start: i64,
    end: i64,
) -> Result<BoundaryAnalysis> {
    let transcripts = collect_transcript_boundaries(bam_path, chrom, start, end)?;
    
    let total = transcripts.len();
    
    // Group by strand
    let plus_strand = transcripts.iter().filter(|t| t.strand == '+').count();
    let minus_strand = total - plus_strand;
    
    // Get 5' end distribution (TSS)
    let mut tss_positions: Vec<i64> = transcripts.iter()
        .map(|t| if t.strand == '+' { t.start } else { t.end })
        .collect();
    tss_positions.sort();
    
    // Get 3' end distribution (TES)
    let mut tes_positions: Vec<i64> = transcripts.iter()
        .map(|t| if t.strand == '+' { t.end } else { t.start })
        .collect();
    tes_positions.sort();
    
    // Calculate TSS density (peaks per 10kb)
    let span = (end - start).max(1) as f64;
    let tss_density = tss_positions.len() as f64 / span * 10000.0;
    
    Ok(BoundaryAnalysis {
        total_transcripts: total,
        plus_strand,
        minus_strand,
        tss_positions,
        tes_positions,
        tss_density,
    })
}

/// Boundary analysis results
#[derive(Debug, Clone)]
pub struct BoundaryAnalysis {
    pub total_transcripts: usize,
    pub plus_strand: usize,
    pub minus_strand: usize,
    pub tss_positions: Vec<i64>,
    pub tes_positions: Vec<i64>,
    pub tss_density: f64,
}
