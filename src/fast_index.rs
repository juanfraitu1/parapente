//! Fast in-memory BAM index for gene family detection
//! 
//! This module provides a pre-built index for fast multi-mapping read queries.
//! It loads all alignments into memory indexed by read name for quick access.

use anyhow::Result;
use bstr::BString;
use fxhash::FxHashSet as HashSet;
use noodles::bam;
use noodles::core::Position;
use noodles::sam::alignment::Record;
use std::collections::HashMap;

/// An alignment record summary for fast access
#[derive(Debug, Clone)]
pub struct AlignmentSummary {
    pub chrom: String,
    pub start: i64,
    pub end: i64,
    pub is_primary: bool,
    pub is_secondary: bool,
    pub is_supplementary: bool,
}

/// In-memory BAM index for fast queries
#[derive(Debug)]
pub struct BamIndex {
    /// Map from read name to all its alignments
    pub index: HashMap<String, Vec<AlignmentSummary>>,
    /// Total number of alignments indexed
    pub total_alignments: usize,
}

impl BamIndex {
    /// Create a new empty index
    pub fn new() -> Self {
        Self {
            index: HashMap::default(),
            total_alignments: 0,
        }
    }

    /// Get all alignments for a set of read names, grouped by chromosome
    pub fn get_alignments_by_chrom(
        &self,
        read_names: &HashSet<String>,
    ) -> HashMap<String, Vec<(i64, i64, String)>> {
        let mut by_chrom: HashMap<String, Vec<(i64, i64, String)>> = HashMap::new();

        for name in read_names {
            if let Some(alignments) = self.index.get(name) {
                for aln in alignments {
                    by_chrom
                        .entry(aln.chrom.clone())
                        .or_default()
                        .push((aln.start, aln.end, name.clone()));
                }
            }
        }

        by_chrom
    }

    /// Get all alignments for a specific read
    pub fn get_alignments(&self, read_name: &str) -> Option<&Vec<AlignmentSummary>> {
        self.index.get(read_name)
    }
}

impl Default for BamIndex {
    fn default() -> Self {
        Self::new()
    }
}

/// Build a BAM index from a BAM file
pub fn build_bam_index(bam_path: &str, _include_supplementary: bool) -> Result<BamIndex> {
    let mut reader = bam::io::indexed_reader::Builder::default()
        .build_from_path(bam_path)?;
    let header = reader.read_header()?;

    let mut index: HashMap<String, Vec<AlignmentSummary>> = HashMap::default();
    let mut total = 0usize;

    // Query each chromosome
    for (name, _seq) in header.reference_sequences().iter() {
        let chrom = String::from_utf8_lossy(name).to_string();
        let ref_seqs = header.reference_sequences();
        let chrom_name: BString = chrom.as_bytes().into();
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

            let name_bytes = match record.name() {
                Some(n) => n,
                None => continue,
            };

            let read_name = match std::str::from_utf8(name_bytes) {
                Ok(s) => s.to_string(),
                Err(_) => continue,
            };

            let aln_start = record.alignment_start()
                .transpose()?
                .map(|p: Position| p.get() as i64)
                .unwrap_or(0);
            let aln_end = record.alignment_end()
                .transpose()?
                .map(|p: Position| p.get() as i64)
                .unwrap_or(aln_start + 1);

            let flags = record.flags();

            let summary = AlignmentSummary {
                chrom: chrom.clone(),
                start: aln_start,
                end: aln_end,
                is_primary: !flags.is_secondary() && !flags.is_supplementary(),
                is_secondary: flags.is_secondary(),
                is_supplementary: flags.is_supplementary(),
            };

            index.entry(read_name).or_default().push(summary);
            total += 1;
        }
    }

    println!("Built BAM index: {} unique reads, {} total alignments", 
             index.len(), total);

    Ok(BamIndex {
        index,
        total_alignments: total,
    })
}

/// Cluster alignments into loci
pub fn cluster_loci(
    by_chrom: HashMap<String, Vec<(i64, i64, String)>>,
    cluster_distance: i64,
    min_reads: usize,
) -> Vec<(String, i64, i64, HashSet<String>)> {
    let mut loci: Vec<(String, i64, i64, HashSet<String>)> = vec![];

    for (chrom, mut aligns) in by_chrom {
        if aligns.len() < min_reads {
            continue;
        }

        aligns.sort_by_key(|a| a.0);

        let mut current_cluster = vec![aligns[0].clone()];
        let mut current_end = aligns[0].1;

        for (start, end, name) in aligns.into_iter().skip(1) {
            if start - current_end <= cluster_distance {
                current_cluster.push((start, end, name));
                if end > current_end {
                    current_end = end;
                }
            } else {
                if current_cluster.len() >= min_reads {
                    let reads: HashSet<String> = current_cluster
                        .iter()
                        .map(|(_, _, name)| name.clone())
                        .collect();
                    loci.push((chrom.clone(), current_cluster[0].0, current_end, reads));
                }
                current_cluster = vec![(start, end, name)];
                current_end = end;
            }
        }

        if current_cluster.len() >= min_reads {
            let reads: HashSet<String> = current_cluster
                .iter()
                .map(|(_, _, name)| name.clone())
                .collect();
            loci.push((chrom, current_cluster[0].0, current_end, reads));
        }
    }

    loci
}
