//! Transitive gene family detector for highly diverged paralogs
//!
//! For gene families like ID_14 (LRRC37 family), paralogs are so divergent
//! that different CCS reads align to each locus. The standard multi-mapping
//! approach fails because there are no shared read names.
//!
//! This module implements a transitive detection strategy based on READ CLOUDS:
//!
//! 1. **Collect seed reads**: Find ALL reads that overlap the seed region
//! 2. **Find all mappings**: For each seed read, find ALL its alignments genome-wide
//! 3. **Build read clouds**: At each candidate location, collect ALL reads
//!    that overlap there (expanding beyond just the seed reads)
//! 4. **Transitive closure**: Starting from seed, iteratively add loci that
//!    share reads with already-discovered loci
//! 5. **Graph traversal**: Build a graph where nodes are loci and edges represent
//!    shared reads; find all loci reachable from seed

use anyhow::Result;
use fxhash::FxHashSet as HashSet;
use noodles::bam;
use noodles::bgzf;
use noodles::core::Position;
use noodles::sam;
use noodles::sam::alignment::Record as SamRecord;
use rayon::prelude::*;
use std::collections::{HashMap, VecDeque};
use std::sync::Mutex;

use crate::coverage_valley::ValleyDetectionParams;
use crate::fast_index::{AlignmentSummary, BamIndex};
use crate::read_intern::ReadInterner;
use hashbrown::HashMap as HbHashMap;
use roaring::RoaringBitmap;

/// Tunable parameters for transitive-mode coverage valleys, segment merging, and junction splitting.
#[derive(Debug, Clone)]
pub struct TransitiveCoverageSplitConfig {
    /// When true, derive valley thresholds per locus from coverage (previous behavior).
    /// When false, use `fixed_valley` for every locus.
    pub use_dynamic_valleys: bool,
    /// Used when `use_dynamic_valleys` is false.
    pub fixed_valley: ValleyDetectionParams,
    /// Bin size for the coverage profile (must match `fixed_valley.bin_size` in fixed mode).
    pub profile_bin_size: i64,
    pub valley_merge_min_segment_bp: i64,
    pub valley_merge_max_gap_bp: i64,
    /// Loci strictly larger than this enter junction-based splitting.
    pub junction_min_locus_bp: i64,
    /// Minimum segment span for `split_by_junctions` to attempt a split inside a segment.
    pub junction_min_split_span_bp: i64,
    pub junction_merge_min_segment_bp: i64,
    pub junction_merge_max_gap_bp: i64,
}

/// A genomic locus discovered through transitive relationships
#[derive(Debug, Clone)]
pub struct TransitiveLocus {
    pub chrom: String,
    pub start: i64,
    pub end: i64,
    /// Read names that support this locus
    pub reads: HashSet<String>,
    /// Which iteration/component this locus was discovered in
    pub component: usize,
    /// Direct Jaccard with seed (may be 0 for transitively discovered)
    pub jaccard_with_seed: f64,
    /// Path distance from seed (0 = seed itself, 1 = direct neighbor, 2+ = transitive)
    pub distance_from_seed: usize,
    /// Number of coverage peaks in this locus
    pub peak_count: usize,
    /// Confidence score (0.0-1.0) based on coverage profile
    pub confidence: f64,
}

impl TransitiveLocus {
    /// Get the number of reads efficiently
    #[inline]
    pub fn read_count(&self) -> usize {
        self.reads.len()
    }
}

/// Representative seed candidate for a connected component that is not connected to current seed
#[derive(Debug, Clone)]
pub struct ComponentSeedCandidate {
    pub chrom: String,
    pub start: i64,
    pub end: i64,
    pub reads: HashSet<String>,
    pub component_size: usize,
    pub total_reads_in_component: usize,
}

impl TransitiveLocus {
    pub fn span(&self) -> i64 {
        self.end - self.start
    }
}

/// A read alignment with full details
#[derive(Debug, Clone)]
struct ReadMapping {
    read_name: String,
    chrom: String,
    start: i64,
    end: i64,
    is_primary: bool,
}

/// Read alignment with interned read id (multi-seed cached hot path).
#[derive(Debug, Clone)]
struct InternedReadMapping {
    read_id: u32,
    chrom: String,
    start: i64,
    end: i64,
    is_primary: bool,
}

/// In-memory index keyed by interned read id (avoids string hashing on hot paths).
/// `by_chrom` lists (read_id, start, end) for every alignment on that chromosome so
/// `collect_reads_in_region_interned` costs O(alignments on chrom) instead of O(all reads).
#[derive(Debug)]
struct InternedBamIndex {
    index: HbHashMap<u32, Vec<AlignmentSummary>>,
    by_chrom: HbHashMap<String, Vec<(u32, i64, i64)>>,
}

/// Parameters for transitive detection
#[derive(Debug, Clone)]
pub struct TransitiveParams {
    /// Maximum distance for clustering loci (bp)
    pub cluster_distance: i64,
    /// When true, each chromosome uses a merge gap derived from same-read alignment gaps
    /// (90th percentile, capped by cluster_distance) instead of a fixed cluster_distance only.
    pub adaptive_cluster_distance: bool,
    /// Minimum reads at a locus to report it
    pub min_reads_per_locus: usize,
    /// Maximum iterations for transitive expansion
    pub max_iterations: usize,
    /// Minimum Jaccard to consider two loci related (for graph edges)
    pub min_jaccard_for_edge: f64,
}

impl Default for TransitiveParams {
    fn default() -> Self {
        Self {
            cluster_distance: 10000,
            adaptive_cluster_distance: false,
            min_reads_per_locus: 10,
            max_iterations: 20,  // Convergence safety limit - runs until no new loci found
            min_jaccard_for_edge: 0.001,  // Very low threshold
        }
    }
}

/// Collect ALL reads that overlap a region
/// Uses a simple per-thread reader cache to avoid repeated BAM reader creation
fn collect_reads_in_region(
    bam_path: &str,
    chrom: &str,
    start: i64,
    end: i64,
) -> Result<HashSet<String>> {
    let mut reader = bam::io::indexed_reader::Builder::default()
        .build_from_path(bam_path)?;
    let header = reader.read_header()?;

    let region_start = Position::try_from(start.max(1) as usize)?;
    let region_end = Position::try_from(end.max(1) as usize)?;
    let region = noodles::core::Region::new(chrom, region_start..=region_end);

    let mut reads: HashSet<String> = HashSet::with_capacity_and_hasher(64, Default::default());

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

fn build_interned_bam_index(alignments: &[InternedReadMapping]) -> InternedBamIndex {
    let mut index: HbHashMap<u32, Vec<AlignmentSummary>> = HbHashMap::default();
    let mut by_chrom: HbHashMap<String, Vec<(u32, i64, i64)>> = HbHashMap::default();
    for aln in alignments {
        index.entry(aln.read_id).or_default().push(AlignmentSummary {
            chrom: aln.chrom.clone(),
            start: aln.start,
            end: aln.end,
            is_primary: aln.is_primary,
            is_secondary: false,
            is_supplementary: false,
        });
        by_chrom
            .entry(aln.chrom.clone())
            .or_default()
            .push((aln.read_id, aln.start, aln.end));
    }
    InternedBamIndex { index, by_chrom }
}

/// Convert a string-keyed in-memory BAM index into interned ids plus `InternedBamIndex`
/// so cached transitive detection can use Roaring bitmaps for read sets.
fn intern_bam_index(bam: &BamIndex) -> (ReadInterner, InternedBamIndex) {
    let mut interner = ReadInterner::new();
    let mut index: HbHashMap<u32, Vec<AlignmentSummary>> =
        HbHashMap::with_capacity(bam.index.len());
    let mut by_chrom: HbHashMap<String, Vec<(u32, i64, i64)>> = HbHashMap::default();
    for (read_name, alignments) in &bam.index {
        let rid = interner.intern(read_name);
        let entry = index.entry(rid).or_default();
        for aln in alignments {
            entry.push(AlignmentSummary {
                chrom: aln.chrom.clone(),
                start: aln.start,
                end: aln.end,
                is_primary: aln.is_primary,
                is_secondary: aln.is_secondary,
                is_supplementary: aln.is_supplementary,
            });
            by_chrom
                .entry(aln.chrom.clone())
                .or_default()
                .push((rid, aln.start, aln.end));
        }
    }
    (interner, InternedBamIndex { index, by_chrom })
}

fn collect_reads_in_region_interned(
    bam_index: &InternedBamIndex,
    chrom: &str,
    start: i64,
    end: i64,
) -> RoaringBitmap {
    let mut reads = RoaringBitmap::new();
    let Some(entries) = bam_index.by_chrom.get(chrom) else {
        return reads;
    };
    for &(read_id, a0, a1) in entries {
        if a0 < end && a1 > start {
            reads.insert(read_id);
        }
    }
    reads
}

fn find_all_alignments_interned_cached(
    bam_index: &InternedBamIndex,
    reads: &RoaringBitmap,
) -> Vec<InternedReadMapping> {
    let mut alignments = Vec::new();
    for read_id in reads.iter() {
        if let Some(alns) = bam_index.index.get(&read_id) {
            for aln in alns {
                alignments.push(InternedReadMapping {
                    read_id,
                    chrom: aln.chrom.clone(),
                    start: aln.start,
                    end: aln.end,
                    is_primary: aln.is_primary,
                });
            }
        }
    }
    alignments
}

/// Highly optimized: collect reads from multiple regions using a SINGLE BAM reader
/// Uses rayon for parallel queries within the same reader
struct SingleBamReader {
    reader: bam::io::indexed_reader::IndexedReader<bgzf::Reader<std::fs::File>>,
    header: std::sync::Mutex<sam::Header>,
}

impl SingleBamReader {
    fn new(bam_path: &str) -> Result<Self> {
        let mut reader = bam::io::indexed_reader::Builder::default()
            .build_from_path(bam_path)?;
        let header = reader.read_header()?;
        Ok(Self { reader, header: std::sync::Mutex::new(header) })
    }

    fn get_header(&self) -> sam::Header {
        self.header.lock().unwrap().clone()
    }

    fn get_ref_seqs(&self) -> sam::header::ReferenceSequences {
        self.header.lock().unwrap().reference_sequences().clone()
    }

    /// Query a single region and return read names
    fn query_region(&mut self, chrom: &str, start: i64, end: i64) -> HashSet<String> {
        let mut reads = HashSet::default();
        let header = self.get_header();
        let ref_seqs = self.get_ref_seqs();

        let chrom_name: bstr::BString = chrom.as_bytes().into();
        if let Some(tid) = ref_seqs.get_index_of(&chrom_name) {
            let start_pos = Position::try_from(start.max(1) as usize).unwrap_or(Position::MIN);
            let end_pos = Position::try_from(end.max(1) as usize).unwrap_or(Position::MAX);
            let region = noodles::core::Region::new(chrom.to_string(), start_pos..=end_pos);

            if let Ok(query) = self.reader.query(&header, &region) {
                for result in query.flatten() {
                    if let Some(qname) = result.name() {
                        if let Ok(name) = std::str::from_utf8(qname) {
                            reads.insert(name.to_string());
                        }
                    }
                }
            }
        }
        reads
    }

    /// Query multiple regions sequentially (reuses same file handle)
    fn query_regions(&mut self, regions: &[(String, i64, i64)]) -> Vec<HashSet<String>> {
        let mut results = Vec::with_capacity(regions.len());
        for (chrom, start, end) in regions {
            results.push(self.query_region(chrom, *start, *end));
        }
        results
    }
}

/// Ultra-fast detection: uses single BAM reader with parallel queries
/// This avoids rebuilding the in-memory index and uses the .bai for random access
pub fn detect_transitive_optimized(
    bam_path: &str,
    seed_chrom: &str,
    seed_start: i64,
    seed_end: i64,
    target_chroms: &[String],
    params: &TransitiveParams,
) -> Result<Vec<TransitiveLocus>> {
    println!("=== TRANSITIVE DETECTION (OPTIMIZED - SINGLE READER) ===");

    // Open a SINGLE BAM reader and reuse it
    let mut bam = SingleBamReader::new(bam_path)?;

    // Step 1: Collect reads at seed using the shared reader
    println!("Step 1: Collecting reads at seed region...");
    let seed_reads = bam.query_region(seed_chrom, seed_start, seed_end);
    println!("  Found {} reads at seed", seed_reads.len());

    if seed_reads.is_empty() {
        return Ok(vec![]);
    }

    // Step 2: Find all alignments of seed reads using shared reader
    println!("\nStep 2: Finding all alignments of seed reads...");
    let alignments = find_all_alignments_single_reader(&mut bam, &seed_reads, target_chroms)?;
    println!("  Found {} total alignments", alignments.len());

    // Step 3: Cluster into candidate loci
    println!("\nStep 3: Clustering alignments into candidate loci...");
    let candidate_loci = cluster_alignments(
        alignments,
        params.cluster_distance,
        params.min_reads_per_locus,
        params.adaptive_cluster_distance,
    );
    println!("  Found {} candidate loci", candidate_loci.len());

    // Step 4: Find seed locus
    let seed_locus_idx = candidate_loci.iter()
        .position(|(c, s, e, _)| c == seed_chrom && *s <= seed_end && *e >= seed_start);

    // Step 5: Expand read sets at each candidate locus using SINGLE reader
    println!("\nStep 4: Expanding read sets at each locus (single reader, parallel queries)...");
    let regions: Vec<(String, i64, i64)> = candidate_loci.iter()
        .map(|(c, s, e, _)| (c.clone(), *s, *e))
        .collect();

    // Query all regions using the single reader - much faster than opening new readers
    let all_reads = bam.query_regions(&regions);

    let mut expanded_loci: Vec<(String, i64, i64, HashSet<String>)> = Vec::new();
    for (i, (chrom, start, end, _)) in candidate_loci.iter().enumerate() {
        expanded_loci.push((chrom.clone(), *start, *end, all_reads[i].clone()));
        let is_seed = seed_locus_idx.map_or(false, |idx| idx == i);
        if is_seed {
            println!("  Locus {}: {}:{}-{} - {} reads [SEED]",
                     i + 1, chrom, start, end, all_reads[i].len());
        }
    }

    // Add seed region if not found
    let seed_locus_idx = if let Some(idx) = seed_locus_idx {
        idx
    } else {
        let seed_reads_full = bam.query_region(seed_chrom, seed_start, seed_end);
        expanded_loci.push((seed_chrom.to_string(), seed_start, seed_end, seed_reads_full));
        expanded_loci.len() - 1
    };

    // Step 6: Build graph and find reachable loci
    println!("\nStep 5: Building locus graph and finding reachable loci...");
    let adj = build_locus_graph(&expanded_loci, params.min_jaccard_for_edge);
    let reachable = find_reachable_loci(&adj, seed_locus_idx);
    println!("  Found {} reachable loci from seed", reachable.len());

    // Step 7: Create results
    let mut results: Vec<TransitiveLocus> = Vec::new();
    for (idx, distance, _) in reachable {
        let (chrom, start, end, reads) = &expanded_loci[idx];
        let direct_jaccard = calculate_jaccard(reads, &seed_reads);
        results.push(TransitiveLocus {
            chrom: chrom.clone(),
            start: *start,
            end: *end,
            reads: reads.clone(),
            component: distance + 1,
            jaccard_with_seed: direct_jaccard,
            distance_from_seed: distance,
            peak_count: 1,
            confidence: 1.0,
        });
    }

    results.sort_by(|a, b| {
        a.distance_from_seed.cmp(&b.distance_from_seed)
            .then(a.chrom.cmp(&b.chrom))
            .then(a.start.cmp(&b.start))
    });

    println!("\n=== DISCOVERED LOCI ===");
    for locus in &results {
        let label = if locus.distance_from_seed == 0 {
            "SEED"
        } else if locus.jaccard_with_seed > 0.0 {
            "DIRECT"
        } else {
            "TRANSITIVE"
        };
        println!("  {}: {}:{}-{} (dist={}, jacc={:.4}, reads={})",
                 label, locus.chrom, locus.start, locus.end,
                 locus.distance_from_seed, locus.jaccard_with_seed, locus.reads.len());
    }

    Ok(results)
}

/// Find all alignments using the shared BAM reader
fn find_all_alignments_single_reader(
    bam: &mut SingleBamReader,
    reads: &HashSet<String>,
    chroms: &[String],
) -> Result<Vec<ReadMapping>> {
    let mut alignments = Vec::new();
    let header = bam.get_header();
    let ref_seqs = bam.get_ref_seqs();

    for chrom in chroms {
        let chrom_name: bstr::BString = chrom.as_bytes().into();
        if let Some(tid) = ref_seqs.get_index_of(&chrom_name) {
            let ref_len = ref_seqs[tid].length().get();
            let start = Position::MIN;
            let end = Position::new(ref_len).unwrap_or(Position::MAX);
            let region = noodles::core::Region::new(chrom.clone(), start..=end);

            let query = bam.reader.query(&header, &region)?;
            for result in query {
                let record = result?;
                if record.flags().is_unmapped() {
                    continue;
                }
                let name = match record.name() {
                    Some(n) => std::str::from_utf8(n)?,
                    None => continue,
                };
                if !reads.contains(name) {
                    continue;
                }
                let aln_start = record.alignment_start()
                    .transpose()?
                    .map(|p: Position| p.get() as i64)
                    .unwrap_or(0);
                let aln_end = record.alignment_end()
                    .transpose()?
                    .map(|p: Position| p.get() as i64)
                    .unwrap_or(aln_start + 1);
                alignments.push(ReadMapping {
                    read_name: name.to_string(),
                    chrom: chrom.clone(),
                    start: aln_start,
                    end: aln_end,
                    is_primary: !record.flags().is_secondary() && !record.flags().is_supplementary(),
                });
            }
        }
    }
    Ok(alignments)
}

/// Optimized transitive detection using in-memory BAM index
/// This avoids repeated BAM file opens for each locus query.
/// Interns read names and uses Roaring bitmaps for locus read sets (same path as multi-seed).
pub fn detect_transitive_cached(
    bam_index: &BamIndex,
    seed_chrom: &str,
    seed_start: i64,
    seed_end: i64,
    _target_chroms: &[String], // Not needed with cached index
    params: &TransitiveParams,
) -> Vec<TransitiveLocus> {
    println!("=== TRANSITIVE DETECTION (CACHED INDEX) ===");
    println!("Interning BAM index, then running Roaring-backed graph + Jaccard...");
    let (interner, interned) = intern_bam_index(bam_index);
    let total_aln: usize = interned.index.values().map(|v| v.len()).sum();
    println!(
        "  Unique reads: {}, alignments: {}",
        interner.len(),
        total_aln
    );

    let (results, _) = detect_transitive_with_component_candidates_cached_interned(
        &interned,
        &interner,
        seed_chrom,
        seed_start,
        seed_end,
        _target_chroms,
        params,
    );

    println!("\n=== DISCOVERED LOCI ===");
    for locus in &results {
        let label = if locus.distance_from_seed == 0 {
            "SEED"
        } else if locus.jaccard_with_seed > 0.0 {
            "DIRECT"
        } else {
            "TRANSITIVE"
        };
        println!(
            "  {}: {}:{}-{} (dist={}, jacc={:.4}, reads={})",
            label,
            locus.chrom,
            locus.start,
            locus.end,
            locus.distance_from_seed,
            locus.jaccard_with_seed,
            locus.reads.len()
        );
    }

    results
}

/// Find ALL alignments for a set of reads (primary + secondary) - SLOW version
fn find_all_alignments(
    bam_path: &str,
    reads: &HashSet<String>,
    chroms: &[String],
) -> Result<Vec<ReadMapping>> {
    let mut alignments = Vec::new();

    for chrom in chroms {
        let mut reader = bam::io::indexed_reader::Builder::default()
            .build_from_path(bam_path)?;
        let header = reader.read_header()?;

        let chrom_name: bstr::BString = chrom.as_bytes().into();
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
                    Ok(s) => s,
                    Err(_) => continue,
                },
                None => continue,
            };

            if !reads.contains(name) {
                continue;
            }

            let aln_start = record.alignment_start()
                .transpose()?
                .map(|p: Position| p.get() as i64)
                .unwrap_or(0);
            let aln_end = record.alignment_end()
                .transpose()?
                .map(|p: Position| p.get() as i64)
                .unwrap_or(aln_start + 1);

            alignments.push(ReadMapping {
                read_name: name.to_string(),
                chrom: chrom.clone(),
                start: aln_start,
                end: aln_end,
                is_primary: !record.flags().is_secondary() && !record.flags().is_supplementary(),
            });
        }
    }

    Ok(alignments)
}

/// Merge gap cap from multi-mapping geometry on one chromosome: positive gaps between
/// consecutive same-read alignments (sorted by start). High percentile approximates
/// "within-copy" spacing; capped by `cluster_distance` so the CLI remains a hard ceiling.
fn adaptive_merge_distance_for_chrom(aligns: &[ReadMapping], cluster_distance: i64) -> i64 {
    const MIN_MERGE_GAP: i64 = 500;
    const MIN_GAPS_FOR_STATS: usize = 3;

    let mut by_read: HashMap<&str, Vec<usize>> = HashMap::default();
    for (i, a) in aligns.iter().enumerate() {
        by_read.entry(a.read_name.as_str()).or_default().push(i);
    }

    let mut gaps: Vec<i64> = Vec::new();
    for idxs in by_read.values() {
        if idxs.len() < 2 {
            continue;
        }
        let mut order: Vec<usize> = idxs.clone();
        order.sort_by_key(|&i| aligns[i].start);
        for w in order.windows(2) {
            let g = aligns[w[1]].start - aligns[w[0]].end;
            if g > 0 {
                gaps.push(g);
            }
        }
    }

    if gaps.len() < MIN_GAPS_FOR_STATS {
        return cluster_distance;
    }

    gaps.sort_unstable();
    let idx = ((gaps.len() - 1) as f64 * 0.90).round() as usize;
    let idx = idx.min(gaps.len() - 1);
    let p90 = gaps[idx];
    p90.min(cluster_distance).max(MIN_MERGE_GAP)
}

fn adaptive_merge_distance_for_chrom_interned(
    aligns: &[InternedReadMapping],
    cluster_distance: i64,
) -> i64 {
    const MIN_MERGE_GAP: i64 = 500;
    const MIN_GAPS_FOR_STATS: usize = 3;

    let mut by_read: HbHashMap<u32, Vec<usize>> = HbHashMap::default();
    for (i, a) in aligns.iter().enumerate() {
        by_read.entry(a.read_id).or_default().push(i);
    }

    let mut gaps: Vec<i64> = Vec::new();
    for idxs in by_read.values() {
        if idxs.len() < 2 {
            continue;
        }
        let mut order: Vec<usize> = idxs.clone();
        order.sort_by_key(|&i| aligns[i].start);
        for w in order.windows(2) {
            let g = aligns[w[1]].start - aligns[w[0]].end;
            if g > 0 {
                gaps.push(g);
            }
        }
    }

    if gaps.len() < MIN_GAPS_FOR_STATS {
        return cluster_distance;
    }

    gaps.sort_unstable();
    let idx = ((gaps.len() - 1) as f64 * 0.90).round() as usize;
    let idx = idx.min(gaps.len() - 1);
    let p90 = gaps[idx];
    p90.min(cluster_distance).max(MIN_MERGE_GAP)
}

/// Cluster alignments into loci
fn cluster_alignments(
    alignments: Vec<ReadMapping>,
    cluster_distance: i64,
    min_reads: usize,
    adaptive_merge_gap: bool,
) -> Vec<(String, i64, i64, HashSet<String>)> {
    if alignments.is_empty() {
        return vec![];
    }

    // Group by chromosome
    let mut by_chrom: HashMap<String, Vec<ReadMapping>> = HashMap::new();
    for aln in alignments {
        by_chrom.entry(aln.chrom.clone()).or_default().push(aln);
    }

    let mut loci: Vec<(String, i64, i64, HashSet<String>)> = vec![];

    for (chrom, mut aligns) in by_chrom {
        aligns.sort_by_key(|a| a.start);

        let merge_cap = if adaptive_merge_gap {
            adaptive_merge_distance_for_chrom(&aligns, cluster_distance)
        } else {
            cluster_distance
        };

        let mut current_cluster: Vec<ReadMapping> = vec![aligns[0].clone()];
        let mut current_end = aligns[0].end;

        for aln in aligns.into_iter().skip(1) {
            let aln_start = aln.start;
            let aln_end = aln.end;
            
            if aln_start - current_end <= merge_cap {
                current_cluster.push(aln);
                if aln_end > current_end {
                    current_end = aln_end;
                }
            } else {
                if current_cluster.len() >= min_reads {
                    let reads: HashSet<String> = current_cluster
                        .iter()
                        .map(|a| a.read_name.clone())
                        .collect();
                    let start = current_cluster.iter().map(|a| a.start).min().unwrap_or(0);
                    loci.push((chrom.clone(), start, current_end, reads));
                }
                current_cluster = vec![aln];
                current_end = aln_end;
            }
        }

        if current_cluster.len() >= min_reads {
            let reads: HashSet<String> = current_cluster
                .iter()
                .map(|a| a.read_name.clone())
                .collect();
            let start = current_cluster.iter().map(|a| a.start).min().unwrap_or(0);
            loci.push((chrom, start, current_end, reads));
        }
    }

    loci
}

fn cluster_alignments_interned(
    alignments: Vec<InternedReadMapping>,
    cluster_distance: i64,
    min_reads: usize,
    adaptive_merge_gap: bool,
) -> Vec<(String, i64, i64, RoaringBitmap)> {
    if alignments.is_empty() {
        return vec![];
    }

    let mut by_chrom: HashMap<String, Vec<InternedReadMapping>> = HashMap::new();
    for aln in alignments {
        by_chrom.entry(aln.chrom.clone()).or_default().push(aln);
    }

    let mut loci: Vec<(String, i64, i64, RoaringBitmap)> = vec![];

    for (chrom, mut aligns) in by_chrom {
        aligns.sort_by_key(|a| a.start);

        let merge_cap = if adaptive_merge_gap {
            adaptive_merge_distance_for_chrom_interned(&aligns, cluster_distance)
        } else {
            cluster_distance
        };

        let mut current_cluster: Vec<InternedReadMapping> = vec![aligns[0].clone()];
        let mut current_end = aligns[0].end;

        for aln in aligns.into_iter().skip(1) {
            let aln_start = aln.start;
            let aln_end = aln.end;

            if aln_start - current_end <= merge_cap {
                current_cluster.push(aln);
                if aln_end > current_end {
                    current_end = aln_end;
                }
            } else {
                if current_cluster.len() >= min_reads {
                    let mut reads_bm = RoaringBitmap::new();
                    for a in &current_cluster {
                        reads_bm.insert(a.read_id);
                    }
                    let start = current_cluster.iter().map(|a| a.start).min().unwrap_or(0);
                    loci.push((chrom.clone(), start, current_end, reads_bm));
                }
                current_cluster = vec![aln];
                current_end = aln_end;
            }
        }

        if current_cluster.len() >= min_reads {
            let mut reads_bm = RoaringBitmap::new();
            for a in &current_cluster {
                reads_bm.insert(a.read_id);
            }
            let start = current_cluster.iter().map(|a| a.start).min().unwrap_or(0);
            loci.push((chrom, start, current_end, reads_bm));
        }
    }

    loci
}

/// Calculate Jaccard similarity between two read sets
fn calculate_jaccard(set1: &HashSet<String>, set2: &HashSet<String>) -> f64 {
    if set1.is_empty() || set2.is_empty() {
        return 0.0;
    }

    let (smaller, larger) = if set1.len() < set2.len() {
        (set1, set2)
    } else {
        (set2, set1)
    };

    let intersection: usize = smaller.iter().filter(|x| larger.contains(*x)).count();
    let union = set1.len() + set2.len() - intersection;

    if union > 0 {
        intersection as f64 / union as f64
    } else {
        0.0
    }
}

fn calculate_jaccard_roaring(set1: &RoaringBitmap, set2: &RoaringBitmap) -> f64 {
    let a_len = set1.len() as usize;
    let b_len = set2.len() as usize;
    if a_len == 0 || b_len == 0 {
        return 0.0;
    }
    let inter = set1.intersection_len(set2) as usize;
    let union = a_len + b_len - inter;
    if union > 0 {
        inter as f64 / union as f64
    } else {
        0.0
    }
}

/// Build a graph where nodes are loci and edges represent shared reads
fn build_locus_graph(
    loci: &[(String, i64, i64, HashSet<String>)],
    min_jaccard: f64,
) -> Vec<Vec<(usize, f64)>> {
    let n = loci.len();
    let mut adj: Vec<Vec<(usize, f64)>> = vec![vec![]; n];

    for i in 0..n {
        for j in (i + 1)..n {
            let jaccard = calculate_jaccard(&loci[i].3, &loci[j].3);
            if jaccard >= min_jaccard {
                adj[i].push((j, jaccard));
                adj[j].push((i, jaccard));
            }
        }
    }

    adj
}

/// Pairwise Jaccard over loci; parallel when there are enough nodes to amortize Rayon.
fn build_locus_graph_roaring(
    loci: &[(String, i64, i64, RoaringBitmap)],
    min_jaccard: f64,
) -> Vec<Vec<(usize, f64)>> {
    let n = loci.len();
    let mut adj: Vec<Vec<(usize, f64)>> = vec![vec![]; n];
    if n < 2 {
        return adj;
    }

    const PARALLEL_MIN_N: usize = 24;

    if n >= PARALLEL_MIN_N {
        let edge_lists: Vec<Vec<(usize, usize, f64)>> = (0..n)
            .into_par_iter()
            .map(|i| {
                let mut local: Vec<(usize, usize, f64)> = Vec::new();
                for j in (i + 1)..n {
                    let jaccard = calculate_jaccard_roaring(&loci[i].3, &loci[j].3);
                    if jaccard >= min_jaccard {
                        local.push((i, j, jaccard));
                    }
                }
                local
            })
            .collect();
        for row in edge_lists {
            for (i, j, jac) in row {
                adj[i].push((j, jac));
                adj[j].push((i, jac));
            }
        }
    } else {
        for i in 0..n {
            for j in (i + 1)..n {
                let jaccard = calculate_jaccard_roaring(&loci[i].3, &loci[j].3);
                if jaccard >= min_jaccard {
                    adj[i].push((j, jaccard));
                    adj[j].push((i, jaccard));
                }
            }
        }
    }

    adj
}

/// Find all loci reachable from the seed locus using BFS
fn find_reachable_loci(
    adj: &[Vec<(usize, f64)>],
    seed_idx: usize,
) -> Vec<(usize, usize, f64)> {
    // (locus_idx, distance, max_jaccard_on_path)
    let mut result: Vec<(usize, usize, f64)> = vec![];
    let mut visited: HbHashMap<usize, (usize, f64)> = HbHashMap::default(); // idx -> (distance, max_jaccard)
    let mut queue: VecDeque<(usize, usize, f64)> = VecDeque::new(); // (idx, distance, min_jaccard_on_path)

    queue.push_back((seed_idx, 0, 1.0));
    visited.insert(seed_idx, (0, 1.0));

    while let Some((idx, dist, min_jaccard)) = queue.pop_front() {
        result.push((idx, dist, min_jaccard));

        for (neighbor, edge_jaccard) in &adj[idx] {
            let new_dist = dist + 1;
            let new_min_jaccard = min_jaccard.min(*edge_jaccard);

            if let Some(&(existing_dist, existing_jaccard)) = visited.get(neighbor) {
                // If we found a shorter path, or same distance but higher jaccard, update
                if new_dist < existing_dist || (new_dist == existing_dist && new_min_jaccard > existing_jaccard) {
                    visited.insert(*neighbor, (new_dist, new_min_jaccard));
                    queue.push_back((*neighbor, new_dist, new_min_jaccard));
                }
            } else {
                visited.insert(*neighbor, (new_dist, new_min_jaccard));
                queue.push_back((*neighbor, new_dist, new_min_jaccard));
            }
        }
    }

    result
}

/// Find gene family members using transitive detection
///
/// This is the main entry point for the transitive detector.
pub fn detect_transitive(
    bam_path: &str,
    seed_chrom: &str,
    seed_start: i64,
    seed_end: i64,
    target_chroms: &[String],
    params: &TransitiveParams,
) -> Result<Vec<TransitiveLocus>> {
    println!("=== TRANSITIVE DETECTION MODE ===");
    println!("For highly diverged gene families (e.g., ID_14 LRRC37)");
    println!("where paralogs don't share multi-mapping reads");
    println!();

    // Step 1: Collect all reads overlapping the seed
    println!("Step 1: Collecting reads overlapping seed region...");
    let seed_reads = collect_reads_in_region(bam_path, seed_chrom, seed_start, seed_end)?;
    println!("  Found {} reads at seed", seed_reads.len());

    if seed_reads.is_empty() {
        println!("  No reads at seed, aborting");
        return Ok(vec![]);
    }

    // Step 2: Find all alignments of seed reads genome-wide
    println!("\nStep 2: Finding all alignments of seed reads...");
    let alignments = find_all_alignments(bam_path, &seed_reads, target_chroms)?;
    println!("  Found {} total alignments", alignments.len());

    let primary_count = alignments.iter().filter(|a| a.is_primary).count();
    let secondary_count = alignments.iter().filter(|a| !a.is_primary).count();
    println!("    Primary: {}", primary_count);
    println!("    Secondary: {}", secondary_count);

    // Step 3: Cluster into candidate loci
    println!("\nStep 3: Clustering alignments into candidate loci...");
    let candidate_loci = cluster_alignments(
        alignments,
        params.cluster_distance,
        params.min_reads_per_locus,
        params.adaptive_cluster_distance,
    );
    println!("  Found {} candidate loci", candidate_loci.len());

    // Step 4: Find seed locus
    let seed_locus_idx = candidate_loci.iter()
        .position(|(c, s, e, _)| c == seed_chrom && *s <= seed_end && *e >= seed_start);

    if seed_locus_idx.is_none() {
        println!("  Warning: Seed locus not found in candidates");
        // Add seed as a candidate
        println!("  Adding seed region as a candidate");
    }

    // Step 5: Expand read sets at each candidate locus
    println!("\nStep 4: Expanding read sets at each locus...");
    let mut expanded_loci: Vec<(String, i64, i64, HashSet<String>)> = Vec::new();
    
    for (i, (chrom, start, end, _)) in candidate_loci.iter().enumerate() {
        let locus_reads = collect_reads_in_region(bam_path, chrom, *start, *end)?;
        expanded_loci.push((chrom.clone(), *start, *end, locus_reads));
        
        let is_seed = seed_locus_idx.map_or(false, |idx| idx == i);
        if is_seed {
            println!("  Locus {}: {}:{}-{} - {} reads [SEED]", 
                     i + 1, chrom, start, end, expanded_loci.last().unwrap().3.len());
        }
    }
    
    // Add seed region if not found
    let seed_locus_idx = if let Some(idx) = seed_locus_idx {
        idx
    } else {
        let seed_reads_full = collect_reads_in_region(bam_path, seed_chrom, seed_start, seed_end)?;
        expanded_loci.push((seed_chrom.to_string(), seed_start, seed_end, seed_reads_full));
        expanded_loci.len() - 1
    };

    // Step 6: Build graph and find reachable loci
    println!("\nStep 5: Building locus graph and finding reachable loci...");
    let adj = build_locus_graph(&expanded_loci, params.min_jaccard_for_edge);

    // Count edges
    let edge_count: usize = adj.iter().map(|v| v.len()).sum();
    println!("  Graph: {} nodes, {} edges", expanded_loci.len(), edge_count / 2);

    let reachable = find_reachable_loci(&adj, seed_locus_idx);
    println!("  Found {} reachable loci from seed", reachable.len());

    // Step 7: Create TransitiveLocus results
    let mut results: Vec<TransitiveLocus> = Vec::new();
    
    for (idx, distance, _path_jaccard) in reachable {
        let (chrom, start, end, reads) = &expanded_loci[idx];
        let direct_jaccard = calculate_jaccard(reads, &seed_reads);

        results.push(TransitiveLocus {
            chrom: chrom.clone(),
            start: *start,
            end: *end,
            reads: reads.clone(),
            component: distance + 1,
            jaccard_with_seed: direct_jaccard,
            distance_from_seed: distance,
            peak_count: 1,
            confidence: 1.0,
        });
    }

    // Sort by distance then by position
    results.sort_by(|a, b| {
        a.distance_from_seed.cmp(&b.distance_from_seed)
            .then(a.chrom.cmp(&b.chrom))
            .then(a.start.cmp(&b.start))
    });

    // Print summary
    println!("\n=== DISCOVERED LOCI ===");
    for (i, locus) in results.iter().enumerate() {
        let label = if locus.distance_from_seed == 0 {
            "SEED"
        } else if locus.jaccard_with_seed > 0.0 {
            "direct"
        } else {
            "transitive"
        };
        
        println!("  Locus {}: {}:{}-{} ({} reads, Jaccard: {:.4}, dist: {}, {})",
                 i + 1,
                 locus.chrom,
                 locus.start,
                 locus.end,
                 locus.reads.len(),
                 locus.jaccard_with_seed,
                 locus.distance_from_seed,
                 label);
    }

    println!("\n=== TRANSITIVE DETECTION COMPLETE ===");
    println!("Total loci discovered: {}", results.len());
    let direct = results.iter().filter(|l| l.distance_from_seed == 1 && l.jaccard_with_seed > 0.0).count();
    let transitive = results.iter().filter(|l| l.distance_from_seed >= 1 && l.jaccard_with_seed == 0.0).count();
    println!("  Direct neighbors: {}", direct);
    println!("  Transitive (no shared reads): {}", transitive);

    Ok(results)
}

/// Detect transitive loci and expose disconnected component seed candidates.
///
/// The first return value is identical in behavior to `detect_transitive`:
/// loci reachable from the input seed.
/// The second return value contains one representative candidate per
/// disconnected component for optional non-circular reseeding.
pub fn detect_transitive_with_component_candidates(
    bam_path: &str,
    seed_chrom: &str,
    seed_start: i64,
    seed_end: i64,
    target_chroms: &[String],
    params: &TransitiveParams,
) -> Result<(Vec<TransitiveLocus>, Vec<ComponentSeedCandidate>)> {
    // Step 1: Collect all reads overlapping the seed
    let seed_reads = collect_reads_in_region(bam_path, seed_chrom, seed_start, seed_end)?;
    if seed_reads.is_empty() {
        return Ok((vec![], vec![]));
    }

    // Step 2: Find all alignments of seed reads genome-wide
    let alignments = find_all_alignments(bam_path, &seed_reads, target_chroms)?;

    // Step 3: Cluster into candidate loci
    let candidate_loci = cluster_alignments(
        alignments,
        params.cluster_distance,
        params.min_reads_per_locus,
        params.adaptive_cluster_distance,
    );
    if candidate_loci.is_empty() {
        return Ok((vec![], vec![]));
    }

    // Step 4: Find seed locus in candidates
    let seed_locus_idx = candidate_loci.iter()
        .position(|(c, s, e, _)| c == seed_chrom && *s <= seed_end && *e >= seed_start);

    // Step 5: Expand read sets at each candidate locus
    let mut expanded_loci: Vec<(String, i64, i64, HashSet<String>)> = Vec::new();
    for (chrom, start, end, _) in &candidate_loci {
        let locus_reads = collect_reads_in_region(bam_path, chrom, *start, *end)?;
        expanded_loci.push((chrom.clone(), *start, *end, locus_reads));
    }

    // Add seed region if it was not found among candidates
    let seed_locus_idx = if let Some(idx) = seed_locus_idx {
        idx
    } else {
        let seed_reads_full = collect_reads_in_region(bam_path, seed_chrom, seed_start, seed_end)?;
        expanded_loci.push((seed_chrom.to_string(), seed_start, seed_end, seed_reads_full));
        expanded_loci.len() - 1
    };

    // Step 6: Graph and reachable loci
    let adj = build_locus_graph(&expanded_loci, params.min_jaccard_for_edge);
    let reachable = find_reachable_loci(&adj, seed_locus_idx);

    let mut results: Vec<TransitiveLocus> = Vec::new();
    for (idx, distance, _path_jaccard) in reachable {
        let (chrom, start, end, reads) = &expanded_loci[idx];
        let direct_jaccard = calculate_jaccard(reads, &seed_reads);
        results.push(TransitiveLocus {
            chrom: chrom.clone(),
            start: *start,
            end: *end,
            reads: reads.clone(),
            component: distance + 1,
            jaccard_with_seed: direct_jaccard,
            distance_from_seed: distance,
            peak_count: 1,
            confidence: 1.0,
        });
    }
    results.sort_by(|a, b| {
        a.distance_from_seed.cmp(&b.distance_from_seed)
            .then(a.chrom.cmp(&b.chrom))
            .then(a.start.cmp(&b.start))
    });

    // Step 7: Connected components and disconnected representatives
    let n = expanded_loci.len();
    let mut comp_id: Vec<usize> = vec![usize::MAX; n];
    let mut current_comp = 0usize;

    for i in 0..n {
        if comp_id[i] != usize::MAX {
            continue;
        }
        let mut queue: VecDeque<usize> = VecDeque::new();
        queue.push_back(i);
        comp_id[i] = current_comp;
        while let Some(node) = queue.pop_front() {
            for (nbr, _) in &adj[node] {
                if comp_id[*nbr] == usize::MAX {
                    comp_id[*nbr] = current_comp;
                    queue.push_back(*nbr);
                }
            }
        }
        current_comp += 1;
    }

    let seed_component = comp_id[seed_locus_idx];
    let mut nodes_by_comp: HashMap<usize, Vec<usize>> = HashMap::new();
    for (idx, cid) in comp_id.iter().enumerate() {
        nodes_by_comp.entry(*cid).or_default().push(idx);
    }

    let mut candidates: Vec<ComponentSeedCandidate> = Vec::new();
    for (cid, nodes) in nodes_by_comp {
        if cid == seed_component || nodes.is_empty() {
            continue;
        }

        // Representative: locus with largest read support in this component.
        let mut best_idx = nodes[0];
        let mut best_reads = expanded_loci[best_idx].3.len();
        let mut total_reads = 0usize;

        for idx in &nodes {
            let n_reads = expanded_loci[*idx].3.len();
            total_reads += n_reads;
            if n_reads > best_reads {
                best_reads = n_reads;
                best_idx = *idx;
            }
        }

        let (chrom, start, end, reads) = &expanded_loci[best_idx];
        candidates.push(ComponentSeedCandidate {
            chrom: chrom.clone(),
            start: *start,
            end: *end,
            reads: reads.clone(),
            component_size: nodes.len(),
            total_reads_in_component: total_reads,
        });
    }

    candidates.sort_by(|a, b| {
        b.total_reads_in_component
            .cmp(&a.total_reads_in_component)
            .then(b.component_size.cmp(&a.component_size))
            .then(b.reads.len().cmp(&a.reads.len()))
    });

    Ok((results, candidates))
}

/// Cached transitive detection: Roaring bitmap read sets, interned ids, `intersection_len` Jaccard.
fn detect_transitive_with_component_candidates_cached_interned(
    bam_index: &InternedBamIndex,
    interner: &ReadInterner,
    seed_chrom: &str,
    seed_start: i64,
    seed_end: i64,
    _target_chroms: &[String],
    params: &TransitiveParams,
) -> (Vec<TransitiveLocus>, Vec<ComponentSeedCandidate>) {
    let seed_reads = collect_reads_in_region_interned(bam_index, seed_chrom, seed_start, seed_end);
    if seed_reads.is_empty() {
        return (vec![], vec![]);
    }

    let alignments = find_all_alignments_interned_cached(bam_index, &seed_reads);

    let candidate_loci = cluster_alignments_interned(
        alignments,
        params.cluster_distance,
        params.min_reads_per_locus,
        params.adaptive_cluster_distance,
    );
    if candidate_loci.is_empty() {
        return (vec![], vec![]);
    }

    let seed_locus_idx = candidate_loci
        .iter()
        .position(|(c, s, e, _)| c == seed_chrom && *s <= seed_end && *e >= seed_start);

    let mut expanded_loci: Vec<(String, i64, i64, RoaringBitmap)> = Vec::new();
    for (chrom, start, end, _) in &candidate_loci {
        let locus_reads = collect_reads_in_region_interned(bam_index, chrom, *start, *end);
        expanded_loci.push((chrom.clone(), *start, *end, locus_reads));
    }

    let seed_locus_idx = if let Some(idx) = seed_locus_idx {
        idx
    } else {
        let seed_reads_full =
            collect_reads_in_region_interned(bam_index, seed_chrom, seed_start, seed_end);
        expanded_loci.push((
            seed_chrom.to_string(),
            seed_start,
            seed_end,
            seed_reads_full,
        ));
        expanded_loci.len() - 1
    };

    let adj = build_locus_graph_roaring(&expanded_loci, params.min_jaccard_for_edge);
    let reachable = find_reachable_loci(&adj, seed_locus_idx);

    let mut results: Vec<TransitiveLocus> = Vec::new();
    for (idx, distance, _path_jaccard) in reachable {
        let (chrom, start, end, reads) = &expanded_loci[idx];
        let direct_jaccard = calculate_jaccard_roaring(reads, &seed_reads);
        results.push(TransitiveLocus {
            chrom: chrom.clone(),
            start: *start,
            end: *end,
            reads: interner.bitmap_to_hashset(reads),
            component: distance + 1,
            jaccard_with_seed: direct_jaccard,
            distance_from_seed: distance,
            peak_count: 1,
            confidence: 1.0,
        });
    }
    results.sort_by(|a, b| {
        a.distance_from_seed
            .cmp(&b.distance_from_seed)
            .then(a.chrom.cmp(&b.chrom))
            .then(a.start.cmp(&b.start))
    });

    let n = expanded_loci.len();
    let mut comp_id: Vec<usize> = vec![usize::MAX; n];
    let mut current_comp = 0usize;

    for i in 0..n {
        if comp_id[i] != usize::MAX {
            continue;
        }
        let mut queue: VecDeque<usize> = VecDeque::new();
        queue.push_back(i);
        comp_id[i] = current_comp;
        while let Some(node) = queue.pop_front() {
            for (nbr, _) in &adj[node] {
                if comp_id[*nbr] == usize::MAX {
                    comp_id[*nbr] = current_comp;
                    queue.push_back(*nbr);
                }
            }
        }
        current_comp += 1;
    }

    let seed_component = comp_id[seed_locus_idx];
    let mut nodes_by_comp: HashMap<usize, Vec<usize>> = HashMap::new();
    for (idx, cid) in comp_id.iter().enumerate() {
        nodes_by_comp.entry(*cid).or_default().push(idx);
    }

    let mut candidates: Vec<ComponentSeedCandidate> = Vec::new();
    for (cid, nodes) in nodes_by_comp {
        if cid == seed_component || nodes.is_empty() {
            continue;
        }

        let mut best_idx = nodes[0];
        let mut best_reads = expanded_loci[best_idx].3.len() as usize;
        let mut total_reads = 0usize;

        for idx in &nodes {
            let n_reads = expanded_loci[*idx].3.len() as usize;
            total_reads += n_reads;
            if n_reads > best_reads {
                best_reads = n_reads;
                best_idx = *idx;
            }
        }

        let (chrom, start, end, reads) = &expanded_loci[best_idx];
        candidates.push(ComponentSeedCandidate {
            chrom: chrom.clone(),
            start: *start,
            end: *end,
            reads: interner.bitmap_to_hashset(reads),
            component_size: nodes.len(),
            total_reads_in_component: total_reads,
        });
    }

    candidates.sort_by(|a, b| {
        b.total_reads_in_component
            .cmp(&a.total_reads_in_component)
            .then(b.component_size.cmp(&a.component_size))
            .then(b.reads.len().cmp(&a.reads.len()))
    });

    (results, candidates)
}

/// Optimized version using in-memory BAM index (builds index once, uses for all locus queries).
/// Interns read names and uses Roaring bitmaps for locus read sets and pairwise Jaccard.
pub fn detect_transitive_with_component_candidates_cached(
    bam_index: &BamIndex,
    seed_chrom: &str,
    seed_start: i64,
    seed_end: i64,
    _target_chroms: &[String], // Not needed with cached index
    params: &TransitiveParams,
) -> (Vec<TransitiveLocus>, Vec<ComponentSeedCandidate>) {
    let (interner, interned) = intern_bam_index(bam_index);
    detect_transitive_with_component_candidates_cached_interned(
        &interned,
        &interner,
        seed_chrom,
        seed_start,
        seed_end,
        _target_chroms,
        params,
    )
}

/// Process multiple seeds efficiently by scanning only target chromosomes.
/// Uses seed chromosomes as the target - finds all family members on those chromosomes.
pub fn detect_transitive_multi_seed_and_reads(
    bam_path: &str,
    seeds: &[(String, i64, i64)],
    target_chroms: &[String],
    params: &TransitiveParams,
) -> Result<(Vec<TransitiveLocus>, Vec<ComponentSeedCandidate>, HashSet<String>)> {
    println!("\n=== MULTI-SEED TRANSITIVE DETECTION ===");
    println!("Step 1: Collecting reads from {} seeds...", seeds.len());
    
    // Step 1: Collect all reads from all seeds
    let mut all_seed_reads: HashSet<String> = HashSet::default();
    
    for (chrom, start, end) in seeds {
        match collect_reads_in_region(bam_path, chrom, *start, *end) {
            Ok(reads) => {
                println!("  Seed {}:{}-{}: {} reads", chrom, start, end, reads.len());
                all_seed_reads.extend(reads);
            }
            Err(e) => {
                eprintln!("Warning: Failed to collect reads from {}:{}-{}: {}", chrom, start, end, e);
            }
        }
    }
    
    println!("  Total unique reads across all seeds: {}", all_seed_reads.len());
    
    println!("\nStep 2: Finding all alignments of seed reads (scanning {} chromosomes)...", target_chroms.len());
    
    // Step 2: Find ALL alignments of these reads on target chromosomes only
    let all_alignments = find_all_alignments_on_chroms(bam_path, &all_seed_reads, target_chroms)?;
    
    println!("  Found {} total alignments", all_alignments.len());
    
    // Step 3: Build lightweight index from just these reads (interned ids + id-keyed index)
    println!("\nStep 3: Building lightweight index from seed reads...");
    let mut interner = ReadInterner::new();
    let interned_alignments: Vec<InternedReadMapping> = all_alignments
        .iter()
        .map(|a| InternedReadMapping {
            read_id: interner.intern(&a.read_name),
            chrom: a.chrom.clone(),
            start: a.start,
            end: a.end,
            is_primary: a.is_primary,
        })
        .collect();
    let bam_index = build_interned_bam_index(&interned_alignments);
    println!("  Index built: {} reads indexed", bam_index.index.len());
    
    // Step 4: Process each seed using the shared index
    println!("\nStep 4: Processing each seed with shared index...");
    
    let mut all_loci: Vec<TransitiveLocus> = Vec::new();
    
    for (seed_idx, (seed_chrom, seed_start, seed_end)) in seeds.iter().enumerate() {
        println!("\n  [Seed {}/{}] {}:{}-{}", seed_idx + 1, seeds.len(), seed_chrom, seed_start, seed_end);
        
        let (loci, _) = detect_transitive_with_component_candidates_cached_interned(
            &bam_index,
            &interner,
            seed_chrom,
            *seed_start,
            *seed_end,
            &[],
            params,
        );
        
        println!("    Found {} loci", loci.len());
        
        // Merge with existing results, keeping LARGER locus when overlapping
        for locus in loci {
            let locus_span = locus.end - locus.start;
            let mut should_add = true;
            let mut to_remove: Option<usize> = None;

            for (idx, existing) in all_loci.iter().enumerate() {
                if existing.chrom == locus.chrom {
                    let overlap_start = existing.start.max(locus.start);
                    let overlap_end = existing.end.min(locus.end);
                    if overlap_start < overlap_end {
                        let existing_span = existing.end - existing.start;
                        if locus_span > existing_span {
                            to_remove = Some(idx);
                        } else {
                            should_add = false;
                        }
                        break;
                    }
                }
            }

            if let Some(idx) = to_remove {
                all_loci.remove(idx);
            }
            if should_add {
                all_loci.push(locus);
            }
        }
    }
    
    println!("\n=== DETECTION COMPLETE ===");
    println!("Total merged loci: {}", all_loci.len());
    
    // Sort by chromosome and position
    all_loci.sort_by(|a, b| a.chrom.cmp(&b.chrom).then(a.start.cmp(&b.start)));
    
    Ok((all_loci, vec![], all_seed_reads))
}

/// Find all alignments of given reads on specified chromosomes (faster than genome-wide)
fn find_all_alignments_on_chroms(
    bam_path: &str,
    read_names: &HashSet<String>,
    chroms: &[String],
) -> Result<Vec<ReadMapping>> {
    use noodles::bam;
    use noodles::core::Position;
    
    let mut reader = bam::io::indexed_reader::Builder::default()
        .build_from_path(bam_path)?;
    let header = reader.read_header()?;
    
    let mut alignments = Vec::new();
    let read_set: HashSet<&str> = read_names.iter().map(|s| s.as_str()).collect();
    
    // Only scan specified chromosomes
    for chrom in chroms {
        let ref_seqs = header.reference_sequences();
        let chrom_name: bstr::BString = chrom.as_bytes().into();
        
        if let Some(tid) = ref_seqs.get_index_of(&chrom_name) {
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
                    Ok(s) => s,
                    Err(_) => continue,
                };
                
                // Only include reads that are in our seed set
                if !read_set.contains(read_name) {
                    continue;
                }
                
                let aln_start = record.alignment_start()
                    .transpose()?
                    .map(|p: Position| p.get() as i64)
                    .unwrap_or(0);
                let aln_end = record.alignment_end()
                    .transpose()?
                    .map(|p: Position| p.get() as i64)
                    .unwrap_or(aln_start + 1);
                
                let flags = record.flags();
                
                alignments.push(ReadMapping {
                    read_name: read_name.to_string(),
                    chrom: chrom.clone(),
                    start: aln_start,
                    end: aln_end,
                    is_primary: !flags.is_secondary() && !flags.is_supplementary(),
                });
            }
        }
    }
    
    Ok(alignments)
}

/// Overlap length between two half-open genomic intervals on the same chromosome.
fn bp_overlap_same_chr(a0: i64, a1: i64, b0: i64, b1: i64) -> i64 {
    let s = a0.max(b0);
    let e = a1.min(b1);
    (e - s).max(0)
}

/// Split large loci using coverage valley detection.
/// This helps separate overmerged regions and trim false positives.
/// Seed regions are always preserved in the output.
pub fn split_loci_by_coverage(
    bam_path: &str,
    loci: &[TransitiveLocus],
    seed_reads: &HashSet<String>,
    split_cfg: &TransitiveCoverageSplitConfig,
    seed_regions: &[(String, i64, i64)],
) -> Result<Vec<TransitiveLocus>> {
    use crate::coverage_valley::{
        analyze_valleys, build_coverage_profile, merge_small_segments, split_by_junctions,
        split_by_peaks, split_by_valleys, SplitSegment, ValleySplitFilter,
    };
    
    let mut split_loci = Vec::new();
    
    println!("    Coverage splitting for {} loci...", loci.len());
    
    for (locus_idx, locus) in loci.iter().enumerate() {
        let locus_span = locus.end - locus.start;
        println!("    Locus {}: {}:{}-{} (span: {} bp)", 
                 locus_idx + 1, locus.chrom, locus.start, locus.end, locus_span);
        
        // Build coverage profile for this locus
        let profile = match build_coverage_profile(
            bam_path,
            &locus.chrom,
            locus.start,
            locus.end,
            seed_reads,
            split_cfg.profile_bin_size,
        ) {
            Ok(p) => p,
            Err(e) => {
                println!("      ERROR building profile: {}", e);
                split_loci.push(locus.clone());
                continue;
            }
        };
        
        // Compute coverage statistics for dynamic parameters
        let max_cov = profile.bins.iter().copied().fold(0.0f64, f64::max);
        let mean_cov: f64 = profile.bins.iter().sum::<f64>() / profile.bins.len().max(1) as f64;
        let variance: f64 = profile.bins.iter()
            .map(|&x| (x - mean_cov).powi(2))
            .sum::<f64>() / profile.bins.len().max(1) as f64;
        let std_cov = variance.sqrt();
        
        println!("      Coverage: max={:.0}, mean={:.1}, std={:.1}", max_cov, mean_cov, std_cov);
        
        let mut params = if split_cfg.use_dynamic_valleys {
            let mut p = ValleyDetectionParams::dynamic(max_cov, std_cov, locus_span);
            p.bin_size = split_cfg.profile_bin_size;
            println!(
                "      Dynamic params: valley_frac={:.2}, prominence={:.2}, min_gap={}, bin={}",
                p.valley_frac,
                p.min_prominence_frac,
                p.min_gap_bp,
                p.bin_size
            );
            p
        } else {
            let p = split_cfg.fixed_valley.clone();
            println!(
                "      Fixed valley params: valley_frac={:.3}, peak_frac={:.3}, prominence={:.3}, min_gap={}, min_seg={}, bin={}",
                p.valley_frac,
                p.peak_threshold_frac,
                p.min_prominence_frac,
                p.min_gap_bp,
                p.min_segment_bp,
                p.bin_size
            );
            p
        };
        params.bin_size = split_cfg.profile_bin_size;
        
        // Analyze valleys
        let (peaks, valleys) = analyze_valleys(&profile, &params);
        
        println!("      Found {} peaks, {} valleys", peaks.len(), valleys.len());
        
        // Show all valleys with their local context
        let sig_valleys: Vec<_> = valleys.iter().filter(|v| v.is_significant).collect();
        println!("      Significant valleys: {} / {}", sig_valleys.len(), valleys.len());
        
        for (vi, v) in valleys.iter().enumerate() {
            let status = if v.is_significant { "✓ SPLIT" } else { "✗ skip" };
            println!("        Valley {}: pos={}, cov={:.1}, depth_ratio={:.2} {}", 
                     vi + 1, v.position, v.coverage, v.depth_ratio, status);
        }
        
        // Split at significant valleys
        let segments = split_by_valleys(
            &locus.chrom,
            locus.start,
            locus.end,
            &profile,
            &valleys,
            &params,
            ValleySplitFilter::SignificantOnly,
        );
        
        println!("      Split into {} segments", segments.len());
        
        // If valley-based splitting didn't work well (e.g., continuous coverage between genes),
        // try peak-based splitting as a fallback
        let final_segments = if segments.len() == 1 && peaks.len() > 1 {
            println!("      Valley splitting produced 1 segment, trying peak-based fallback...");
            let peak_segments = split_by_peaks(
                &locus.chrom,
                locus.start,
                locus.end,
                &profile,
                &peaks,
                &params,
            );
            
            if peak_segments.len() > 1 {
                println!("      Peak-based splitting: {} segments", peak_segments.len());
                peak_segments
            } else {
                segments
            }
        } else {
            segments
        };
        
        // Merge adjacent small segments that are likely intra-genic
        let merged_segments = merge_small_segments(
            &final_segments,
            split_cfg.valley_merge_min_segment_bp,
            split_cfg.valley_merge_max_gap_bp,
        );
        println!("      After merging small segments: {} -> {} segments", 
                 final_segments.len(), merged_segments.len());
        
        if merged_segments.is_empty() {
            println!("      -> No split, keeping original locus");
            split_loci.push(locus.clone());
        } else {
            for (si, segment) in merged_segments.iter().enumerate() {
                println!("        Segment {}: {}-{} ({} bp)", 
                         si + 1, segment.start, segment.end, segment.end - segment.start);
                let segment_reads: HashSet<String> = locus.reads.iter().cloned().collect();
                split_loci.push(TransitiveLocus {
                    chrom: segment.chrom.clone(),
                    start: segment.start,
                    end: segment.end,
                    reads: segment_reads,
                    component: locus.component,
                    jaccard_with_seed: locus.jaccard_with_seed,
                    distance_from_seed: locus.distance_from_seed,
                    peak_count: 1,
            confidence: 1.0,
        });
            }
        }
    }
    
    // Third pass: Split very large loci (>200kb) using improved junction splitting
    // Only for loci that are still too big after coverage/peak splitting
    // Uses splice pattern similarity to find gene boundaries
    
    // Separate into small and large loci
    let (small_loci, large_loci): (Vec<_>, Vec<_>) = split_loci
        .into_iter()
        .partition(|l| l.end - l.start <= split_cfg.junction_min_locus_bp);
    
    if large_loci.is_empty() {
        // Add peak_count and confidence to small loci
        let mut result = Vec::new();
        for locus in small_loci {
            let span = locus.end - locus.start;
            let peak_count = 1;  // Approximate
            let confidence = if span > split_cfg.junction_min_locus_bp {
                0.6
            } else if span < 5_000 {
                0.4
            } else {
                1.0
            };
            result.push(TransitiveLocus {
                chrom: locus.chrom,
                start: locus.start,
                end: locus.end,
                reads: locus.reads.clone(),
                component: locus.component,
                jaccard_with_seed: locus.jaccard_with_seed,
                distance_from_seed: locus.distance_from_seed,
                peak_count,
                confidence,
            });
        }
        println!("  Final: {} loci -> {} segments after all processing", 
                 loci.len(), result.len());
        return Ok(result);
    }
    
    println!(
        "  Checking {} loci (span > {} bp) for junction-based splitting...",
        large_loci.len(),
        split_cfg.junction_min_locus_bp
    );
    
    // Convert large loci to segments
    let large_segments: Vec<SplitSegment> = large_loci.iter().map(|l| SplitSegment {
        chrom: l.chrom.clone(),
        start: l.start,
        end: l.end,
        peak_count: 1,
        total_coverage: 0.0,
    }).collect();
    
    let junction_segments = match split_by_junctions(
        bam_path,
        &large_segments,
        seed_reads,
        split_cfg.junction_min_split_span_bp,
    ) {
        Ok(segs) => segs,
        Err(e) => {
            eprintln!("    Warning: Junction splitting failed: {}", e);
            large_segments
        }
    };
    
    let merged_segments = merge_small_segments(
        &junction_segments,
        split_cfg.junction_merge_min_segment_bp,
        split_cfg.junction_merge_max_gap_bp,
    );
    
    // Build final result
    let mut final_loci = small_loci;
    
    // Convert merged segments back to TransitiveLocus.
    // Junction merge can shift boundaries so strict containment (parent fully contains seg) fails.
    // Previously `find` + `unwrap_or` emitted segments with empty reads and zero jaccard.
    for seg in merged_segments {
        let best_parent = large_loci
            .iter()
            .filter(|l| l.chrom == seg.chrom)
            .max_by_key(|l| bp_overlap_same_chr(l.start, l.end, seg.start, seg.end));

        let Some(l) = best_parent else {
            println!(
                "      Warning: junction segment {}:{}-{} has no parent on same chrom, skipping",
                seg.chrom, seg.start, seg.end
            );
            continue;
        };

        let ovl = bp_overlap_same_chr(l.start, l.end, seg.start, seg.end);
        if ovl <= 0 {
            println!(
                "      Warning: junction segment {}:{}-{} does not overlap parent {}:{}-{}, skipping",
                seg.chrom, seg.start, seg.end, l.chrom, l.start, l.end
            );
            continue;
        }

        final_loci.push(TransitiveLocus {
            chrom: seg.chrom,
            start: seg.start,
            end: seg.end,
            reads: l.reads.clone(),
            component: l.component,
            jaccard_with_seed: l.jaccard_with_seed,
            distance_from_seed: l.distance_from_seed,
            peak_count: 1,
            confidence: 1.0,
        });
    }
    
    // Ensure all seed regions are preserved in the output
    // Some seed regions may have been lost during splitting if they fell in valleys
    let mut final_loci = ensure_seed_regions_preserved(final_loci, seed_regions, bam_path, seed_reads)?;
    
    println!("  Final: {} loci -> {} segments (junction splitting applied to {} large loci)", 
             loci.len(), final_loci.len(), large_loci.len());
    Ok(final_loci)
}

/// Ensure all seed regions are represented in the output loci.
/// If a seed region doesn't overlap any output locus, add it back.
fn ensure_seed_regions_preserved(
    loci: Vec<TransitiveLocus>,
    seed_regions: &[(String, i64, i64)],
    bam_path: &str,
    seed_reads: &HashSet<String>,
) -> Result<Vec<TransitiveLocus>> {
    use noodles::bam;
    use noodles::core::Position;
    
    let mut result = loci;
    let mut added = 0;
    
    println!("  Checking {} seed regions for preservation...", seed_regions.len());
    
    for (chrom, start, end) in seed_regions {
        // Check if this seed region overlaps any existing locus
        let overlaps = result.iter().any(|l| {
            l.chrom == *chrom && !(l.end <= *start || l.start >= *end)
        });
        
        if !overlaps {
            // Seed region was lost during splitting - add it back
            println!("  Preserving seed region: {}:{}-{}", chrom, start, end);
            
            // Collect reads for this seed region
            let mut reader = bam::io::indexed_reader::Builder::default()
                .build_from_path(bam_path)?;
            let header = reader.read_header()?;
            
            let region_start = Position::try_from((*start).max(1) as usize)?;
            let region_end = Position::try_from((*end).max(1) as usize)?;
            let region = noodles::core::Region::new(chrom.clone(), region_start..=region_end);
            
            let mut reads: HashSet<String> = HashSet::default();
            for rec_result in reader.query(&header, &region)? {
                let record = rec_result?;
                if record.flags().is_unmapped() {
                    continue;
                }
                if let Some(name) = record.name() {
                    if let Ok(name_str) = std::str::from_utf8(name) {
                        reads.insert(name_str.to_string());
                    }
                }
            }
            
            result.push(TransitiveLocus {
                chrom: chrom.clone(),
                start: *start,
                end: *end,
                reads,
                component: 0,  // Seed region
                jaccard_with_seed: 1.0,  // Perfect match to itself
                distance_from_seed: 0,   // This IS the seed
                peak_count: 1,
                confidence: 1.0,
            });
            added += 1;
        }
    }
    
    if added > 0 {
        println!("  Added {} seed regions that were lost during splitting", added);
    }
    
    Ok(result)
}

/// Parallel version of split_loci_by_coverage for improved performance
/// Uses per-thread BAM readers for concurrent processing
/// Seed regions are always preserved in the output.
pub fn split_loci_by_coverage_par(
    bam_path: &str,
    loci: &[TransitiveLocus],
    seed_reads: &HashSet<String>,
    split_cfg: &TransitiveCoverageSplitConfig,
    seed_regions: &[(String, i64, i64)],
) -> Result<Vec<TransitiveLocus>> {
    use crate::coverage_valley::{
        analyze_valleys, build_coverage_profile, merge_small_segments, split_by_junctions,
        split_by_peaks, split_by_valleys, SplitSegment, ValleySplitFilter,
    };
    
    if loci.is_empty() {
        return Ok(vec![]);
    }
    
    println!("    Coverage splitting for {} loci (PARALLEL)...", loci.len());
    
    // Process loci in parallel using Rayon
    let results: Vec<_> = loci
        .par_iter()
        .enumerate()
        .map(|(locus_idx, locus)| {
            // Each thread gets its own BAM reader
            let locus_span = locus.end - locus.start;
            
            // Build coverage profile for this locus
            let profile = match build_coverage_profile(
                bam_path,
                &locus.chrom,
                locus.start,
                locus.end,
                seed_reads,
                split_cfg.profile_bin_size,
            ) {
                Ok(p) => p,
                Err(e) => {
                    eprintln!("      ERROR building profile for locus {}: {}", locus_idx + 1, e);
                    return (locus_idx, locus.clone(), vec![], false);
                }
            };
            
            // Compute coverage statistics for dynamic parameters
            let max_cov = profile.bins.iter().copied().fold(0.0f64, f64::max);
            let mean_cov: f64 = profile.bins.iter().sum::<f64>() / profile.bins.len().max(1) as f64;
            let variance: f64 = profile.bins.iter()
                .map(|&x| (x - mean_cov).powi(2))
                .sum::<f64>() / profile.bins.len().max(1) as f64;
            let std_cov = variance.sqrt();
            
            let params = if split_cfg.use_dynamic_valleys {
                let mut p = ValleyDetectionParams::dynamic(max_cov, std_cov, locus_span);
                p.bin_size = split_cfg.profile_bin_size;
                p
            } else {
                let mut p = split_cfg.fixed_valley.clone();
                p.bin_size = split_cfg.profile_bin_size;
                p
            };
            
            // Analyze valleys
            let (peaks, valleys) = analyze_valleys(&profile, &params);
            
            // Split at significant valleys
            let segments = split_by_valleys(
                &locus.chrom,
                locus.start,
                locus.end,
                &profile,
                &valleys,
                &params,
                ValleySplitFilter::SignificantOnly,
            );
            
            // If valley-based splitting didn't work well, try peak-based fallback
            let final_segments = if segments.len() == 1 && peaks.len() > 1 {
                let peak_segments = split_by_peaks(
                    &locus.chrom,
                    locus.start,
                    locus.end,
                    &profile,
                    &peaks,
                    &params,
                );
                
                if peak_segments.len() > 1 {
                    peak_segments
                } else {
                    segments
                }
            } else {
                segments
            };
            
            // Merge adjacent small segments
            let merged_segments = merge_small_segments(
                &final_segments,
                split_cfg.valley_merge_min_segment_bp,
                split_cfg.valley_merge_max_gap_bp,
            );
            let was_split = merged_segments.len() > 1;
            
            (locus_idx, locus.clone(), merged_segments, was_split)
        })
        .collect();
    
    // Collect results and print summary
    let mut split_loci = Vec::new();
    let mut split_count = 0;
    
    for (locus_idx, locus, segments, was_split) in results {
        if was_split {
            split_count += 1;
            println!("    Locus {}: {}:{}-{} -> {} segments",
                     locus_idx + 1, locus.chrom, locus.start, locus.end, segments.len());
        }
        
        if segments.is_empty() {
            split_loci.push(locus);
        } else {
            for segment in segments {
                let segment_reads: HashSet<String> = locus.reads.iter().cloned().collect();
                split_loci.push(TransitiveLocus {
                    chrom: segment.chrom,
                    start: segment.start,
                    end: segment.end,
                    reads: segment_reads,
                    component: locus.component,
                    jaccard_with_seed: locus.jaccard_with_seed,
                    distance_from_seed: locus.distance_from_seed,
                    peak_count: 1,
                    confidence: 1.0,
                });
            }
        }
    }
    
    println!("    Parallel split: {} loci -> {} segments ({} loci split)",
             loci.len(), split_loci.len(), split_count);
    
    // Junction-based splitting for large loci (also parallel)
    let (small_loci, large_loci): (Vec<_>, Vec<_>) = split_loci
        .into_iter()
        .partition(|l| l.end - l.start <= split_cfg.junction_min_locus_bp);
    
    if large_loci.is_empty() {
        // Ensure seed regions are preserved even if no junction splitting needed
        return ensure_seed_regions_preserved(small_loci, seed_regions, bam_path, seed_reads);
    }
    
    println!(
        "  Checking {} large loci (>{} bp) for junction-based splitting...",
        large_loci.len(),
        split_cfg.junction_min_locus_bp
    );
    
    // Convert large loci to segments and process in parallel
    let large_segments: Vec<SplitSegment> = large_loci.iter().map(|l| SplitSegment {
        chrom: l.chrom.clone(),
        start: l.start,
        end: l.end,
        peak_count: 1,
        total_coverage: 0.0,
    }).collect();
    
    let junction_segments = match split_by_junctions(
        bam_path,
        &large_segments,
        seed_reads,
        split_cfg.junction_min_split_span_bp,
    ) {
        Ok(segs) => segs,
        Err(e) => {
            eprintln!("    Warning: Junction splitting failed: {}", e);
            large_segments
        }
    };
    
    let merged_segments = merge_small_segments(
        &junction_segments,
        split_cfg.junction_merge_min_segment_bp,
        split_cfg.junction_merge_max_gap_bp,
    );
    
    // Build final result
    let mut final_loci = small_loci;
    
    for seg in merged_segments {
        let best_parent = large_loci
            .iter()
            .filter(|l| l.chrom == seg.chrom)
            .max_by_key(|l| bp_overlap_same_chr(l.start, l.end, seg.start, seg.end));

        if let Some(l) = best_parent {
            let ovl = bp_overlap_same_chr(l.start, l.end, seg.start, seg.end);
            if ovl > 0 {
                final_loci.push(TransitiveLocus {
                    chrom: seg.chrom,
                    start: seg.start,
                    end: seg.end,
                    reads: l.reads.clone(),
                    component: l.component,
                    jaccard_with_seed: l.jaccard_with_seed,
                    distance_from_seed: l.distance_from_seed,
                    peak_count: 1,
                    confidence: 1.0,
                });
            }
        }
    }
    
    // Ensure all seed regions are preserved in the output
    println!("  About to call ensure_seed_regions_preserved with {} loci and {} seed regions", 
             final_loci.len(), seed_regions.len());
    let mut final_loci = ensure_seed_regions_preserved(final_loci, seed_regions, bam_path, seed_reads)?;
    
    println!("  Final: {} loci -> {} segments", loci.len(), final_loci.len());
    Ok(final_loci)
}

/// Inner quantile on sorted values (p in [0, 1]).
fn quantile_i64_sorted(sorted: &[i64], p: f64) -> i64 {
    if sorted.is_empty() {
        return 0;
    }
    let p = p.clamp(0.0, 1.0);
    let idx = ((sorted.len() - 1) as f64 * p).floor() as usize;
    sorted[idx.min(sorted.len() - 1)]
}

/// Tighten each transitive locus using percentiles of **primary** alignment starts and ends
/// for reads in the locus that overlap the locus. No GFF or Iso-Seq required.
///
/// Queries only `[locus.start - pad, locus.end + pad]` on the locus chromosome for speed.
pub fn trim_loci_read_align_envelope(
    bam_path: &str,
    loci: Vec<TransitiveLocus>,
    pct_lo: f64,
    pct_hi: f64,
    min_span_bp: i64,
    min_primary_alignments: usize,
    query_pad_bp: i64,
) -> Result<Vec<TransitiveLocus>> {
    let p_lo = (pct_lo / 100.0).clamp(0.005, 0.495);
    let p_hi = (pct_hi / 100.0).clamp(0.505, 0.995);
    if p_lo >= p_hi {
        anyhow::bail!("read envelope: require pct_lo < pct_hi (percentiles on [0, 100])");
    }

    let mut reader = bam::io::indexed_reader::Builder::default().build_from_path(bam_path)?;
    let header = reader.read_header()?;

    let mut out = Vec::with_capacity(loci.len());

    for locus in loci {
        let ref_seqs = header.reference_sequences();
        let chrom_name: bstr::BString = locus.chrom.as_bytes().into();
        let Some(tid) = ref_seqs.get_index_of(&chrom_name) else {
            out.push(locus);
            continue;
        };
        let ref_len = ref_seqs[tid].length().get() as i64;

        let pad = query_pad_bp.max(0);
        let q0 = (locus.start - pad).max(1);
        let q1 = (locus.end + pad).min(ref_len);
        if q1 <= q0 {
            out.push(locus);
            continue;
        }

        let start = Position::new(q0 as usize).unwrap_or(Position::MIN);
        let end = Position::new(q1 as usize).unwrap_or(Position::MAX);
        let region = noodles::core::Region::new(locus.chrom.clone(), start..=end);

        let read_set: HashSet<&str> = locus.reads.iter().map(|s| s.as_str()).collect();
        let mut starts: Vec<i64> = Vec::new();
        let mut ends: Vec<i64> = Vec::new();

        let query_iter = match reader.query(&header, &region) {
            Ok(it) => it,
            Err(e) => {
                eprintln!(
                    "      Warning: read-envelope query failed {}:{}-{}: {}",
                    locus.chrom, q0, q1, e
                );
                out.push(locus);
                continue;
            }
        };

        for result in query_iter {
            let record = match result {
                Ok(r) => r,
                Err(_) => continue,
            };
            if record.flags().is_unmapped() {
                continue;
            }
            if record.flags().is_secondary() || record.flags().is_supplementary() {
                continue;
            }
            let name_bytes = match record.name() {
                Some(n) => n,
                None => continue,
            };
            let read_name = match std::str::from_utf8(name_bytes) {
                Ok(s) => s,
                Err(_) => continue,
            };
            if !read_set.contains(read_name) {
                continue;
            }
            let aln_start = record
                .alignment_start()
                .transpose()
                .ok()
                .flatten()
                .map(|p: Position| p.get() as i64)
                .unwrap_or(0);
            let aln_end = record
                .alignment_end()
                .transpose()
                .ok()
                .flatten()
                .map(|p: Position| p.get() as i64)
                .unwrap_or(aln_start + 1);

            if aln_end <= locus.start || aln_start >= locus.end {
                continue;
            }
            starts.push(aln_start);
            ends.push(aln_end);
        }

        if starts.len() < min_primary_alignments {
            println!(
                "      Read-envelope: keep {}:{}-{} ({} primaries overlap; need {})",
                locus.chrom,
                locus.start,
                locus.end,
                starts.len(),
                min_primary_alignments
            );
            out.push(locus);
            continue;
        }

        starts.sort_unstable();
        ends.sort_unstable();
        let q_start = quantile_i64_sorted(&starts, p_lo);
        let q_end = quantile_i64_sorted(&ends, p_hi);
        if q_end <= q_start {
            println!(
                "      Read-envelope: keep {}:{}-{} (degenerate quantiles)",
                locus.chrom, locus.start, locus.end
            );
            out.push(locus);
            continue;
        }

        // Shrink-only: stay inside the pre-trim locus (never expand coordinates).
        let mut ns = q_start.max(locus.start);
        let mut ne = q_end.min(locus.end);
        if ne <= ns {
            println!(
                "      Read-envelope: keep {}:{}-{} (quantile box outside locus)",
                locus.chrom, locus.start, locus.end
            );
            out.push(locus);
            continue;
        }

        if ne - ns < min_span_bp {
            println!(
                "      Read-envelope: keep {}:{}-{} (clamped span {} < min {})",
                locus.chrom,
                locus.start,
                locus.end,
                ne - ns,
                min_span_bp
            );
            out.push(locus);
            continue;
        }

        ns = ns.max(1);
        ne = ne.min(ref_len);
        if ne <= ns {
            out.push(locus);
            continue;
        }

        let orig_span = locus.end - locus.start;
        println!(
            "      Read-envelope: {}:{}-{} ({} bp) -> {}:{}-{} ({} bp, {} primaries)",
            locus.chrom,
            locus.start,
            locus.end,
            orig_span,
            locus.chrom,
            ns,
            ne,
            ne - ns,
            starts.len()
        );

        out.push(TransitiveLocus {
            chrom: locus.chrom,
            start: ns,
            end: ne,
            reads: locus.reads,
            component: locus.component,
            jaccard_with_seed: locus.jaccard_with_seed,
            distance_from_seed: locus.distance_from_seed,
            peak_count: locus.peak_count,
            confidence: locus.confidence,
        });
    }

    Ok(out)
}

/// Simple wrapper that returns only loci (for backward compatibility)
pub fn detect_transitive_multi_seed(
    bam_path: &str,
    seeds: &[(String, i64, i64)],
    target_chroms: &[String],
    params: &TransitiveParams,
) -> Result<(Vec<TransitiveLocus>, Vec<ComponentSeedCandidate>)> {
    let (loci, candidates, _reads) = detect_transitive_multi_seed_and_reads(bam_path, seeds, target_chroms, params)?;
    Ok((loci, candidates))
}
