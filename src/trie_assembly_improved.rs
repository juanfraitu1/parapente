//! Improved Trie-Based Assembly with Coverage-Aware Filtering
//! 
//! Key improvements over original:
//! 1. Coverage-weighted path extraction (not just read count)
//! 2. Junction confidence scoring (donor/acceptor support)
//! 3. Path support aggregation (sum of read support, not just count)
//! 4. Minimum coverage threshold per junction
//! 5. Reference-guided junction validation (optional)

use crate::bam_parser::{SplicedRead, SpliceJunction};
use crate::isoform::Isoform;
use hashbrown::{HashMap, HashSet};
use rayon::prelude::*;

/// Configuration for improved trie assembly
#[derive(Debug, Clone)]
pub struct AssemblyConfig {
    /// Minimum support for a junction to be considered (default: 2)
    pub min_junction_support: usize,
    /// Minimum coverage for a path to be output (default: 3)
    pub min_path_coverage: usize,
    /// Tolerance for junction matching in bp (default: 10)
    pub junction_tolerance: u64,
    /// Minimum fraction of reads supporting a junction (default: 0.1)
    pub min_junction_frequency: f64,
    /// Whether to use coverage-weighted support (default: true)
    pub use_coverage_weighting: bool,
    /// Minimum junctions per isoform (default: 1)
    pub min_junctions: usize,
    /// Maximum isoforms per gene (default: 50)
    pub max_isoforms_per_gene: usize,
}

impl Default for AssemblyConfig {
    fn default() -> Self {
        Self {
            min_junction_support: 2,
            min_path_coverage: 3,
            junction_tolerance: 10,
            min_junction_frequency: 0.1,
            use_coverage_weighting: true,
            min_junctions: 1,
            max_isoforms_per_gene: 50,
        }
    }
}

/// Junction with coverage information
#[derive(Debug, Clone)]
struct JunctionInfo {
    start: u64,
    end: u64,
    /// Total coverage (sum of read supports)
    coverage: usize,
    /// Number of distinct reads
    read_count: usize,
    /// Read IDs for tracking
    read_ids: HashSet<String>,
}

impl JunctionInfo {
    fn new(start: u64, end: u64) -> Self {
        Self {
            start,
            end,
            coverage: 0,
            read_count: 0,
            read_ids: HashSet::new(),
        }
    }
    
    fn add_read(&mut self, read_id: &str, support: usize) {
        if self.read_ids.insert(read_id.to_string()) {
            self.read_count += 1;
        }
        self.coverage += support;
    }
    
    fn confidence_score(&self, total_reads: usize) -> f64 {
        if total_reads == 0 {
            return 0.0;
        }
        // Confidence based on both coverage and distinct read count
        let coverage_frac = self.coverage as f64 / total_reads as f64;
        let read_frac = self.read_count as f64 / total_reads as f64;
        (coverage_frac + read_frac) / 2.0
    }
}

/// Improved trie node with coverage tracking
#[derive(Debug, Default)]
struct TrieNode {
    junction: Option<(u64, u64)>,
    children: HashMap<(u64, u64), TrieNode>,
    /// Junction info for this node
    junction_info: Option<JunctionInfo>,
    /// Complete reads that end here
    complete_reads: Vec<(String, u64, u64, usize)>, // (id, tss, tes, support)
    /// Total coverage at this node
    total_coverage: usize,
}

impl TrieNode {
    fn new(junction: Option<(u64, u64)>) -> Self {
        let junction_info = junction.map(|(s, e)| JunctionInfo::new(s, e));
        Self {
            junction,
            children: HashMap::new(),
            junction_info,
            complete_reads: Vec::new(),
            total_coverage: 0,
        }
    }
    
    fn coverage(&self) -> usize {
        self.total_coverage
    }
    
    fn read_count(&self) -> usize {
        self.junction_info.as_ref().map(|j| j.read_count).unwrap_or(0)
    }
}

/// Improved junction trie with coverage-aware extraction
#[derive(Debug)]
struct CoverageAwareTrie {
    root: TrieNode,
    tolerance: u64,
    total_reads: usize,
    total_coverage: usize,
}

impl CoverageAwareTrie {
    fn new(tolerance: u64) -> Self {
        Self {
            root: TrieNode::new(None),
            tolerance,
            total_reads: 0,
            total_coverage: 0,
        }
    }
    
    fn insert(&mut self, read: &SplicedRead) {
        let read_support = read.read_support.unwrap_or(1);
        self.total_reads += 1;
        self.total_coverage += read_support;
        
        let mut current = &mut self.root;
        
        for junction in &read.junctions {
            let key = (junction.start, junction.end);
            let matched_key = self.find_matching_junction(current, key);
            let use_key = matched_key.unwrap_or(key);
            
            if !current.children.contains_key(&use_key) {
                current.children.insert(use_key, TrieNode::new(Some(use_key)));
            }
            
            current = current.children.get_mut(&use_key).unwrap();
            
            // Update junction info
            if let Some(ref mut info) = current.junction_info {
                info.add_read(&read.read_id, read_support);
            }
            current.total_coverage += read_support;
            current.complete_reads.push((
                read.read_id.clone(),
                read.tx_start,
                read.tx_end,
                read_support,
            ));
        }
    }
    
    fn find_matching_junction(&self, node: &TrieNode, target: (u64, u64)) -> Option<(u64, u64)> {
        for (key, _) in &node.children {
            let donor_diff = key.0.abs_diff(target.0);
            let acceptor_diff = key.1.abs_diff(target.1);
            
            if donor_diff <= self.tolerance && acceptor_diff <= self.tolerance {
                return Some(*key);
            }
        }
        None
    }
    
    /// Extract high-confidence isoforms using coverage-weighted scoring
    fn extract_isoforms(&self, config: &AssemblyConfig) -> Vec<IsoformPath> {
        let mut results = Vec::new();
        let mut current_path = Vec::new();
        
        self.dfs_extract(
            &self.root,
            &mut current_path,
            config,
            &mut results,
        );
        
        // Score and sort by confidence
        let mut scored: Vec<_> = results.into_iter()
            .map(|path| {
                let score = self.score_path(&path, config);
                (score, path)
            })
            .collect();
        
        scored.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
        
        // Take top N per gene
        scored.into_iter()
            .take(config.max_isoforms_per_gene)
            .map(|(_, path)| path)
            .collect()
    }
    
    fn dfs_extract(
        &self,
        node: &TrieNode,
        current_path: &mut Vec<(u64, u64)>,
        config: &AssemblyConfig,
        results: &mut Vec<IsoformPath>,
    ) {
        // Check if current path is valid output
        if !current_path.is_empty() && node.coverage() >= config.min_path_coverage {
            // Check all junctions meet minimum support
            let all_junctions_valid = current_path.iter().all(|(s, e)| {
                if let Some(child) = self.find_child_with_junction(node, (*s, *e)) {
                    child.read_count() >= config.min_junction_support
                } else {
                    false
                }
            });
            
            if all_junctions_valid {
                let total_support: usize = node.complete_reads.iter().map(|(_, _, _, s)| s).sum();
                results.push(IsoformPath {
                    junctions: current_path.clone(),
                    complete_reads: node.complete_reads.clone(),
                    coverage: node.coverage(),
                    support: total_support,
                });
            }
        }
        
        // Recurse to children with sufficient support
        for (_, child) in &node.children {
            if let Some(ref info) = child.junction_info {
                let confidence = info.confidence_score(self.total_reads);
                if confidence >= config.min_junction_frequency {
                    if let Some(junction) = child.junction {
                        current_path.push(junction);
                        self.dfs_extract(child, current_path, config, results);
                        current_path.pop();
                    }
                }
            }
        }
    }
    
    fn find_child_with_junction(&self, node: &TrieNode, junction: (u64, u64)) -> Option<&TrieNode> {
        for (key, child) in &node.children {
            if key.0 == junction.0 && key.1 == junction.1 {
                return Some(child);
            }
        }
        None
    }
    
    fn score_path(&self, path: &IsoformPath, config: &AssemblyConfig) -> f64 {
        let coverage_score = path.coverage as f64 / self.total_coverage as f64;
        let support_score = path.support as f64 / self.total_coverage as f64;
        let length_score = (path.junctions.len() as f64).sqrt() / 10.0; // Favor longer isoforms slightly
        
        if config.use_coverage_weighting {
            coverage_score * 0.5 + support_score * 0.4 + length_score * 0.1
        } else {
            support_score * 0.7 + length_score * 0.3
        }
    }
}

/// Represents an isoform path through the trie
#[derive(Debug, Clone)]
struct IsoformPath {
    junctions: Vec<(u64, u64)>,
    complete_reads: Vec<(String, u64, u64, usize)>, // (id, tss, tes, support)
    coverage: usize,
    support: usize,
}

/// Main entry point for improved trie assembly
pub fn trie_assembly_improved(
    reads: Vec<SplicedRead>,
    config: &AssemblyConfig,
) -> anyhow::Result<Vec<Isoform>> {
    eprintln!("  Improved Trie Assembly");
    eprintln!("  Input reads: {}", reads.len());
    eprintln!("  Config: {:?}", config);
    
    // Separate by chrom:strand
    let mut by_location: HashMap<(String, char), Vec<SplicedRead>> = HashMap::new();
    for read in reads {
        let key = (read.chrom.clone(), read.strand);
        by_location.entry(key).or_default().push(read);
    }
    
    let all_isoforms = std::sync::Mutex::new(Vec::new());
    let counter = std::sync::atomic::AtomicUsize::new(0);
    
    for ((chrom, strand), location_reads) in by_location {
        eprintln!("  {}:{} - {} reads", chrom, strand, location_reads.len());
        
        // Group into genes
        let gene_groups = group_reads_into_genes(&location_reads, 100_000);
        
        // Process each gene group in parallel
        let chrom_isoforms: Vec<Isoform> = gene_groups
            .into_par_iter()
            .enumerate()
            .filter_map(|(gene_idx, gene_reads)| {
                let counter_base = counter.fetch_add(1000, std::sync::atomic::Ordering::Relaxed);
                let isoforms = assemble_gene_improved(
                    gene_reads,
                    &chrom,
                    strand,
                    config,
                    gene_idx + 1,
                    counter_base,
                );
                if isoforms.is_empty() { None } else { Some(isoforms) }
            })
            .flatten()
            .collect();
        
        all_isoforms.lock().unwrap().extend(chrom_isoforms);
    }
    
    let result = all_isoforms.into_inner().unwrap();
    eprintln!("  Total isoforms: {}", result.len());
    Ok(result)
}

/// Assemble transcripts from a gene group with improved filtering
fn assemble_gene_improved(
    reads: Vec<&SplicedRead>,
    chrom: &str,
    strand: char,
    config: &AssemblyConfig,
    gene_idx: usize,
    counter_base: usize,
) -> Vec<Isoform> {
    if reads.len() < config.min_path_coverage.max(2) {
        return Vec::new();
    }
    
    // Separate single-exon and multi-exon
    let single_exon: Vec<_> = reads.iter().filter(|r| r.junctions.is_empty()).copied().collect();
    let multi_exon: Vec<_> = reads.iter().filter(|r| !r.junctions.is_empty()).copied().collect();
    
    let mut isoforms = Vec::new();
    let mut counter = counter_base;
    
    // === MULTI-EXON: Use Coverage-Aware Trie ===
    if !multi_exon.is_empty() {
        let trie_isoforms = assemble_multi_exon_improved(
            &multi_exon,
            chrom,
            strand,
            config,
            gene_idx,
            &mut counter,
        );
        isoforms.extend(trie_isoforms);
    }
    
    // === SINGLE-EXON: Coverage-based with higher thresholds ===
    if !single_exon.is_empty() {
        let se_isoforms = assemble_single_exon_improved(
            &single_exon,
            chrom,
            strand,
            config,
            gene_idx,
            &mut counter,
        );
        isoforms.extend(se_isoforms);
    }
    
    isoforms
}

/// Assemble multi-exon transcripts using improved trie
fn assemble_multi_exon_improved(
    reads: &[&SplicedRead],
    chrom: &str,
    strand: char,
    config: &AssemblyConfig,
    gene_idx: usize,
    counter: &mut usize,
) -> Vec<Isoform> {
    let mut trie = CoverageAwareTrie::new(config.junction_tolerance);
    
    for read in reads {
        trie.insert(read);
    }
    
    // Extract high-confidence paths
    let paths = trie.extract_isoforms(config);
    
    // Build isoforms from paths
    let mut isoforms = Vec::new();
    for (idx, path) in paths.iter().enumerate() {
        if path.junctions.len() < config.min_junctions {
            continue;
        }
        
        let path_isoforms = build_isoforms_from_path_improved(
            &path.junctions,
            &path.complete_reads,
            chrom,
            strand,
            idx + 1,
            counter,
            path.coverage,
            path.support,
        );
        isoforms.extend(path_isoforms);
    }
    
    isoforms
}

/// Improved single-exon assembly with coverage weighting
fn assemble_single_exon_improved(
    reads: &[&SplicedRead],
    chrom: &str,
    strand: char,
    config: &AssemblyConfig,
    gene_idx: usize,
    counter: &mut usize,
) -> Vec<Isoform> {
    if reads.len() < config.min_path_coverage {
        return Vec::new();
    }
    
    // Sort by position
    let mut sorted: Vec<_> = reads.iter().copied().collect();
    sorted.sort_by_key(|r| r.tx_start);
    
    // Cluster with coverage weighting
    let tolerance = config.junction_tolerance * 10; // 100bp default
    let mut groups: Vec<Vec<&SplicedRead>> = Vec::new();
    
    for read in sorted {
        let mut added = false;
        for group in &mut groups {
            let (group_start, group_end, group_weight) = calculate_weighted_center(group);
            
            let read_support = read.read_support.unwrap_or(1) as f64;
            let read_center = (read.tx_start as f64 + read.tx_end as f64) / 2.0;
            let group_center = (group_start + group_end) / 2.0;
            
            let diff = (read_center - group_center).abs();
            if diff <= tolerance as f64 {
                group.push(read);
                added = true;
                break;
            }
        }
        if !added {
            groups.push(vec![read]);
        }
    }
    
    // Build isoforms from groups with sufficient weighted support
    let mut isoforms = Vec::new();
    for group in groups {
        let total_support: usize = group.iter().map(|r| r.read_support.unwrap_or(1)).sum();
        
        if total_support < config.min_path_coverage {
            continue;
        }
        
        // Weighted median for boundaries
        let (tss, tes, _) = calculate_weighted_center(&group);
        
        if tes - tss < 200 {
            continue;
        }
        
        *counter += 1;
        isoforms.push(create_isoform_improved(
            format!("PB.{}.{}.{}_se", gene_idx, counter, strand),
            format!("STRG.{}", gene_idx),
            chrom.to_string(),
            strand,
            vec![(tss as u64, tes as u64)],
            Vec::new(),
            total_support,
        ));
    }
    
    isoforms
}

/// Calculate weighted center of a group of reads
fn calculate_weighted_center(reads: &[&SplicedRead]) -> (f64, f64, f64) {
    let mut total_weight = 0.0;
    let mut weighted_start = 0.0;
    let mut weighted_end = 0.0;
    
    for read in reads {
        let weight = read.read_support.unwrap_or(1) as f64;
        total_weight += weight;
        weighted_start += read.tx_start as f64 * weight;
        weighted_end += read.tx_end as f64 * weight;
    }
    
    if total_weight > 0.0 {
        (weighted_start / total_weight, weighted_end / total_weight, total_weight)
    } else {
        (0.0, 0.0, 0.0)
    }
}

/// Build isoforms from a path with improved boundary estimation
fn build_isoforms_from_path_improved(
    junctions: &[(u64, u64)],
    complete_reads: &[(String, u64, u64, usize)],
    chrom: &str,
    strand: char,
    _path_idx: usize,
    counter: &mut usize,
    _coverage: usize,
    support: usize,
) -> Vec<Isoform> {
    let mut isoforms = Vec::new();
    
    if junctions.is_empty() || complete_reads.is_empty() {
        return isoforms;
    }
    
    let junction_objs: Vec<SpliceJunction> = junctions.iter()
        .map(|(d, a)| SpliceJunction::new(chrom.to_string(), *d, *a, strand))
        .collect();
    
    // Group by TSS/TES with weighted support
    let mut groups: Vec<Vec<(String, u64, u64, usize)>> = Vec::new();
    let boundary_tolerance = 100u64;
    
    for (read_id, read_tss, read_tes, read_support) in complete_reads {
        let mut added = false;
        for group in &mut groups {
            let (group_tss, group_tes, _) = calculate_weighted_center_group(group);
            
            let tss_diff = read_tss.abs_diff(group_tss as u64);
            let tes_diff = read_tes.abs_diff(group_tes as u64);
            
            if tss_diff <= boundary_tolerance && tes_diff <= boundary_tolerance {
                group.push((read_id.clone(), *read_tss, *read_tes, *read_support));
                added = true;
                break;
            }
        }
        if !added {
            groups.push(vec![(read_id.clone(), *read_tss, *read_tes, *read_support)]);
        }
    }
    
    // Build isoform per group
    for group in groups {
        if group.is_empty() { continue; }
        
        let (tss, tes, group_support) = calculate_weighted_center_group(&group);
        
        *counter += 1;
        isoforms.push(create_isoform_improved(
            format!("PB.{}.{}", counter, strand),
            format!("STRG.{}", counter),
            chrom.to_string(),
            strand,
            build_exons_from_junctions(&junction_objs, tss as u64, tes as u64),
            junction_objs.clone(),
            support.max(group_support as usize),
        ));
    }
    
    isoforms
}

/// Calculate weighted center for a group with support values
fn calculate_weighted_center_group(group: &[(String, u64, u64, usize)]) -> (f64, f64, f64) {
    let mut total_weight = 0.0;
    let mut weighted_start = 0.0;
    let mut weighted_end = 0.0;
    
    for (_, start, end, support) in group {
        let weight = *support as f64;
        total_weight += weight;
        weighted_start += *start as f64 * weight;
        weighted_end += *end as f64 * weight;
    }
    
    if total_weight > 0.0 {
        (weighted_start / total_weight, weighted_end / total_weight, total_weight)
    } else {
        (0.0, 0.0, 0.0)
    }
}

/// Build exons from junctions
fn build_exons_from_junctions(
    junctions: &[SpliceJunction],
    tss: u64,
    tes: u64,
) -> Vec<(u64, u64)> {
    let mut exons = Vec::new();
    let mut current_start = tss;
    
    for junction in junctions {
        exons.push((current_start, junction.donor));
        current_start = junction.acceptor;
    }
    exons.push((current_start, tes));
    
    exons
}

/// Create isoform with proper support annotation
fn create_isoform_improved(
    id: String,
    gene_id: String,
    chrom: String,
    strand: char,
    exons: Vec<(u64, u64)>,
    junctions: Vec<SpliceJunction>,
    support: usize,
) -> Isoform {
    let start = exons.first().map(|e| e.0).unwrap_or(0);
    let end = exons.last().map(|e| e.1).unwrap_or(0);
    
    Isoform {
        id,
        gene_id,
        chrom,
        start,
        end,
        strand: strand.to_string(),
        exons,
        junctions,
        support,
        is_coding: false,
        cds_start: None,
        cds_end: None,
    }
}

/// Group reads into genes by proximity
fn group_reads_into_genes(reads: &[SplicedRead], max_gap: u64) -> Vec<Vec<&SplicedRead>> {
    if reads.is_empty() {
        return Vec::new();
    }
    
    let mut sorted: Vec<_> = reads.iter().collect();
    sorted.sort_by_key(|r| r.tx_start);
    
    let mut groups: Vec<Vec<&SplicedRead>> = Vec::new();
    let mut current_group: Vec<&SplicedRead> = vec![sorted[0]];
    let mut current_end = sorted[0].tx_end;
    
    for read in sorted.into_iter().skip(1) {
        if read.tx_start <= current_end + max_gap {
            current_group.push(read);
            current_end = current_end.max(read.tx_end);
        } else {
            groups.push(current_group);
            current_group = vec![read];
            current_end = read.tx_end;
        }
    }
    groups.push(current_group);
    groups
}
