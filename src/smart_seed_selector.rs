//! Optimized Smart Seed Selection using bitsets for Jaccard computation
//! 
//! Performance optimizations:
//! - Single BAM reader for all queries
//! - HashMap-based read tracking (read_name -> bit_index)
//! - Vec<u128> bitsets for fast Jaccard via bit operations
//! - Rayon parallelization for independent region queries

use anyhow::{Context, Result};
use hashbrown::HashMap;
use noodles::bam;
use noodles::core::Position;
use rayon::prelude::*;
use std::sync::Mutex;

/// Number of reads that fit in one u128 block (128 bits)
const BITS_PER_BLOCK: usize = 128;

/// A bitset for representing read presence in a region
#[derive(Clone)]
struct ReadBitset {
    blocks: Vec<u128>,
}

impl ReadBitset {
    fn new(n_reads: usize) -> Self {
        let n_blocks = (n_reads + BITS_PER_BLOCK - 1) / BITS_PER_BLOCK;
        Self {
            blocks: vec![0u128; n_blocks],
        }
    }
    
    fn set(&mut self, idx: usize) {
        self.blocks[idx / BITS_PER_BLOCK] |= 1u128 << (idx % BITS_PER_BLOCK);
    }
    
    fn contains(&self, idx: usize) -> bool {
        self.blocks[idx / BITS_PER_BLOCK] & (1u128 << (idx % BITS_PER_BLOCK)) != 0
    }
    
    /// Fast union count using popcnt
    fn union_count(&self, other: &ReadBitset) -> usize {
        let mut count = 0;
        for i in 0..self.blocks.len() {
            count += (self.blocks[i] | other.blocks[i]).count_ones() as usize;
        }
        count
    }
    
    /// Fast intersection count using popcnt
    fn intersection_count(&self, other: &ReadBitset) -> usize {
        let mut count = 0;
        for i in 0..self.blocks.len() {
            count += (self.blocks[i] & other.blocks[i]).count_ones() as usize;
        }
        count
    }
}

/// A candidate gene region
#[derive(Debug, Clone)]
pub struct GeneRegion {
    pub chrom: String,
    pub start: i64,
    pub end: i64,
    pub name: Option<String>,
}

/// Parameters for building the gene-gene connectivity graph
#[derive(Debug, Clone)]
pub struct ConnectivityParams {
    /// Never use an edge threshold below this (CLI `--min-jaccard`)
    pub min_jaccard_floor: f64,
    /// If set, edge threshold is max(min_jaccard_floor, q-quantile of observed Jaccards for pairs
    /// with at least min_shared_reads). q in [0, 1], e.g. 0.15 keeps edges down to the 15th
    /// percentile of positive pairwise signal (more permissive than higher quantiles).
    pub adaptive_quantile: Option<f64>,
    /// Require at least this many shared read names to draw an edge (reduces spurious links)
    pub min_shared_reads: usize,
    /// If set, require one-sided hypergeometric p <= this value (enrichment of shared reads vs random
    /// draw null with fixed |A|,|B|,N) in addition to the Jaccard rules.
    pub hypergeom_max_pvalue: Option<f64>,
}

impl Default for ConnectivityParams {
    fn default() -> Self {
        Self {
            min_jaccard_floor: 0.001,
            adaptive_quantile: None,
            min_shared_reads: 1,
            hypergeom_max_pvalue: None,
        }
    }
}

/// Result of connectivity analysis
#[derive(Debug)]
pub struct ConnectivityAnalysis {
    pub n_genes: usize,
    pub n_components: usize,
    pub seed_indices: Vec<usize>,
    pub read_counts: Vec<usize>,
    pub components: Vec<Vec<usize>>,
    /// Jaccard cutoff used after optional adaptive adjustment
    pub effective_jaccard_threshold: f64,
    pub min_shared_reads_used: usize,
    /// Unique read names pooled across all input regions (hypergeometric universe size N)
    pub n_universe_reads: usize,
    pub hypergeom_p_max_used: Option<f64>,
}

impl ConnectivityAnalysis {
    pub fn n_seeds(&self) -> usize {
        self.seed_indices.len()
    }
    
    pub fn reduction_pct(&self) -> f64 {
        if self.n_genes == 0 {
            return 0.0;
        }
        100.0 * (1.0 - self.n_seeds() as f64 / self.n_genes as f64)
    }
}

/// Optimized: Collect reads from multiple regions using a single BAM reader
pub fn collect_reads_parallel(
    bam_path: &str,
    genes: &[GeneRegion],
) -> Result<(Vec<ReadBitset>, Vec<usize>, HashMap<String, usize>)> {
    let n = genes.len();
    if n == 0 {
        return Ok((vec![], vec![], HashMap::new()));
    }
    
    // Phase 1: Collect all reads from all regions in parallel
    println!("  Collecting reads from {} regions...", n);
    
    let results: Vec<Result<Vec<String>>> = genes.par_iter()
        .map(|gene| {
            collect_reads_from_region(bam_path, &gene.chrom, gene.start, gene.end)
        })
        .collect();
    
    // Check for errors
    let all_reads: Vec<Vec<String>> = results.into_iter()
        .map(|r| r.unwrap_or_default())
        .collect();
    
    // Phase 2: Build global read index (shared across all genes)
    let mut read_to_idx: HashMap<String, usize> = HashMap::new();
    for reads in &all_reads {
        for read in reads {
            if !read_to_idx.contains_key(read) {
                read_to_idx.insert(read.clone(), read_to_idx.len());
            }
        }
    }
    
    let n_reads = read_to_idx.len();
    println!("  Total unique reads: {}", n_reads);
    
    // Phase 3: Build bitsets (parallel)
    println!("  Building bitsets...");
    let bitsets: Vec<ReadBitset> = all_reads.par_iter()
        .map(|reads| {
            let mut bitset = ReadBitset::new(n_reads);
            for read in reads {
                if let Some(&idx) = read_to_idx.get(read) {
                    bitset.set(idx);
                }
            }
            bitset
        })
        .collect();
    
    let read_counts: Vec<usize> = bitsets.iter()
        .map(|b| b.blocks.iter().map(|&x| x.count_ones() as usize).sum())
        .collect();
    
    Ok((bitsets, read_counts, read_to_idx))
}

/// Collect read names from a single region
fn collect_reads_from_region(
    bam_path: &str,
    chrom: &str,
    start: i64,
    end: i64,
) -> Result<Vec<String>> {
    let mut reader = bam::io::indexed_reader::Builder::default()
        .build_from_path(bam_path)?;
    
    let header = reader.read_header()?;
    
    let start_pos = Position::try_from(start.max(1) as usize)?;
    let end_pos = Position::try_from(end.max(1) as usize)?;
    let region = noodles::core::Region::new(chrom, start_pos..=end_pos);
    
    let mut reads = Vec::new();
    
    for result in reader.query(&header, &region)? {
        let record = result?;
        if let Some(qname) = record.name() {
            reads.push(qname.to_string());
        }
    }
    
    Ok(reads)
}

/// log C(n,k) for 0 <= k <= n using sum of logs (stable for moderate n).
fn log_binomial(n: usize, k: usize) -> f64 {
    if k > n {
        return f64::NEG_INFINITY;
    }
    let k = k.min(n - k);
    let mut s = 0.0;
    for i in 0..k {
        s += ((n - i) as f64).ln() - ((i + 1) as f64).ln();
    }
    s
}

/// log PMF of Hypergeom(N_total, K_success_in_pop, n_draws) at k.
fn hypergeom_log_pmf(k: usize, k_pop: usize, n_draw: usize, n_total: usize) -> f64 {
    if k_pop > n_total || n_draw > n_total || k > n_draw.min(k_pop) {
        return f64::NEG_INFINITY;
    }
    if n_draw + k_pop > n_total + k {
        return f64::NEG_INFINITY;
    }
    log_binomial(k_pop, k) + log_binomial(n_total - k_pop, n_draw - k) - log_binomial(n_total, n_draw)
}

/// One-sided upper tail P(X >= k_obs) under drawing n_draw items without replacement
/// from a population of n_total with k_pop successes (scipy: hypergeom.sf(k_obs-1, M, n_good, n_draw)).
fn hypergeom_upper_tail(k_obs: usize, k_pop: usize, n_draw: usize, n_total: usize) -> f64 {
    if k_obs == 0 {
        return 1.0;
    }
    if k_pop == 0 || n_draw == 0 {
        return 0.0;
    }
    let max_i = n_draw.min(k_pop);
    if k_obs > max_i {
        return 0.0;
    }
    let n_terms = max_i.saturating_sub(k_obs).saturating_add(1);
    if n_terms > 256 {
        return hypergeom_upper_tail_normal(k_obs, k_pop, n_draw, n_total).clamp(0.0, 1.0);
    }
    let mut sum = 0.0;
    for i in k_obs..=max_i {
        let lp = hypergeom_log_pmf(i, k_pop, n_draw, n_total);
        if lp.is_finite() && lp > -745.0 {
            sum += lp.exp();
        }
    }
    sum.min(1.0)
}

/// Normal approximation with continuity correction when the exact tail has many terms.
fn hypergeom_upper_tail_normal(k_obs: usize, k_pop: usize, n_draw: usize, n_total: usize) -> f64 {
    let nf = n_total as f64;
    let mean = n_draw as f64 * k_pop as f64 / nf;
    let denom = (nf - 1.0).max(1.0);
    let var = mean * (1.0 - k_pop as f64 / nf) * (nf - n_draw as f64) / denom;
    if var <= 0.0 {
        return if (k_obs as f64) > mean { 0.0 } else { 1.0 };
    }
    let sd = var.sqrt();
    let z = ((k_obs as f64) - 0.5 - mean) / sd;
    normal_sf(z)
}

/// Standard normal survival 1 - Phi(z).
fn normal_sf(z: f64) -> f64 {
    0.5 * (1.0 - erf_approx(z / std::f64::consts::SQRT_2))
}

/// erf(x) via Abramowitz and Stegun 7.1.26 (absolute error < 1.5e-7).
fn erf_approx(x: f64) -> f64 {
    let sign = if x < 0.0 { -1.0 } else { 1.0 };
    let x = x.abs();
    let p = 0.3275911_f64;
    let a1 = 0.254829592_f64;
    let a2 = -0.284496736_f64;
    let a3 = 1.421413741_f64;
    let a4 = -1.453152027_f64;
    let a5 = 1.061405429_f64;
    let t = 1.0 / (1.0 + p * x);
    sign * (1.0
        - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * (-x * x).exp())
}

/// Passes hypergeometric enrichment filter: shared read count k vs random draw null (fixed margins).
/// Uses K = n_b (success states), n = n_a (draws), consistent with gene index order i < j.
fn hypergeom_pair_ok(k: usize, n_a: usize, n_b: usize, n_total: usize, p_max: f64) -> bool {
    if n_total == 0 || n_a == 0 || n_b == 0 || k == 0 {
        return false;
    }
    if n_a > n_total || n_b > n_total {
        return false;
    }
    let p = hypergeom_upper_tail(k, n_b, n_a, n_total);
    p.is_finite() && p <= p_max
}

/// Fast Jaccard using bitsets
fn jaccard_bitset(set1: &ReadBitset, set2: &ReadBitset) -> f64 {
    let inter = set1.intersection_count(set2);
    if inter == 0 {
        return 0.0;
    }
    let union = set1.union_count(set2);
    inter as f64 / union as f64
}

/// Analyze connectivity between genes (union-find on a Jaccard graph).
pub fn analyze_connectivity(
    bam_path: &str,
    genes: &[GeneRegion],
    params: &ConnectivityParams,
) -> Result<ConnectivityAnalysis> {
    if genes.is_empty() {
        return Ok(ConnectivityAnalysis {
            n_genes: 0,
            n_components: 0,
            seed_indices: vec![],
            read_counts: vec![],
            components: vec![],
            effective_jaccard_threshold: params.min_jaccard_floor,
            min_shared_reads_used: params.min_shared_reads,
            n_universe_reads: 0,
            hypergeom_p_max_used: params.hypergeom_max_pvalue,
        });
    }

    let n = genes.len();

    // Collect reads and build bitsets
    let (bitsets, read_counts, read_to_idx) = collect_reads_parallel(bam_path, genes)?;
    let n_universe = read_to_idx.len();

    // Compute Jaccard matrix using bitsets
    println!("  Computing Jaccard matrix ({})...", n * (n - 1) / 2);

    let pairs: Vec<(usize, usize, f64, usize)> = (0..n)
        .flat_map(|i| ((i + 1)..n).map(move |j| (i, j)))
        .par_bridge()
        .map(|(i, j)| {
            let inter = bitsets[i].intersection_count(&bitsets[j]);
            let jacc = if inter == 0 {
                0.0
            } else {
                let u = bitsets[i].union_count(&bitsets[j]);
                inter as f64 / u as f64
            };
            (i, j, jacc, inter)
        })
        .collect();

    let mut edge_threshold = params.min_jaccard_floor;
    if let Some(q) = params.adaptive_quantile {
        let q = q.clamp(0.0, 1.0);
        let mut jac_for_quantile: Vec<f64> = pairs
            .iter()
            .filter(|(_, _, _, inter)| *inter >= params.min_shared_reads)
            .map(|(_, _, j, _)| *j)
            .filter(|j| *j > 0.0)
            .collect();
        jac_for_quantile.sort_by(|a, b| a.partial_cmp(b).unwrap());
        if jac_for_quantile.len() >= 3 {
            let idx = ((jac_for_quantile.len() - 1) as f64 * q).round() as usize;
            let idx = idx.min(jac_for_quantile.len() - 1);
            edge_threshold = edge_threshold.max(jac_for_quantile[idx]);
            println!(
                "  Adaptive connectivity: quantile {:.3} over {} qualifying pairs -> Jaccard cutoff {:.6} (floor {:.6}, min_shared {})",
                q,
                jac_for_quantile.len(),
                edge_threshold,
                params.min_jaccard_floor,
                params.min_shared_reads
            );
        } else {
            println!(
                "  Adaptive connectivity: only {} qualifying pairs; using floor {:.6}",
                jac_for_quantile.len(),
                params.min_jaccard_floor
            );
        }
    }

    if let Some(pmax) = params.hypergeom_max_pvalue {
        println!(
            "  Hypergeometric edge filter: p_max={:.6} (universe N={} unique reads)",
            pmax, n_universe
        );
    }

    // Union-find for connected components
    let mut parent: Vec<usize> = (0..n).collect();
    let mut rank: Vec<usize> = vec![0; n];

    fn find(x: usize, parent: &mut [usize]) -> usize {
        if parent[x] != x {
            parent[x] = find(parent[x], parent);
        }
        parent[x]
    }

    fn union(x: usize, y: usize, parent: &mut [usize], rank: &mut [usize]) {
        let px = find(x, parent);
        let py = find(y, parent);
        if px == py {
            return;
        }
        if rank[px] < rank[py] {
            parent[px] = py;
        } else if rank[px] > rank[py] {
            parent[py] = px;
        } else {
            parent[py] = px;
            rank[px] += 1;
        }
    }

    for (i, j, jacc, inter) in &pairs {
        let na = read_counts[*i];
        let nb = read_counts[*j];
        let jacc_ok = *inter >= params.min_shared_reads && *jacc >= edge_threshold;
        let stat_ok = match params.hypergeom_max_pvalue {
            None => true,
            Some(pmax) => hypergeom_pair_ok(*inter, na, nb, n_universe, pmax),
        };
        if jacc_ok && stat_ok {
            union(*i, *j, &mut parent, &mut rank);
        }
    }

    // Group by component
    let mut component_map: HashMap<usize, Vec<usize>> = HashMap::new();
    for i in 0..n {
        let root = find(i, &mut parent);
        component_map.entry(root).or_default().push(i);
    }

    let components: Vec<Vec<usize>> = component_map.into_values().collect();
    let n_components = components.len();

    // Select best seed per component (gene with most reads)
    let seed_indices: Vec<usize> = components
        .iter()
        .map(|comp| {
            *comp
                .iter()
                .max_by_key(|&&idx| read_counts[idx])
                .unwrap()
        })
        .collect();

    println!(
        "  Found {} components, selected {} seeds ({}% reduction)",
        n_components,
        seed_indices.len(),
        100.0 * (1.0 - seed_indices.len() as f64 / n as f64)
    );

    Ok(ConnectivityAnalysis {
        n_genes: n,
        n_components,
        seed_indices,
        read_counts,
        components,
        effective_jaccard_threshold: edge_threshold,
        min_shared_reads_used: params.min_shared_reads,
        n_universe_reads: n_universe,
        hypergeom_p_max_used: params.hypergeom_max_pvalue,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_bitset() {
        let mut bs = ReadBitset::new(300);
        bs.set(0);
        bs.set(128);
        bs.set(256);
        assert_eq!(bs.blocks[0], 1u128);
        assert_eq!(bs.blocks[1], 1u128);
        assert_eq!(bs.blocks[2], 1u128);
    }
    
    #[test]
    fn test_jaccard_bitset() {
        let mut bs1 = ReadBitset::new(10);
        let mut bs2 = ReadBitset::new(10);
        
        // 50% overlap: {0,1,2,3} vs {2,3,4,5}
        for i in 0..4 { bs1.set(i); }
        for i in 2..6 { bs2.set(i); }
        
        let jacc = jaccard_bitset(&bs1, &bs2);
        assert!((jacc - 0.5).abs() < 0.01);
    }

    #[test]
    fn hypergeom_tail_bounds() {
        let p = hypergeom_upper_tail(10, 50, 30, 200);
        assert!(p >= 0.0 && p <= 1.0);
    }

    #[test]
    fn hypergeom_weak_overlap_high_p() {
        // Small overlap vs large margins: upper tail should be large (not significant)
        let p = hypergeom_upper_tail(2, 100, 100, 10000);
        assert!(p > 0.2, "p={}", p);
    }

    #[test]
    fn hypergeom_strong_overlap_low_p() {
        let p = hypergeom_upper_tail(40, 50, 50, 100);
        assert!(p < 0.05, "p={}", p);
    }
}

/// Get the selected seeds as GeneRegion objects
pub fn get_selected_seeds(analysis: &ConnectivityAnalysis, genes: &[GeneRegion]) -> Vec<GeneRegion> {
    analysis.seed_indices.iter()
        .map(|&idx| genes[idx].clone())
        .collect()
}

/// Print a summary of the connectivity analysis
pub fn print_analysis_summary(analysis: &ConnectivityAnalysis, gene_names: &[String]) {
    println!("\n=== Connectivity Analysis ===");
    println!("  Total genes analyzed: {}", analysis.n_genes);
    println!("  Connected components: {}", analysis.n_components);
    println!(
        "  Effective Jaccard threshold: {:.6} (edges require >= {} shared reads)",
        analysis.effective_jaccard_threshold, analysis.min_shared_reads_used
    );
    println!("  Universe size N (unique reads across regions): {}", analysis.n_universe_reads);
    if let Some(p) = analysis.hypergeom_p_max_used {
        println!("  Hypergeometric filter active: p_max={:.6}", p);
    }
    println!("  Seeds selected: {} ({:.1}% reduction)", 
             analysis.n_seeds(), analysis.reduction_pct());
    
    println!("\n  Connected components:");
    for (i, comp) in analysis.components.iter().enumerate() {
        let names: Vec<String> = comp.iter()
            .map(|&idx| gene_names.get(idx).cloned().unwrap_or_else(|| format!("gene_{}", idx)))
            .collect();
        println!("    Component {}: {} genes", i + 1, comp.len());
    }
}
