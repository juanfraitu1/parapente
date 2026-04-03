use anyhow::{Context, Result};
use bstr::BString;
use clap::Parser;
use fxhash::FxHashSet as HashSet;
use noodles::bam;
use noodles::core::Position;
use noodles::sam::alignment::Record as SamRecord;
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
use std::sync::atomic::{AtomicUsize, Ordering};

mod fast_index;
use fast_index::{BamIndex, cluster_loci};

mod isoseq;
use isoseq::{analyze_locus_transcripts, validate_with_isoseq};

mod isoseq_advanced;
use isoseq_advanced::{validate_advanced, filter_by_advanced_features};

mod spliced_read;
use spliced_read::split_region_by_transcript_boundaries;

mod read_cloud;
use read_cloud::{find_read_clouds, filter_clouds, ReadCloud};

mod cloud_splitter;
use cloud_splitter::split_cloud_by_isoseq_signals;

mod splice_clouds;
use splice_clouds::{create_splice_subclouds, refine_sub_clouds};

mod coverage_valley;
use coverage_valley::ValleyDetectionParams;

mod coverage_validation;
use coverage_validation::{
    analyze_coverage_quality, validate_gene_core, CoverageQuality,
    ValidationResult, ValidationParams, ValidationConfidence as CoverageValConfidence
};

mod junction_refinement;

mod coverage_std_filter;
use coverage_std_filter::{
    adaptive_std_threshold, coverage_mean_and_std_in_region, depth_scaled_min_std,
};

mod integrated_detector;
use integrated_detector::{
    detect_gene_family_integrated, IntegratedParams, GeneCore, CoreConfidence
};
use junction_refinement::JunctionRefinementParams;

mod read_intern;

mod jaccard_cache;
use jaccard_cache::{JaccardCache, GLOBAL_JACCARD_CACHE};

mod coverage_cache;
use coverage_cache::{CoverageCache, GLOBAL_COVERAGE_CACHE, CachedCoverage};

mod transitive_detector;
use transitive_detector::{
    detect_transitive, detect_transitive_with_component_candidates,
    detect_transitive_multi_seed,
    detect_transitive_multi_seed_and_reads,
    split_loci_by_coverage, split_loci_by_coverage_par,
    trim_loci_read_align_envelope,
    TransitiveCoverageSplitConfig,
    TransitiveParams as TransDetParams, TransitiveLocus, ComponentSeedCandidate
};

mod boundary_refinement;

mod tight_locus_split;
use tight_locus_split::{tight_resplit_loci, TightSplitParams};

mod soft_clip;

mod clip_tail_trim;
use clip_tail_trim::{trim_loci_clip_tail, ClipTailTrimParams};

mod smart_seed_selector;
use smart_seed_selector::{
    analyze_connectivity, ConnectivityAnalysis, ConnectivityParams, GeneRegion, get_selected_seeds,
    print_analysis_summary,
};

use boundary_refinement::{
    refine_loci, write_refined_bed, IsoBoundaryMode, RefinementParams, RefinedLocus,
};
use pathfinding::matrix::Matrix;
use pathfinding::prelude::kuhn_munkres_min;

/// Default `--std-cv-bypass-mean-floor`: bypass stays off until you set a lower value explicitly.
const STD_CV_BYPASS_DISABLED_MEAN_FLOOR: f64 = 1_000_000_000.0;

/// Reference mean floor used only for `--debug-std-cv-bypass` hints (typical manual enable).
const STD_CV_BYPASS_DEBUG_HINT_MEAN_FLOOR: f64 = 100.0;

#[derive(Parser, Debug)]
#[command(name = "gene_family_detector")]
#[command(about = "Iso-Seq aware gene family detector")]
#[command(version)]
struct Args {
    #[arg(short, long)]
    bam: String,

    /// Seed chromosome (optional if using --seeds-bed or --seeds-list)
    #[arg(long, required = false)]
    seed_chrom: Option<String>,

    /// Seed start position (optional if using --seeds-bed or --seeds-list)
    #[arg(long, required = false)]
    seed_start: Option<usize>,

    /// Seed end position (optional if using --seeds-bed or --seeds-list)
    #[arg(long, required = false)]
    seed_end: Option<usize>,

    #[arg(short, long)]
    output: String,

    #[arg(long)]
    chromosomes: Option<String>,

    #[arg(long, default_value_t = 0.001)]
    min_jaccard: f64,

    /// Maximum iterations for discovery convergence (0 = unlimited until no new loci found)
    #[arg(long, default_value_t = 20)]
    max_iterations: usize,

    #[arg(long, default_value_t = 5)]
    min_reads: usize,

    #[arg(long, default_value_t = 10000)]
    cluster_distance: i64,

    /// Merge seed-read alignments into loci using a per-chromosome gap from same-read
    /// multi-map spacing (90th percentile of positive inter-alignment gaps), never above
    /// --cluster-distance. Off keeps a single fixed merge gap for the whole run.
    #[arg(long, default_value_t = false)]
    adaptive_cluster_distance: bool,

    #[arg(long)]
    include_supplementary: bool,

    #[arg(long, default_value_t = 4)]
    threads: usize,

    /// Minimum core reads required for a candidate locus.
    /// Core reads are those present in ALL seed regions. This filters artifacts
    /// that don't share the conserved multi-mapping reads of true family members.
    #[arg(long, default_value_t = 1)]
    min_core_reads: usize,

    /// Minimum number of seeds a locus must share reads with.
    /// For disjoint families (seeds don't share reads), this filters artifacts
    /// that only connect to a single seed. True family members typically share
    /// reads with multiple seeds through transitive connectivity.
    #[arg(long, default_value_t = 1)]
    min_seeds: usize,

    /// Minimum coverage quality score (0-1) for a locus to be reported.
    /// Coverage quality measures how "gene-like" the coverage profile is:
    /// - Symmetric peak shape
    /// - Steep boundaries (sharp transitions)
    /// - Consistent peak width
    /// Artifacts often have flat, noisy, or asymmetric coverage.
    #[arg(long, default_value_t = 0.0)]
    min_coverage_quality: f64,

    /// Use Iso-Seq validation
    #[arg(long)]
    use_isoseq: bool,

    /// Minimum Iso-Seq similarity score (0-1)
    #[arg(long, default_value_t = 0.6)]
    min_isoseq_sim: f64,

    /// Require read-space overlap (not just shared name)
    #[arg(long)]
    require_read_overlap: bool,

    /// Minimum read-space overlap in bp
    #[arg(long, default_value_t = 100)]
    min_read_overlap_bp: i64,

    /// Use splice-based splitting for tandem arrays
    #[arg(long)]
    split_splicing: bool,

    /// Use fast pre-indexed mode (builds in-memory index, much faster for multiple queries)
    #[arg(long)]
    fast: bool,

    /// Split loci using coverage valleys (separates merged genes)
    #[arg(long)]
    split_loci: bool,

    /// Maximum distance for clustering transcript boundaries (bp)
    #[arg(long, default_value_t = 5000)]
    max_boundary_distance: i64,

    /// Minimum transcripts required for boundary-based splitting
    #[arg(long, default_value_t = 10)]
    min_transcripts_for_split: usize,
    
    /// Expected number of genes in the region (optional)
    #[arg(long)]
    expected_genes: Option<usize>,

    /// Use advanced Iso-Seq validation (quality filters)
    #[arg(long)]
    use_advanced_isoseq: bool,

    /// Minimum MAPQ for Iso-Seq reads (DEPRECATED: Not used - multi-mapping reads have MAPQ=0)
    #[arg(long, default_value_t = 0)]
    min_mapq: u8,

    /// Maximum divergence for Iso-Seq reads
    #[arg(long, default_value_t = 0.05)]
    max_divergence: f64,

    /// Minimum strand consistency (0-1)
    #[arg(long, default_value_t = 0.7)]
    min_strand_consistency: f64,

    /// Use read cloud analysis (dense overlapping reads) instead of simple mapping
    #[arg(long)]
    use_read_clouds: bool,

    /// Use strand-aware cloud detection (separate clouds by strand)
    #[arg(long)]
    strand_aware_clouds: bool,

    /// Minimum number of overlapping reads to form a cloud
    #[arg(long, default_value_t = 5)]
    min_cloud_reads: usize,

    /// Minimum cloud span in bp
    #[arg(long, default_value_t = 1000)]
    min_cloud_span: i64,

    /// Maximum cloud span in bp
    #[arg(long, default_value_t = 500000)]
    max_cloud_span: i64,

    /// Minimum read density (reads per kb)
    #[arg(long, default_value_t = 1.0)]
    min_density: f64,

    /// Require secondary alignments in cloud
    #[arg(long)]
    require_secondary: bool,

    /// Split large clouds using Iso-Seq signals
    #[arg(long)]
    split_clouds: bool,

    /// Use splice-aware sub-clouds for splitting
    #[arg(long)]
    splice_subclouds: bool,

    /// Minimum splice similarity for clustering (0-1)
    #[arg(long, default_value_t = 0.3)]
    min_splice_sim: f64,

    /// Minimum Jaccard similarity for a locus to be reported (filters low-confidence FPs)
    #[arg(long, default_value_t = 0.01)]
    min_jaccard_locus: f64,

    /// Merge adjacent sub-loci within this distance (bp) with high Jaccard
    #[arg(long, default_value_t = 10000)]
    merge_adjacent_bp: i64,

    /// Minimum Jaccard for merging adjacent loci
    #[arg(long, default_value_t = 0.3)]
    merge_min_jaccard: f64,

    /// Use integrated mode: coverage valleys + Iso-Seq validation + read-space overlap
    #[arg(long)]
    integrated: bool,

    /// Disable transitive multi-seed detection (use legacy iterative single-component mode).
    #[arg(long = "no-transitive", default_value_t = false, action = clap::ArgAction::SetTrue)]
    no_transitive: bool,

    /// Ignored (transitive is default). Kept so older scripts that pass `--transitive` still run.
    #[allow(dead_code)]
    #[arg(long = "transitive", hide = true, action = clap::ArgAction::SetTrue)]
    _transitive_deprecated_noop: bool,

    /// Minimum binned coverage standard deviation (legacy Python `--min_std_coverage`).
    /// When adaptive std is on, this baseline is relaxed for windows with mean_bin below std-mean-ref
    /// (sqrt scaling, capped at the baseline for deeper windows), then span/shared adjustment applies.
    /// Default 50 provides good filtering of low-variability artifacts without losing true genes.
    #[arg(long, default_value_t = 50.0)]
    min_std_coverage: f64,

    /// Disable adaptive scaling of the std threshold for weak or short loci (Python `--no-dynamic-std-coverage`).
    #[arg(long, action = clap::ArgAction::SetTrue)]
    no_dynamic_std_coverage: bool,

    /// Reference mean reads per bin for depth scaling (default 200). Used with --min-std-coverage.
    #[arg(long, default_value_t = 200.0)]
    std_mean_ref: f64,

    /// With --min-std-coverage, keep a locus if std/mean_bin is at least this (gene-like variability)
    /// when mean_bin is at least std-cv-bypass-mean-floor. Supplements the absolute std threshold.
    #[arg(long, default_value_t = 0.085)]
    std_cv_bypass_min_cv: f64,

    /// Minimum mean reads per bin required for the coeff-var bypass. Default matches
    /// STD_CV_BYPASS_DISABLED_MEAN_FLOOR (bypass off); use e.g. 100 to opt in.
    #[arg(long, default_value_t = STD_CV_BYPASS_DISABLED_MEAN_FLOOR)]
    std_cv_bypass_mean_floor: f64,

    /// Stderr hints when a reject would have been kept with CV bypass at the debug reference
    /// mean floor (100) and current std-cv-bypass-min-cv. Does not change filtering.
    #[arg(long, default_value_t = false)]
    debug_std_cv_bypass: bool,

    #[arg(long, default_value_t = 8.0)]
    std_min_floor: f64,

    #[arg(long, default_value_t = 500.0)]
    std_shared_ref: f64,

    #[arg(long, default_value_t = 150000.0)]
    std_span_ref: f64,

    /// Bin size (bp) for coverage std (Python default 100).
    #[arg(long, default_value_t = 100)]
    coverage_bin_size: i64,

    /// In transitive mode, use fixed --valley-*, --peak-threshold-frac, --coverage-bin-size for every locus
    /// instead of per-locus dynamic valley thresholds (reduces overmerging when tuned stricter).
    #[arg(long, action = clap::ArgAction::SetTrue)]
    fixed_transitive_valleys: bool,

    /// After valley or peak splits, merge segments shorter than this (bp) with neighbors if gap <= transitive-valley-merge-gap-bp.
    #[arg(long, default_value_t = 5000)]
    transitive_valley_merge_min_segment_bp: i64,

    /// Max gap (bp) between segments for small-segment merging after valley splitting.
    #[arg(long, default_value_t = 1000)]
    transitive_valley_merge_max_gap_bp: i64,

    /// Loci wider than this (bp) are eligible for junction-based splitting after valley splitting.
    #[arg(long, default_value_t = 200000)]
    transitive_junction_min_locus_bp: i64,

    /// Inside junction splitting, only segments at least this wide (bp) are subdivided.
    #[arg(long, default_value_t = 200000)]
    transitive_junction_min_split_span_bp: i64,

    /// After junction splits, merge segments shorter than this (bp) when adjacent gap allows.
    #[arg(long, default_value_t = 10000)]
    transitive_junction_merge_min_segment_bp: i64,

    #[arg(long, default_value_t = 5000)]
    transitive_junction_merge_max_gap_bp: i64,

    /// Valley fraction: coverage minimum must be below (max * this) to count as a valley (lower = more sensitive)
    #[arg(long, default_value_t = 0.08)]
    valley_frac: f64,

    /// Peak height as fraction of global max (lower = more local peaks detected)
    #[arg(long, default_value_t = 0.12)]
    peak_threshold_frac: f64,

    /// Minimum separation between coverage peaks (bp); lower can split tandem genes more often
    #[arg(long, default_value_t = 2200)]
    valley_min_gap_bp: i64,

    /// Minimum segment length after a valley split (bp)
    #[arg(long, default_value_t = 1500)]
    valley_min_segment_bp: i64,

    /// Minimum prominence (peak minus valley) relative to taller peak; lower accepts more valleys
    #[arg(long, default_value_t = 0.32)]
    valley_min_prominence: f64,

    /// Auto-split threshold (valley depth ratio below this = auto split)
    #[arg(long, default_value_t = 0.1)]
    auto_split_threshold: f64,

    /// Validation threshold (valley depth ratio below this = use Iso-Seq validation)
    #[arg(long, default_value_t = 0.3)]
    validate_threshold: f64,

    /// Disable per-strand adaptive valley merging in integrated mode (default: enabled)
    #[arg(long, action = clap::ArgAction::SetTrue)]
    no_integrated_strand_split: bool,

    /// Minimum locus span (bp) to always run per-strand coverage graphs in integrated mode
    #[arg(long, default_value_t = 50000)]
    strand_adaptive_min_span_bp: i64,

    /// Integrated mode: tighten each core using CIGAR exon intervals (Skip-aware), no annotation
    #[arg(long, action = clap::ArgAction::SetTrue)]
    integrated_junction_refine: bool,

    /// Per-side cap as a fraction of span for junction refinement (exon union vs original bounds)
    #[arg(long, default_value_t = 0.35)]
    junction_refine_max_trim_frac: f64,

    /// Abort junction trim if the refined span would drop below this many bp
    #[arg(long, default_value_t = 500)]
    junction_refine_min_span: i64,

    /// Minimum reads per clustered splice acceptor anchor (0 disables acceptor anchor)
    #[arg(long, default_value_t = 3)]
    junction_refine_min_junction_reads: usize,

    #[arg(long, default_value_t = 12)]
    junction_refine_cluster_slack_bp: i64,

    /// Anti-overtrim: clamp to median alignment start/end plus this slack bp (0 disables)
    #[arg(long, default_value_t = 2000)]
    junction_refine_read_end_slack_bp: i64,

    #[arg(long, default_value_t = 2500)]
    junction_refine_acceptor_win_left_bp: i64,

    #[arg(long, default_value_t = 1200)]
    junction_refine_acceptor_win_right_bp: i64,

    /// BED file with multiple seed regions (one per line: chrom start end [name])
    #[arg(long)]
    seeds_bed: Option<String>,

    /// Comma-separated list of seeds (format: chrom:start-end,chrom:start-end)
    #[arg(long)]
    seeds_list: Option<String>,

    /// After coverage splitting and std filtering, tighten transitive locus bounds using
    /// percentiles of primary alignment start/end for reads in the locus (BAM-only).
    #[arg(long, default_value_t = false)]
    trim_read_envelope: bool,

    /// Low percentile for alignment starts (0 to 100), e.g. 5 = trim aggressive 5' tails.
    #[arg(long, default_value_t = 5.0)]
    read_envelope_pct_lo: f64,

    /// High percentile for alignment ends (0 to 100), e.g. 95 = trim 3' tails.
    #[arg(long, default_value_t = 95.0)]
    read_envelope_pct_hi: f64,

    /// If quantile span is below this, keep the pre-trim interval.
    #[arg(long, default_value_t = 2000)]
    read_envelope_min_span: i64,

    /// Minimum overlapping primary alignments required to apply the envelope.
    #[arg(long, default_value_t = 8)]
    read_envelope_min_alignments: usize,

    /// Bases to extend each side of the locus for the indexed BAM query window.
    #[arg(long, default_value_t = 500000)]
    read_envelope_query_pad_bp: i64,

    /// Optional final pass: stricter coverage valleys (and peak fallback) to split wide loci.
    /// BAM and seed-read coverage only; default off. Can add fragments; capped by tight-max-segments-per-locus.
    #[arg(long, default_value_t = false)]
    tight_coverage_resplit: bool,

    #[arg(long, default_value_t = 100)]
    tight_coverage_bin_size: i64,

    #[arg(long, default_value_t = 0.10)]
    tight_valley_frac: f64,

    #[arg(long, default_value_t = 0.12)]
    tight_peak_threshold_frac: f64,

    #[arg(long, default_value_t = 1200)]
    tight_min_gap_bp: i64,

    #[arg(long, default_value_t = 800)]
    tight_min_segment_bp: i64,

    #[arg(long, default_value_t = 0.20)]
    tight_min_prominence_frac: f64,

    #[arg(long, default_value_t = 2500)]
    tight_merge_min_segment_bp: i64,

    #[arg(long, default_value_t = 400)]
    tight_merge_max_gap_bp: i64,

    /// Only intervals at least this wide are considered for tight resplit.
    #[arg(long, default_value_t = 8000)]
    tight_min_input_span_bp: i64,

    /// If splitting would produce more segments, the original interval is kept unchanged.
    #[arg(long, default_value_t = 12)]
    tight_max_segments_per_locus: usize,

    /// Shrink loci using A-rich 3' soft clips (polyA proxy) on primary alignments. Default off.
    #[arg(long, default_value_t = false)]
    clip_tail_trim: bool,

    #[arg(long, default_value_t = 12)]
    clip_tail_min_len: usize,

    #[arg(long, default_value_t = 0.62)]
    clip_tail_min_a_frac: f64,

    #[arg(long, default_value_t = 6)]
    clip_tail_min_reads: usize,

    #[arg(long, default_value_t = 400)]
    clip_tail_end_slack_bp: i64,

    #[arg(long, default_value_t = 400)]
    clip_tail_start_slack_bp: i64,

    #[arg(long, default_value_t = 0.65)]
    clip_tail_strand_dom_frac: f64,

    #[arg(long, default_value_t = 500000)]
    clip_tail_query_pad_bp: i64,

    #[arg(long, default_value_t = 2000)]
    clip_tail_min_result_span: i64,

    /// Per-locus diagnostic lines for clip-tail trim (overlap, no 3' clip, short, low A, skip reason).
    #[arg(long, default_value_t = false)]
    clip_tail_verbose: bool,

    /// Refine boundaries using Iso-Seq transcripts and Jaccard trimming
    /// Note: For highly clustered families like ID_14, refinement may create gaps
    #[arg(long, default_value_t = false)]
    refine_boundaries: bool,

    /// Maximum fraction of region to trim during refinement (0-1)
    #[arg(long, default_value_t = 0.95)]
    max_trim_frac: f64,

    /// Minimum transcript density for Iso-Seq boundaries (transcripts per kb)
    #[arg(long, default_value_t = 0.05)]
    min_tx_density: f64,

    /// Iso-Seq coordinate model for --refine-boundaries (strand + CIGAR exon shell vs legacy percentiles)
    #[arg(long, value_enum, default_value_t = IsoBoundaryMode::StrandExon)]
    iso_boundary_mode: IsoBoundaryMode,

    /// Keep only the best segment (default: keep all significant segments)
    /// Use this for single-seed runs where you want maximum precision
    #[arg(long)]
    single_best_segment: bool,
    
    /// Minimum Jaccard for a segment to be kept (default: 0.01 to keep more regions)
    #[arg(long, default_value_t = 0.01)]
    min_segment_jaccard: f64,

    /// Ground truth BED file for evaluation (optional)
    /// If provided, will calculate sensitivity, precision, and overlap statistics
    #[arg(long)]
    ground_truth: Option<String>,

    /// Family ID for ground truth evaluation (e.g., ID_14)
    #[arg(long)]
    family_id: Option<String>,

    /// Use minimum seed selection algorithm
    /// Iteratively adds seeds until all ground truth genes are found
    #[arg(long)]
    min_seed: bool,

    /// Alternate non-circular mode: discover disconnected components via graph connected components,
    /// then add representative component loci as new seeds.
    #[arg(long)]
    connected_components_mode: bool,

    /// Improved non-circular adaptive discovery mode with stricter reseed guardrails
    #[arg(long)]
    adaptive_discovery_v2: bool,

    /// Use connectivity-based seed selection (Jaccard graph connected components)
    /// Requires --seeds-bed or --seeds-list to provide initial gene regions
    #[arg(long)]
    connectivity_seed: bool,

    /// With --connectivity-seed, raise the Jaccard edge cutoff to the q-quantile of observed
    /// pairwise Jaccards (pairs with at least --connectivity-min-shared-reads), never below --min-jaccard
    #[arg(long, default_value_t = false)]
    adaptive_connectivity: bool,

    /// Quantile q in [0,1] for --adaptive-connectivity (default 0.15: cutoff near lower tail of positive edges)
    #[arg(long, default_value_t = 0.15)]
    connectivity_adaptive_quantile: f64,

    /// Minimum shared read names required for a graph edge (and to enter the adaptive quantile sample)
    #[arg(long, default_value_t = 1)]
    connectivity_min_shared_reads: usize,

    /// With --connectivity-seed, require one-sided hypergeometric p <= this (read overlap vs random
    /// margins with universe N = unique reads across all regions). Used together with Jaccard rules.
    #[arg(long)]
    connectivity_hypergeom_p_max: Option<f64>,

    /// Exit after seed selection (for batch evaluation)
    #[arg(long)]
    seeds_only: bool,

    /// Maximum number of new component-derived seeds to add per iteration
    #[arg(long, default_value_t = 1)]
    cc_max_new_seeds_per_iter: usize,

    /// Minimum total reads in a disconnected component to be considered for reseeding
    #[arg(long, default_value_t = 20)]
    cc_min_component_reads: usize,

    /// Minimum exploratory candidate Jaccard to the pooled discovered-family reads
    #[arg(long, default_value_t = 0.02)]
    adv2_min_explore_jaccard: f64,

    /// Maximum exploratory seed span (bp)
    #[arg(long, default_value_t = 200000)]
    adv2_max_seed_span: i64,

    /// Minimum fraction of candidate locus reads that are shared with discovered family reads
    #[arg(long, default_value_t = 0.20)]
    adv2_min_shared_frac: f64,

    /// Minimum composite score for exploratory candidate acceptance
    #[arg(long, default_value_t = 0.12)]
    adv2_min_seed_score: f64,

    /// Minimum distance (bp) from known seeds/discovered loci for exploratory reseeds
    #[arg(long, default_value_t = 5000000)]
    adv2_min_novel_distance_bp: i64,

    /// Enable RNA skeleton scoring (multi-mapping + Iso-Seq features) for exploratory reseeds
    #[arg(long)]
    skeleton_mode: bool,

    /// Minimum skeleton score for exploratory reseed acceptance
    #[arg(long, default_value_t = 0.35)]
    adv2_skeleton_min_score: f64,

    /// Maximum number of exploratory candidates to evaluate with Iso-Seq skeleton scoring per iteration
    #[arg(long, default_value_t = 25)]
    adv2_skeleton_top_k: usize,

    /// Enable secondary RNA-only candidate source from transcript-density hotspots
    #[arg(long)]
    adv2_hotspot_source: bool,

    /// Bin size (bp) for hotspot scan
    #[arg(long, default_value_t = 50000)]
    adv2_hotspot_bin_bp: i64,

    /// Minimum alignments in a hotspot bin to propose it
    #[arg(long, default_value_t = 200)]
    adv2_hotspot_min_bin_reads: usize,

    /// Max number of hotspot seeds considered per iteration
    #[arg(long, default_value_t = 3)]
    adv2_hotspot_top_k: usize,

    /// Minimum fraction of hotspot-bin reads shared with current discovered family read pool
    #[arg(long, default_value_t = 0.05)]
    adv2_hotspot_min_shared_frac: f64,

    /// Minimum number of current seeds whose discovered loci must support a hotspot region
    #[arg(long, default_value_t = 2)]
    adv2_hotspot_min_support_seeds: usize,

    /// Max gap (bp) between hotspot and per-seed discovered loci to count as support
    #[arg(long, default_value_t = 500000)]
    adv2_hotspot_support_gap_bp: i64,

    /// Early-stop if loci count grows too quickly with little signal control
    #[arg(long, default_value_t = 2.0)]
    adv2_max_loci_growth_factor: f64,

    /// Maximum iterations for minimum seed selection
    #[arg(long, default_value_t = 5)]
    max_seed_iterations: usize,

    /// Minimum overlap fraction to count as "found" (default: 0.5 = 50%)
    #[arg(long, default_value_t = 0.5)]
    min_overlap_frac: f64,

    /// Ground-truth eval: legacy mode uses a single best detected interval per gene and
    /// requires overlap_bp/gene_span >= min-overlap-frac. Default is bipartite IoU matching.
    #[arg(long, default_value_t = false)]
    eval_legacy_best_interval: bool,

    /// Minimum genomic IoU (intersection over bounding span) for counting a GT gene as found
    /// when using default bipartite evaluation.
    #[arg(long, default_value_t = 0.3)]
    eval_min_iou: f64,

    /// With bipartite eval: count a gene as found if overlap_bp/gene_span is at least this,
    /// even when IoU is below eval-min-iou (e.g. extended detected loci that fully cover the gene).
    #[arg(long, default_value_t = 0.5)]
    eval_min_ref_frac: f64,

    /// If overlap_bp/gene_span vs any provided seed is at least this, the gene counts as found
    /// (sanity: all GT regions given as seeds yields 100% recall). Set to 1.01 to disable.
    #[arg(long, default_value_t = 0.95)]
    eval_seed_confirm_ref_frac: f64,
}

/// A seed region for gene family detection
#[derive(Debug, Clone)]
struct Seed {
    chrom: String,
    start: i64,
    end: i64,
    name: Option<String>,
}

impl Seed {
    fn span(&self) -> i64 {
        self.end - self.start
    }
}

#[derive(Debug, Clone)]
struct Locus {
    chrom: String,
    start: i64,
    end: i64,
    reads: HashSet<String>,
    isoseq_sim: Option<f64>,
    seed_spliced: Option<usize>,
    target_spliced: Option<usize>,
}

impl Locus {
    fn span(&self) -> i64 {
        self.end - self.start
    }

    fn is_valid(&self, min_sim: f64) -> bool {
        match self.isoseq_sim {
            Some(sim) => sim >= min_sim,
            None => true, // If no validation performed, accept
        }
    }
}

#[derive(Debug, Clone)]
struct SkeletonModel {
    mean: [f64; 5],
    std: [f64; 5],
}

/// Returns `Some(reason)` if the locus should be rejected (legacy Python std gate).
fn min_std_coverage_reject_reason(
    args: &Args,
    bam: &str,
    chrom: &str,
    start: i64,
    end: i64,
    shared_with_seed: usize,
) -> Result<Option<String>> {
    let base = args.min_std_coverage;
    if base <= 0.0 {
        return Ok(None); // Disabled
    }
    let (mean_bin, std_dev) = coverage_mean_and_std_in_region(
        bam,
        chrom,
        start,
        end,
        args.coverage_bin_size,
        args.include_supplementary,
    )?;
    let span = end - start;
    let base_scaled = depth_scaled_min_std(base, mean_bin, args.std_mean_ref);
    let thresh = if args.no_dynamic_std_coverage {
        base
    } else {
        adaptive_std_threshold(
            base_scaled,
            span,
            shared_with_seed,
            args.std_min_floor,
            args.std_shared_ref,
            args.std_span_ref,
        )
    };
    let cv = std_dev / mean_bin.max(1.0);
    let cv_bypass = mean_bin >= args.std_cv_bypass_mean_floor
        && cv >= args.std_cv_bypass_min_cv;
    println!(
        "    Coverage mean_bin={:.1} std={:.1} cv={:.3} (threshold {:.1}; depth_scaled_base {:.1}, min_std {}){}",
        mean_bin,
        std_dev,
        cv,
        thresh,
        base_scaled,
        base,
        if cv_bypass { " [cv bypass]" } else { "" }
    );
    if std_dev < thresh && !cv_bypass {
        if args.debug_std_cv_bypass {
            let hint_would_keep = mean_bin >= STD_CV_BYPASS_DEBUG_HINT_MEAN_FLOOR
                && cv >= args.std_cv_bypass_min_cv;
            if hint_would_keep {
                eprintln!(
                    "[debug std-cv-bypass] rejected on absolute std only; CV bypass would apply at mean floor >= {:.0} with std-cv-bypass-min-cv {:.3}: {}:{}-{} mean_bin={:.1} std={:.1} cv={:.3} thresh={:.1} (active floor {:.1})",
                    STD_CV_BYPASS_DEBUG_HINT_MEAN_FLOOR,
                    args.std_cv_bypass_min_cv,
                    chrom,
                    start,
                    end,
                    mean_bin,
                    std_dev,
                    cv,
                    thresh,
                    args.std_cv_bypass_mean_floor,
                );
                eprintln!(
                    "[debug std-cv-bypass] try: --std-cv-bypass-mean-floor {:.0}",
                    STD_CV_BYPASS_DEBUG_HINT_MEAN_FLOOR
                );
            }
        }
        Ok(Some(format!(
            "std_coverage {:.1} < {:.1}",
            std_dev, thresh
        )))
    } else {
        Ok(None)
    }
}

fn locus_overlaps_seed_region(
    chrom: &str,
    start: i64,
    end: i64,
    seeds: &[(String, i64, i64)],
) -> bool {
    seeds.iter().any(|(sc, ss, se)| {
        sc == chrom && !(end <= *ss || start >= *se)
    })
}

/// Compute coverage quality score for a locus.
/// Returns a score from 0.0-1.0 indicating how "gene-like" the coverage profile is.
fn compute_locus_coverage_quality(
    bam_path: &str,
    chrom: &str,
    start: i64,
    end: i64,
    bin_size: i64,
    _include_supplementary: bool,
) -> Result<f64> {
    use crate::coverage_valley::build_coverage_profile;
    
    // Collect all reads in this region (we want total coverage, not just seed reads)
    let all_reads = collect_seed_reads_fast(
        bam_path,
        chrom,
        start as usize,
        end as usize,
        true, // Include supplementary to get complete picture
    ).unwrap_or_default();
    
    // Build coverage profile for this locus
    let profile = build_coverage_profile(
        bam_path,
        chrom,
        start,
        end,
        &all_reads,
        bin_size,
    )?;
    
    if profile.bins.is_empty() {
        return Ok(0.0);
    }
    
    // Analyze coverage quality using the existing module
    let quality = analyze_coverage_quality(&profile.bins, bin_size);
    
    Ok(quality.overall_score)
}

/// Returns `Some(reason)` if the locus should be rejected due to low coverage quality.
fn coverage_quality_reject_reason(
    args: &Args,
    bam: &str,
    chrom: &str,
    start: i64,
    end: i64,
) -> Result<Option<String>> {
    if args.min_coverage_quality <= 0.0 {
        return Ok(None);
    }
    
    let score = compute_locus_coverage_quality(
        bam,
        chrom,
        start,
        end,
        args.coverage_bin_size,
        args.include_supplementary,
    )?;
    
    if score < args.min_coverage_quality {
        Ok(Some(format!(
            "coverage_quality {:.3} < {:.3}",
            score, args.min_coverage_quality
        )))
    } else {
        Ok(None)
    }
}

fn interval_gap(a_start: i64, a_end: i64, b_start: i64, b_end: i64) -> i64 {
    if a_end < b_start {
        b_start - a_end
    } else if b_end < a_start {
        a_start - b_end
    } else {
        0
    }
}

fn is_novel_seed_region(
    chrom: &str,
    start: i64,
    end: i64,
    existing_seeds: &[Seed],
    discovered_loci: &[TransitiveLocus],
    min_novel_distance_bp: i64,
) -> bool {
    let on_new_chrom = !existing_seeds.iter().any(|s| s.chrom == chrom);
    if on_new_chrom {
        return true;
    }

    let far_from_seeds = existing_seeds.iter().all(|s| {
        if s.chrom != chrom {
            return true;
        }
        interval_gap(start, end, s.start, s.end) >= min_novel_distance_bp
    });
    let far_from_discovered = discovered_loci.iter().all(|d| {
        if d.chrom != chrom {
            return true;
        }
        interval_gap(start, end, d.start, d.end) >= min_novel_distance_bp
    });
    far_from_seeds && far_from_discovered
}

fn propose_hotspot_seeds_from_transcript_density(
    args: &Args,
    target_chroms: &[String],
    discovered_loci: &[TransitiveLocus],
    existing_seeds: &[Seed],
    per_seed_loci: &[Vec<TransitiveLocus>],
) -> Result<Vec<Seed>> {
    let pooled_reads: HashSet<String> = discovered_loci
        .iter()
        .flat_map(|l| l.reads.iter().cloned())
        .collect();
    let mut hotspots: Vec<(String, i64, i64, usize, usize)> = Vec::new(); // chrom,start,end,total,shared

    for chrom in target_chroms {
        let mut reader = bam::io::indexed_reader::Builder::default()
            .build_from_path(&args.bam)?;
        let header = reader.read_header()?;

        let chrom_name: bstr::BString = chrom.as_bytes().into();
        let ref_seqs = header.reference_sequences();
        let tid = match ref_seqs.get_index_of(&chrom_name) {
            Some(t) => t,
            None => continue,
        };
        let ref_len = ref_seqs[tid].length().get() as i64;

        let bin_bp = args.adv2_hotspot_bin_bp.max(1000);
        let n_bins = ((ref_len + bin_bp - 1) / bin_bp) as usize;
        let mut bins_total = vec![0usize; n_bins];
        let mut bins_shared = vec![0usize; n_bins];

        let start = Position::MIN;
        let end = Position::new(ref_len as usize).unwrap_or(Position::MAX);
        let region = noodles::core::Region::new(chrom.clone(), start..=end);

        for result in reader.query(&header, &region)? {
            let record = result?;
            if record.flags().is_unmapped() {
                continue;
            }
            let read_name = match record.name() {
                Some(n) => match std::str::from_utf8(n) {
                    Ok(s) => s,
                    Err(_) => continue,
                },
                None => continue,
            };
            let aln_start = record.alignment_start()
                .transpose()?
                .map(|p: Position| p.get() as i64)
                .unwrap_or(1);
            let idx = ((aln_start - 1) / bin_bp).clamp(0, n_bins as i64 - 1) as usize;
            bins_total[idx] += 1;
            if pooled_reads.contains(read_name) {
                bins_shared[idx] += 1;
            }
        }

        for i in 0..n_bins {
            let c = bins_total[i];
            let shared = bins_shared[i];
            let shared_frac = if c > 0 { shared as f64 / c as f64 } else { 0.0 };
            if c >= args.adv2_hotspot_min_bin_reads && shared_frac >= args.adv2_hotspot_min_shared_frac {
                let s = i as i64 * bin_bp + 1;
                let e = ((i as i64 + 1) * bin_bp).min(ref_len);
                hotspots.push((chrom.clone(), s, e, c, shared));
            }
        }
    }

    hotspots.sort_by(|a, b| {
        b.4.cmp(&a.4) // prioritize shared-support
            .then(b.3.cmp(&a.3))
    });

    let mut seeds: Vec<Seed> = Vec::new();
    let mut used_chroms: HashSet<String> = HashSet::default();
    for (chrom, start, end, count, shared_count) in hotspots {
        let span = end - start;
        if span > args.adv2_max_seed_span || span < 1000 {
            continue;
        }
        if used_chroms.contains(&chrom) {
            continue;
        }
        if !is_novel_seed_region(
            &chrom,
            start,
            end,
            existing_seeds,
            discovered_loci,
            args.adv2_min_novel_distance_bp,
        ) {
            continue;
        }

        let support_seed_count = per_seed_loci
            .iter()
            .filter(|loci| {
                loci.iter().any(|l| {
                    l.chrom == chrom
                        && interval_gap(start, end, l.start, l.end) <= args.adv2_hotspot_support_gap_bp
                })
            })
            .count();
        if support_seed_count < args.adv2_hotspot_min_support_seeds {
            continue;
        }

        seeds.push(Seed {
            chrom: chrom.clone(),
            start,
            end,
            name: Some(format!("hotspot_seed_{}_{}_s{}", count, shared_count, support_seed_count)),
        });
        used_chroms.insert(chrom);

        if seeds.len() >= args.adv2_hotspot_top_k.min(args.cc_max_new_seeds_per_iter) {
            break;
        }
    }

    Ok(seeds)
}

fn locus_feature_vector(
    reads_len: usize,
    start: i64,
    end: i64,
    tx_count: usize,
    coverage_std: f64,
) -> [f64; 5] {
    let span = (end - start).max(1) as f64;
    let density = reads_len as f64 / (span / 1000.0).max(1.0);
    [
        (reads_len as f64 + 1.0).ln(),
        (span + 1.0).ln(),
        (density + 1.0).ln(),
        (tx_count as f64 + 1.0).ln(),
        (coverage_std + 1.0).ln(),
    ]
}

fn build_skeleton_model(args: &Args, loci: &[TransitiveLocus]) -> Result<Option<SkeletonModel>> {
    let trusted: Vec<&TransitiveLocus> = loci
        .iter()
        .filter(|l| l.distance_from_seed <= 1)
        .collect();
    if trusted.len() < 2 {
        return Ok(None);
    }

    let mut rows: Vec<[f64; 5]> = Vec::new();
    for l in trusted {
        let tx = analyze_locus_transcripts(&args.bam, &l.chrom, l.start, l.end)?;
        rows.push(locus_feature_vector(
            l.reads.len(),
            l.start,
            l.end,
            tx.total_transcripts,
            tx.coverage_std,
        ));
    }
    if rows.len() < 2 {
        return Ok(None);
    }

    let mut mean = [0.0; 5];
    for r in &rows {
        for i in 0..5 {
            mean[i] += r[i];
        }
    }
    for m in &mut mean {
        *m /= rows.len() as f64;
    }

    let mut std = [0.0; 5];
    for r in &rows {
        for i in 0..5 {
            std[i] += (r[i] - mean[i]).powi(2);
        }
    }
    for s in &mut std {
        *s = (*s / rows.len() as f64).sqrt().max(0.20);
    }

    Ok(Some(SkeletonModel { mean, std }))
}

fn skeleton_score(model: &SkeletonModel, fv: &[f64; 5]) -> f64 {
    let mut z_sum = 0.0;
    for (i, val) in fv.iter().enumerate().take(5) {
        z_sum += (*val - model.mean[i]).abs() / model.std[i];
    }
    let z_avg = z_sum / 5.0;
    1.0 / (1.0 + z_avg)
}

fn propose_exploratory_seeds_from_discovered(
    args: &Args,
    target_chroms: &[String],
    discovered_loci: &[TransitiveLocus],
    existing_seeds: &[Seed],
    strict_v2: bool,
    skeleton_model: Option<&SkeletonModel>,
) -> Result<Vec<Seed>> {
    let mut pooled_reads: HashSet<String> = HashSet::default();
    for locus in discovered_loci {
        pooled_reads.extend(locus.reads.iter().cloned());
    }
    if pooled_reads.is_empty() {
        return Ok(vec![]);
    }

    let mut exclude_regions: Vec<(String, i64, i64)> = discovered_loci
        .iter()
        .map(|l| (l.chrom.clone(), l.start, l.end))
        .collect();
    exclude_regions.extend(
        existing_seeds
            .iter()
            .map(|s| (s.chrom.clone(), s.start, s.end)),
    );

    // Mine additional unlabeled loci from discovered read pools.
    let exploratory_loci = find_where_reads_map_parallel(
        &args.bam,
        &pooled_reads,
        target_chroms,
        &exclude_regions,
        args.include_supplementary,
        args.cluster_distance,
        args.min_reads,
        args.threads,
    )?;

    if exploratory_loci.is_empty() {
        return Ok(vec![]);
    }

    let mut ranked = exploratory_loci;
    ranked.sort_by(|a, b| {
        b.reads
            .len()
            .cmp(&a.reads.len())
            .then((b.end - b.start).cmp(&(a.end - a.start)))
    });

    let mut seeds: Vec<Seed> = Vec::new();
    let mut used_chroms: HashSet<String> = HashSet::default();
    let mut skeleton_evaluated = 0usize;
    for locus in ranked {
        let span = locus.end - locus.start;
        if locus.reads.len() < args.cc_min_component_reads {
            continue;
        }
        if strict_v2 && (span > args.adv2_max_seed_span || span < 1000) {
            continue;
        }
        if strict_v2 {
            let sim = calculate_jaccard(&locus.reads, &pooled_reads);
            if sim < args.adv2_min_explore_jaccard {
                continue;
            }
            let intersection = locus.reads.iter().filter(|r| pooled_reads.contains(*r)).count();
            let shared_frac = if !locus.reads.is_empty() {
                intersection as f64 / locus.reads.len() as f64
            } else {
                0.0
            };
            if shared_frac < args.adv2_min_shared_frac {
                continue;
            }
            let density = if span > 0 {
                locus.reads.len() as f64 / ((span as f64) / 1000.0).max(1.0)
            } else {
                0.0
            };
            let density_score = (density / 20.0).min(1.0);
            let seed_score = (0.50 * sim) + (0.35 * shared_frac) + (0.15 * density_score);
            if seed_score < args.adv2_min_seed_score {
                continue;
            }
            if let Some(model) = skeleton_model {
                if skeleton_evaluated >= args.adv2_skeleton_top_k {
                    continue;
                }
                let tx = analyze_locus_transcripts(&args.bam, &locus.chrom, locus.start, locus.end)?;
                let fv = locus_feature_vector(
                    locus.reads.len(),
                    locus.start,
                    locus.end,
                    tx.total_transcripts,
                    tx.coverage_std,
                );
                let sk_score = skeleton_score(model, &fv);
                skeleton_evaluated += 1;
                if sk_score < args.adv2_skeleton_min_score {
                    continue;
                }
            }
            if used_chroms.contains(&locus.chrom) {
                continue;
            }
        }

        let near_existing = existing_seeds.iter().any(|s| {
            s.chrom == locus.chrom
                && (s.start - locus.start).abs() < 20000
                && (s.end - locus.end).abs() < 20000
        });
        if near_existing {
            continue;
        }
        if strict_v2 {
            if !is_novel_seed_region(
                &locus.chrom,
                locus.start,
                locus.end,
                existing_seeds,
                discovered_loci,
                args.adv2_min_novel_distance_bp,
            ) {
                continue;
            }
        }

        seeds.push(Seed {
            chrom: locus.chrom.clone(),
            start: locus.start,
            end: locus.end,
            name: Some(format!("explore_seed_{}", locus.reads.len())),
        });
        if strict_v2 {
            used_chroms.insert(locus.chrom.clone());
        }

        if seeds.len() >= args.cc_max_new_seeds_per_iter {
            break;
        }
    }

    Ok(seeds)
}

/// Ground truth gene annotation
#[derive(Debug, Clone)]
struct Gene {
    chrom: String,
    start: i64,
    end: i64,
    name: String,
    span: i64,
}

/// Evaluation result for a single gene
#[derive(Debug, Clone)]
struct GeneEvaluation {
    gene_name: String,
    found: bool,
    overlap_bp: i64,
    overlap_pct: f64,
    detected_region: Option<(String, i64, i64)>,
    /// Genomic IoU for the bipartite-assigned locus, if applicable.
    match_iou: Option<f64>,
    /// Matched by provided seed interval (eval-seed-confirm-ref-frac), not only by detector output.
    confirmed_by_seed: bool,
}

/// Per-locus support for ranking novel (unmatched) candidates: len must match detected loci order.
#[derive(Debug, Clone)]
struct NovelLocusSupport {
    n_reads: Vec<usize>,
    jaccard: Vec<f64>,
    distance_from_seed: Vec<usize>,
}

#[derive(Debug, Clone)]
struct NovelLocusScore {
    chrom: String,
    start: i64,
    end: i64,
    n_reads: usize,
    jaccard: f64,
    distance_from_seed: usize,
    rank_score: f64,
}

/// Comprehensive evaluation statistics
#[derive(Debug, Clone)]
struct EvaluationStats {
    total_genes: usize,
    genes_found: usize,
    genes_missed: usize,
    /// Detected loci not assigned as a qualifying match to any GT (novel candidates).
    extra_regions: usize,
    sensitivity: f64,
    precision: f64,
    avg_overlap_pct: f64,
    total_gt_span: i64,
    total_detected_span: i64,
    per_gene: Vec<GeneEvaluation>,
    /// True when stats used bipartite IoU matching.
    eval_bipartite: bool,
    /// Novel detections ranked by read support and seed Jaccard (strongest first).
    novel_candidates_ranked: Vec<NovelLocusScore>,
}

fn print_separator() {
    println!("{}", "=".repeat(80));
}

/// Fast detection using pre-built BAM index
fn detect_with_index(
    bam_index: &BamIndex,
    seed_chrom: &str,
    seed_start: usize,
    seed_end: usize,
    include_supplementary: bool,
    cluster_distance: i64,
    min_reads: usize,
) -> Result<Vec<Locus>> {
    use fxhash::FxHashSet as HashSet;
    
    // Collect seed reads from the index
    let mut seed_reads: HashSet<String> = HashSet::default();
    
    // We need to get reads from the seed region
    // For now, query the BAM for seed reads (this is fast since it's a small region)
    let mut reader = bam::io::indexed_reader::Builder::default()
        .build_from_path("")?; // Placeholder - we need the actual path
    
    // Actually, let's use the index directly
    // Find reads that align to the seed region
    for (read_name, alignments) in bam_index.index.iter() {
        let read_name: &String = read_name;
        for aln in alignments {
            if aln.chrom == seed_chrom 
                && aln.start < seed_end as i64 
                && aln.end > seed_start as i64 {
                seed_reads.insert(read_name.clone());
                break;
            }
        }
    }
    
    println!("  Seed reads from index: {}", seed_reads.len());
    
    if seed_reads.is_empty() {
        return Ok(vec![]);
    }
    
    // Get all alignments for seed reads
    let by_chrom = bam_index.get_alignments_by_chrom(&seed_reads);
    
    // Cluster into loci
    let loci = cluster_loci(by_chrom, cluster_distance, min_reads);
    
    // Convert to Locus structs
    let result: Vec<Locus> = loci
        .into_iter()
        .map(|(chrom, start, end, reads)| Locus {
            chrom,
            start,
            end,
            reads,
            isoseq_sim: None,
            seed_spliced: None,
            target_spliced: None,
        })
        .collect();
    
    Ok(result)
}

fn collect_seed_reads_fast(
    bam_path: &str,
    chrom: &str,
    start: usize,
    end: usize,
    include_supplementary: bool,
) -> Result<HashSet<String>> {
    let mut reader = bam::io::indexed_reader::Builder::default()
        .build_from_path(bam_path)
        .with_context(|| format!("Failed to open BAM: {}", bam_path))?;

    let header = reader.read_header()?;
    let reference_sequences = header.reference_sequences();

    let chrom_name: BString = chrom.into();
    let _chrom_id = reference_sequences
        .get_index_of(&chrom_name)
        .ok_or_else(|| anyhow::anyhow!("Chromosome {} not found", chrom))?;

    let start_pos = Position::try_from(start)?;
    let end_pos = Position::try_from(end)?;
    let region = noodles::core::Region::new(chrom, start_pos..=end_pos);

    let mut reads = HashSet::with_capacity_and_hasher(4096, Default::default());

    for result in reader.query(&header, &region)? {
        let record = result?;

        if record.flags().is_unmapped() {
            continue;
        }
        if !include_supplementary && record.flags().is_supplementary() {
            continue;
        }
        // Include secondary alignments (0x100 flag) - they are vital for multi-mapping detection
        // No filter here - we want ALL alignments including secondary

        if let Some(name) = record.name() {
            if let Ok(s) = std::str::from_utf8(name) {
                reads.insert(s.to_string());
            }
        }
    }

    Ok(reads)
}

fn scan_chromosome(
    bam_path: &str,
    chrom: &str,
    target_reads: &HashSet<String>,
    exclude_regions: &[(String, i64, i64)],
    include_supplementary: bool,
) -> Result<Vec<(String, i64, i64, String)>> {
    let mut reader = bam::io::indexed_reader::Builder::default()
        .build_from_path(bam_path)
        .with_context(|| format!("Failed to open BAM: {}", bam_path))?;

    let header = reader.read_header()?;
    let reference_sequences = header.reference_sequences();

    let chrom_name: BString = chrom.into();
    let tid = reference_sequences
        .get_index_of(&chrom_name)
        .ok_or_else(|| anyhow::anyhow!("Chromosome {} not found", chrom))?;

    let start = Position::MIN;
    let ref_seq = &reference_sequences[tid];
    let end = Position::new(ref_seq.length().get()).unwrap_or(Position::MAX);
    let region = noodles::core::Region::new(chrom, start..=end);

    let mut matches = Vec::with_capacity(1024);

    for result in reader.query(&header, &region)? {
        let record = result?;

        if record.flags().is_unmapped() {
            continue;
        }
        if !include_supplementary && record.flags().is_supplementary() {
            continue;
        }
        // Include secondary alignments - vital for multi-mapping detection

        let name_bytes = match record.name() {
            Some(n) => n,
            None => continue,
        };

        let name = match std::str::from_utf8(name_bytes) {
            Ok(s) => s,
            Err(_) => continue,
        };

        if !target_reads.contains(name) {
            continue;
        }

        let start = record.alignment_start()
            .transpose()?
            .map(|p| usize::from(p) as i64)
            .unwrap_or(0);
        let end = record.alignment_end()
            .transpose()?
            .map(|p| usize::from(p) as i64)
            .unwrap_or(start + 1);

        let mut skip = false;
        for (exc_chrom, exc_start, exc_end) in exclude_regions {
            if chrom == exc_chrom && start < *exc_end && end > *exc_start {
                skip = true;
                break;
            }
        }

        if !skip {
            matches.push((chrom.to_string(), start, end, name.to_string()));
        }
    }

    Ok(matches)
}

fn find_where_reads_map_parallel(
    bam_path: &str,
    target_reads: &HashSet<String>,
    target_chroms: &[String],
    exclude_regions: &[(String, i64, i64)],
    include_supplementary: bool,
    cluster_distance: i64,
    min_reads: usize,
    num_threads: usize,
) -> Result<Vec<Locus>> {
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();

    println!(
        "  Scanning {} chromosomes using {} threads...",
        target_chroms.len(),
        num_threads
    );

    let scanned_count = AtomicUsize::new(0);
    let matched_count = AtomicUsize::new(0);

    let all_matches: Vec<(String, i64, i64, String)> = pool.install(|| {
        target_chroms
            .par_iter()
            .map(|chrom| {
                let result = scan_chromosome(
                    bam_path,
                    chrom,
                    target_reads,
                    exclude_regions,
                    include_supplementary,
                );
                
                match result {
                    Ok(matches) => {
                        scanned_count.fetch_add(1, Ordering::Relaxed);
                        matched_count.fetch_add(matches.len(), Ordering::Relaxed);
                        matches
                    }
                    Err(e) => {
                        eprintln!("  Error scanning {}: {}", chrom, e);
                        vec![]
                    }
                }
            })
            .flatten()
            .collect()
    });

    println!(
        "  Scanned {} chromosomes, found {} matches",
        scanned_count.load(Ordering::Relaxed),
        all_matches.len()
    );

    let mut by_chrom: HashMap<String, Vec<(i64, i64, String)>> = HashMap::new();
    for (chrom, start, end, name) in all_matches {
        by_chrom.entry(chrom).or_default().push((start, end, name));
    }

    let mut raw_loci: Vec<Locus> = Vec::new();

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
                    raw_loci.push(Locus {
                        chrom: chrom.clone(),
                        start: current_cluster[0].0,
                        end: current_end,
                        reads,
                        isoseq_sim: None,
                        seed_spliced: None,
                        target_spliced: None,
                    });
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
            raw_loci.push(Locus {
                chrom,
                start: current_cluster[0].0,
                end: current_end,
                reads,
                isoseq_sim: None,
                seed_spliced: None,
                target_spliced: None,
            });
        }
    }

    println!("  Found {} candidate loci after clustering", raw_loci.len());

    Ok(raw_loci)
}

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

    intersection as f64 / union as f64
}

/// Filter loci by minimum Jaccard similarity to seed
fn filter_by_jaccard(
    loci: Vec<(Locus, f64, usize)>,
    min_jaccard: f64,
) -> Vec<(Locus, f64, usize)> {
    let before = loci.len();
    let filtered: Vec<_> = loci.into_iter()
        .filter(|(_, jaccard, _)| *jaccard >= min_jaccard)
        .collect();
    let after = filtered.len();
    
    if before > after {
        println!("  Filtered by Jaccard (≥{:.3}): {} → {} loci", min_jaccard, before, after);
    }
    filtered
}

/// Merge adjacent sub-loci that are close together and have high Jaccard
fn merge_adjacent_loci(
    loci: Vec<(Locus, f64, usize)>,
    max_gap_bp: i64,
    min_jaccard: f64,
) -> Vec<(Locus, f64, usize)> {
    if loci.len() < 2 {
        return loci;
    }
    
    // Sort by chrom, start
    let mut sorted = loci;
    sorted.sort_by(|a, b| {
        a.0.chrom.cmp(&b.0.chrom)
            .then(a.0.start.cmp(&b.0.start))
    });
    
    let before = sorted.len();
    let mut merged: Vec<(Locus, f64, usize)> = Vec::new();
    let mut current = sorted[0].clone();
    
    for (locus, jaccard, component) in sorted.into_iter().skip(1) {
        // Check if adjacent on same chrom with small gap and both have high Jaccard
        let gap = locus.start - current.0.end;
        let same_chrom = locus.chrom == current.0.chrom;
        let close = gap <= max_gap_bp && gap >= 0;
        let both_high_jaccard = jaccard >= min_jaccard && current.1 >= min_jaccard;
        
        if same_chrom && close && both_high_jaccard {
            // Merge: extend end, union reads, use max Jaccard
            current.0.end = current.0.end.max(locus.end);
            current.0.reads.extend(locus.reads);
            current.1 = current.1.max(jaccard); // Use higher Jaccard
            // Keep component from first locus
        } else {
            merged.push(current);
            current = (locus, jaccard, component);
        }
    }
    merged.push(current);
    
    let after = merged.len();
    if before > after {
        println!("  Merged adjacent loci: {} → {} loci", before, after);
    }
    
    merged
}

/// Parse seeds from a BED file
/// Format: chrom\tstart\tend\t[name]
fn parse_seeds_bed(bed_path: &str) -> Result<Vec<Seed>> {
    let content = std::fs::read_to_string(bed_path)
        .with_context(|| format!("Failed to read seeds BED file: {}", bed_path))?;
    
    let mut seeds = Vec::new();
    
    for (line_num, line) in content.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 3 {
            eprintln!("  Warning: Skipping line {} (expected at least 3 fields)", line_num + 1);
            continue;
        }
        
        let chrom = fields[0].to_string();
        let start: i64 = fields[1].parse()
            .with_context(|| format!("Invalid start position on line {}: {}", line_num + 1, fields[1]))?;
        let end: i64 = fields[2].parse()
            .with_context(|| format!("Invalid end position on line {}: {}", line_num + 1, fields[2]))?;
        
        let name = if fields.len() >= 4 {
            Some(fields[3].to_string())
        } else {
            None
        };
        
        seeds.push(Seed { chrom, start, end, name });
    }
    
    println!("  Parsed {} seeds from BED file", seeds.len());
    Ok(seeds)
}

/// Parse seeds from a comma-separated list
/// Format: chrom:start-end,chrom:start-end
fn parse_seeds_list(list: &str) -> Result<Vec<Seed>> {
    let mut seeds = Vec::new();
    
    for seed_str in list.split(',') {
        let seed_str = seed_str.trim();
        if seed_str.is_empty() {
            continue;
        }
        
        // Parse format: chrom:start-end
        let parts: Vec<&str> = seed_str.split(':').collect();
        if parts.len() != 2 {
            anyhow::bail!("Invalid seed format: {} (expected chrom:start-end)", seed_str);
        }
        
        let chrom = parts[0].to_string();
        
        let range_parts: Vec<&str> = parts[1].split('-').collect();
        if range_parts.len() != 2 {
            anyhow::bail!("Invalid range format: {} (expected start-end)", parts[1]);
        }
        
        let start: i64 = range_parts[0].parse()
            .with_context(|| format!("Invalid start position: {}", range_parts[0]))?;
        let end: i64 = range_parts[1].parse()
            .with_context(|| format!("Invalid end position: {}", range_parts[1]))?;
        
        seeds.push(Seed { chrom, start, end, name: None });
    }
    
    println!("  Parsed {} seeds from list", seeds.len());
    Ok(seeds)
}

/// Collect seeds from all sources (CLI single seed, BED file, or list)
/// At least one seed source must be provided (CLI args, --seeds-bed, or --seeds-list)
fn collect_seeds(args: &Args) -> Result<Vec<Seed>> {
    let mut all_seeds = Vec::new();
    
    // Add CLI seed if provided
    if let (Some(chrom), Some(start), Some(end)) = (&args.seed_chrom, args.seed_start, args.seed_end) {
        all_seeds.push(Seed {
            chrom: chrom.clone(),
            start: start as i64,
            end: end as i64,
            name: Some("seed".to_string()),
        });
    }
    
    // Add seeds from BED file if provided
    if let Some(ref bed_path) = args.seeds_bed {
        let bed_seeds = parse_seeds_bed(bed_path)?;
        all_seeds.extend(bed_seeds);
    }
    
    // Add seeds from list if provided
    if let Some(ref list) = args.seeds_list {
        let list_seeds = parse_seeds_list(list)?;
        all_seeds.extend(list_seeds);
    }
    
    // Ensure at least one seed was provided
    if all_seeds.is_empty() {
        anyhow::bail!("At least one seed must be provided via --seed-chrom/--seed-start/--seed-end, --seeds-bed, or --seeds-list");
    }
    
    // Remove duplicate seeds (same position)
    let mut unique_seeds: Vec<Seed> = Vec::new();
    for seed in all_seeds {
        let is_duplicate = unique_seeds.iter().any(|s| {
            s.chrom == seed.chrom && 
            (s.start - seed.start).abs() < 100 && 
            (s.end - seed.end).abs() < 100
        });
        if !is_duplicate {
            unique_seeds.push(seed);
        }
    }
    
    Ok(unique_seeds)
}

/// Compute core reads: intersection of reads from ALL seed regions
/// These are the multi-mapping reads present in every true family member
/// Compute core reads per connected component of seeds.
/// Seeds form a graph where edges exist if seeds share reads.
/// For each connected component, we compute the intersection of reads
/// across all seeds in that component. This handles families with
/// multiple sub-clusters that share different sets of reads.
fn compute_core_reads_by_component(
    bam_path: &str,
    seeds: &[Seed],
    include_supplementary: bool,
) -> Result<Vec<HashSet<String>>> {
    use noodles::bam;
    use noodles::core::Position;
    use noodles::sam::alignment::Record;
    
    if seeds.is_empty() {
        return Ok(vec![]);
    }
    
    // Collect reads from each seed
    let mut seed_read_sets: Vec<HashSet<String>> = Vec::new();
    
    for seed in seeds {
        let mut reader = bam::io::indexed_reader::Builder::default()
            .build_from_path(bam_path)?;
        let header = reader.read_header()?;
        
        let region_start = Position::try_from(seed.start.max(1) as usize)?;
        let region_end = Position::try_from(seed.end.max(1) as usize)?;
        let region = noodles::core::Region::new(seed.chrom.clone(), region_start..=region_end);
        
        let mut reads: HashSet<String> = HashSet::default();
        
        for result in reader.query(&header, &region)? {
            let record = result?;
            
            if record.flags().is_unmapped() {
                continue;
            }
            if !include_supplementary && record.flags().is_supplementary() {
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
        
        seed_read_sets.push(reads);
    }
    
    // Build adjacency graph: two seeds are connected if they share reads
    let n = seeds.len();
    let mut adjacency: Vec<Vec<usize>> = vec![vec![]; n];
    
    for i in 0..n {
        for j in (i + 1)..n {
            let shared: HashSet<_> = seed_read_sets[i].intersection(&seed_read_sets[j]).collect();
            if !shared.is_empty() {
                adjacency[i].push(j);
                adjacency[j].push(i);
            }
        }
    }
    
    // Find connected components using BFS
    let mut visited = vec![false; n];
    let mut components: Vec<Vec<usize>> = vec![];
    
    for start in 0..n {
        if visited[start] {
            continue;
        }
        
        let mut component = vec![];
        let mut queue = vec![start];
        visited[start] = true;
        
        while let Some(node) = queue.pop() {
            component.push(node);
            for &neighbor in &adjacency[node] {
                if !visited[neighbor] {
                    visited[neighbor] = true;
                    queue.push(neighbor);
                }
            }
        }
        
        components.push(component);
    }
    
    // Compute core reads for each component (intersection of all seeds in component)
    let mut component_cores: Vec<HashSet<String>> = vec![];
    
    for component in &components {
        if component.is_empty() {
            continue;
        }
        
        // Start with reads from first seed in component
        let mut core_reads = seed_read_sets[component[0]].clone();
        
        // Intersect with all other seeds in component
        for &seed_idx in &component[1..] {
            core_reads.retain(|r| seed_read_sets[seed_idx].contains(r));
        }
        
        if !core_reads.is_empty() {
            component_cores.push(core_reads);
        }
    }
    
    Ok(component_cores)
}

/// Merge gene cores from multiple seeds, removing duplicates
fn merge_cores_from_multiple_seeds(all_cores: Vec<Vec<GeneCore>>) -> Vec<GeneCore> {
    if all_cores.is_empty() {
        return vec![];
    }
    
    // Flatten all cores
    let mut merged: Vec<GeneCore> = all_cores.into_iter().flatten().collect();
    
    // Sort by chrom, then by span DESCENDING (larger loci first)
    // This ensures we keep the larger locus when deduplicating
    merged.sort_by(|a, b| {
        a.chrom.cmp(&b.chrom)
            .then((b.end - b.start).cmp(&(a.end - a.start)))
    });
    
    // Remove duplicates: same segment reported from multiple seeds must overlap heavily.
    // Since we sorted by span descending, we keep the LARGER locus.
    let mut unique: Vec<GeneCore> = vec![];

    for core in merged {
        let core_span = core.end - core.start;
        let is_duplicate = unique.iter().any(|u| {
            if u.chrom != core.chrom {
                return false;
            }
            let o0 = core.start.max(u.start);
            let o1 = core.end.min(u.end);
            let overlap = (o1 - o0).max(0);
            if overlap <= 0 {
                return false;
            }
            let uspan = u.end - u.start;
            let smaller = core_span.min(uspan).max(1);
            overlap * 10 >= 8 * smaller
        });

        if !is_duplicate {
            unique.push(core);
        }
    }
    
    // Sort by position for consistent output
    unique.sort_by(|a, b| {
        a.chrom.cmp(&b.chrom)
            .then(a.start.cmp(&b.start))
    });
    
    unique
}

/// Parse ground truth BED file and extract genes for a specific family
fn parse_ground_truth(bed_path: &str, family_id: &str) -> Result<Vec<Gene>> {
    let content = std::fs::read_to_string(bed_path)
        .with_context(|| format!("Failed to read ground truth BED file: {}", bed_path))?;
    
    let mut genes = Vec::new();
    
    for line in content.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 4 {
            continue;
        }
        
        let chrom = parts[0].to_string();
        let start: i64 = parts[1].parse().unwrap_or(0);
        let end: i64 = parts[2].parse().unwrap_or(0);
        let name = parts[3].to_string();
        
        // Check if this gene belongs to the target family
        // Match exact family ID (e.g., "|ID_14\b" to avoid matching ID_141)
        let pattern = format!("|{}\t", family_id);
        let pattern2 = format!("|{}|", family_id);
        if name.contains(&pattern) || name.ends_with(&format!("|{}", family_id)) || 
           name.contains(&pattern2) || name == family_id {
            genes.push(Gene {
                chrom,
                start,
                end,
                name: name.clone(),
                span: end - start,
            });
        }
    }
    
    // Sort by size (largest first for seed selection)
    genes.sort_by(|a, b| b.span.cmp(&a.span));
    
    println!("  Parsed {} genes for {} from ground truth", genes.len(), family_id);
    Ok(genes)
}

/// Calculate overlap between two regions
fn calculate_overlap(a: &Gene, b: &(String, i64, i64)) -> i64 {
    if a.chrom != b.0 {
        return 0;
    }
    let overlap_start = a.start.max(b.1);
    let overlap_end = a.end.min(b.2);
    (overlap_end - overlap_start).max(0)
}

/// Genomic IoU: intersection length over span of the bounding interval (matches batch Python eval).
fn genomic_iou(s1: i64, e1: i64, s2: i64, e2: i64) -> f64 {
    let inter = (e1.min(e2) - s1.max(s2)).max(0);
    if inter <= 0 {
        return 0.0;
    }
    let union = (e1.max(e2) - s1.min(s2)).max(1);
    inter as f64 / union as f64
}

fn iou_to_assignment_cost(iou: f64) -> i64 {
    let c = ((1.0 - iou.clamp(0.0, 1.0)) * 1_000_000.0).round() as i64;
    c.clamp(0, 1_000_000)
}

/// Best seed support for a GT gene: ref coverage, IoU vs that seed, overlap bp, seed interval.
fn best_seed_support_for_gene(gene: &Gene, seeds: &[Seed]) -> (f64, f64, i64, Option<(String, i64, i64)>) {
    let mut best_ref = 0.0_f64;
    let mut best_iou = 0.0_f64;
    let mut best_ob = 0_i64;
    let mut best_reg: Option<(String, i64, i64)> = None;
    for s in seeds {
        if s.chrom != gene.chrom {
            continue;
        }
        let ob = calculate_overlap(gene, &(s.chrom.clone(), s.start, s.end));
        let ref_cov = if gene.span > 0 {
            ob as f64 / gene.span as f64
        } else {
            0.0
        };
        let iou = genomic_iou(gene.start, gene.end, s.start, s.end);
        if ref_cov > best_ref + 1e-9 || (ref_cov >= best_ref - 1e-9 && iou > best_iou) {
            best_ref = ref_cov;
            best_iou = iou;
            best_ob = ob;
            best_reg = Some((s.chrom.clone(), s.start, s.end));
        }
    }
    (best_ref, best_iou, best_ob, best_reg)
}

fn build_novel_ranked(
    detected: &[(String, i64, i64, HashSet<String>)],
    matched_det: &HashSet<usize>,
    novel_support: Option<&NovelLocusSupport>,
) -> Vec<NovelLocusScore> {
    let mut out: Vec<NovelLocusScore> = Vec::new();
    for j in 0..detected.len() {
        if matched_det.contains(&j) {
            continue;
        }
        let (chr, s, e, reads) = &detected[j];
        let nr = novel_support
            .and_then(|ns| ns.n_reads.get(j).copied())
            .unwrap_or_else(|| reads.len());
        let jac = novel_support
            .and_then(|ns| ns.jaccard.get(j).copied())
            .unwrap_or(1.0);
        let dist = novel_support
            .and_then(|ns| ns.distance_from_seed.get(j).copied())
            .unwrap_or(0);
        let rank_score = (nr as f64 + 1.0).ln() * jac.max(1e-6) / (1.0 + dist as f64);
        out.push(NovelLocusScore {
            chrom: chr.clone(),
            start: *s,
            end: *e,
            n_reads: nr,
            jaccard: jac,
            distance_from_seed: dist,
            rank_score,
        });
    }
    out.sort_by(|a, b| {
        b.rank_score
            .partial_cmp(&a.rank_score)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    out
}

fn novel_support_from_transitive(loci: &[TransitiveLocus]) -> NovelLocusSupport {
    NovelLocusSupport {
        n_reads: loci.iter().map(|l| l.reads.len()).collect(),
        jaccard: loci.iter().map(|l| l.jaccard_with_seed).collect(),
        distance_from_seed: loci.iter().map(|l| l.distance_from_seed).collect(),
    }
}

fn novel_support_from_refined(loci: &[RefinedLocus]) -> NovelLocusSupport {
    NovelLocusSupport {
        n_reads: loci.iter().map(|l| l.reads.len()).collect(),
        jaccard: loci.iter().map(|l| l.jaccard_with_seed).collect(),
        distance_from_seed: vec![0; loci.len()],
    }
}

/// Default ground-truth eval: seed confirmation plus bipartite matching to detector output.
fn evaluate_ground_truth(
    args: &Args,
    ground_truth: &[Gene],
    detected: &[(String, i64, i64, HashSet<String>)],
    seeds: &[Seed],
    novel_support: Option<&NovelLocusSupport>,
) -> EvaluationStats {
    if args.eval_legacy_best_interval {
        let mut s = evaluate_detection(ground_truth, detected, args.min_overlap_frac);
        let matched_legacy = legacy_matched_det_indices(ground_truth, detected, args.min_overlap_frac);
        s.novel_candidates_ranked = build_novel_ranked(detected, &matched_legacy, novel_support);
        s
    } else {
        evaluate_detection_bipartite(
            ground_truth,
            detected,
            args.eval_min_iou,
            args.eval_min_ref_frac,
            seeds,
            args.eval_seed_confirm_ref_frac,
            novel_support,
        )
    }
}

fn legacy_matched_det_indices(
    ground_truth: &[Gene],
    detected: &[(String, i64, i64, HashSet<String>)],
    min_overlap_frac: f64,
) -> HashSet<usize> {
    let mut matched = HashSet::default();
    for gene in ground_truth {
        let mut best_overlap = 0i64;
        let mut best_idx = None;
        for (idx, (chrom, start, end, _)) in detected.iter().enumerate() {
            let overlap = calculate_overlap(gene, &(chrom.clone(), *start, *end));
            if overlap > best_overlap {
                best_overlap = overlap;
                best_idx = Some(idx);
            }
        }
        let overlap_pct = if gene.span > 0 {
            (best_overlap as f64 / gene.span as f64) * 100.0
        } else {
            0.0
        };
        if overlap_pct >= (min_overlap_frac * 100.0) {
            if let Some(idx) = best_idx {
                matched.insert(idx);
            }
        }
    }
    matched
}

fn evaluate_detection_bipartite(
    ground_truth: &[Gene],
    detected: &[(String, i64, i64, HashSet<String>)],
    min_iou: f64,
    min_ref_frac: f64,
    seeds: &[Seed],
    seed_confirm_ref_frac: f64,
    novel_support: Option<&NovelLocusSupport>,
) -> EvaluationStats {
    const NO_ASSIGN: i64 = 1_000_000;
    let n_gt = ground_truth.len();
    let n_det = detected.len();
    let total_gt_span: i64 = ground_truth.iter().map(|g| g.span).sum();
    let total_detected_span: i64 = detected.iter().map(|(_, s, e, _)| e - s).sum();

    let novel_support_effective: Option<&NovelLocusSupport> =
        if let Some(ns) = novel_support {
            if ns.n_reads.len() == n_det
                && ns.jaccard.len() == n_det
                && ns.distance_from_seed.len() == n_det
            {
                Some(ns)
            } else {
                eprintln!(
                    "Warning: novel support lengths {} / {} / {} != n_det {}; ignoring support",
                    ns.n_reads.len(),
                    ns.jaccard.len(),
                    ns.distance_from_seed.len(),
                    n_det
                );
                None
            }
        } else {
            None
        };

    if n_gt == 0 {
        let matched = HashSet::default();
        return EvaluationStats {
            total_genes: 0,
            genes_found: 0,
            genes_missed: 0,
            extra_regions: n_det,
            sensitivity: 0.0,
            precision: if n_det == 0 { 100.0 } else { 0.0 },
            avg_overlap_pct: 0.0,
            total_gt_span: 0,
            total_detected_span,
            per_gene: Vec::new(),
            eval_bipartite: true,
            novel_candidates_ranked: build_novel_ranked(detected, &matched, novel_support_effective),
        };
    }

    let seed_disabled = seed_confirm_ref_frac > 1.0;
    let mut seed_confirmed = vec![false; n_gt];
    let mut seed_ref_cov = vec![0.0_f64; n_gt];
    let mut seed_iou = vec![0.0_f64; n_gt];
    let mut seed_ob = vec![0_i64; n_gt];
    let mut seed_reg: Vec<Option<(String, i64, i64)>> = vec![None; n_gt];
    if !seed_disabled {
        for i in 0..n_gt {
            let g = &ground_truth[i];
            let (rc, iou, ob, reg) = best_seed_support_for_gene(g, seeds);
            seed_ref_cov[i] = rc;
            seed_iou[i] = iou;
            seed_ob[i] = ob;
            seed_reg[i] = reg;
            seed_confirmed[i] = rc >= seed_confirm_ref_frac;
        }
    }

    if n_det == 0 {
        let mut per_gene = Vec::with_capacity(n_gt);
        let mut genes_found = 0usize;
        let mut total_overlap_pct = 0.0_f64;
        for i in 0..n_gt {
            let g = &ground_truth[i];
            let sc = seed_confirmed[i];
            let overlap_pct = if g.span > 0 {
                (seed_ob[i] as f64 / g.span as f64) * 100.0
            } else {
                0.0
            };
            if sc {
                genes_found += 1;
                total_overlap_pct += overlap_pct;
            }
            per_gene.push(GeneEvaluation {
                gene_name: g.name.clone(),
                found: sc,
                overlap_bp: seed_ob[i],
                overlap_pct,
                detected_region: seed_reg[i].clone(),
                match_iou: Some(seed_iou[i]),
                confirmed_by_seed: sc,
            });
        }
        let matched = HashSet::default();
        return EvaluationStats {
            total_genes: n_gt,
            genes_found,
            genes_missed: n_gt - genes_found,
            extra_regions: 0,
            sensitivity: (genes_found as f64 / n_gt as f64) * 100.0,
            precision: if genes_found > 0 { 100.0 } else { 0.0 },
            avg_overlap_pct: if genes_found > 0 {
                total_overlap_pct / genes_found as f64
            } else {
                0.0
            },
            total_gt_span,
            total_detected_span: 0,
            per_gene,
            eval_bipartite: true,
            novel_candidates_ranked: build_novel_ranked(detected, &matched, novel_support_effective),
        };
    }

    let need_idx: Vec<usize> = (0..n_gt)
        .filter(|&i| !seed_confirmed[i])
        .collect();
    let n_need = need_idx.len();

    let mut matched_det: HashSet<usize> = HashSet::default();

    let assign_need: Vec<usize> = if n_need == 0 {
        Vec::new()
    } else {
        let n = n_need.max(n_det);
        let mut cost = vec![vec![0_i64; n]; n];
        for rr in 0..n_need {
            let gi = need_idx[rr];
            let g = &ground_truth[gi];
            for j in 0..n_det {
                let (dc, ds, de, _) = &detected[j];
                let iou = if g.chrom == *dc {
                    genomic_iou(g.start, g.end, *ds, *de)
                } else {
                    0.0
                };
                cost[rr][j] = iou_to_assignment_cost(iou);
            }
            for j in n_det..n {
                cost[rr][j] = NO_ASSIGN;
            }
        }
        for i in n_need..n {
            for j in 0..n_det {
                cost[i][j] = 0;
            }
            for j in n_det..n {
                cost[i][j] = 0;
            }
        }
        let matrix = Matrix::from_rows(cost).expect("square cost matrix for bipartite eval");
        let (_tot, assign) = kuhn_munkres_min(&matrix);
        assign
    };

    let mut det_result_for_gene: Vec<Option<(usize, bool)>> = vec![None; n_gt];
    for rr in 0..n_need {
        let gi = need_idx[rr];
        let g = &ground_truth[gi];
        let j = assign_need[rr];
        if j >= n_det {
            det_result_for_gene[gi] = None;
            continue;
        }
        let (dc, ds, de, _) = &detected[j];
        if g.chrom != *dc {
            det_result_for_gene[gi] = Some((j, false));
            continue;
        }
        let iou = genomic_iou(g.start, g.end, *ds, *de);
        let ob = calculate_overlap(g, &(dc.clone(), *ds, *de));
        let ref_cov = if g.span > 0 {
            ob as f64 / g.span as f64
        } else {
            0.0
        };
        let ok = iou >= min_iou || ref_cov >= min_ref_frac;
        if ok {
            matched_det.insert(j);
        }
        det_result_for_gene[gi] = Some((j, ok));
    }

    let mut per_gene = Vec::with_capacity(n_gt);
    let mut genes_found = 0usize;
    let mut total_overlap_pct = 0.0_f64;

    for i in 0..n_gt {
        let g = &ground_truth[i];
        if seed_confirmed[i] {
            let overlap_pct = if g.span > 0 {
                (seed_ob[i] as f64 / g.span as f64) * 100.0
            } else {
                0.0
            };
            genes_found += 1;
            total_overlap_pct += overlap_pct;
            per_gene.push(GeneEvaluation {
                gene_name: g.name.clone(),
                found: true,
                overlap_bp: seed_ob[i],
                overlap_pct,
                detected_region: seed_reg[i].clone(),
                match_iou: Some(seed_iou[i]),
                confirmed_by_seed: true,
            });
            continue;
        }
        match det_result_for_gene[i] {
            Some((j, true)) => {
                let (dc, ds, de, _) = &detected[j];
                let iou = genomic_iou(g.start, g.end, *ds, *de);
                let ob = calculate_overlap(g, &(dc.clone(), *ds, *de));
                let overlap_pct = if g.span > 0 {
                    (ob as f64 / g.span as f64) * 100.0
                } else {
                    0.0
                };
                genes_found += 1;
                total_overlap_pct += overlap_pct;
                per_gene.push(GeneEvaluation {
                    gene_name: g.name.clone(),
                    found: true,
                    overlap_bp: ob,
                    overlap_pct,
                    detected_region: Some((dc.clone(), *ds, *de)),
                    match_iou: Some(iou),
                    confirmed_by_seed: false,
                });
            }
            Some((_, false)) | None => {
                let j_opt = det_result_for_gene[i].map(|(j, _)| j);
                let (overlap_bp, overlap_pct, region, match_iou) = if let Some(j) = j_opt {
                    if j < n_det {
                        let (dc, ds, de, _) = &detected[j];
                        if g.chrom == *dc {
                            let iou = genomic_iou(g.start, g.end, *ds, *de);
                            let ob = calculate_overlap(g, &(dc.clone(), *ds, *de));
                            let pct = if g.span > 0 {
                                (ob as f64 / g.span as f64) * 100.0
                            } else {
                                0.0
                            };
                            (
                                ob,
                                pct,
                                Some((dc.clone(), *ds, *de)),
                                Some(iou),
                            )
                        } else {
                            (0_i64, 0.0, None, None)
                        }
                    } else {
                        (0, 0.0, None, None)
                    }
                } else {
                    (0, 0.0, None, None)
                };
                per_gene.push(GeneEvaluation {
                    gene_name: g.name.clone(),
                    found: false,
                    overlap_bp,
                    overlap_pct,
                    detected_region: region,
                    match_iou,
                    confirmed_by_seed: false,
                });
            }
        }
    }

    let novel = (0..n_det).filter(|j| !matched_det.contains(j)).count();
    let genes_missed = n_gt - genes_found;
    let sensitivity = (genes_found as f64 / n_gt as f64) * 100.0;
    let tp = genes_found as f64;
    let fp = novel as f64;
    let precision = if tp + fp > 0.0 {
        (tp / (tp + fp)) * 100.0
    } else {
        100.0
    };
    let avg_overlap_pct = if genes_found > 0 {
        total_overlap_pct / genes_found as f64
    } else {
        0.0
    };

    let novel_ranked = build_novel_ranked(detected, &matched_det, novel_support_effective);

    EvaluationStats {
        total_genes: n_gt,
        genes_found,
        genes_missed,
        extra_regions: novel,
        sensitivity,
        precision,
        avg_overlap_pct,
        total_gt_span,
        total_detected_span,
        per_gene,
        eval_bipartite: true,
        novel_candidates_ranked: novel_ranked,
    }
}

/// Legacy: each GT gene uses its single best overlapping detected interval and min overlap fraction of span.
fn evaluate_detection(
    ground_truth: &[Gene],
    detected: &[(String, i64, i64, HashSet<String>)],
    min_overlap_frac: f64,
) -> EvaluationStats {
    let mut per_gene = Vec::new();
    let mut genes_found = 0;
    let mut total_overlap_pct = 0.0;
    let mut total_gt_span: i64 = ground_truth.iter().map(|g| g.span).sum();
    
    // Track which detected regions match ground truth
    let mut matched_detected = vec![false; detected.len()];
    
    for gene in ground_truth {
        let mut best_overlap = 0i64;
        let mut best_idx = None;
        
        for (idx, (chrom, start, end, _)) in detected.iter().enumerate() {
            let overlap = calculate_overlap(gene, &(chrom.clone(), *start, *end));
            if overlap > best_overlap {
                best_overlap = overlap;
                best_idx = Some(idx);
            }
        }
        
        let overlap_pct = if gene.span > 0 {
            (best_overlap as f64 / gene.span as f64) * 100.0
        } else {
            0.0
        };
        
        let found = overlap_pct >= (min_overlap_frac * 100.0);
        if found {
            genes_found += 1;
            total_overlap_pct += overlap_pct;
            if let Some(idx) = best_idx {
                matched_detected[idx] = true;
            }
        }
        
        let detected_region = best_idx.map(|idx| {
            (detected[idx].0.clone(), detected[idx].1, detected[idx].2)
        });
        
        per_gene.push(GeneEvaluation {
            gene_name: gene.name.clone(),
            found,
            overlap_bp: best_overlap,
            overlap_pct,
            detected_region,
            match_iou: None,
            confirmed_by_seed: false,
        });
    }
    
    let genes_missed = ground_truth.len() - genes_found;
    let extra_regions = matched_detected.iter().filter(|&&m| !m).count();
    
    let sensitivity = if !ground_truth.is_empty() {
        (genes_found as f64 / ground_truth.len() as f64) * 100.0
    } else {
        0.0
    };
    
    let precision = if !detected.is_empty() {
        let tp = genes_found as f64;
        let fp = extra_regions as f64;
        (tp / (tp + fp)) * 100.0
    } else {
        0.0
    };
    
    let avg_overlap_pct = if genes_found > 0 {
        total_overlap_pct / genes_found as f64
    } else {
        0.0
    };
    
    let total_detected_span: i64 = detected.iter().map(|(_, s, e, _)| e - s).sum();
    
    EvaluationStats {
        total_genes: ground_truth.len(),
        genes_found,
        genes_missed,
        extra_regions,
        sensitivity,
        precision,
        avg_overlap_pct,
        total_gt_span,
        total_detected_span,
        per_gene,
        eval_bipartite: false,
        novel_candidates_ranked: Vec::new(),
    }
}

/// Print evaluation results in a formatted table
fn print_evaluation_stats(stats: &EvaluationStats) {
    print_separator();
    println!("EVALUATION RESULTS");
    print_separator();
    println!();
    
    println!("Summary Statistics:");
    println!("  Total ground truth genes: {}", stats.total_genes);
    if stats.eval_bipartite {
        let n_by_seed = stats.per_gene.iter().filter(|e| e.confirmed_by_seed).count();
        println!(
            "  Genes matched (seed + detector; seed-confirmed: {}): {}",
            n_by_seed, stats.genes_found
        );
        println!("  Genes missed: {}", stats.genes_missed);
        println!(
            "  Novel candidate regions (no qualifying GT match): {}",
            stats.extra_regions
        );
        println!(
            "  (Rules: --eval-seed-confirm-ref-frac, --eval-min-iou, --eval-min-ref-frac; legacy: --eval-legacy-best-interval)"
        );
    } else {
        println!(
            "  Genes found (legacy best-interval vs min-overlap-frac): {}",
            stats.genes_found
        );
        println!("  Genes missed: {}", stats.genes_missed);
        println!("  Extra regions detected: {}", stats.extra_regions);
    }
    println!();
    println!("  Sensitivity (Recall): {:.1}%", stats.sensitivity);
    println!("  Precision: {:.1}%", stats.precision);
    println!("  Average ref overlap (found genes): {:.1}%", stats.avg_overlap_pct);
    println!();
    println!("  Total ground truth span: {} bp", stats.total_gt_span);
    println!("  Total detected span: {} bp", stats.total_detected_span);
    println!();
    
    println!("Per-Gene Results:");
    if stats.eval_bipartite {
        println!(
            "{:<36} {:<8} {:<6} {:<8} {:<12} {:<15}",
            "Gene", "Status", "Seed", "IoU", "Ref ov %", "Overlap (bp)"
        );
        println!("{}", "-".repeat(90));
        for eval in &stats.per_gene {
            let status = if eval.found { "FOUND" } else { "MISSING" };
            let seed_m = if eval.confirmed_by_seed { "yes" } else { "no" };
            let iou_s = eval
                .match_iou
                .map(|x| format!("{:.3}", x))
                .unwrap_or_else(|| "-".to_string());
            println!(
                "{:<36} {:<8} {:<6} {:>8} {:>10.1}% {:>13}",
                eval.gene_name, status, seed_m, iou_s, eval.overlap_pct, eval.overlap_bp
            );
        }
    } else {
        println!(
            "{:<40} {:<10} {:<12} {:<15}",
            "Gene", "Status", "Overlap %", "Overlap (bp)"
        );
        println!("{}", "-".repeat(80));
        for eval in &stats.per_gene {
            let status = if eval.found { "FOUND" } else { "MISSING" };
            println!(
                "{:<40} {:<10} {:>10.1}% {:>13}",
                eval.gene_name, status, eval.overlap_pct, eval.overlap_bp
            );
        }
    }

    if !stats.novel_candidates_ranked.is_empty() {
        println!();
        println!(
            "Top novel candidates (rank_score = ln(reads+1)*jaccard/(1+dist); strongest first):"
        );
        println!(
            "{:<4} {:<28} {:>7} {:>8} {:>5} {:>8} {:>10}",
            "Rank", "Locus", "Reads", "Jaccard", "Dist", "Score", "Span_bp"
        );
        let n_show = stats.novel_candidates_ranked.len().min(25);
        for (i, nv) in stats.novel_candidates_ranked.iter().take(n_show).enumerate() {
            let span = nv.end - nv.start;
            println!(
                "{:<4} {:<28} {:>7} {:>8.4} {:>5} {:>8.3} {:>10}",
                i + 1,
                format!("{}:{}-{}", nv.chrom, nv.start, nv.end),
                nv.n_reads,
                nv.jaccard,
                nv.distance_from_seed,
                nv.rank_score,
                span
            );
        }
        if stats.novel_candidates_ranked.len() > n_show {
            println!(
                "  ... {} more (see companion *_novel_candidates.csv)",
                stats.novel_candidates_ranked.len() - n_show
            );
        }
    }
    
    print_separator();
}

/// Write evaluation results to CSV
fn write_evaluation_csv(stats: &EvaluationStats, output_path: &str) -> Result<()> {
    let mut file = File::create(output_path)?;
    
    // Write summary
    writeln!(file, "metric,value")?;
    writeln!(file, "total_genes,{}", stats.total_genes)?;
    writeln!(file, "genes_found,{}", stats.genes_found)?;
    writeln!(file, "genes_missed,{}", stats.genes_missed)?;
    writeln!(file, "extra_regions,{}", stats.extra_regions)?;
    writeln!(file, "sensitivity,{:.2}", stats.sensitivity)?;
    writeln!(file, "precision,{:.2}", stats.precision)?;
    writeln!(file, "avg_overlap_pct,{:.2}", stats.avg_overlap_pct)?;
    writeln!(file, "total_gt_span,{}", stats.total_gt_span)?;
    writeln!(file, "total_detected_span,{}", stats.total_detected_span)?;
    writeln!(
        file,
        "eval_bipartite,{}",
        if stats.eval_bipartite { "1" } else { "0" }
    )?;
    writeln!(file)?;
    
    // Write per-gene results
    writeln!(
        file,
        "gene_name,found,confirmed_by_seed,overlap_bp,overlap_pct,match_iou"
    )?;
    for eval in &stats.per_gene {
        let iou_cell = eval
            .match_iou
            .map(|x| format!("{:.4}", x))
            .unwrap_or_default();
        writeln!(
            file,
            "{},{},{},{},{:.2},{}",
            eval.gene_name,
            if eval.found { "1" } else { "0" },
            if eval.confirmed_by_seed { "1" } else { "0" },
            eval.overlap_bp,
            eval.overlap_pct,
            iou_cell
        )?;
    }
    
    println!("Evaluation results written to: {}", output_path);
    write_novel_candidates_csv(stats, output_path)?;
    Ok(())
}

fn write_novel_candidates_csv(stats: &EvaluationStats, eval_csv_path: &str) -> Result<()> {
    if stats.novel_candidates_ranked.is_empty() {
        return Ok(());
    }
    let novel_path = if let Some(base) = eval_csv_path.strip_suffix("_evaluation.csv") {
        format!("{}_novel_candidates.csv", base)
    } else if let Some(base) = eval_csv_path.strip_suffix(".csv") {
        format!("{}_novel_candidates.csv", base)
    } else {
        format!("{}.novel_candidates.csv", eval_csv_path)
    };
    let mut f = File::create(&novel_path)?;
    writeln!(
        f,
        "rank,chrom,start,end,span_bp,n_reads,jaccard,distance_from_seed,rank_score"
    )?;
    for (i, nv) in stats.novel_candidates_ranked.iter().enumerate() {
        writeln!(
            f,
            "{},{},{},{},{},{},{:.6},{},{:.6}",
            i + 1,
            nv.chrom,
            nv.start,
            nv.end,
            nv.end - nv.start,
            nv.n_reads,
            nv.jaccard,
            nv.distance_from_seed,
            nv.rank_score
        )?;
    }
    println!("Novel candidate table written to: {}", novel_path);
    Ok(())
}

/// Convert Gene to Seed
fn gene_to_seed(gene: &Gene) -> Seed {
    Seed {
        chrom: gene.chrom.clone(),
        start: gene.start,
        end: gene.end,
        name: Some(gene.name.clone()),
    }
}

/// Run minimum seed selection algorithm
/// Iteratively adds seeds until all ground truth genes are found
fn run_minimum_seed_selection(
    bam_path: &str,
    ground_truth: &[Gene],
    target_chroms: &[String],
    args: &Args,
) -> Result<(Vec<Seed>, Vec<(String, i64, i64, HashSet<String>)>)> {
    println!("Running minimum seed selection algorithm...");
    println!("  Max iterations: {}", args.max_seed_iterations);
    println!();
    
    let mut seeds: Vec<Seed> = Vec::new();
    let mut all_detected: Vec<(String, i64, i64, HashSet<String>)> = Vec::new();
    
    for iteration in 1..=args.max_seed_iterations {
        // If no seeds yet, start with the largest gene
        if seeds.is_empty() {
            if let Some(first_gene) = ground_truth.first() {
                seeds.push(gene_to_seed(first_gene));
                println!("Iteration {}: Starting with largest gene {}",
                    iteration, first_gene.name);
            } else {
                anyhow::bail!("No ground truth genes provided");
            }
        } else {
            // Find missing genes and add the largest as new seed
            let detected_set: HashSet<String> = all_detected.iter()
                .flat_map(|(_, _, _, reads)| reads.iter().cloned())
                .collect();
            
            let missing: Vec<&Gene> = ground_truth.iter()
                .filter(|g| {
                    // Check if this gene is covered by any detected region
                    !all_detected.iter().any(|(chrom, s, e, _)| {
                        g.chrom == *chrom && g.start < *e && g.end > *s
                    })
                })
                .collect();
            
            if missing.is_empty() {
                println!("\n✓ All genes found with {} seed(s)!", seeds.len());
                break;
            }
            
            // Add largest missing gene as new seed
            let new_seed_gene = missing.iter().max_by_key(|g| g.span).unwrap();
            seeds.push(gene_to_seed(new_seed_gene));
            println!("Iteration {}: Added seed {} ({} missing genes remaining)",
                iteration, new_seed_gene.name, missing.len());
        }
        
        println!("  Running with {} seed(s)...", seeds.len());
        
        // Run detection with current seeds
        // ... (reuse existing detection logic)
        // For now, return the seeds - actual detection happens in main
        
        if iteration == args.max_seed_iterations {
            println!("\n⚠ Reached max iterations ({}) with {} seed(s)",
                args.max_seed_iterations, seeds.len());
        }
    }
    
    Ok((seeds, all_detected))
}

fn run_transitive_detection_collect(
    args: &Args,
    seeds: &[Seed],
    target_chroms: &[String],
    refine_boundaries: bool,
    max_trim_frac: f64,
    min_tx_density: f64,
    min_segment_jaccard: f64,
) -> Result<(Vec<TransitiveLocus>, Option<Vec<RefinedLocus>>)> {
    if seeds.is_empty() {
        anyhow::bail!("No seeds provided to transitive detection");
    }

    println!("Running transitive detection for {} seed(s)...", seeds.len());

    let params = TransDetParams {
        cluster_distance: args.cluster_distance,
        adaptive_cluster_distance: args.adaptive_cluster_distance,
        min_reads_per_locus: args.min_reads,
        max_iterations: args.max_iterations,
        min_jaccard_for_edge: args.min_jaccard,
    };

    let mut all_loci: Vec<TransitiveLocus> = Vec::new();

    for (seed_idx, seed) in seeds.iter().enumerate() {
        if seeds.len() > 1 {
            print_separator();
            println!("Adaptive seed {}/{}: {}:{}-{}", seed_idx + 1, seeds.len(), seed.chrom, seed.start, seed.end);
            print_separator();
        }

        let loci = detect_transitive(
            &args.bam,
            &seed.chrom,
            seed.start,
            seed.end,
            target_chroms,
            &params,
        )?;

        for locus in loci {
            // Keep the LARGER locus when overlapping
            let locus_span = locus.end - locus.start;
            let mut should_add = true;
            let mut to_remove: Option<usize> = None;

            for (idx, l) in all_loci.iter().enumerate() {
                if l.chrom == locus.chrom {
                    let overlap_start = l.start.max(locus.start);
                    let overlap_end = l.end.min(locus.end);
                    if overlap_start < overlap_end {
                        let existing_span = l.end - l.start;
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

    if !refine_boundaries {
        return Ok((all_loci, None));
    }

    // Boundary refinement uses the seed reads for Jaccard trimming.
    // We follow the existing behavior by using the first seed in this batch.
    let first_seed = &seeds[0];
    let seed_reads = collect_seed_reads_fast(
        &args.bam,
        &first_seed.chrom,
        first_seed.start as usize,
        first_seed.end as usize,
        args.include_supplementary,
    )?;

    let refine_params = RefinementParams {
        max_trim_fraction: max_trim_frac,
        min_transcript_density: min_tx_density,
        keep_all_segments: !args.single_best_segment,
        min_segment_jaccard,
        iso_boundary_mode: args.iso_boundary_mode,
        ..Default::default()
    };

    let loci_for_refinement: Vec<(String, i64, i64, HashSet<String>)> = all_loci
        .iter()
        .map(|l| (l.chrom.clone(), l.start, l.end, l.reads.clone()))
        .collect();

    let refined = refine_loci(
        &args.bam,
        &loci_for_refinement,
        &seed_reads,
        &refine_params,
    )?;

    Ok((all_loci, Some(refined)))
}

fn main() -> Result<()> {
    let args = Args::parse();

    print_separator();
    println!("GENE FAMILY DETECTOR (Iso-Seq Aware)");
    print_separator();
    println!();
    if args.integrated {
        println!("Strategy: INTEGRATED MODE");
        println!("  Multi-mapping reads + Coverage valleys + Iso-Seq validation");
        println!("  Valley fraction: {:.2}", args.valley_frac);
        println!("  Peak threshold frac: {:.2}", args.peak_threshold_frac);
        println!("  Min peak gap: {} bp", args.valley_min_gap_bp);
        println!("  Min segment: {} bp", args.valley_min_segment_bp);
        println!("  Min prominence: {:.2}", args.valley_min_prominence);
        println!("  Auto-split threshold: {:.2}", args.auto_split_threshold);
        println!("  Validation threshold: {:.2}", args.validate_threshold);
        println!(
            "  Strand-adaptive split: {}",
            if args.no_integrated_strand_split {
                "off"
            } else {
                "on"
            }
        );
        println!(
            "  Strand-adaptive min span: {} bp",
            args.strand_adaptive_min_span_bp
        );
        println!(
            "  Junction (CIGAR exon) refine: {}",
            if args.integrated_junction_refine {
                "on"
            } else {
                "off"
            }
        );
        if args.integrated_junction_refine {
            println!(
                "  Junction refine max trim/side: {:.2} frac, min span {} bp",
                args.junction_refine_max_trim_frac,
                args.junction_refine_min_span
            );
            println!(
                "  Junction RNA cues: acceptor min reads {} (slack {} bp), read-end slack {} bp",
                args.junction_refine_min_junction_reads,
                args.junction_refine_cluster_slack_bp,
                args.junction_refine_read_end_slack_bp
            );
        }
    } else if args.use_read_clouds {
        println!("Strategy: Read cloud analysis (dense overlapping multi-mapping reads)");
        println!("  Min cloud reads: {}", args.min_cloud_reads);
        println!("  Min density: {:.1} reads/kb", args.min_density);
        if args.strand_aware_clouds {
            println!("  Strand-aware: YES (+/- separated)");
        }
        println!("  Improvements: soft-clip extension, full-length weighting");
        if args.require_secondary {
            println!("  Require secondary: YES");
        }
    } else {
        println!("Strategy: Multi-mapping reads + Iso-Seq transcript validation");
    }
    println!("Threads: {}", args.threads);
    if args.use_isoseq {
        println!("Iso-Seq validation: ENABLED (min similarity: {})", args.min_isoseq_sim);
    }
    if args.require_read_overlap {
        println!("Read-space overlap: REQUIRED (min: {} bp)", args.min_read_overlap_bp);
    }
    println!();

    let target_chroms: Vec<String> = if let Some(ref chroms) = args.chromosomes {
        chroms.split(',').map(|x| x.to_string()).collect()
    } else {
        let mut reader = bam::io::indexed_reader::Builder::default()
            .build_from_path(&args.bam)?;
        let header = reader.read_header()?;
        header
            .reference_sequences()
            .iter()
            .filter_map(|(name, _)| std::str::from_utf8(name).ok().map(|s| s.to_string()))
            .collect()
    };

    // Collect all seeds
    let seeds = collect_seeds(&args)?;
    let is_multi_seed = seeds.len() > 1;
    
    if is_multi_seed {
        println!("=== MULTI-SEED MODE ===");
        println!("Processing {} seeds", seeds.len());
        for (i, seed) in seeds.iter().enumerate() {
            let _name_str = seed.name.as_ref().map(|n| format!(" ({})", n)).unwrap_or_default();
            println!("  Seed {}: {}:{}-{}", i + 1, seed.chrom, seed.start, seed.end);
        }
    } else {
        println!("=== SINGLE SEED MODE ===");
        println!("Seed: {}:{}-{} ({})", 
                 seeds[0].chrom, seeds[0].start, seeds[0].end,
                 seeds[0].name.as_ref().unwrap_or(&"unnamed".to_string()));
    }

    let start_time = std::time::Instant::now();

    if args.min_std_coverage > 0.0 {
        println!(
            "Coverage std filtering: ENABLED (baseline min std >= {:.1}, dynamic: {})",
            args.min_std_coverage,
            !args.no_dynamic_std_coverage
        );
    } else {
        println!("Coverage std filtering: DISABLED");
    }

    if args.connected_components_mode || args.adaptive_discovery_v2 {
        if args.no_transitive {
            anyhow::bail!(
                "--connected-components-mode/--adaptive-discovery-v2 is incompatible with --no-transitive"
            );
        }
        if args.integrated {
            anyhow::bail!("--connected-components-mode/--adaptive-discovery-v2 is not supported with --integrated");
        }
        if args.connected_components_mode && args.adaptive_discovery_v2 {
            anyhow::bail!("Choose only one of --connected-components-mode or --adaptive-discovery-v2");
        }

        let strict_v2 = args.adaptive_discovery_v2;

        println!();
        if strict_v2 {
            println!("=== ADAPTIVE DISCOVERY V2 MODE ===");
            println!("Non-circular reseeding with strict guardrails");
        } else {
            println!("=== CONNECTED COMPONENTS RESEED MODE ===");
            println!("Non-circular alternate mode using disconnected graph components");
        }
        println!("Max seed iterations: {}", args.max_seed_iterations);
        println!("Max new seeds per iteration: {}", args.cc_max_new_seeds_per_iter);
        println!("Min component reads: {}", args.cc_min_component_reads);
        if strict_v2 {
            println!("Min exploratory Jaccard: {:.4}", args.adv2_min_explore_jaccard);
            println!("Max exploratory seed span: {}", args.adv2_max_seed_span);
            println!("Min shared fraction: {:.3}", args.adv2_min_shared_frac);
            println!("Min exploratory seed score: {:.3}", args.adv2_min_seed_score);
            println!("Min exploratory novelty distance: {}", args.adv2_min_novel_distance_bp);
            println!("Skeleton mode: {}", args.skeleton_mode);
            if args.skeleton_mode {
                println!("Skeleton min score: {:.3}", args.adv2_skeleton_min_score);
                println!("Skeleton top-k eval: {}", args.adv2_skeleton_top_k);
            }
            println!("Hotspot source: {}", args.adv2_hotspot_source);
            if args.adv2_hotspot_source {
                println!("Hotspot bin bp: {}", args.adv2_hotspot_bin_bp);
                println!("Hotspot min bin reads: {}", args.adv2_hotspot_min_bin_reads);
                println!("Hotspot min shared frac: {:.3}", args.adv2_hotspot_min_shared_frac);
                println!("Hotspot min support seeds: {}", args.adv2_hotspot_min_support_seeds);
                println!("Hotspot support gap bp: {}", args.adv2_hotspot_support_gap_bp);
                println!("Hotspot top-k: {}", args.adv2_hotspot_top_k);
            }
            println!("Max loci growth factor: {:.2}", args.adv2_max_loci_growth_factor);
        }
        println!();

        let params = TransDetParams {
            cluster_distance: args.cluster_distance,
            adaptive_cluster_distance: args.adaptive_cluster_distance,
            min_reads_per_locus: args.min_reads,
            max_iterations: args.max_iterations,
            min_jaccard_for_edge: args.min_jaccard,
        };

        let mut seeds_current: Vec<Seed> = vec![seeds[0].clone()];
        let mut last_all_loci: Vec<TransitiveLocus> = Vec::new();
        let mut prev_loci_count: usize = 0;
        let mut last_added_chrom: Option<String> = None;

        for iter in 1..=args.max_seed_iterations {
            println!("--- CC iteration {}/{} ---", iter, args.max_seed_iterations);
            println!("Seeds in this round: {}", seeds_current.len());
            print_separator();

            let mut merged_loci: Vec<TransitiveLocus> = Vec::new();
            let mut all_candidates: Vec<ComponentSeedCandidate> = Vec::new();
            let mut per_seed_loci: Vec<Vec<TransitiveLocus>> = Vec::new();

            for seed in &seeds_current {
                let (loci, candidates) = detect_transitive_with_component_candidates(
                    &args.bam,
                    &seed.chrom,
                    seed.start,
                    seed.end,
                    &target_chroms,
                    &params,
                )?;
                per_seed_loci.push(loci.clone());

                for locus in loci {
                    let is_duplicate = merged_loci.iter().any(|l| {
                        l.chrom == locus.chrom
                            && ((l.start - locus.start).abs() < 5000 || (l.end - locus.end).abs() < 5000)
                    });
                    if !is_duplicate {
                        merged_loci.push(locus);
                    }
                }

                for c in candidates {
                    if c.total_reads_in_component >= args.cc_min_component_reads {
                        all_candidates.push(c);
                    }
                }
            }

            last_all_loci = merged_loci;
            println!("Discovered loci so far: {}", last_all_loci.len());

            if strict_v2 && prev_loci_count > 0 {
                let growth = (last_all_loci.len() as f64) / (prev_loci_count as f64);
                if growth > args.adv2_max_loci_growth_factor {
                    println!(
                        "Early stop: loci growth factor {:.2} exceeds {:.2}",
                        growth, args.adv2_max_loci_growth_factor
                    );
                    break;
                }
            }
            prev_loci_count = last_all_loci.len();

            if iter == args.max_seed_iterations {
                break;
            }

            // Deduplicate and rank component candidates
            all_candidates.sort_by(|a, b| {
                b.total_reads_in_component
                    .cmp(&a.total_reads_in_component)
                    .then(b.component_size.cmp(&a.component_size))
                    .then(b.reads.len().cmp(&a.reads.len()))
            });

            let mut new_seeds: Vec<Seed> = Vec::new();
            for c in all_candidates {
                let near_existing_seed = seeds_current.iter().any(|s| {
                    s.chrom == c.chrom
                        && (s.start - c.start).abs() < 20000
                        && (s.end - c.end).abs() < 20000
                });
                if near_existing_seed {
                    continue;
                }

                let near_discovered = last_all_loci.iter().any(|l| {
                    l.chrom == c.chrom && l.start <= c.end && l.end >= c.start
                });
                if near_discovered {
                    continue;
                }
                if strict_v2 {
                    let span = c.end - c.start;
                    if span > args.adv2_max_seed_span || span < 1000 {
                        continue;
                    }
                    if let Some(ref last_chrom) = last_added_chrom {
                        if &c.chrom == last_chrom {
                            continue;
                        }
                    }
                }

                new_seeds.push(Seed {
                    chrom: c.chrom.clone(),
                    start: c.start,
                    end: c.end,
                    name: Some(format!("cc_seed_{}_{}", c.component_size, c.total_reads_in_component)),
                });

                if new_seeds.len() >= args.cc_max_new_seeds_per_iter {
                    break;
                }
            }

            if new_seeds.is_empty() {
                println!("No eligible disconnected components. Mining exploratory unlabeled seeds...");
                let skeleton_model = if strict_v2 && args.skeleton_mode {
                    build_skeleton_model(&args, &last_all_loci)?
                } else {
                    None
                };
                new_seeds = propose_exploratory_seeds_from_discovered(
                    &args,
                    &target_chroms,
                    &last_all_loci,
                    &seeds_current,
                    strict_v2,
                    skeleton_model.as_ref(),
                )?;
                if new_seeds.is_empty() && strict_v2 && args.adv2_hotspot_source {
                    println!("No exploratory seeds passed. Trying transcript-density hotspots...");
                    new_seeds = propose_hotspot_seeds_from_transcript_density(
                        &args,
                        &target_chroms,
                        &last_all_loci,
                        &seeds_current,
                        &per_seed_loci,
                    )?;
                }
            }

            if new_seeds.is_empty() {
                println!("No additional exploratory seeds found. Stopping.");
                break;
            }

            println!("Adding {} new component-derived seed(s)", new_seeds.len());
            for s in new_seeds {
                println!("  + {}:{}-{}", s.chrom, s.start, s.end);
                if strict_v2 {
                    last_added_chrom = Some(s.chrom.clone());
                }
                seeds_current.push(s);
            }
        }

        let mut file = File::create(&args.output)?;
        writeln!(file, "#chrom\tstart\tend\tjaccard\tdistance\tcomponent\treads")?;
        let mut sorted_loci = last_all_loci.clone();
        sorted_loci.sort_by(|a, b| a.chrom.cmp(&b.chrom).then(a.start.cmp(&b.start)));

        for locus in &sorted_loci {
            writeln!(
                file,
                "{}\t{}\t{}\t{:.4}\t{}\t{}\t{}",
                locus.chrom,
                locus.start,
                locus.end,
                locus.jaccard_with_seed,
                locus.distance_from_seed,
                locus.component,
                locus.reads.len()
            )?;
        }

        // Optional evaluation if user provides ground truth (for benchmarking only)
        if let (Some(ref gt_path), Some(ref family_id)) = (&args.ground_truth, &args.family_id) {
            match parse_ground_truth(gt_path, family_id) {
                Ok(ground_truth) if !ground_truth.is_empty() => {
                    let detected: Vec<(String, i64, i64, HashSet<String>)> = last_all_loci
                        .iter()
                        .map(|l| (l.chrom.clone(), l.start, l.end, HashSet::default()))
                        .collect();
                    let novel_sup = novel_support_from_transitive(&last_all_loci);
                    let stats = evaluate_ground_truth(
                        &args,
                        &ground_truth,
                        &detected,
                        &seeds_current,
                        Some(&novel_sup),
                    );
                    print_evaluation_stats(&stats);
                    let csv_path = args.output.replace(".bed", "_evaluation.csv");
                    if let Err(e) = write_evaluation_csv(&stats, &csv_path) {
                        eprintln!("Warning: Failed to write evaluation CSV: {}", e);
                    }
                }
                Ok(_) => println!("Warning: No ground truth genes found for {}", family_id),
                Err(e) => eprintln!("Warning: Failed to parse ground truth: {}", e),
            }
        }

        print_separator();
        if strict_v2 {
            println!("SUMMARY (Adaptive Discovery V2 Mode)");
        } else {
            println!("SUMMARY (Connected Components Reseed Mode)");
        }
        print_separator();
        println!("Final seeds used: {}", seeds_current.len());
        println!("Final loci found: {}", last_all_loci.len());
        println!("Total time: {:?}", start_time.elapsed());
        println!("Output BED: {}", args.output);
        return Ok(());
    }

    // === CONNECTIVITY-BASED SEED SELECTION MODE ===
    // Uses Jaccard similarity graph connected components to select minimum seeds
    // This avoids circular reasoning: we analyze read-sharing patterns to determine
    // which genes need separate seeds, then run transitive detection with those seeds
    if args.connectivity_seed {
        if args.no_transitive {
            anyhow::bail!("--connectivity-seed is incompatible with --no-transitive");
        }
        if args.integrated {
            anyhow::bail!("--connectivity-seed is not supported with --integrated");
        }
        // Require seeds BED or list to provide initial gene regions
        if args.seeds_bed.is_none() && args.seeds_list.is_none() {
            anyhow::bail!("--connectivity-seed requires --seeds-bed or --seeds-list to provide initial gene regions");
        }
        
        // Parse all gene regions from BED/list (ignore the single seed arg)
        let mut gene_regions: Vec<GeneRegion> = Vec::new();
        
        if let Some(ref bed_path) = args.seeds_bed {
            let seeds_from_bed = parse_seeds_bed(bed_path)?;
            for seed in seeds_from_bed {
                gene_regions.push(GeneRegion {
                    chrom: seed.chrom,
                    start: seed.start,
                    end: seed.end,
                    name: seed.name,
                });
            }
        }
        
        if let Some(ref list) = args.seeds_list {
            let seeds_from_list = parse_seeds_list(list)?;
            for seed in seeds_from_list {
                gene_regions.push(GeneRegion {
                    chrom: seed.chrom,
                    start: seed.start,
                    end: seed.end,
                    name: seed.name,
                });
            }
        }
        
        if gene_regions.is_empty() {
            anyhow::bail!("No gene regions parsed from seeds BED/list");
        }
        
        // Also include the CLI seed if provided and different from parsed seeds
        if let (Some(chrom), Some(start), Some(end)) = (&args.seed_chrom, args.seed_start, args.seed_end) {
            let cli_seed = Seed {
                chrom: chrom.clone(),
                start: start as i64,
                end: end as i64,
                name: Some("cli_seed".to_string()),
            };
            
            let cli_is_new = !gene_regions.iter().any(|g| {
                g.chrom == cli_seed.chrom && 
                (g.start - cli_seed.start).abs() < 100 && 
                (g.end - cli_seed.end).abs() < 100
            });
            
            if cli_is_new {
                gene_regions.push(GeneRegion {
                    chrom: cli_seed.chrom.clone(),
                    start: cli_seed.start,
                    end: cli_seed.end,
                    name: cli_seed.name,
                });
            }
        }
        
        println!();
        println!("=== CONNECTIVITY-BASED SEED SELECTION ===");
        println!("Total gene regions to analyze: {}", gene_regions.len());
        println!("Jaccard floor (--min-jaccard): {}", args.min_jaccard);
        if args.adaptive_connectivity {
            println!(
                "Adaptive graph: quantile {:.3}, min shared reads per edge {}",
                args.connectivity_adaptive_quantile, args.connectivity_min_shared_reads
            );
        }
        if let Some(p) = args.connectivity_hypergeom_p_max {
            println!("Hypergeometric edge filter: p_max={:.6}", p);
        }
        println!();

        // Run connectivity analysis to find connected components
        let gene_names: Vec<String> = gene_regions.iter()
            .enumerate()
            .map(|(i, g)| {
                g.name.clone().unwrap_or_else(|| format!("gene_{}", i + 1))
            })
            .collect();

        let conn_params = ConnectivityParams {
            min_jaccard_floor: args.min_jaccard,
            adaptive_quantile: if args.adaptive_connectivity {
                Some(args.connectivity_adaptive_quantile)
            } else {
                None
            },
            min_shared_reads: args.connectivity_min_shared_reads,
            hypergeom_max_pvalue: args.connectivity_hypergeom_p_max,
        };

        let connectivity = analyze_connectivity(&args.bam, &gene_regions, &conn_params)?;
        
        print_analysis_summary(&connectivity, &gene_names);
        
        let mut seeds_current: Vec<Seed>;
        
        if connectivity.seed_indices.is_empty() {
            println!("Warning: No seeds selected (all regions may be isolated)");
            println!("Falling back to CLI seed only");
            seeds_current = vec![seeds[0].clone()];
        } else {
            // Get selected seeds from connectivity analysis
            let selected_seed_regions = get_selected_seeds(&connectivity, &gene_regions);
            
            // Convert to Seed objects
            seeds_current = selected_seed_regions.into_iter()
                .map(|g| Seed {
                    chrom: g.chrom,
                    start: g.start,
                    end: g.end,
                    name: g.name,
                })
                .collect();
        }
        
        println!();
        println!("Selected {} seeds from connectivity analysis", seeds_current.len());
        for (i, seed) in seeds_current.iter().enumerate() {
            let name = seed.name.as_deref().unwrap_or("unknown");
            println!("  Seed {}: {}:{}-{} ({})", i + 1, seed.chrom, seed.start, seed.end, name);
        }
        println!();
        
        // Exit early if --seeds-only is set (for batch evaluation)
        if args.seeds_only {
            println!("Exiting after seed selection (--seeds-only mode)");
            return Ok(());
        }
        
        // Run transitive detection with selected seeds
        
        // Run transitive detection with selected seeds
        let params = TransDetParams {
            cluster_distance: args.cluster_distance,
            adaptive_cluster_distance: args.adaptive_cluster_distance,
            min_reads_per_locus: args.min_reads,
            max_iterations: args.max_iterations,
            min_jaccard_for_edge: args.min_jaccard,
        };
        
        let mut all_loci: Vec<TransitiveLocus> = Vec::new();
        
        for (seed_idx, seed) in seeds_current.iter().enumerate() {
            println!("Running detection for seed {}: {}:{}-{}", 
                     seed_idx + 1, seed.chrom, seed.start, seed.end);
            
            let (loci, _) = detect_transitive_with_component_candidates(
                &args.bam,
                &seed.chrom,
                seed.start,
                seed.end,
                &target_chroms,
                &params,
            )?;
            
            for locus in loci {
                // Keep the LARGER locus when overlapping
                let locus_span = locus.end - locus.start;
                let mut should_add = true;
                let mut to_remove: Option<usize> = None;

                for (idx, l) in all_loci.iter().enumerate() {
                    if l.chrom == locus.chrom {
                        let overlap_start = l.start.max(locus.start);
                        let overlap_end = l.end.min(locus.end);
                        if overlap_start < overlap_end {
                            let existing_span = l.end - l.start;
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
        
        println!();
        println!("Total loci discovered: {}", all_loci.len());
        
        // Write output
        let mut file = File::create(&args.output)?;
        writeln!(file, "#chrom\tstart\tend\tjaccard\tdistance\tcomponent\treads")?;
        let mut sorted_loci = all_loci.clone();
        sorted_loci.sort_by(|a, b| a.chrom.cmp(&b.chrom).then(a.start.cmp(&b.start)));
        
        for locus in &sorted_loci {
            writeln!(
                file,
                "{}\t{}\t{}\t{:.4}\t{}\t{}\t{}",
                locus.chrom,
                locus.start,
                locus.end,
                locus.jaccard_with_seed,
                locus.distance_from_seed,
                locus.component,
                locus.reads.len()
            )?;
        }
        
        // Optional evaluation
        if let (Some(ref gt_path), Some(ref family_id)) = (&args.ground_truth, &args.family_id) {
            match parse_ground_truth(gt_path, family_id) {
                Ok(ground_truth) if !ground_truth.is_empty() => {
                    let detected: Vec<(String, i64, i64, HashSet<String>)> = all_loci
                        .iter()
                        .map(|l| (l.chrom.clone(), l.start, l.end, HashSet::default()))
                        .collect();
                    let novel_sup = novel_support_from_transitive(&all_loci);
                    let stats = evaluate_ground_truth(
                        &args,
                        &ground_truth,
                        &detected,
                        &seeds_current,
                        Some(&novel_sup),
                    );
                    print_evaluation_stats(&stats);
                }
                Ok(_) => println!("Warning: No ground truth genes found for {}", family_id),
                Err(e) => eprintln!("Warning: Failed to parse ground truth: {}", e),
            }
        }
        
        println!();
        println!("SUMMARY (Connectivity-Based Seed Selection)");
        println!("Seeds used: {}", seeds_current.len());
        println!("Loci found: {}", all_loci.len());
        println!("Output: {}", args.output);
        return Ok(());
    }

    // Adaptive minimum seed selection:
    // - After each detection run, identify missed genes with zero overlap (disconnected components)
    // - Add seeds for those disconnected genes, without lowering Jaccard thresholds
    // - For partial-overlap misses, try boundary refinement adjustments (no further Jaccard relaxation)
    if args.min_seed {
        if args.no_transitive {
            anyhow::bail!("--min-seed is incompatible with --no-transitive");
        }
        let (gt_path, family_id) = match (&args.ground_truth, &args.family_id) {
            (Some(gt_path), Some(family_id)) => (gt_path.clone(), family_id.clone()),
            _ => anyhow::bail!("--min-seed requires --ground-truth and --family-id"),
        };

        if args.integrated {
            anyhow::bail!("--min-seed is not supported with --integrated");
        }

        let ground_truth = parse_ground_truth(&gt_path, &family_id)?;
        if ground_truth.is_empty() {
            anyhow::bail!("No ground truth genes parsed for {}", family_id);
        }

        println!();
        println!("=== ADAPTIVE MINIMUM SEED SELECTION ===");
        println!("Family: {}", family_id);
        println!("Target genes: {}", ground_truth.len());
        println!("Max seed iterations: {}", args.max_seed_iterations);
        println!("Min overlap fraction: {}", args.min_overlap_frac);
        println!();

        // Start from the required CLI seed only (ignore any additional seeds provided via seeds-bed/list).
        // The goal is to learn the missing disconnected components.
        let mut seeds_current: Vec<Seed> = vec![seeds[0].clone()];

        let mut refine_boundaries_current = args.refine_boundaries;
        let mut max_trim_frac_current = args.max_trim_frac;
        let mut min_tx_density_current = args.min_tx_density;
        let mut min_segment_jaccard_current = args.min_segment_jaccard;

        let mut last_all_loci: Vec<TransitiveLocus> = Vec::new();
        let mut last_refined: Option<Vec<RefinedLocus>> = None;
        let mut last_stats: Option<EvaluationStats> = None;

        for iter in 1..=args.max_seed_iterations {
            println!("--- Adaptive iteration {}/{} ---", iter, args.max_seed_iterations);
            println!("Seeds in this round: {}", seeds_current.len());
            println!(
                "Refinement: {} (max_trim_frac={:.3}, min_tx_density={:.3}, min_segment_jaccard={:.3})",
                refine_boundaries_current, max_trim_frac_current, min_tx_density_current, min_segment_jaccard_current
            );
            print_separator();

            let (all_loci, refined_opt) = run_transitive_detection_collect(
                &args,
                &seeds_current,
                &target_chroms,
                refine_boundaries_current,
                max_trim_frac_current,
                min_tx_density_current,
                min_segment_jaccard_current,
            )?;

            // Choose what to evaluate:
            // - If refinement ran, evaluate refined coordinates
            // - Otherwise, evaluate raw transitive loci
            let detected_for_eval: Vec<(String, i64, i64, HashSet<String>)> = if let Some(refined) = refined_opt.as_ref() {
                refined
                    .iter()
                    .map(|r| (r.chrom.clone(), r.start, r.end, HashSet::default()))
                    .collect()
            } else {
                all_loci
                    .iter()
                    .map(|l| (l.chrom.clone(), l.start, l.end, HashSet::default()))
                    .collect()
            };

            let novel_sup = if let Some(r) = refined_opt.as_ref() {
                novel_support_from_refined(r)
            } else {
                novel_support_from_transitive(&all_loci)
            };
            let stats = evaluate_ground_truth(
                &args,
                &ground_truth,
                &detected_for_eval,
                &seeds_current,
                Some(&novel_sup),
            );
            print_evaluation_stats(&stats);

            last_all_loci = all_loci;
            last_refined = refined_opt;
            last_stats = Some(stats.clone());

            if stats.genes_found == ground_truth.len() {
                println!("Adaptive selection complete: all genes found.");
                break;
            }

            let gene_by_name: HashMap<String, Gene> = ground_truth
                .iter()
                .map(|g| (g.name.clone(), g.clone()))
                .collect();

            let mut disconnected_missing: Vec<Gene> = Vec::new();
            let mut partial_overlap_missing: Vec<Gene> = Vec::new();

            for per_gene in &stats.per_gene {
                if per_gene.found {
                    continue;
                }
                let gene = match gene_by_name.get(&per_gene.gene_name) {
                    Some(g) => g,
                    None => continue,
                };
                if per_gene.overlap_bp == 0 {
                    disconnected_missing.push(gene.clone());
                } else {
                    partial_overlap_missing.push(gene.clone());
                }
            }

            println!(
                "Missing genes: {} (disconnected={}, partial-overlap={})",
                stats.genes_missed,
                disconnected_missing.len(),
                partial_overlap_missing.len()
            );

            if iter == args.max_seed_iterations {
                break;
            }

            let disconnected_empty = disconnected_missing.is_empty();
            let mut added_seed = false;
            for gene in &disconnected_missing {
                let new_seed = gene_to_seed(gene);
                let is_duplicate = seeds_current.iter().any(|s| {
                    s.chrom == new_seed.chrom
                        && (s.start - new_seed.start).abs() < 100
                        && (s.end - new_seed.end).abs() < 100
                });
                if !is_duplicate {
                    seeds_current.push(new_seed);
                    added_seed = true;
                }
            }

            // Boundary/overextension case:
            // enable refinement (if not already) and make trimming slightly less aggressive.
            // This avoids changing Jaccard thresholds and only tunes boundary trimming.
            let mut refined_adjusted = false;
            if !partial_overlap_missing.is_empty() {
                if !refine_boundaries_current {
                    refine_boundaries_current = true;
                    refined_adjusted = true;
                }
                // Conservative adjustment: reduce allowed trimming so detected loci keep more signal.
                let new_max_trim = (max_trim_frac_current - 0.05).max(0.30);
                if new_max_trim < max_trim_frac_current {
                    max_trim_frac_current = new_max_trim;
                    refined_adjusted = true;
                }
                let new_min_seg = (min_segment_jaccard_current * 0.7).max(0.001);
                if new_min_seg < min_segment_jaccard_current {
                    min_segment_jaccard_current = new_min_seg;
                    refined_adjusted = true;
                }
            }

            if disconnected_empty && !refined_adjusted && !added_seed {
                println!("No further progress possible under current rules; stopping adaptive selection.");
                break;
            }
        }

        // Write final outputs for the last evaluated set.
        // Output BED always uses the raw transitive loci coordinates (consistent with existing behavior).
        // If refinement was enabled in the last round, also write the refined BED.
        let mut file = File::create(&args.output)?;
        writeln!(file, "#chrom\tstart\tend\tjaccard\tdistance\tcomponent\treads")?;

        let mut sorted_loci = last_all_loci.clone();
        sorted_loci.sort_by(|a, b| a.chrom.cmp(&b.chrom).then(a.start.cmp(&b.start)));

        for locus in &sorted_loci {
            writeln!(
                file,
                "{}\t{}\t{}\t{:.4}\t{}\t{}\t{}",
                locus.chrom,
                locus.start,
                locus.end,
                locus.jaccard_with_seed,
                locus.distance_from_seed,
                locus.component,
                locus.reads.len()
            )?;
        }

        if refine_boundaries_current {
            if let Some(refined) = last_refined.as_ref() {
                let first_seed = &seeds_current[0];
                let seed_reads = collect_seed_reads_fast(
                    &args.bam,
                    &first_seed.chrom,
                    first_seed.start as usize,
                    first_seed.end as usize,
                    args.include_supplementary,
                )?;

                let refined_output = args.output.replace(".bed", "_refined.bed");
                write_refined_bed(&refined_output, refined, &seed_reads)?;
                println!("Refined output BED: {}", refined_output);
            }
        }

        if let Some(stats) = last_stats.as_ref() {
            let csv_path = args.output.replace(".bed", "_evaluation.csv");
            if let Err(e) = write_evaluation_csv(stats, &csv_path) {
                eprintln!("Warning: Failed to write evaluation CSV: {}", e);
            }
        }

        print_separator();
        println!("SUMMARY (Adaptive Minimum Seed Selection)");
        print_separator();
        println!("Final seeds used: {}", seeds_current.len());
        println!("Final loci found: {}", last_all_loci.len());
        println!("Total time: {:?}", start_time.elapsed());
        println!("Output BED: {}", args.output);
        print_separator();

        return Ok(());
    }

    // === INTEGRATED MODE WITH MULTI-SEED SUPPORT ===
    if args.integrated {
        println!("\n=== INTEGRATED DETECTION MODE ===");
        println!("Combining: multi-mapping reads + coverage valleys + Iso-Seq validation");
        println!();

        let integrated_params = IntegratedParams {
            valley_params: ValleyDetectionParams {
                bin_size: 100,
                valley_frac: args.valley_frac,
                peak_threshold_frac: args.peak_threshold_frac,
                min_gap_bp: args.valley_min_gap_bp,
                min_segment_bp: args.valley_min_segment_bp,
                min_prominence_frac: args.valley_min_prominence,
            },
            auto_split_threshold: args.auto_split_threshold,
            validate_threshold: args.validate_threshold,
            min_read_overlap_bp: args.min_read_overlap_bp,
            require_isoseq_for_ambiguous: true,
            strand_adaptive_split: !args.no_integrated_strand_split,
            strand_adaptive_min_span_bp: args.strand_adaptive_min_span_bp,
            junction_refinement: if args.integrated_junction_refine {
                Some(JunctionRefinementParams {
                    max_trim_frac_per_side: args.junction_refine_max_trim_frac,
                    min_segment_span: args.junction_refine_min_span,
                    include_supplementary: args.include_supplementary,
                    min_junction_read_support: args.junction_refine_min_junction_reads,
                    junction_cluster_slack_bp: args.junction_refine_cluster_slack_bp,
                    read_end_guard_slack_bp: args.junction_refine_read_end_slack_bp,
                    acceptor_window_left_bp: args.junction_refine_acceptor_win_left_bp,
                    acceptor_window_right_bp: args.junction_refine_acceptor_win_right_bp,
                })
            } else {
                None
            },
        };

        let seed_regions: Vec<(String, i64, i64)> = seeds
            .iter()
            .map(|s| (s.chrom.clone(), s.start, s.end))
            .collect();
        let mut union_seed_reads: HashSet<String> = HashSet::default();

        // Process each seed
        let mut all_cores_per_seed: Vec<Vec<GeneCore>> = vec![];
        
        for (seed_idx, seed) in seeds.iter().enumerate() {
            if is_multi_seed {
                print_separator();
                println!("Processing seed {}/{}: {}:{}-{}",
                         seed_idx + 1, seeds.len(), seed.chrom, seed.start, seed.end);
                print_separator();
                println!();
            }
            
            // Collect seed reads
            let seed_reads = collect_seed_reads_fast(
                &args.bam,
                &seed.chrom,
                seed.start as usize,
                seed.end as usize,
                args.include_supplementary,
            )?;
            
            println!("  Seed reads: {}", seed_reads.len());
            union_seed_reads.extend(seed_reads.iter().cloned());
            
            if seed_reads.is_empty() {
                println!("  Warning: No reads found for this seed, skipping...");
                continue;
            }
            
            // Detect gene family for this seed
            let cores = detect_gene_family_integrated(
                &args.bam,
                &seed.chrom,
                seed.start,
                seed.end,
                &seed_reads,
                &target_chroms,
                &integrated_params,
                args.cluster_distance,
                args.min_reads,
            )?;
            
            all_cores_per_seed.push(cores);
        }
        
        // Merge cores from all seeds
        let mut merged_cores = merge_cores_from_multiple_seeds(all_cores_per_seed);

        if args.min_std_coverage > 0.0 {
            let before_ic = merged_cores.len();
            let mut filtered_cores: Vec<GeneCore> = Vec::new();
            for c in merged_cores {
                if locus_overlaps_seed_region(&c.chrom, c.start, c.end, &seed_regions) {
                    filtered_cores.push(c);
                    continue;
                }
                let shared = c.reads.intersection(&union_seed_reads).count();
                match min_std_coverage_reject_reason(
                    &args,
                    &args.bam,
                    &c.chrom,
                    c.start,
                    c.end,
                    shared,
                )? {
                    None => filtered_cores.push(c),
                    Some(reason) => println!(
                        "  REJECTED integrated core {}:{}-{}: {}",
                        c.chrom, c.start, c.end, reason
                    ),
                }
            }
            println!(
                "  Coverage std filter (integrated): {} -> {} cores",
                before_ic,
                filtered_cores.len()
            );
            merged_cores = filtered_cores;
        }

        // Write integrated results to BED
        let mut file = File::create(&args.output)?;
        writeln!(
            file,
            "#chrom\tstart\tend\tjaccard\tconfidence\tvalleys\tsplits\tisoseq_validated"
        )?;

        for core in &merged_cores {
            writeln!(
                file,
                "{}\t{}\t{}\t{:.4}\t{}\t{}\t{}\t{}",
                core.chrom,
                core.start,
                core.end,
                core.jaccard_with_seed,
                core.confidence,
                core.coverage_valleys,
                core.significant_splits,
                core.isoseq_validated
            )?;
        }

        print_separator();
        if is_multi_seed {
            println!("SUMMARY (Integrated Multi-Seed Mode)");
            print_separator();
            println!("Seeds processed: {}", seeds.len());
        } else {
            println!("SUMMARY (Integrated Mode)");
            print_separator();
        }
        println!("Gene cores found: {}", merged_cores.len());
        
        let high_conf = merged_cores.iter().filter(|c| c.confidence == CoreConfidence::High).count();
        let med_conf = merged_cores.iter().filter(|c| c.confidence == CoreConfidence::Medium).count();
        let low_conf = merged_cores.iter().filter(|c| c.confidence == CoreConfidence::Low).count();
        
        println!("  High confidence: {}", high_conf);
        println!("  Medium confidence: {}", med_conf);
        println!("  Low confidence: {}", low_conf);
        println!("Total time: {:?}", start_time.elapsed());
        println!("Output BED: {}", args.output);
        print_separator();

        return Ok(());
    }

    // === TRANSITIVE MODE ===
    if !args.no_transitive {
        println!("\n=== TRANSITIVE DETECTION MODE ===");
        println!("For highly diverged gene families (e.g., ID_14 LRRC37)");
        println!("Using optimized multi-seed processing (single genome scan)");
        println!();
        if args.adaptive_cluster_distance {
            println!(
                "Adaptive cluster merge gap: ON (per-chrom same-read gap 90p, cap {} bp)",
                args.cluster_distance
            );
        }

        let params = TransDetParams {
            cluster_distance: args.cluster_distance,
            adaptive_cluster_distance: args.adaptive_cluster_distance,
            min_reads_per_locus: args.min_reads,
            max_iterations: args.max_iterations,
            min_jaccard_for_edge: args.min_jaccard,
        };

        // Convert seeds to the format expected by detect_transitive_multi_seed
        let seed_tuples: Vec<(String, i64, i64)> = seeds.iter()
            .map(|s| (s.chrom.clone(), s.start, s.end))
            .collect();
        
        // Extract target chromosomes from seeds (family members are on these chromosomes)
        let target_chroms: Vec<String> = seeds.iter()
            .map(|s| s.chrom.clone())
            .collect::<std::collections::HashSet<_>>()
            .into_iter()
            .collect();
        
        // Use the optimized multi-seed function that scans only target chromosomes
        let detection_start = std::time::Instant::now();
        let (mut all_loci, _, all_seed_reads) = detect_transitive_multi_seed_and_reads(
            &args.bam,
            &seed_tuples,
            &target_chroms,
            &params,
        )?;
        println!("  Multi-seed detection: {:?}", detection_start.elapsed());
        
        // === POST-HOC CLASSIFICATION: Identify single-copy vs multi-copy ===
        // Check which seeds only found themselves (single-copy) vs found paralogs (multi-copy)
        println!("\n=== COPY NUMBER CLASSIFICATION ===");
        let mut single_copy_seeds = Vec::new();
        let mut multi_copy_count = 0;
        
        for (idx, seed) in seeds.iter().enumerate() {
            // Count how many loci overlap this seed
            let overlapping_loci: Vec<_> = all_loci.iter()
                .filter(|l| {
                    l.chrom == seed.chrom &&
                    l.start < seed.end &&
                    l.end > seed.start
                })
                .collect();
            
            // Check if this seed connects to other loci via shared reads
            let seed_reads = collect_seed_reads_fast(
                &args.bam,
                &seed.chrom,
                seed.start as usize,
                seed.end as usize,
                args.include_supplementary
            ).unwrap_or_default();
            
            let connected_loci: Vec<_> = all_loci.iter()
                .filter(|l| {
                    if l.chrom == seed.chrom && l.start == seed.start && l.end == seed.end {
                        return false; // Skip self
                    }
                    // Check if locus shares reads with seed
                    l.reads.iter().any(|r| seed_reads.contains(r))
                })
                .collect();
            
            if connected_loci.is_empty() {
                single_copy_seeds.push(idx);
                println!("  Seed {} ({}:{}-{}): SINGLE-COPY (no paralogs found)",
                         seed.name.as_deref().unwrap_or("?"),
                         seed.chrom, seed.start, seed.end);
            } else {
                multi_copy_count += 1;
                println!("  Seed {} ({}:{}-{}): MULTI-COPY ({} connected loci)",
                         seed.name.as_deref().unwrap_or("?"),
                         seed.chrom, seed.start, seed.end,
                         connected_loci.len());
            }
        }
        
        println!("  Summary: {} multi-copy families, {} single-copy genes",
                 multi_copy_count, single_copy_seeds.len());
        
        // Compute core reads per connected component of seeds
        // This handles families with multiple sub-clusters that share different reads
        let core_read_components = compute_core_reads_by_component(&args.bam, &seeds, args.include_supplementary)?;
        let total_core_reads: usize = core_read_components.iter().map(|c| c.len()).sum();
        println!("  Core read components: {} (total {} core reads)", 
                 core_read_components.len(), total_core_reads);
        
        // Multi-seed support filter: require loci to share reads with at least N seeds
        // This helps filter artifacts in disjoint families where seeds don't share reads directly
        if args.min_seeds > 1 && seeds.len() > 1 {
            let filter_start = std::time::Instant::now();
            let before_filter = all_loci.len();
            
            // Collect per-seed read sets
            let mut seed_read_sets: Vec<HashSet<String>> = Vec::new();
            for seed in &seeds {
                let reads = collect_seed_reads_fast(
                    &args.bam, 
                    &seed.chrom, 
                    seed.start as usize, 
                    seed.end as usize, 
                    args.include_supplementary
                ).unwrap_or_default();
                seed_read_sets.push(reads);
            }
            
            // Count how many seeds each locus shares reads with
            let mut kept = Vec::new();
            let mut rejected = 0;
            
            for locus in all_loci {
                // Skip seed regions themselves
                if locus_overlaps_seed_region(&locus.chrom, locus.start, locus.end, &seed_tuples) {
                    kept.push(locus);
                    continue;
                }
                
                let seed_count = seed_read_sets.iter()
                    .filter(|seed_reads| {
                        // Check if locus shares at least one read with this seed
                        locus.reads.iter().any(|r| seed_reads.contains(r))
                    })
                    .count();
                
                if seed_count >= args.min_seeds {
                    kept.push(locus);
                } else {
                    rejected += 1;
                    if rejected <= 3 {
                        println!("  REJECTED (insufficient seeds) {}:{}-{}: shares reads with {} seeds (min: {})",
                                 locus.chrom, locus.start, locus.end, seed_count, args.min_seeds);
                    }
                }
            }
            
            if rejected > 3 {
                println!("  ... and {} more rejected (insufficient seeds)", rejected - 3);
            }
            
            println!("  Multi-seed filter ({}+ seeds): {} -> {} loci ({} rejected) in {:?}",
                     args.min_seeds, before_filter, kept.len(), rejected, filter_start.elapsed());
            all_loci = kept;
        }
        
        // Always split loci by coverage valleys (mandatory)
        println!("\n=== SPLITTING LOCI BY COVERAGE ===");
        let split_start = std::time::Instant::now();
        let fixed_valley = ValleyDetectionParams {
            bin_size: args.coverage_bin_size,
            valley_frac: args.valley_frac,
            peak_threshold_frac: args.peak_threshold_frac,
            min_gap_bp: args.valley_min_gap_bp,
            min_segment_bp: args.valley_min_segment_bp,
            min_prominence_frac: args.valley_min_prominence,
        };
        let split_cfg = TransitiveCoverageSplitConfig {
            use_dynamic_valleys: !args.fixed_transitive_valleys,
            fixed_valley,
            profile_bin_size: args.coverage_bin_size,
            valley_merge_min_segment_bp: args.transitive_valley_merge_min_segment_bp,
            valley_merge_max_gap_bp: args.transitive_valley_merge_max_gap_bp,
            junction_min_locus_bp: args.transitive_junction_min_locus_bp,
            junction_min_split_span_bp: args.transitive_junction_min_split_span_bp,
            junction_merge_min_segment_bp: args.transitive_junction_merge_min_segment_bp,
            junction_merge_max_gap_bp: args.transitive_junction_merge_max_gap_bp,
        };
        if args.fixed_transitive_valleys {
            println!(
                "  Transitive valley mode: FIXED (valley_frac={}, peak_frac={}, min_gap={}, min_seg={}, prominence={})",
                args.valley_frac,
                args.peak_threshold_frac,
                args.valley_min_gap_bp,
                args.valley_min_segment_bp,
                args.valley_min_prominence
            );
        } else {
            println!("  Transitive valley mode: DYNAMIC (per-locus; profile bin size {} bp)", args.coverage_bin_size);
        }
        println!(
            "  Valley merge: min_seg={} bp max_gap={} bp | Junction: locus_min={} split_min={} merge_min_seg={} merge_max_gap={}",
            split_cfg.valley_merge_min_segment_bp,
            split_cfg.valley_merge_max_gap_bp,
            split_cfg.junction_min_locus_bp,
            split_cfg.junction_min_split_span_bp,
            split_cfg.junction_merge_min_segment_bp,
            split_cfg.junction_merge_max_gap_bp
        );
        // Use parallel version for multi-threaded runs
        // Pass seed_tuples to ensure seed regions are preserved during splitting
        let mut all_loci = if args.threads > 1 {
            println!("  Using parallel coverage splitting ({} threads)", args.threads);
            match split_loci_by_coverage_par(&args.bam, &all_loci, &all_seed_reads, &split_cfg, &seed_tuples) {
                Ok(split) => {
                    println!("  Split {} loci into {} segments in {:?}", all_loci.len(), split.len(), split_start.elapsed());
                    split
                }
                Err(e) => {
                    eprintln!("Warning: Parallel coverage splitting failed: {}", e);
                    all_loci
                }
            }
        } else {
            match split_loci_by_coverage(&args.bam, &all_loci, &all_seed_reads, &split_cfg, &seed_tuples) {
                Ok(split) => {
                    println!("  Split {} loci into {} segments in {:?}", all_loci.len(), split.len(), split_start.elapsed());
                    split
                }
                Err(e) => {
                    eprintln!("Warning: Coverage splitting failed: {}", e);
                    all_loci
                }
            }
        };

        if args.min_std_coverage > 0.0 {
            let filter_start = std::time::Instant::now();
            let before_tc = all_loci.len();
            
            // Parallel std coverage filtering for significant speedup
            let filter_results: Vec<_> = if args.threads > 1 {
                println!("  Using parallel std coverage filtering ({} threads)", args.threads);
                all_loci
                    .into_par_iter()
                    .map(|l| {
                        if locus_overlaps_seed_region(&l.chrom, l.start, l.end, &seed_tuples) {
                            return (l, None); // Keep seed-overlapping loci
                        }
                        let shared = l.reads.intersection(&all_seed_reads).count();
                        match min_std_coverage_reject_reason(
                            &args,
                            &args.bam,
                            &l.chrom,
                            l.start,
                            l.end,
                            shared,
                        ) {
                            Ok(None) => (l, None),
                            Ok(Some(reason)) => (l, Some(reason)),
                            Err(_) => (l, None), // Keep on error
                        }
                    })
                    .collect()
            } else {
                all_loci
                    .into_iter()
                    .map(|l| {
                        if locus_overlaps_seed_region(&l.chrom, l.start, l.end, &seed_tuples) {
                            return (l, None);
                        }
                        let shared = l.reads.intersection(&all_seed_reads).count();
                        match min_std_coverage_reject_reason(
                            &args,
                            &args.bam,
                            &l.chrom,
                            l.start,
                            l.end,
                            shared,
                        ) {
                            Ok(None) => (l, None),
                            Ok(Some(reason)) => {
                                println!(
                                    "  REJECTED (transitive) {}:{}-{}: {}",
                                    l.chrom, l.start, l.end, reason
                                );
                                (l, Some(reason))
                            }
                            Err(_) => (l, None),
                        }
                    })
                    .collect()
            };
            
            // Collect kept loci and rejected count
            let mut kept = Vec::new();
            let mut rejected_count = 0;
            for (locus, reason) in filter_results {
                if reason.is_none() {
                    kept.push(locus);
                } else {
                    rejected_count += 1;
                }
            }
            
            println!(
                "  Coverage std filter (transitive): {} -> {} loci ({} rejected) in {:?}",
                before_tc,
                kept.len(),
                rejected_count,
                filter_start.elapsed()
            );
            all_loci = kept;
        }
        
        // Post-splitting multi-seed filter: recompute reads for each segment and filter
        // This catches artifacts created during coverage/junction splitting
        if args.min_seeds > 1 && seeds.len() > 1 && !all_loci.is_empty() {
            let filter_start = std::time::Instant::now();
            let before_filter = all_loci.len();
            
            // Collect per-seed read sets (if not already done)
            let seed_read_sets: Vec<HashSet<String>> = seeds.iter()
                .map(|seed| {
                    collect_seed_reads_fast(
                        &args.bam, 
                        &seed.chrom, 
                        seed.start as usize, 
                        seed.end as usize, 
                        args.include_supplementary
                    ).unwrap_or_default()
                })
                .collect();
            
            // Re-compute reads for each locus based on actual overlap
            // (Split segments may have inherited parent's reads)
            let mut kept = Vec::new();
            let mut rejected = 0;
            
            for mut locus in all_loci {
                // Skip seed regions themselves
                if locus_overlaps_seed_region(&locus.chrom, locus.start, locus.end, &seed_tuples) {
                    kept.push(locus);
                    continue;
                }
                
                // Re-compute reads for this specific segment
                let segment_reads = collect_seed_reads_fast(
                    &args.bam,
                    &locus.chrom,
                    locus.start as usize,
                    locus.end as usize,
                    args.include_supplementary
                ).unwrap_or_default();
                
                // Update locus with accurate reads
                locus.reads = segment_reads;
                
                // Count how many seeds this segment shares reads with
                let seed_count = seed_read_sets.iter()
                    .filter(|seed_reads| {
                        locus.reads.iter().any(|r| seed_reads.contains(r))
                    })
                    .count();
                
                if seed_count >= args.min_seeds {
                    kept.push(locus);
                } else {
                    rejected += 1;
                    if rejected <= 3 {
                        println!("  REJECTED (post-split, insufficient seeds) {}:{}-{}: shares reads with {} seeds (min: {})",
                                 locus.chrom, locus.start, locus.end, seed_count, args.min_seeds);
                    }
                }
            }
            
            if rejected > 3 {
                println!("  ... and {} more rejected (post-split insufficient seeds)", rejected - 3);
            }
            
            println!("  Post-split multi-seed filter ({}+ seeds): {} -> {} loci ({} rejected) in {:?}",
                     args.min_seeds, before_filter, kept.len(), rejected, filter_start.elapsed());
            all_loci = kept;
        }
        
        // Core read filter: remove candidates that don't share core multi-mapping reads
        // A locus is valid if it shares at least min_core_reads with ANY component's core set
        if args.min_core_reads > 0 && !core_read_components.is_empty() {
            let before_filter = all_loci.len();
            let mut kept = Vec::new();
            let mut rejected = 0;
            
            for locus in all_loci {
                // Skip seed regions themselves
                if locus_overlaps_seed_region(&locus.chrom, locus.start, locus.end, &seed_tuples) {
                    kept.push(locus);
                    continue;
                }
                
                // Check if locus shares reads with ANY component's core set
                let max_core_count = core_read_components.iter()
                    .map(|component| locus.reads.intersection(component).count())
                    .max()
                    .unwrap_or(0);
                
                if max_core_count >= args.min_core_reads {
                    kept.push(locus);
                } else {
                    rejected += 1;
                    // Only print a sample of rejections to avoid spam
                    if rejected <= 5 {
                        println!(
                            "  REJECTED (no core reads) {}:{}-{}: max {} core reads across {} components (min: {})",
                            locus.chrom, locus.start, locus.end, max_core_count, 
                            core_read_components.len(), args.min_core_reads
                        );
                    }
                }
            }
            
            if rejected > 5 {
                println!("  ... and {} more rejected (no core reads)", rejected - 5);
            }
            
            println!(
                "  Core read filter: {} -> {} loci ({} rejected, min {} core reads, {} components)",
                before_filter, kept.len(), rejected, args.min_core_reads, core_read_components.len()
            );
            all_loci = kept;
        }
        
        // Jaccard filter for transitive path: remove low-Jaccard false positives
        // This catches artifacts that share many reads but have low similarity to seeds
        if args.min_jaccard_locus > 0.0 {
            let before_filter = all_loci.len();
            let mut kept = Vec::new();
            let mut rejected = 0;
            
            for locus in all_loci {
                // Skip seed regions themselves
                if locus_overlaps_seed_region(&locus.chrom, locus.start, locus.end, &seed_tuples) {
                    kept.push(locus);
                    continue;
                }
                
                if locus.jaccard_with_seed >= args.min_jaccard_locus {
                    kept.push(locus);
                } else {
                    rejected += 1;
                    if rejected <= 3 {
                        println!("  REJECTED (low Jaccard) {}:{}-{}: jaccard={:.4} (min: {:.4})",
                                 locus.chrom, locus.start, locus.end, 
                                 locus.jaccard_with_seed, args.min_jaccard_locus);
                    }
                }
            }
            
            if rejected > 3 {
                println!("  ... and {} more rejected (low Jaccard)", rejected - 3);
            }
            
            println!("  Jaccard filter (transitive, ≥{:.3}): {} -> {} loci ({} rejected)",
                     args.min_jaccard_locus, before_filter, kept.len(), rejected);
            all_loci = kept;
        }
        
        // Coverage quality filter: remove artifacts with non-gene-like coverage profiles
        // Note: Disabled by default (min_coverage_quality=0.0) to preserve pseudogene detection
        if args.min_coverage_quality > 0.0 {
            let filter_start = std::time::Instant::now();
            let before_filter = all_loci.len();
            
            let filter_results: Vec<_> = if args.threads > 1 {
                println!("  Using parallel coverage quality filtering ({} threads)", args.threads);
                all_loci
                    .into_par_iter()
                    .map(|l| {
                        if locus_overlaps_seed_region(&l.chrom, l.start, l.end, &seed_tuples) {
                            return (l, None); // Keep seed-overlapping loci
                        }
                        match coverage_quality_reject_reason(&args, &args.bam, &l.chrom, l.start, l.end) {
                            Ok(None) => (l, None),
                            Ok(Some(reason)) => (l, Some(reason)),
                            Err(_) => (l, None), // Keep on error
                        }
                    })
                    .collect()
            } else {
                all_loci
                    .into_iter()
                    .map(|l| {
                        if locus_overlaps_seed_region(&l.chrom, l.start, l.end, &seed_tuples) {
                            return (l, None);
                        }
                        match coverage_quality_reject_reason(&args, &args.bam, &l.chrom, l.start, l.end) {
                            Ok(None) => (l, None),
                            Ok(Some(reason)) => {
                                println!("  REJECTED (coverage quality) {}:{}-{}: {}",
                                         l.chrom, l.start, l.end, reason);
                                (l, Some(reason))
                            }
                            Err(_) => (l, None),
                        }
                    })
                    .collect()
            };
            
            let mut kept = Vec::new();
            let mut rejected_count = 0;
            for (locus, reason) in filter_results {
                if reason.is_none() {
                    kept.push(locus);
                } else {
                    rejected_count += 1;
                }
            }
            
            println!(
                "  Coverage quality filter (transitive): {} -> {} loci ({} rejected, min: {:.2}) in {:?}",
                before_filter, kept.len(), rejected_count, args.min_coverage_quality, filter_start.elapsed()
            );
            all_loci = kept;
        }

        if args.trim_read_envelope {
            println!("\n=== READ-ENVELOPE TRIM (primary alignments) ===");
            all_loci = trim_loci_read_align_envelope(
                &args.bam,
                all_loci,
                args.read_envelope_pct_lo,
                args.read_envelope_pct_hi,
                args.read_envelope_min_span,
                args.read_envelope_min_alignments,
                args.read_envelope_query_pad_bp,
            )?;
        }

        if args.tight_coverage_resplit {
            println!("\n=== TIGHT COVERAGE RESPLIT (optional, BAM-only) ===");
            let tp = TightSplitParams {
                profile_bin_size: args.tight_coverage_bin_size,
                valley: ValleyDetectionParams {
                    bin_size: args.tight_coverage_bin_size,
                    valley_frac: args.tight_valley_frac,
                    peak_threshold_frac: args.tight_peak_threshold_frac,
                    min_gap_bp: args.tight_min_gap_bp,
                    min_segment_bp: args.tight_min_segment_bp,
                    min_prominence_frac: args.tight_min_prominence_frac,
                },
                merge_min_segment_bp: args.tight_merge_min_segment_bp,
                merge_max_gap_bp: args.tight_merge_max_gap_bp,
                min_input_span_bp: args.tight_min_input_span_bp,
                max_segments_per_locus: args.tight_max_segments_per_locus,
            };
            println!(
                "  tight valley_frac={} peak_frac={} min_gap={} min_seg={} prominence={} merge_min_seg={} merge_max_gap={}",
                tp.valley.valley_frac,
                tp.valley.peak_threshold_frac,
                tp.valley.min_gap_bp,
                tp.valley.min_segment_bp,
                tp.valley.min_prominence_frac,
                tp.merge_min_segment_bp,
                tp.merge_max_gap_bp
            );
            all_loci = tight_resplit_loci(&args.bam, all_loci, &all_seed_reads, &tp)?;
        }

        if args.clip_tail_trim {
            println!("\n=== CLIP-TAIL TRIM (A-rich 3' soft clips, optional) ===");
            let ctp = ClipTailTrimParams {
                min_clip_len: args.clip_tail_min_len,
                min_a_frac: args.clip_tail_min_a_frac,
                min_qualifying_reads: args.clip_tail_min_reads,
                end_slack_bp: args.clip_tail_end_slack_bp,
                start_slack_bp: args.clip_tail_start_slack_bp,
                strand_dom_frac: args.clip_tail_strand_dom_frac,
                query_pad_bp: args.clip_tail_query_pad_bp,
                min_result_span: args.clip_tail_min_result_span,
                verbose: args.clip_tail_verbose,
            };
            all_loci = trim_loci_clip_tail(&args.bam, all_loci, &ctp)?;
        }
        
        let mut refined_for_eval: Option<Vec<RefinedLocus>> = None;
        
        // Print summary
        println!("\nProcessing complete. Total loci found: {}", all_loci.len());

        // Write results to BED
        let mut file = File::create(&args.output)?;
        writeln!(
            file,
            "#chrom\tstart\tend\tjaccard\tdistance\tcomponent\treads"
        )?;

        // Sort by chrom, start
        let mut sorted_loci = all_loci.clone();
        sorted_loci.sort_by(|a, b| {
            a.chrom.cmp(&b.chrom)
                .then(a.start.cmp(&b.start))
        });

        // Refine boundaries if requested
        if args.refine_boundaries {
            println!("\n=== REFINING BOUNDARIES ===");
            println!("Using Iso-Seq transcripts and Jaccard trimming...");
            
            // Collect seed reads for Jaccard calculation
            let first_seed = &seeds[0];
            let seed_reads = collect_seed_reads_fast(
                &args.bam,
                &first_seed.chrom,
                first_seed.start as usize,
                first_seed.end as usize,
                args.include_supplementary,
            )?;
            
            let refine_params = RefinementParams {
                max_trim_fraction: args.max_trim_frac,
                min_transcript_density: args.min_tx_density,
                keep_all_segments: !args.single_best_segment,  // Invert: default is keep-all
                min_segment_jaccard: args.min_segment_jaccard,
                iso_boundary_mode: args.iso_boundary_mode,
                ..Default::default()
            };
            
            // Convert TransitiveLocus to format for refinement
            let loci_for_refinement: Vec<(String, i64, i64, HashSet<String>)> = sorted_loci
                .iter()
                .map(|l| (l.chrom.clone(), l.start, l.end, l.reads.clone()))
                .collect();
            
            let refined = refine_loci(
                &args.bam,
                &loci_for_refinement,
                &seed_reads,
                &refine_params,
            )?;
            
            // Write refined results
            let refined_output = args.output.replace(".bed", "_refined.bed");
            write_refined_bed(&refined_output, &refined, &seed_reads)?;
            refined_for_eval = Some(refined);
            
            println!("  Refined output: {}", refined_output);
        }

        for locus in &sorted_loci {
            writeln!(
                file,
                "{}\t{}\t{}\t{:.4}\t{}\t{}\t{}",
                locus.chrom,
                locus.start,
                locus.end,
                locus.jaccard_with_seed,
                locus.distance_from_seed,
                locus.component,
                locus.reads.len()
            )?;
        }

        print_separator();
        if is_multi_seed {
            println!("SUMMARY (Transitive Multi-Seed Mode)");
            print_separator();
            println!("Seeds processed: {}", seeds.len());
        } else {
            println!("SUMMARY (Transitive Mode)");
            print_separator();
        }
        println!("Total loci found: {}", all_loci.len());
        let direct = all_loci.iter().filter(|l| l.distance_from_seed == 1).count();
        let transitive = all_loci.iter().filter(|l| l.distance_from_seed >= 2).count();
        println!("  Direct neighbors: {}", direct);
        println!("  Transitive: {}", transitive);
        println!("Total time: {:?}", start_time.elapsed());
        println!("Output BED: {}", args.output);
        print_separator();

        // Evaluate against ground truth if provided
        if let (Some(ref gt_path), Some(ref family_id)) = (&args.ground_truth, &args.family_id) {
            println!();
            println!("Evaluating against ground truth...");
            
            match parse_ground_truth(gt_path, family_id) {
                Ok(ground_truth) => {
                    if !ground_truth.is_empty() {
                        // Convert loci to evaluation format
                        // If boundary refinement ran, evaluate the refined coordinates.
                        let detected: Vec<(String, i64, i64, HashSet<String>)> = if let Some(ref refined) = refined_for_eval.as_ref() {
                            refined
                                .iter()
                                .map(|r| (r.chrom.clone(), r.start, r.end, HashSet::default()))
                                .collect()
                        } else {
                            all_loci
                                .iter()
                                .map(|l| (l.chrom.clone(), l.start, l.end, HashSet::default()))
                                .collect()
                        };

                        let novel_sup = if let Some(ref r) = refined_for_eval.as_ref() {
                            novel_support_from_refined(r)
                        } else {
                            novel_support_from_transitive(&all_loci)
                        };
                        let stats = evaluate_ground_truth(
                            &args,
                            &ground_truth,
                            &detected,
                            &seeds,
                            Some(&novel_sup),
                        );
                        print_evaluation_stats(&stats);
                        
                        // Write CSV
                        let csv_path = args.output.replace(".bed", "_evaluation.csv");
                        if let Err(e) = write_evaluation_csv(&stats, &csv_path) {
                            eprintln!("Warning: Failed to write evaluation CSV: {}", e);
                        }
                    } else {
                        println!("Warning: No ground truth genes found for {}", family_id);
                    }
                }
                Err(e) => {
                    eprintln!("Warning: Failed to parse ground truth: {}", e);
                }
            }
        }

        return Ok(());
    }

    // Non-integrated mode: Use first seed for backward compatibility
    // (Multi-seed is only fully supported in integrated mode)
    if is_multi_seed {
        println!("WARNING: Multi-seed mode is only fully supported with --integrated");
        println!("Using first seed only for non-integrated mode: {}:{}-{}",
                 seeds[0].chrom, seeds[0].start, seeds[0].end);
        println!();
    }
    
    let first_seed = &seeds[0];
    
    // Analyze seed with Iso-Seq (non-integrated mode)
    let seed_analysis = if args.use_isoseq {
        println!("  Analyzing seed transcript structure...");
        Some(analyze_locus_transcripts(
            &args.bam,
            &first_seed.chrom,
            first_seed.start,
            first_seed.end,
        )?)
    } else {
        None
    };

    if let Some(ref analysis) = seed_analysis {
        println!("  Seed transcripts: {} total", analysis.total_transcripts);
        println!("  Mean length: {:.0} bp", analysis.mean_length);
    }
    
    let original_seed_reads = collect_seed_reads_fast(
        &args.bam,
        &first_seed.chrom,
        first_seed.start as usize,
        first_seed.end as usize,
        args.include_supplementary,
    )?;

    println!("Seed reads: {}", original_seed_reads.len());
    println!("Setup time: {:?}", start_time.elapsed());
    println!();

    let mut all_discovered_loci: Vec<(Locus, f64, usize)> = Vec::new();
    let mut discovered_regions: Vec<(String, i64, i64)> = Vec::new();
    let mut current_component_reads = original_seed_reads.clone();

    // Use read cloud analysis if requested
    if args.use_read_clouds {
        println!("=== READ CLOUD ANALYSIS ===");
        println!("Looking for dense clusters of overlapping multi-mapping reads...");
        println!();
        
        let clouds = find_read_clouds(
            &args.bam,
            &original_seed_reads,
            &target_chroms,
            args.min_cloud_reads,
            args.min_cloud_span,
            args.max_cloud_span,
            args.strand_aware_clouds,
        )?;
        
        let filtered_clouds = filter_clouds(
            clouds,
            args.min_cloud_reads,
            args.min_density,
            args.require_secondary,
        );
        
        println!("  {} clouds passed quality filters", filtered_clouds.len());
        println!();
        
        for (i, cloud) in filtered_clouds.iter().enumerate() {
            let jaccard_seed = calculate_jaccard(&cloud.read_names, &original_seed_reads);
            
            // Check if this is the seed region
            let seed_region = (first_seed.chrom.clone(), first_seed.start, first_seed.end);
            let is_seed = cloud.chrom == seed_region.0 
                && cloud.start < seed_region.2 
                && cloud.end > seed_region.1;
            
            let seed_label = if is_seed { " [SEED]" } else { "" };
            println!(
                "  Cloud {}: {}:{}-{} | {} unique reads | Jaccard: {:.4}{}",
                i + 1,
                cloud.chrom,
                cloud.start,
                cloud.end,
                cloud.n_unique_reads(),
                jaccard_seed,
                seed_label
            );
            
            // Split large clouds if requested
            // This works even for seed-overlapping clouds - we'll split then filter
            let should_split = (args.split_clouds || args.splice_subclouds) && cloud.span() > 50000;
            
            if should_split {
                if args.splice_subclouds {
                    // Use splice-aware sub-clouds
                    println!("    Using splice-aware sub-cloud clustering...");
                    match create_splice_subclouds(
                        &args.bam,
                        &cloud.chrom,
                        cloud.start,
                        cloud.end,
                        Some(&cloud.read_names),
                        args.min_splice_sim,
                        10, // min_cluster_size
                    ) {
                        Ok(subclouds) => {
                            let refined = refine_sub_clouds(subclouds, &cloud.read_names, 5);
                            println!("    Created {} splice-based sub-clouds", refined.len());
                            
                            for (j, sub) in refined.iter().enumerate() {
                                // Skip sub-loci that overlap the seed region
                                let sub_overlaps_seed = sub.chrom == seed_region.0 
                                    && sub.start < seed_region.2 
                                    && sub.end > seed_region.1;
                                
                                if sub_overlaps_seed {
                                    println!("      Skipping sub-cloud {}.{} - overlaps seed", i + 1, j + 1);
                                    continue;
                                }
                                
                                let sub_reads = collect_seed_reads_fast(
                                    &args.bam,
                                    &sub.chrom,
                                    sub.start.max(1) as usize,
                                    sub.end.max(1) as usize,
                                    args.include_supplementary,
                                )?;
                                let sub_jaccard = calculate_jaccard(&sub_reads, &original_seed_reads);
                                
                                println!(
                                    "      Sub-cloud {}.{}: {}:{}-{} | {} reads | Jaccard: {:.4}",
                                    i + 1,
                                    j + 1,
                                    sub.chrom,
                                    sub.start,
                                    sub.end,
                                    sub.read_names.len(),
                                    sub_jaccard,
                                );
                                
                                let locus = Locus {
                                    chrom: sub.chrom.clone(),
                                    start: sub.start,
                                    end: sub.end,
                                    reads: sub_reads,
                                    isoseq_sim: None,
                                    seed_spliced: Some(sub.n_spliced_reads),
                                    target_spliced: None,
                                };
                                
                                all_discovered_loci.push((locus, sub_jaccard, 1));
                            }
                        }
                        Err(e) => {
                            eprintln!("    Warning: Splice sub-cloud creation failed: {}", e);
                        }
                    }
                } else {
                    // Use Iso-Seq boundary-based splitting (seed-aware)
                    let seed_tuple = (first_seed.chrom.clone(), first_seed.start, first_seed.end);
                    match split_cloud_by_isoseq_signals(
                        &args.bam,
                        &cloud.chrom,
                        cloud.start,
                        cloud.end,
                        Some(&cloud.read_names),
                        args.expected_genes,
                        Some(seed_tuple),
                    ) {
                        Ok(subloci) => {
                            for (j, sub) in subloci.iter().enumerate() {
                                // Skip sub-loci that overlap the seed region
                                let sub_overlaps_seed = sub.chrom == seed_region.0 
                                    && sub.start < seed_region.2 
                                    && sub.end > seed_region.1;
                                
                                if sub_overlaps_seed {
                                    println!("      Skipping sub-locus {}.{} - overlaps seed", i + 1, j + 1);
                                    continue;
                                }
                                
                                let sub_reads = collect_seed_reads_fast(
                                    &args.bam,
                                    &sub.chrom,
                                    sub.start.max(1) as usize,
                                    sub.end.max(1) as usize,
                                    args.include_supplementary,
                                )?;
                                let sub_jaccard = calculate_jaccard(&sub_reads, &original_seed_reads);
                                
                                println!(
                                    "      Sub-locus {}.{}: {}:{}-{} | Jaccard: {:.4} | {}",
                                    i + 1,
                                    j + 1,
                                    sub.chrom,
                                    sub.start,
                                    sub.end,
                                    sub_jaccard,
                                    sub.split_evidence.join(", ")
                                );
                                
                                let locus = Locus {
                                    chrom: sub.chrom.clone(),
                                    start: sub.start,
                                    end: sub.end,
                                    reads: sub_reads,
                                    isoseq_sim: None,
                                    seed_spliced: None,
                                    target_spliced: None,
                                };
                                
                                all_discovered_loci.push((locus, sub_jaccard, 1));
                            }
                        }
                        Err(e) => {
                            eprintln!("    Warning: Cloud splitting failed: {}", e);
                            if !is_seed {
                                // Fall back to using whole cloud
                                let locus = Locus {
                                    chrom: cloud.chrom.clone(),
                                    start: cloud.start,
                                    end: cloud.end,
                                    reads: cloud.read_names.clone(),
                                    isoseq_sim: None,
                                    seed_spliced: None,
                                    target_spliced: None,
                                };
                                all_discovered_loci.push((locus, jaccard_seed, 1));
                            }
                        }
                    }
                }
            } else if !is_seed {
                // Use whole cloud (not splitting and not seed)
                let locus = Locus {
                    chrom: cloud.chrom.clone(),
                    start: cloud.start,
                    end: cloud.end,
                    reads: cloud.read_names.clone(),
                    isoseq_sim: None,
                    seed_spliced: None,
                    target_spliced: None,
                };
                
                all_discovered_loci.push((locus, jaccard_seed, 1));
            } else {
                // Skip seed cloud when not splitting
                println!("      Skipping seed cloud (use --split-clouds to extract paralogs)");
            }
        }
        
        // Apply quick-win filters
        all_discovered_loci = filter_by_jaccard(all_discovered_loci, args.min_jaccard_locus);
        all_discovered_loci = merge_adjacent_loci(
            all_discovered_loci, 
            args.merge_adjacent_bp, 
            args.merge_min_jaccard
        );
        
        write_results(
            &args.output,
            &all_discovered_loci,
            &original_seed_reads,
            &first_seed.chrom,
            first_seed.start as usize,
            first_seed.end as usize,
            args.use_isoseq,
        )?;
        
        print_separator();
        println!("SUMMARY (Read Cloud Mode)");
        print_separator();
        println!("Seed reads: {}", original_seed_reads.len());
        println!("Clouds found: {}", all_discovered_loci.len());
        println!("Total time: {:?}", start_time.elapsed());
        println!("Output BED: {}", args.output);
        print_separator();
        
        return Ok(());
    }

    let mut last_read_count = current_component_reads.len();
    let mut converged = false;
    
    for iteration in 1..=args.max_iterations {
        let iter_start = std::time::Instant::now();
        
        println!("=== ITERATION {} (max: {}) ===", iteration, args.max_iterations);
        println!("Searching with {} reads", current_component_reads.len());
        println!();

        let mut new_loci = find_where_reads_map_parallel(
            &args.bam,
            &current_component_reads,
            &target_chroms,
            &discovered_regions,
            args.include_supplementary,
            args.cluster_distance,
            args.min_reads,
            args.threads,
        )?;

        if new_loci.is_empty() {
            println!("  No new loci found. Done!");
            break;
        }

        println!("  Found {} candidate loci", new_loci.len());
        
        // Validate with Iso-Seq
        if args.use_isoseq {
            println!("  Validating with Iso-Seq features...");
            let val_start = std::time::Instant::now();
            
            for locus in &mut new_loci {
                match validate_with_isoseq(
                    &args.bam,
                    &first_seed.chrom,
                    first_seed.start as usize,
                    first_seed.end as usize,
                    &locus.chrom,
                    locus.start,
                    locus.end,
                ) {
                    Ok((sim, seed_tx, target_tx)) => {
                        locus.isoseq_sim = Some(sim);
                        locus.seed_spliced = Some(seed_tx.total_transcripts);
                        locus.target_spliced = Some(target_tx.total_transcripts);
                    }
                    Err(e) => {
                        eprintln!("    Warning: Iso-Seq validation failed: {}", e);
                    }
                }
            }
            
            // Filter by Iso-Seq similarity
            let before_filter = new_loci.len();
            new_loci.retain(|l| l.is_valid(args.min_isoseq_sim));
            let after_filter = new_loci.len();
            
            println!("  Iso-Seq validation: {} -> {} loci (filtered {})",
                before_filter, after_filter, before_filter - after_filter);
            println!("  Validation time: {:?}", val_start.elapsed());
        }
        
        println!();

        let mut accepted_loci: Vec<Locus> = Vec::new();
        let mut bridge_reads: HashSet<String> = HashSet::default();

        for (i, locus) in new_loci.iter().enumerate() {
            let locus_start = locus.start.max(1) as usize;
            let locus_end = locus.end.max(1) as usize;

            let locus_reads = collect_seed_reads_fast(
                &args.bam,
                &locus.chrom,
                locus_start,
                locus_end,
                args.include_supplementary,
            )?;

            let jaccard_seed = calculate_jaccard(&locus_reads, &original_seed_reads);
            let jaccard_component = calculate_jaccard(&locus_reads, &current_component_reads);

            let isoseq_info = if let Some(sim) = locus.isoseq_sim {
                format!(" | IsoSeq: {:.2}", sim)
            } else {
                String::new()
            };

            println!(
                "  Locus {}: {}:{}-{} (span: {}){}",
                i + 1,
                locus.chrom,
                locus.start,
                locus.end,
                locus.span(),
                isoseq_info
            );
            println!("    Jaccard with seed: {:.4}", jaccard_seed);
            println!(
                "    Reads: {} total, {} shared with seed",
                locus_reads.len(),
                locus_reads.intersection(&original_seed_reads).count()
            );

            if jaccard_component < args.min_jaccard {
                println!("    REJECTED: Jaccard too low");
                continue;
            }

            let use_splice_split = args.split_splicing && locus.span() > 100_000;
            let shared_seed = locus_reads.intersection(&original_seed_reads).count();

            if !use_splice_split {
                if let Some(reason) = min_std_coverage_reject_reason(
                    &args,
                    &args.bam,
                    &locus.chrom,
                    locus.start,
                    locus.end,
                    shared_seed,
                )? {
                    println!("    REJECTED: {}", reason);
                    continue;
                }
            }

            // Splice-based splitting for large merged loci
            let loci_to_add = if use_splice_split {
                println!("    Large locus detected ({} bp), attempting splice-based splitting...", 
                         locus.span());
                match split_region_by_transcript_boundaries(
                    &args.bam,
                    &locus.chrom,
                    locus.start,
                    locus.end,
                    args.max_boundary_distance,
                    args.min_transcripts_for_split,
                    args.expected_genes,
                ) {
                    Ok(gene_loci) => {
                        if gene_loci.len() > 1 {
                            println!("    Split into {} gene loci", gene_loci.len());
                            gene_loci.iter().enumerate().for_each(|(j, g)| {
                                println!("      Gene {}: {}:{}-{} ({} transcripts)",
                                         j + 1, g.chrom, g.start, g.end, 
                                         g.n_transcripts);
                            });
                            gene_loci.into_iter().map(|g| Locus {
                                chrom: g.chrom,
                                start: g.start,
                                end: g.end,
                                reads: HashSet::default(), // Will be filled later
                                isoseq_sim: locus.isoseq_sim,
                                seed_spliced: locus.seed_spliced,
                                target_spliced: locus.target_spliced,
                            }).collect()
                        } else {
                            vec![Locus {
                                chrom: locus.chrom.clone(),
                                start: locus.start,
                                end: locus.end,
                                reads: locus_reads.clone(),
                                isoseq_sim: locus.isoseq_sim,
                                seed_spliced: locus.seed_spliced,
                                target_spliced: locus.target_spliced,
                            }]
                        }
                    }
                    Err(e) => {
                        eprintln!("    Warning: Splice-based splitting failed: {}", e);
                        vec![Locus {
                            chrom: locus.chrom.clone(),
                            start: locus.start,
                            end: locus.end,
                            reads: locus_reads.clone(),
                            isoseq_sim: locus.isoseq_sim,
                            seed_spliced: locus.seed_spliced,
                            target_spliced: locus.target_spliced,
                        }]
                    }
                }
            } else {
                vec![Locus {
                    chrom: locus.chrom.clone(),
                    start: locus.start,
                    end: locus.end,
                    reads: locus_reads.clone(),
                    isoseq_sim: locus.isoseq_sim,
                    seed_spliced: locus.seed_spliced,
                    target_spliced: locus.target_spliced,
                }]
            };

            for mut split_locus in loci_to_add {
                if split_locus.reads.is_empty() {
                    let ls = split_locus.start.max(1) as usize;
                    let le = split_locus.end.max(1) as usize;
                    split_locus.reads = collect_seed_reads_fast(
                        &args.bam,
                        &split_locus.chrom,
                        ls,
                        le,
                        args.include_supplementary,
                    )?;
                }

                if use_splice_split {
                    let sh = split_locus
                        .reads
                        .intersection(&original_seed_reads)
                        .count();
                    if let Some(reason) = min_std_coverage_reject_reason(
                        &args,
                        &args.bam,
                        &split_locus.chrom,
                        split_locus.start,
                        split_locus.end,
                        sh,
                    )? {
                        println!(
                            "    REJECTED segment {}:{}-{}: {}",
                            split_locus.chrom, split_locus.start, split_locus.end, reason
                        );
                        continue;
                    }
                }

                println!("    ACCEPTED");

                accepted_loci.push(split_locus.clone());
                discovered_regions.push((split_locus.chrom.clone(), split_locus.start, split_locus.end));
                bridge_reads.extend(split_locus.reads.clone());
                all_discovered_loci.push((split_locus, jaccard_seed, iteration));
            }
        }

        if accepted_loci.is_empty() {
            println!("\n  No loci accepted. Converged!");
            converged = true;
            break;
        }

        current_component_reads = bridge_reads;
        
        // Check for read-level convergence
        let new_read_count = current_component_reads.len();
        let reads_added = new_read_count.saturating_sub(last_read_count);
        last_read_count = new_read_count;

        println!();
        println!("  Iteration {}: {} loci accepted", iteration, accepted_loci.len());
        println!("  Reads added: {} (total: {})", reads_added, new_read_count);
        println!("  Iteration time: {:?}", iter_start.elapsed());
        println!("  Total loci: {}", all_discovered_loci.len());
        println!();
    }

    // Report convergence status
    if !converged && args.max_iterations > 0 {
        println!("\n  WARNING: Hit max iterations ({}) before convergence.", args.max_iterations);
        println!("  Consider increasing --max-iterations if more genes are expected.");
    } else if converged {
        println!("\n  Converged after discovering all connected components.");
    }

    write_results(
        &args.output,
        &all_discovered_loci,
        &original_seed_reads,
        &first_seed.chrom,
        first_seed.start as usize,
        first_seed.end as usize,
        args.use_isoseq,
    )?;

    // Evaluate against ground truth if provided
    if let (Some(ref gt_path), Some(ref family_id)) = (&args.ground_truth, &args.family_id) {
        println!();
        println!("Evaluating against ground truth...");
        
        match parse_ground_truth(gt_path, family_id) {
            Ok(ground_truth) => {
                if !ground_truth.is_empty() {
                    // Convert discovered loci to evaluation format
                    let detected: Vec<(String, i64, i64, HashSet<String>)> = all_discovered_loci
                        .iter()
                        .map(|(locus, _, _)| {
                            (locus.chrom.clone(), locus.start, locus.end, locus.reads.clone())
                        })
                        .collect();
                    
                    let novel_sup = NovelLocusSupport {
                        n_reads: all_discovered_loci
                            .iter()
                            .map(|(l, _, _)| l.reads.len())
                            .collect(),
                        jaccard: all_discovered_loci.iter().map(|(_, j, _)| *j).collect(),
                        distance_from_seed: all_discovered_loci
                            .iter()
                            .map(|(_, _, it)| *it)
                            .collect(),
                    };
                    let stats = evaluate_ground_truth(
                        &args,
                        &ground_truth,
                        &detected,
                        &seeds,
                        Some(&novel_sup),
                    );
                    print_evaluation_stats(&stats);
                    
                    // Write CSV if output path specified
                    let csv_path = args.output.replace(".bed", "_evaluation.csv");
                    if let Err(e) = write_evaluation_csv(&stats, &csv_path) {
                        eprintln!("Warning: Failed to write evaluation CSV: {}", e);
                    }
                } else {
                    println!("Warning: No ground truth genes found for {}", family_id);
                }
            }
            Err(e) => {
                eprintln!("Warning: Failed to parse ground truth: {}", e);
            }
        }
    }

    print_separator();
    println!("SUMMARY");
    print_separator();
    println!("Seed reads: {}", original_seed_reads.len());
    println!("Loci found: {}", all_discovered_loci.len());
    println!("Total time: {:?}", start_time.elapsed());
    println!("Output BED: {}", args.output);
    print_separator();

    Ok(())
}

fn write_results(
    output_path: &str,
    loci: &[(Locus, f64, usize)],
    seed_reads: &HashSet<String>,
    seed_chrom: &str,
    seed_start: usize,
    seed_end: usize,
    use_isoseq: bool,
) -> Result<()> {
    let mut file = File::create(output_path)
        .with_context(|| format!("Failed to create output file: {}", output_path))?;

    writeln!(file, "# Gene Family Detection Results (Iso-Seq Aware)")?;
    writeln!(
        file,
        "# Seed: {}:{}-{}",
        seed_chrom, seed_start, seed_end
    )?;
    writeln!(file, "# Seed reads: {}", seed_reads.len())?;
    writeln!(file, "# Loci found: {}", loci.len())?;
    writeln!(file, "#")?;
    
    if use_isoseq {
        writeln!(
            file,
            "# Columns: chrom, start, end, name, jaccard_seed, component, shared_with_seed, total_reads, span, isoseq_sim"
        )?;
    } else {
        writeln!(
            file,
            "# Columns: chrom, start, end, name, jaccard_seed, component, shared_with_seed, total_reads, span"
        )?;
    }
    writeln!(file, "#")?;

    let mut sorted_loci = loci.to_vec();
    sorted_loci.sort_by(|a, b| {
        a.0.chrom
            .cmp(&b.0.chrom)
            .then(a.0.start.cmp(&b.0.start))
    });

    for (i, (locus, jaccard_seed, component)) in sorted_loci.iter().enumerate() {
        let name = format!("locus_{}", i + 1);
        let shared_with_seed = locus.reads.intersection(seed_reads).count();

        if use_isoseq {
            let isoseq_sim = locus.isoseq_sim.unwrap_or(0.0);
            writeln!(
                file,
                "{}\t{}\t{}\t{}\t{:.6}\t{}\t{}\t{}\t{}\t{:.3}",
                locus.chrom,
                locus.start,
                locus.end,
                name,
                jaccard_seed,
                component,
                shared_with_seed,
                locus.reads.len(),
                locus.span(),
                isoseq_sim
            )?;
        } else {
            writeln!(
                file,
                "{}\t{}\t{}\t{}\t{:.6}\t{}\t{}\t{}\t{}",
                locus.chrom,
                locus.start,
                locus.end,
                name,
                jaccard_seed,
                component,
                shared_with_seed,
                locus.reads.len(),
                locus.span()
            )?;
        }
    }

    println!("\nWrote {} loci to {}", loci.len(), output_path);

    Ok(())
}
