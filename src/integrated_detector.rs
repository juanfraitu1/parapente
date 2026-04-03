//! Integrated gene family detector combining three methods:
//! 
//! 1. **Multi-mapping read discovery**: Find where seed reads map genome-wide
//! 2. **Coverage valley splitting**: Split overextended loci at coverage minima
//! 3. **Iso-Seq validation**: Resolve ambiguous valleys using transcript structure
//! 4. **Strand-adaptive augmentation** (optional): If combined-strand coverage hides dips,
//!    or the locus is wide, repeat peak/valley detection on plus and minus strand profiles
//!    (same query reads) and merge non-redundant split points.
//!
//! Architecture:
//! ```
//! Seed reads → Multi-mapping discovery → Candidate loci
//!                                         ↓
//!                            Coverage profile building (combined; optional per-strand)
//!                                         ↓
//!                         Peak/valley detection
//!                                         ↓
//!                    ┌────────────────────┼────────────────────┐
//!                    ↓                    ↓                    ↓
//!              Deep valley          Shallow valley      Near splice site
//!              (< 10% max)          (10-30% max)        (spanning reads)
//!                    ↓                    ↓                    ↓
//!              AUTO-SPLIT         ISO-SEQ VALIDATE    ISO-SEQ VALIDATE
//!                    ↓                    ↓                    ↓
//!                                    Split or Keep
//!                                         ↓
//!                              Read-space overlap refinement
//!                                         ↓
//!                                  Final gene cores
//! ```

use anyhow::Result;
use fxhash::FxHashSet as HashSet;
use noodles::bam;
use noodles::core::Position;
use noodles::sam::alignment::Record as SamRecord;
use std::collections::HashMap;

use crate::coverage_valley::{
    analyze_valleys, build_coverage_profile, build_strand_coverage_profiles, split_by_valleys,
    CoverageValley, ValleyDetectionParams, ValleySplitFilter,
};
use crate::isoseq::analyze_locus_transcripts;
use crate::junction_refinement::{refine_core_bounds_from_spliced_exons, JunctionRefinementParams};

/// Iso-Seq validation result for a valley
#[derive(Debug, Clone)]
pub struct ValleyValidation {
    pub has_spanning_reads: bool,    // Reads spanning across the valley
    pub has_splice_junctions: bool,  // Splice sites at valley position
    pub mean_transcript_length: f64, // Average transcript length
    pub n_transcripts: usize,        // Number of transcripts
}

/// A detected gene core after all processing
#[derive(Debug, Clone)]
pub struct GeneCore {
    pub chrom: String,
    pub start: i64,
    pub end: i64,
    pub reads: HashSet<String>,
    pub jaccard_with_seed: f64,
    pub coverage_valleys: usize,     // Number of valleys detected
    pub significant_splits: usize,   // Number of actual splits made
    pub isoseq_validated: bool,      // Whether Iso-Seq validated boundaries
    pub confidence: CoreConfidence,
}

/// Confidence level for a gene core
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CoreConfidence {
    High,    // Deep valley + Iso-Seq support OR clear read-space overlap
    Medium,  // Shallow valley validated by Iso-Seq
    Low,     // Uncertain boundaries, needs manual review
}

impl std::fmt::Display for CoreConfidence {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            CoreConfidence::High => write!(f, "HIGH"),
            CoreConfidence::Medium => write!(f, "MEDIUM"),
            CoreConfidence::Low => write!(f, "LOW"),
        }
    }
}

/// Parameters for integrated detection
#[derive(Debug, Clone)]
pub struct IntegratedParams {
    pub valley_params: ValleyDetectionParams,
    pub auto_split_threshold: f64,    // Valley depth ratio for auto-split (default: 0.1)
    pub validate_threshold: f64,      // Valley depth ratio for validation (default: 0.3)
    pub min_read_overlap_bp: i64,     // Minimum read-space overlap for core refinement
    pub require_isoseq_for_ambiguous: bool, // Always use Iso-Seq for shallow valleys
    /// When true, merge split candidates from plus and minus strand coverage (same query reads).
    pub strand_adaptive_split: bool,
    /// Minimum locus span (bp) to always run per-strand peak/valley analysis in addition to combined.
    pub strand_adaptive_min_span_bp: i64,
    /// Optional CIGAR exon-union refinement (no gene annotation).
    pub junction_refinement: Option<JunctionRefinementParams>,
}

impl Default for IntegratedParams {
    fn default() -> Self {
        Self {
            valley_params: ValleyDetectionParams::default(),
            auto_split_threshold: 0.1,   // < 10% of max = auto split
            validate_threshold: 0.3,     // 10-30% = validate with Iso-Seq
            min_read_overlap_bp: 1000,
            require_isoseq_for_ambiguous: true,
            strand_adaptive_split: true,
            strand_adaptive_min_span_bp: 50_000,
            junction_refinement: None,
        }
    }
}

fn valley_near_existing(position: i64, existing: &[CoverageValley], sep_bp: i64) -> bool {
    existing
        .iter()
        .any(|v| (v.position - position).abs() < sep_bp)
}

/// Add valleys found on each strand's coverage graph when they are not redundant with `validated`.
fn augment_validated_with_strand_valleys(
    bam_path: &str,
    chrom: &str,
    start: i64,
    end: i64,
    query_reads: &HashSet<String>,
    params: &IntegratedParams,
    validated: &mut Vec<CoverageValley>,
    validations_used: &mut usize,
) -> Result<()> {
    let (plus_p, minus_p) = build_strand_coverage_profiles(
        bam_path,
        chrom,
        start,
        end,
        query_reads,
        params.valley_params.bin_size,
    )?;
    let sep = params.valley_params.bin_size.max(100) * 2;
    for prof in [plus_p, minus_p] {
        let (peaks, valleys) = analyze_valleys(&prof, &params.valley_params);
        if peaks.len() < 2 {
            continue;
        }
        for valley in &valleys {
            if valley_near_existing(valley.position, validated, sep) {
                continue;
            }
            let (should_split, was_validated, _) =
                decide_split(valley, params, bam_path, chrom)?;
            if was_validated {
                *validations_used += 1;
            }
            if should_split {
                validated.push(valley.clone());
            }
        }
    }
    Ok(())
}

/// Validate a valley using Iso-Seq signals
/// 
/// Returns true if valley should be used as a split point
fn validate_valley_with_isoseq(
    bam_path: &str,
    chrom: &str,
    valley_pos: i64,
    window_size: i64,  // Window around valley to check
) -> Result<ValleyValidation> {
    let start = (valley_pos - window_size).max(1);
    let end = valley_pos + window_size;

    // Analyze transcripts in window around valley
    let analysis = analyze_locus_transcripts(bam_path, chrom, start, end)?;

    // Check for spanning reads and splice patterns
    // This is a simplified version - full implementation would:
    // 1. Check for reads that span across the valley position
    // 2. Look for splice junctions near the valley
    // 3. Compare transcript continuity
    
    let has_spanning = analysis.total_transcripts > 0 && analysis.coverage_std < 20.0;
    let has_splice = analysis.coverage_std > 30.0;  // High variability = introns

    Ok(ValleyValidation {
        has_spanning_reads: has_spanning,
        has_splice_junctions: has_splice,
        mean_transcript_length: analysis.mean_length,
        n_transcripts: analysis.total_transcripts,
    })
}

/// Decide whether to split at a valley based on depth and Iso-Seq validation
/// 
/// Returns: (should_split, was_validated, confidence)
fn decide_split(
    valley: &CoverageValley,
    params: &IntegratedParams,
    bam_path: &str,
    chrom: &str,
) -> Result<(bool, bool, CoreConfidence)> {
    // Deep valley: auto-split
    if valley.depth_ratio < params.auto_split_threshold {
        return Ok((true, false, CoreConfidence::High));
    }

    // Shallow valley: needs validation
    if params.require_isoseq_for_ambiguous && valley.depth_ratio < params.validate_threshold {
        let validation = validate_valley_with_isoseq(bam_path, chrom, valley.position, 5000)?;
        
        // If spanning reads or splice junctions cross the valley, don't split (it's intragenic)
        if validation.has_spanning_reads || validation.has_splice_junctions {
            return Ok((false, true, CoreConfidence::Medium));
        }
        
        // No spanning reads: safe to split
        return Ok((true, true, CoreConfidence::Medium));
    }

    // Too shallow, don't split
    Ok((false, false, CoreConfidence::Low))
}

/// Process a candidate locus with integrated detection
/// 
/// 1. Build coverage profile from query reads
/// 2. Detect peaks and valleys
/// 3. Validate ambiguous valleys with Iso-Seq
/// 4. Split into gene cores
/// 5. Refine boundaries using read-space overlap
pub fn process_locus_integrated(
    bam_path: &str,
    chrom: &str,
    start: i64,
    end: i64,
    query_reads: &HashSet<String>,
    seed_reads: &HashSet<String>,
    params: &IntegratedParams,
) -> Result<Vec<GeneCore>> {
    // Step 1: Build coverage profile
    let profile = build_coverage_profile(
        bam_path,
        chrom,
        start,
        end,
        query_reads,
        params.valley_params.bin_size,
    )?;

    // Step 2: Detect peaks and valleys (combined-strand signal)
    let (peaks, valleys) = analyze_valleys(&profile, &params.valley_params);
    let span = end - start;

    // Step 3: Validate valleys and decide splits (combined graph)
    let mut validated_valleys: Vec<CoverageValley> = vec![];
    let mut validations_used = 0;

    if peaks.len() >= 2 {
        for valley in &valleys {
            let (should_split, was_validated, _confidence) =
                decide_split(valley, params, bam_path, chrom)?;

            if was_validated {
                validations_used += 1;
            }

            if should_split {
                validated_valleys.push(valley.clone());
            }
        }
    }

    let n_combined_validated = validated_valleys.len();
    let run_strand = params.strand_adaptive_split
        && (span >= params.strand_adaptive_min_span_bp
            || (validated_valleys.is_empty() && !valleys.is_empty())
            || (peaks.len() < 2 && span >= params.strand_adaptive_min_span_bp));

    if run_strand {
        augment_validated_with_strand_valleys(
            bam_path,
            chrom,
            start,
            end,
            query_reads,
            params,
            &mut validated_valleys,
            &mut validations_used,
        )?;
    }

    if validated_valleys.is_empty() {
        let mut rstart = start;
        let mut rend = end;
        let mut reads = query_reads.clone();
        if let Some(jr) = params.junction_refinement.as_ref() {
            let (ns, ne) = refine_core_bounds_from_spliced_exons(
                bam_path,
                chrom,
                rstart,
                rend,
                &reads,
                jr,
            )?;
            if ns != rstart || ne != rend {
                println!(
                    "      Junction refinement (no valley split): {}:{}-{} -> {}:{}-{}",
                    chrom, rstart, rend, chrom, ns, ne
                );
                rstart = ns;
                rend = ne;
                reads = collect_reads_in_segment(bam_path, chrom, rstart, rend, query_reads)?;
            }
        }
        let jaccard = calculate_jaccard(&reads, seed_reads);
        return Ok(vec![GeneCore {
            chrom: chrom.to_string(),
            start: rstart,
            end: rend,
            reads,
            jaccard_with_seed: jaccard,
            coverage_valleys: valleys.len(),
            significant_splits: 0,
            isoseq_validated: false,
            confidence: CoreConfidence::High,
        }]);
    }

    let n_strand_extra = validated_valleys.len().saturating_sub(n_combined_validated);

    // Step 4: Split locus
    let segments = split_by_valleys(
        chrom,
        start,
        end,
        &profile,
        &validated_valleys,
        &params.valley_params,
        ValleySplitFilter::CallerValidated,
    );

    // Step 5: Build GeneCores from segments
    let mut cores = vec![];
    for (_i, seg) in segments.iter().enumerate() {
        // Collect reads for this segment
        let mut seg_start = seg.start;
        let mut seg_end = seg.end;
        let mut seg_reads = collect_reads_in_segment(bam_path, chrom, seg_start, seg_end, query_reads)?;

        if let Some(jr) = params.junction_refinement.as_ref() {
            let (ns, ne) = refine_core_bounds_from_spliced_exons(
                bam_path,
                chrom,
                seg_start,
                seg_end,
                &seg_reads,
                jr,
            )?;
            if ns != seg_start || ne != seg_end {
                println!(
                    "      Junction refinement: {}:{}-{} -> {}:{}-{}",
                    chrom, seg_start, seg_end, chrom, ns, ne
                );
                seg_start = ns;
                seg_end = ne;
                seg_reads =
                    collect_reads_in_segment(bam_path, chrom, seg_start, seg_end, query_reads)?;
            }
        }

        let jaccard = calculate_jaccard(&seg_reads, seed_reads);

        // Determine confidence
        let confidence = if validated_valleys.len() == valleys.iter().filter(|v| v.is_significant).count() {
            CoreConfidence::High
        } else if validations_used > 0 {
            CoreConfidence::Medium
        } else {
            CoreConfidence::High // All valleys were deep (auto-split)
        };

        cores.push(GeneCore {
            chrom: seg.chrom.clone(),
            start: seg_start,
            end: seg_end,
            reads: seg_reads,
            jaccard_with_seed: jaccard,
            coverage_valleys: valleys.len(),
            significant_splits: validated_valleys.len(),
            isoseq_validated: validations_used > 0,
            confidence,
        });
    }

    // Print summary
    println!("    Integrated analysis for {}:{}-{}", chrom, start, end);
    println!("      Coverage peaks: {}", peaks.len());
    println!("      Valleys detected: {} ({} significant)", 
             valleys.len(), 
             valleys.iter().filter(|v| v.is_significant).count());
    if run_strand {
        println!(
            "      Strand-adaptive: yes (extra split points vs combined: {})",
            n_strand_extra
        );
    } else {
        println!("      Strand-adaptive: no");
    }
    println!("      Iso-Seq validations: {}", validations_used);
    println!("      Segments created: {}", cores.len());
    for (i, core) in cores.iter().enumerate() {
        println!("        Core {}: {}:{}-{} (span: {} bp, Jaccard: {:.3}, confidence: {})",
                 i + 1,
                 core.chrom,
                 core.start,
                 core.end,
                 core.end - core.start,
                 core.jaccard_with_seed,
                 core.confidence
        );
    }

    Ok(cores)
}

/// Collect reads that map within a specific segment
fn collect_reads_in_segment(
    bam_path: &str,
    chrom: &str,
    start: i64,
    end: i64,
    candidate_reads: &HashSet<String>,
) -> Result<HashSet<String>> {
    let mut reader = bam::io::indexed_reader::Builder::default()
        .build_from_path(bam_path)?;
    let header = reader.read_header()?;

    let region_start = Position::try_from(start.max(1) as usize)?;
    let region_end = Position::try_from(end.max(1) as usize)?;
    let region = noodles::core::Region::new(chrom, region_start..=region_end);

    let mut found_reads = HashSet::default();

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

        if candidate_reads.contains(name) {
            found_reads.insert(name.to_string());
        }
    }

    Ok(found_reads)
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

/// Find gene family members using the integrated approach
///
/// This is the main entry point for the integrated detector.
/// It combines multi-mapping discovery with coverage valley splitting
/// and Iso-Seq validation.
pub fn detect_gene_family_integrated(
    bam_path: &str,
    _seed_chrom: &str,
    _seed_start: i64,
    _seed_end: i64,
    seed_reads: &HashSet<String>,
    target_chroms: &[String],
    params: &IntegratedParams,
    cluster_distance: i64,
    min_reads: usize,
) -> Result<Vec<GeneCore>> {
    println!("=== INTEGRATED DETECTION ===");
    println!("Phase 1: Multi-mapping read discovery");
    println!("  Seed reads: {}", seed_reads.len());
    
    // Discover candidate loci using multi-mapping reads
    let candidate_loci = discover_candidate_loci(
        bam_path,
        seed_reads,
        target_chroms,
        cluster_distance,
        min_reads,
    )?;
    
    println!("  Discovered {} candidate loci", candidate_loci.len());
    println!();

    // Process each locus with integrated detection
    println!("Phase 2: Coverage valley splitting + Iso-Seq validation");
    let mut all_cores: Vec<GeneCore> = vec![];

    for (i, (chrom, start, end, locus_reads)) in candidate_loci.iter().enumerate() {
        println!("  Locus {}: {}:{}-{} ({} reads)", 
                 i + 1, chrom, start, end, locus_reads.len());

        let cores = process_locus_integrated(
            bam_path,
            chrom,
            *start,
            *end,
            locus_reads,
            seed_reads,
            params,
        )?;

        all_cores.extend(cores);
    }

    println!();
    println!("=== INTEGRATED DETECTION COMPLETE ===");
    println!("Total gene cores found: {}", all_cores.len());
    
    let high_conf = all_cores.iter().filter(|c| c.confidence == CoreConfidence::High).count();
    let med_conf = all_cores.iter().filter(|c| c.confidence == CoreConfidence::Medium).count();
    let low_conf = all_cores.iter().filter(|c| c.confidence == CoreConfidence::Low).count();
    
    println!("  High confidence: {}", high_conf);
    println!("  Medium confidence: {}", med_conf);
    println!("  Low confidence: {}", low_conf);

    Ok(all_cores)
}

/// Discover candidate loci using multi-mapping reads
/// 
/// Scans chromosomes for regions where seed reads map
fn discover_candidate_loci(
    bam_path: &str,
    seed_reads: &HashSet<String>,
    target_chroms: &[String],
    cluster_distance: i64,
    min_reads: usize,
) -> Result<Vec<(String, i64, i64, HashSet<String>)>> {
    let mut all_matches: Vec<(String, i64, i64, String)> = vec![];

    for chrom in target_chroms {
        let mut reader = bam::io::indexed_reader::Builder::default()
            .build_from_path(bam_path)?;
        let header = reader.read_header()?;

        let ref_seqs = header.reference_sequences();
        let chrom_name: bstr::BString = chrom.as_bytes().into();
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

            if !seed_reads.contains(name) {
                continue;
            }

            let start = record.alignment_start()
                .transpose()?
                .map(|p: Position| p.get() as i64)
                .unwrap_or(0);
            let end = record.alignment_end()
                .transpose()?
                .map(|p: Position| p.get() as i64)
                .unwrap_or(start + 1);

            all_matches.push((chrom.clone(), start, end, name.to_string()));
        }
    }

    // Cluster into loci
    let mut by_chrom: HashMap<String, Vec<(i64, i64, String)>> = HashMap::new();
    for (chrom, start, end, name) in all_matches {
        by_chrom.entry(chrom).or_default().push((start, end, name));
    }

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
            loci.push((chrom.clone(), current_cluster[0].0, current_end, reads));
        }
    }

    Ok(loci)
}
