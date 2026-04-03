//! Optional late pass: subdivide wide transitive loci using stricter coverage valleys
//! (same machinery as the main splitter, fixed parameters, no annotation).

use anyhow::Result;
use fxhash::FxHashSet as HashSet;

use crate::coverage_valley::{
    analyze_valleys, build_coverage_profile, merge_small_segments, split_by_peaks, split_by_valleys,
    SplitSegment, ValleyDetectionParams, ValleySplitFilter,
};
use crate::transitive_detector::TransitiveLocus;

#[derive(Debug, Clone)]
pub struct TightSplitParams {
    pub profile_bin_size: i64,
    pub valley: ValleyDetectionParams,
    pub merge_min_segment_bp: i64,
    pub merge_max_gap_bp: i64,
    /// Skip tight split for intervals narrower than this.
    pub min_input_span_bp: i64,
    /// If a split would yield more segments than this, keep the original locus unchanged.
    pub max_segments_per_locus: usize,
}

/// Apply a stricter coverage valley (and optional peak) split per locus. Intended to run
/// after envelope trim. Uses the same seed-read set as the main transitive pipeline.
pub fn tight_resplit_loci(
    bam_path: &str,
    loci: Vec<TransitiveLocus>,
    seed_reads: &HashSet<String>,
    params: &TightSplitParams,
) -> Result<Vec<TransitiveLocus>> {
    let n_in = loci.len();
    let mut out: Vec<TransitiveLocus> = Vec::new();
    println!(
        "  Input {} loci (min span {} bp, max {} segs per locus)",
        n_in,
        params.min_input_span_bp,
        params.max_segments_per_locus
    );

    for (idx, locus) in loci.into_iter().enumerate() {
        let span = locus.end - locus.start;
        if span < params.min_input_span_bp {
            out.push(locus);
            continue;
        }

        let profile = match build_coverage_profile(
            bam_path,
            &locus.chrom,
            locus.start,
            locus.end,
            seed_reads,
            params.profile_bin_size,
        ) {
            Ok(p) => p,
            Err(e) => {
                eprintln!(
                    "    Tight-split [{}] {}:{}-{}: profile failed: {}",
                    idx + 1,
                    locus.chrom,
                    locus.start,
                    locus.end,
                    e
                );
                out.push(locus);
                continue;
            }
        };

        let mut vp = params.valley.clone();
        vp.bin_size = params.profile_bin_size;

        let (peaks, valleys) = analyze_valleys(&profile, &vp);

        let mut segments = split_by_valleys(
            &locus.chrom,
            locus.start,
            locus.end,
            &profile,
            &valleys,
            &vp,
            ValleySplitFilter::SignificantOnly,
        );

        if segments.len() == 1 && peaks.len() > 1 {
            let peak_segments = split_by_peaks(
                &locus.chrom,
                locus.start,
                locus.end,
                &profile,
                &peaks,
                &vp,
            );
            if peak_segments.len() > 1 {
                segments = peak_segments;
            }
        }

        let merged: Vec<SplitSegment> = merge_small_segments(
            &segments,
            params.merge_min_segment_bp,
            params.merge_max_gap_bp,
        );

        if merged.is_empty() {
            out.push(locus);
            continue;
        }

        if merged.len() > params.max_segments_per_locus {
            println!(
                "    Tight-split: keep {}:{}-{} ({} segments > max {})",
                locus.chrom,
                locus.start,
                locus.end,
                merged.len(),
                params.max_segments_per_locus
            );
            out.push(locus);
            continue;
        }

        if merged.len() == 1 {
            let s = &merged[0];
            if s.start == locus.start && s.end == locus.end {
                out.push(locus);
                continue;
            }
        }

        println!(
            "    Tight-split: {}:{}-{} ({} bp) -> {} segments",
            locus.chrom,
            locus.start,
            locus.end,
            span,
            merged.len()
        );

        for seg in merged {
            out.push(TransitiveLocus {
                chrom: seg.chrom,
                start: seg.start,
                end: seg.end,
                reads: locus.reads.clone(),
                component: locus.component,
                jaccard_with_seed: locus.jaccard_with_seed,
                distance_from_seed: locus.distance_from_seed,
                peak_count: locus.peak_count,
                confidence: locus.confidence,
            });
        }
    }

    println!(
        "  Tight resplit done: {} input -> {} output loci",
        n_in,
        out.len()
    );
    Ok(out)
}
