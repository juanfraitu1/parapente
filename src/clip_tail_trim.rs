//! Optional shrink-only trim using A-rich 3' soft clips (polyA tail proxy) on primary alignments.
//! Strand-dominant: plus uses alignment ends; minus uses alignment starts (TES in ref coords).

use anyhow::Result;
use fxhash::FxHashSet as HashSet;
use noodles::bam;
use noodles::core::Position;
use noodles::sam::alignment::Record as SamRecord;
use crate::soft_clip;
use crate::transitive_detector::TransitiveLocus;

#[derive(Debug, Clone)]
pub struct ClipTailTrimParams {
    pub min_clip_len: usize,
    pub min_a_frac: f64,
    pub min_qualifying_reads: usize,
    pub end_slack_bp: i64,
    pub start_slack_bp: i64,
    /// Fraction of qualifying reads on one strand required to call that strand dominant.
    pub strand_dom_frac: f64,
    pub query_pad_bp: i64,
    pub min_result_span: i64,
    /// Per-locus counts and skip reasons (stdout).
    pub verbose: bool,
}

#[derive(Debug, Clone, Copy)]
enum ClipOutcome {
    NoThreePrime,
    TooShort,
    LowA,
    OkPlus,
    OkMinus,
}

#[derive(Default)]
struct ClipCounters {
    overlap: usize,
    no_three_prime: usize,
    too_short: usize,
    low_a: usize,
    ok_plus: usize,
    ok_minus: usize,
}

impl ClipCounters {
    fn record(&mut self, o: ClipOutcome) {
        match o {
            ClipOutcome::NoThreePrime => self.no_three_prime += 1,
            ClipOutcome::TooShort => self.too_short += 1,
            ClipOutcome::LowA => self.low_a += 1,
            ClipOutcome::OkPlus => self.ok_plus += 1,
            ClipOutcome::OkMinus => self.ok_minus += 1,
        }
    }
}

fn median_i64(xs: &mut Vec<i64>) -> Option<i64> {
    if xs.is_empty() {
        return None;
    }
    xs.sort_unstable();
    let n = xs.len();
    let mid = n / 2;
    if n % 2 == 1 {
        Some(xs[mid])
    } else {
        Some((xs[mid - 1] + xs[mid]) / 2)
    }
}

fn frac_a_bases(slice: &[u8]) -> f64 {
    if slice.is_empty() {
        return 0.0;
    }
    let mut c = 0usize;
    for &b in slice {
        if b == b'A' || b == b'a' {
            c += 1;
        }
    }
    c as f64 / slice.len() as f64
}

/// Classify 3' soft clip for one overlapping primary read (polyA proxy).
fn clip_outcome(
    record: &bam::Record,
    params: &ClipTailTrimParams,
    seq_buf: &mut Vec<u8>,
) -> ClipOutcome {
    let Some((q0, q1)) = soft_clip::three_prime_soft_clip_query_range(record) else {
        return ClipOutcome::NoThreePrime;
    };
    let clip_len = q1.saturating_sub(q0);
    if clip_len < params.min_clip_len {
        return ClipOutcome::TooShort;
    }
    seq_buf.clear();
    let seq = record.sequence();
    let slen = seq.len();
    if q1 > slen {
        return ClipOutcome::TooShort;
    }
    for i in q0..q1 {
        let Some(b) = seq.get(i) else {
            return ClipOutcome::TooShort;
        };
        seq_buf.push(b);
    }
    if frac_a_bases(seq_buf) < params.min_a_frac {
        return ClipOutcome::LowA;
    }
    if record.flags().is_reverse_complemented() {
        ClipOutcome::OkMinus
    } else {
        ClipOutcome::OkPlus
    }
}

/// Shrink loci using median TES anchors from reads with A-rich 3' soft clips. Never expands bounds.
pub fn trim_loci_clip_tail(
    bam_path: &str,
    loci: Vec<TransitiveLocus>,
    params: &ClipTailTrimParams,
) -> Result<Vec<TransitiveLocus>> {
    let mut reader = bam::io::indexed_reader::Builder::default().build_from_path(bam_path)?;
    let header = reader.read_header()?;
    let mut out = Vec::with_capacity(loci.len());
    let mut seq_buf: Vec<u8> = Vec::new();

    println!(
        "  Clip-tail trim: min_clip={} min_a_frac={} min_reads={} dom_frac={} verbose={}",
        params.min_clip_len,
        params.min_a_frac,
        params.min_qualifying_reads,
        params.strand_dom_frac,
        params.verbose
    );

    for locus in loci {
        let read_set: HashSet<&str> = locus.reads.iter().map(|s| s.as_str()).collect();
        let ref_seqs = header.reference_sequences();
        let chrom_name: bstr::BString = locus.chrom.as_bytes().into();
        let Some(tid) = ref_seqs.get_index_of(&chrom_name) else {
            out.push(locus);
            continue;
        };
        let ref_len = ref_seqs[tid].length().get() as i64;

        let pad = params.query_pad_bp.max(0);
        let q0 = (locus.start - pad).max(1);
        let q1 = (locus.end + pad).min(ref_len);
        if q1 <= q0 {
            out.push(locus);
            continue;
        }

        let start = Position::new(q0 as usize).unwrap_or(Position::MIN);
        let end = Position::new(q1 as usize).unwrap_or(Position::MAX);
        let region = noodles::core::Region::new(locus.chrom.clone(), start..=end);

        let mut plus_anchors: Vec<i64> = Vec::new();
        let mut minus_anchors: Vec<i64> = Vec::new();
        let mut n_plus = 0usize;
        let mut n_minus = 0usize;
        let mut ctr = ClipCounters::default();

        let query_iter = match reader.query(&header, &region) {
            Ok(it) => it,
            Err(e) => {
                eprintln!(
                    "    clip-tail: query failed {}:{}-{}: {}",
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
            if record.flags().is_unmapped()
                || record.flags().is_secondary()
                || record.flags().is_supplementary()
            {
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

            ctr.overlap += 1;
            let o = clip_outcome(&record, params, &mut seq_buf);
            if params.verbose {
                ctr.record(o);
            }
            match o {
                ClipOutcome::OkPlus => {
                    plus_anchors.push(aln_end);
                    n_plus += 1;
                }
                ClipOutcome::OkMinus => {
                    minus_anchors.push(aln_start);
                    n_minus += 1;
                }
                _ => {}
            }
        }

        let total = n_plus + n_minus;
        let mut ns = locus.start;
        let mut ne = locus.end;
        let mut changed = false;
        let mut skip_reason: Option<&'static str> = None;

        if total > 0 {
            let f_plus = n_plus as f64 / total as f64;
            let f_minus = n_minus as f64 / total as f64;

            if f_plus >= params.strand_dom_frac && plus_anchors.len() >= params.min_qualifying_reads {
                if let Some(med) = median_i64(&mut plus_anchors) {
                    let cand = (med + params.end_slack_bp).min(ref_len);
                    let new_end = cand.min(ne);
                    if new_end > ns && new_end - ns >= params.min_result_span && new_end < ne {
                        ne = new_end;
                        changed = true;
                    } else if params.verbose {
                        skip_reason = Some("plus_median_does_not_shrink_or_span_floor");
                    }
                }
            } else if f_minus >= params.strand_dom_frac
                && minus_anchors.len() >= params.min_qualifying_reads
            {
                if let Some(med) = median_i64(&mut minus_anchors) {
                    let cand = (med - params.start_slack_bp).max(1);
                    let new_start = cand.max(ns);
                    if new_start < ne && ne - new_start >= params.min_result_span && new_start > ns {
                        ns = new_start;
                        changed = true;
                    } else if params.verbose {
                        skip_reason = Some("minus_median_does_not_shrink_or_span_floor");
                    }
                }
            } else if params.verbose {
                if f_plus < params.strand_dom_frac && f_minus < params.strand_dom_frac {
                    skip_reason = Some("no_strand_dominance");
                } else if f_plus >= params.strand_dom_frac
                    && plus_anchors.len() < params.min_qualifying_reads
                {
                    skip_reason = Some("plus_dominant_but_insufficient_qualifying_reads");
                } else if f_minus >= params.strand_dom_frac
                    && minus_anchors.len() < params.min_qualifying_reads
                {
                    skip_reason = Some("minus_dominant_but_insufficient_qualifying_reads");
                } else {
                    skip_reason = Some("dominance_branch_not_applicable");
                }
            }
        } else if params.verbose {
            skip_reason = Some(if ctr.overlap == 0 {
                "no_primary_overlap_in_readset"
            } else {
                "no_qualifying_reads"
            });
        }

        if params.verbose {
            let sr = skip_reason.unwrap_or(if changed { "applied_trim" } else { "no_change" });
            println!(
                "    clip-tail diag: {}:{}-{} overlap={} no_3p={} short={} low_a={} ok+={} ok-={} q_total={} skip={}",
                locus.chrom,
                locus.start,
                locus.end,
                ctr.overlap,
                ctr.no_three_prime,
                ctr.too_short,
                ctr.low_a,
                ctr.ok_plus,
                ctr.ok_minus,
                total,
                sr
            );
        }

        if changed {
            println!(
                "    clip-tail: {}:{}-{} -> {}:{}-{} (plus_q={} minus_q={})",
                locus.chrom,
                locus.start,
                locus.end,
                locus.chrom,
                ns,
                ne,
                plus_anchors.len(),
                minus_anchors.len()
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
        } else {
            out.push(locus);
        }
    }

    Ok(out)
}
