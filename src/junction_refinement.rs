//! Tighten predicted locus bounds using RNA alignment structure (no GFF).
//! Signals: CIGAR exon blocks with N introns, clustered splice acceptors, read alignment extent.

use anyhow::Result;
use fxhash::FxHashSet as HashSet;
use noodles::bam;
use noodles::core::Position;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::Record as SamRecord;
use std::cmp::{max, min};

/// Bounds refinement from BAM-only cues.
#[derive(Debug, Clone)]
pub struct JunctionRefinementParams {
    pub max_trim_frac_per_side: f64,
    pub min_segment_span: i64,
    pub include_supplementary: bool,
    /// Minimum distinct reads supporting a clustered splice site for acceptor anchor (0 = off).
    pub min_junction_read_support: usize,
    pub junction_cluster_slack_bp: i64,
    /// If > 0, clamp to median alignment start/end plus/minus this slack (0 = off).
    pub read_end_guard_slack_bp: i64,
    /// Acceptors (intron 3' end / downstream exon start) within this window of exon-union left edge qualify.
    pub acceptor_window_left_bp: i64,
    pub acceptor_window_right_bp: i64,
}

impl Default for JunctionRefinementParams {
    fn default() -> Self {
        Self {
            max_trim_frac_per_side: 0.35,
            min_segment_span: 500,
            include_supplementary: false,
            min_junction_read_support: 3,
            junction_cluster_slack_bp: 12,
            read_end_guard_slack_bp: 2000,
            acceptor_window_left_bp: 2500,
            acceptor_window_right_bp: 1200,
        }
    }
}

fn exon_blocks_from_record(record: &bam::Record) -> Option<Vec<(i64, i64)>> {
    if record.flags().is_unmapped() {
        return None;
    }

    let aln_start = record
        .alignment_start()
        .transpose()
        .ok()??
        .get() as i64;

    let mut ref_pos = aln_start;
    let mut block_start: Option<i64> = None;
    let mut blocks: Vec<(i64, i64)> = Vec::new();

    for result in record.cigar().iter() {
        let op = result.ok()?;
        let len = op.len() as i64;
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch | Kind::Deletion => {
                if block_start.is_none() {
                    block_start = Some(ref_pos);
                }
                ref_pos += len;
            }
            Kind::Skip => {
                if let Some(bs) = block_start.take() {
                    if ref_pos > bs {
                        blocks.push((bs, ref_pos));
                    }
                }
                ref_pos += len;
            }
            _ => {}
        }
    }

    if let Some(bs) = block_start {
        if ref_pos > bs {
            blocks.push((bs, ref_pos));
        }
    }

    Some(blocks)
}

/// Acceptors: first reference base of each exon after an intron (CIGAR N), same as splice_clouds `intron_end`.
fn acceptor_sites_from_record(record: &bam::Record) -> Option<Vec<i64>> {
    if record.flags().is_unmapped() {
        return None;
    }
    let aln_start = record
        .alignment_start()
        .transpose()
        .ok()??
        .get() as i64;
    let mut ref_pos = aln_start;
    let mut sites = Vec::new();

    for result in record.cigar().iter() {
        let op = result.ok()?;
        let len = op.len() as i64;
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch | Kind::Deletion => {
                ref_pos += len;
            }
            Kind::Skip => {
                ref_pos += len;
                sites.push(ref_pos);
            }
            _ => {}
        }
    }
    if sites.is_empty() {
        None
    } else {
        Some(sites)
    }
}

fn clip_block_to_region(b: (i64, i64), lo: i64, hi: i64) -> Option<(i64, i64)> {
    let s = max(b.0, lo);
    let e = min(b.1, hi);
    if e > s {
        Some((s, e))
    } else {
        None
    }
}

fn merge_intervals(mut ivs: Vec<(i64, i64)>) -> Vec<(i64, i64)> {
    if ivs.is_empty() {
        return vec![];
    }
    ivs.sort_by_key(|x| x.0);
    let mut out = vec![ivs[0]];
    for (s, e) in ivs.into_iter().skip(1) {
        let last = out.len() - 1;
        if out[last].1 >= s {
            out[last].1 = out[last].1.max(e);
        } else {
            out.push((s, e));
        }
    }
    out
}

fn median_sorted(sorted: &mut [i64]) -> i64 {
    if sorted.is_empty() {
        return 0;
    }
    let n = sorted.len();
    sorted.sort_unstable();
    sorted[n / 2]
}

/// Cluster genomic positions; each cluster records unique read count.
fn cluster_site_support(
    mut hits: Vec<(i64, String)>,
    slack: i64,
) -> Vec<(i64, usize)> {
    if hits.is_empty() {
        return vec![];
    }
    hits.sort_by_key(|x| x.0);
    let mut out: Vec<(i64, usize)> = Vec::new();
    let mut i = 0;
    while i < hits.len() {
        let mut j = i + 1;
        let mut names = HashSet::default();
        names.insert(hits[i].1.clone());
        let mut sum = hits[i].0;
        let mut count = 1_i64;
        while j < hits.len() && hits[j].0 - hits[i].0 <= slack {
            names.insert(hits[j].1.clone());
            sum += hits[j].0;
            count += 1;
            j += 1;
        }
        let center = sum / count;
        out.push((center, names.len()));
        i = j;
    }
    out
}

/// Cap left over-trim (do not move start right of bulk 5' plus slack) and
/// cap right over-trim (do not move end left of bulk 3' minus slack).
fn apply_read_end_guards(ns: i64, ne: i64, aln_starts: &[i64], aln_ends: &[i64], slack: i64) -> (i64, i64) {
    if aln_starts.is_empty() || aln_ends.is_empty() || slack <= 0 {
        return (ns, ne);
    }
    let mut s = aln_starts.to_vec();
    let mut e = aln_ends.to_vec();
    let ms = median_sorted(&mut s);
    let me = median_sorted(&mut e);
    let a = min(ns, ms + slack);
    let b = max(ne, me - slack);
    (a, b)
}

fn apply_acceptor_guard(
    ns: i64,
    u_start: i64,
    end: i64,
    acceptor_hits: Vec<(i64, String)>,
    min_reads: usize,
    slack: i64,
    win_left: i64,
    win_right: i64,
) -> i64 {
    if min_reads == 0 || acceptor_hits.is_empty() {
        return ns;
    }
    let clustered = cluster_site_support(acceptor_hits, slack);
    let mut best: Option<i64> = None;
    for (pos, sup) in clustered {
        if sup < min_reads {
            continue;
        }
        if pos < u_start - win_left || pos > u_start + win_right {
            continue;
        }
        if pos > end || pos < 0 {
            continue;
        }
        match best {
            None => best = Some(pos),
            Some(p0) => {
                if pos < p0 {
                    best = Some(pos);
                }
            }
        }
    }
    match best {
        Some(cap) => min(ns, cap),
        None => ns,
    }
}

/// Refine core bounds: exon-union trim, then median alignment extent, then acceptor anchor.
pub fn refine_core_bounds_from_spliced_exons(
    bam_path: &str,
    chrom: &str,
    start: i64,
    end: i64,
    read_names: &HashSet<String>,
    params: &JunctionRefinementParams,
) -> Result<(i64, i64)> {
    if read_names.is_empty() || end <= start {
        return Ok((start, end));
    }

    let mut reader = bam::io::indexed_reader::Builder::default()
        .build_from_path(bam_path)?;
    let header = reader.read_header()?;

    let region_start = Position::try_from(start.max(1) as usize)?;
    let region_end = Position::try_from(end.max(1) as usize)?;
    let region = noodles::core::Region::new(chrom, region_start..=region_end);

    let mut clipped_blocks: Vec<(i64, i64)> = Vec::new();
    let mut aln_starts: Vec<i64> = Vec::new();
    let mut aln_ends: Vec<i64> = Vec::new();
    let mut acceptor_hits: Vec<(i64, String)> = Vec::new();

    for result in reader.query(&header, &region)? {
        let record = result?;
        if record.flags().is_unmapped() {
            continue;
        }
        if record.flags().is_secondary() {
            continue;
        }
        if !params.include_supplementary && record.flags().is_supplementary() {
            continue;
        }

        let name = match record.name() {
            Some(n) => match std::str::from_utf8(n) {
                Ok(s) => s,
                Err(_) => continue,
            },
            None => continue,
        };

        if !read_names.contains(name) {
            continue;
        }

        if let Some(rs) = record.alignment_start().transpose().ok().flatten() {
            aln_starts.push(rs.get() as i64);
        }
        if let Some(re) = record.alignment_end().transpose().ok().flatten() {
            aln_ends.push(re.get() as i64);
        }

        if let Some(blocks) = exon_blocks_from_record(&record) {
            for b in blocks {
                if let Some(c) = clip_block_to_region(b, start, end) {
                    clipped_blocks.push(c);
                }
            }
        }

        if let Some(acc) = acceptor_sites_from_record(&record) {
            for p in acc {
                if p >= start && p <= end {
                    acceptor_hits.push((p, name.to_string()));
                }
            }
        }
    }

    let merged = merge_intervals(clipped_blocks);
    if merged.is_empty() {
        return Ok((start, end));
    }

    let u_start = merged.iter().map(|x| x.0).min().unwrap_or(start);
    let u_end = merged.iter().map(|x| x.1).max().unwrap_or(end);

    let span = end - start;
    let max_trim = max(
        1,
        ((span as f64) * params.max_trim_frac_per_side).ceil() as i64,
    );

    let left_delta = (u_start - start).max(0).min(max_trim);
    let right_delta = (end - u_end).max(0).min(max_trim);

    let mut new_start = start + left_delta;
    let mut new_end = end - right_delta;

    if new_end - new_start < params.min_segment_span {
        return Ok((start, end));
    }
    if new_start >= new_end {
        return Ok((start, end));
    }

    let (ns0, ne0) = (new_start, new_end);

    if params.read_end_guard_slack_bp > 0 && !aln_starts.is_empty() && !aln_ends.is_empty() {
        let (a, b) = apply_read_end_guards(
            new_start,
            new_end,
            &aln_starts,
            &aln_ends,
            params.read_end_guard_slack_bp,
        );
        new_start = a;
        new_end = b;
    }

    new_start = apply_acceptor_guard(
        new_start,
        u_start,
        end,
        acceptor_hits,
        params.min_junction_read_support,
        params.junction_cluster_slack_bp,
        params.acceptor_window_left_bp,
        params.acceptor_window_right_bp,
    );

    if new_end - new_start < params.min_segment_span || new_start >= new_end {
        return Ok((ns0, ne0));
    }

    Ok((new_start, new_end))
}
