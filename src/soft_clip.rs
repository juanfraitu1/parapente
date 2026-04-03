//! Shared 3' / 5' soft-clip handling (CIGAR + query coordinates).
//! Plus strand: 3' clip is usually the last `S` op; if the last op is not `S`, we fall back to the
//! rightmost soft-clip block in query order (handles `100M 50S 10N 20M` style intronic tails).
//! Minus strand: 3' clip is only the leading `S` (first op), matching read_cloud / Iso-Seq convention.

use noodles::bam;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::cigar::Op;

fn query_bases_consumed(op: &Op) -> usize {
    if op.kind().consumes_read() {
        op.len()
    } else {
        0
    }
}

fn query_len_before_op_index(ops: &[Op], idx: usize) -> usize {
    let mut q = 0usize;
    for op in ops.iter().take(idx) {
        q += query_bases_consumed(op);
    }
    q
}

/// All soft-clip intervals in query order (5' to 3' along SEQ).
fn all_soft_clip_segments(ops: &[Op]) -> Vec<(usize, usize)> {
    let mut out = Vec::new();
    let mut q = 0usize;
    for op in ops.iter() {
        if op.kind() == Kind::SoftClip {
            let len = op.len();
            out.push((q, q + len));
        }
        q += query_bases_consumed(op);
    }
    out
}

/// Query interval `[start, end)` for 3' soft clip (polyA proxy). `reverse` is SAM reverse flag.
pub fn three_prime_soft_clip_from_ops(ops: &[Op], reverse: bool) -> Option<(usize, usize)> {
    if ops.is_empty() {
        return None;
    }
    if reverse {
        if ops.first().map(|o| o.kind()) == Some(Kind::SoftClip) {
            let len = ops.first().unwrap().len();
            return Some((0, len));
        }
        return None;
    }

    if ops.last().map(|o| o.kind()) == Some(Kind::SoftClip) {
        let i = ops.len() - 1;
        let q = query_len_before_op_index(ops, i);
        let len = ops[i].len();
        return Some((q, q + len));
    }

    let segs = all_soft_clip_segments(ops);
    if segs.is_empty() {
        return None;
    }
    let best = *segs.iter().max_by_key(|(s, _)| *s)?;
    if segs.len() == 1 && best.0 == 0 {
        return None;
    }
    Some(best)
}

pub fn three_prime_soft_clip_query_range(record: &bam::Record) -> Option<(usize, usize)> {
    let ops: Vec<_> = record.cigar().iter().filter_map(|r| r.ok()).collect();
    if ops.is_empty() {
        return None;
    }
    let reverse = record.flags().is_reverse_complemented();
    three_prime_soft_clip_from_ops(&ops, reverse)
}

/// Returns (5' soft-clip length, 3' soft-clip length) in read orientation (matches read_cloud mapping).
pub fn extract_soft_clip_5p_3p_lens(record: &bam::Record) -> (i64, i64) {
    let ops: Vec<_> = record.cigar().iter().filter_map(|r| r.ok()).collect();
    if ops.is_empty() {
        return (0, 0);
    }
    let reverse = record.flags().is_reverse_complemented();
    if !reverse {
        let five = if ops.first().map(|o| o.kind()) == Some(Kind::SoftClip) {
            ops.first().unwrap().len() as i64
        } else {
            0
        };
        let three = three_prime_soft_clip_from_ops(&ops, false)
            .map(|(a, b)| (b - a) as i64)
            .unwrap_or(0);
        (five, three)
    } else {
        let five = if ops.last().map(|o| o.kind()) == Some(Kind::SoftClip) {
            ops.last().unwrap().len() as i64
        } else {
            0
        };
        let three = if ops.first().map(|o| o.kind()) == Some(Kind::SoftClip) {
            ops.first().unwrap().len() as i64
        } else {
            0
        };
        (five, three)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn plus_5s_100m_rejects_lone_5prime_clip() {
        let ops = vec![Op::new(Kind::SoftClip, 5), Op::new(Kind::Match, 100)];
        assert!(three_prime_soft_clip_from_ops(&ops, false).is_none());
    }

    #[test]
    fn plus_100m_5s_trailing_strict() {
        let ops = vec![Op::new(Kind::Match, 100), Op::new(Kind::SoftClip, 5)];
        let r = three_prime_soft_clip_from_ops(&ops, false).unwrap();
        assert_eq!(r, (100, 105));
    }

    #[test]
    fn plus_100m_50s_skip_20m_uses_rightmost_s_before_tail() {
        let ops = vec![
            Op::new(Kind::Match, 100),
            Op::new(Kind::SoftClip, 50),
            Op::new(Kind::Skip, 10),
            Op::new(Kind::Match, 20),
        ];
        let r = three_prime_soft_clip_from_ops(&ops, false).unwrap();
        assert_eq!(r, (100, 150));
    }

    #[test]
    fn plus_5s_100m_3s_last_op_is_3prime() {
        let ops = vec![
            Op::new(Kind::SoftClip, 5),
            Op::new(Kind::Match, 100),
            Op::new(Kind::SoftClip, 3),
        ];
        let r = three_prime_soft_clip_from_ops(&ops, false).unwrap();
        assert_eq!(r, (105, 108));
    }

    #[test]
    fn minus_leading_s_only() {
        let ops = vec![Op::new(Kind::SoftClip, 12), Op::new(Kind::Match, 80)];
        let r = three_prime_soft_clip_from_ops(&ops, true).unwrap();
        assert_eq!(r, (0, 12));
    }
}
