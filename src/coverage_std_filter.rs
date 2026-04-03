//! Per-locus binned coverage standard deviation (legacy Python `calculate_coverage_profile_stats` / `std`).
//! Used to filter segmental-duplication false positives for dispersed families.

use anyhow::Result;
use noodles::bam;
use noodles::core::Position;
use noodles::sam::alignment::Record;

/// Mean and population standard deviation of per-bin read counts (matches NumPy `std`, ddof=0).
pub fn coverage_mean_and_std_in_region(
    bam_path: &str,
    chrom: &str,
    start: i64,
    end: i64,
    bin_size: i64,
    include_supplementary: bool,
) -> Result<(f64, f64)> {
    let bs = bin_size.max(1);
    let span = (end - start).max(0);
    let num_bins = ((span / bs) + 1).max(1) as usize;

    let mut reader = bam::io::indexed_reader::Builder::default().build_from_path(bam_path)?;
    let header = reader.read_header()?;

    let region_start = Position::try_from(start.max(1) as usize)?;
    let region_end = Position::try_from(end.max(1) as usize)?;
    let region = noodles::core::Region::new(chrom, region_start..=region_end);

    let mut bins = vec![0i64; num_bins];

    for result in reader.query(&header, &region)? {
        let record = result?;
        if record.flags().is_unmapped() {
            continue;
        }
        if !include_supplementary && record.flags().is_supplementary() {
            continue;
        }
        let r_start = record
            .alignment_start()
            .transpose()?
            .map(|p: Position| p.get() as i64)
            .unwrap_or(1)
            - 1;
        let r_end = record
            .alignment_end()
            .transpose()?
            .map(|p: Position| p.get() as i64)
            .unwrap_or(r_start + 1);

        let read_start = r_start.max(start);
        let read_end = r_end.min(end);
        if read_end <= read_start {
            continue;
        }
        let start_bin = ((read_start - start) / bs).clamp(0, num_bins as i64 - 1) as usize;
        let end_bin = ((read_end - start) / bs).clamp(0, num_bins as i64 - 1) as usize;
        for b in start_bin..=end_bin.min(num_bins.saturating_sub(1)) {
            bins[b] += 1;
        }
    }

    let n = num_bins as f64;
    let sum: f64 = bins.iter().map(|&c| c as f64).sum();
    let mean = sum / n;
    let var = bins
        .iter()
        .map(|&c| {
            let d = c as f64 - mean;
            d * d
        })
        .sum::<f64>()
        / n;
    Ok((mean, var.sqrt()))
}

/// Population standard deviation of per-bin read counts (matches NumPy `std` with default ddof=0).
#[allow(dead_code)]
pub fn coverage_std_in_region(
    bam_path: &str,
    chrom: &str,
    start: i64,
    end: i64,
    bin_size: i64,
    include_supplementary: bool,
) -> Result<f64> {
    Ok(coverage_mean_and_std_in_region(
        bam_path,
        chrom,
        start,
        end,
        bin_size,
        include_supplementary,
    )?
    .1)
}

/// Scale a legacy minimum std (e.g. 50) down when mean bin depth is below the reference.
/// Capped at the baseline when mean is at or above reference so high-coverage windows are not
/// penalized with a higher bar than the legacy absolute value.
pub fn depth_scaled_min_std(base_min_std: f64, mean_bin: f64, mean_ref: f64) -> f64 {
    let r = mean_ref.max(1.0);
    let scale = ((mean_bin + 1.0) / r).sqrt().min(1.0);
    base_min_std * scale
}

/// Adaptive minimum std threshold (legacy Python `adaptive_std_threshold`).
pub fn adaptive_std_threshold(
    base_min_std: f64,
    locus_span: i64,
    shared_with_seed: usize,
    min_floor: f64,
    shared_ref: f64,
    span_ref: f64,
) -> f64 {
    let sr_scale = ((shared_with_seed.max(1) as f64) / shared_ref).sqrt();
    let span_scale = (locus_span.max(1000) as f64 / span_ref).sqrt();
    let mut scaled = base_min_std * (sr_scale * span_scale).min(1.0);
    const TINY_SHARED_MAX: f64 = 25.0;
    const TINY_SHARED_RELAX: f64 = 0.7;
    if (shared_with_seed as f64) <= TINY_SHARED_MAX {
        scaled *= TINY_SHARED_RELAX;
    }
    scaled.max(min_floor)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn adaptive_matches_python_example() {
        let t = adaptive_std_threshold(50.0, 150_000, 500, 8.0, 500.0, 150_000.0);
        assert!((t - 50.0).abs() < 1e-6);
        let t2 = adaptive_std_threshold(50.0, 50_000, 50, 8.0, 500.0, 150_000.0);
        assert!(t2 >= 8.0 && t2 < 50.0);
    }

    #[test]
    fn depth_scaled_matches_reference_mean() {
        let s = depth_scaled_min_std(50.0, 200.0, 200.0);
        assert!((s - 50.0).abs() < 1e-6);
        let s_hi = depth_scaled_min_std(50.0, 400.0, 200.0);
        assert!((s_hi - 50.0).abs() < 1e-6);
        let s2 = depth_scaled_min_std(50.0, 50.0, 200.0);
        assert!((s2 - 50.0 * (51.0f64 / 200.0).sqrt()).abs() < 1e-6);
    }
}
