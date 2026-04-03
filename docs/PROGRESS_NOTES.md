# Gene Family Detection Progress Notes

**Date:** 2026-04-01
**Session Summary:** Testing and improving gene family detection for 83 families

## Current Results (27/83 families processed)

| Metric | Value |
|--------|-------|
| Total GT copies | 138 |
| Found (≥50% coverage) | 128 |
| Average Sensitivity | 90.9% (IoU) / 84.9% (50% cov) |
| Average Specificity | 11.5% |
| Main Issue | Over-merged loci (many false positives) |

## Key Discovery: Read Identity Distinguishes Adjacent Genes

**Problem:** Coverage-based splitting doesn't work because multi-mapping reads create continuous coverage across adjacent genes.

**Solution:** Read identity changes sharply between genes:

| Gene Pair | Distance | Shared Reads | Jaccard | Interpretation |
|-----------|----------|--------------|---------|----------------|
| LRRC37A vs LRRC37A2 | 174 kb | 4.1% | 0.014 | **Distinct genes** |
| LRRC37B vs AC090616.5 | 35 kb | 0.6% | 0.006 | **Distinct genes** |
| LRRC37A vs LRRC37A3 | 18.5 Mb | 0% | 0.000 | **Distant genes** |
| AC005562.2 vs LRRC37BP1 | 12 kb | 99% | 0.616 | **Same locus** |

**Key insight:** Adjacent gene family members share FEW reads (1-13%), while distant genes share 0%. This can be used for locus splitting.

## Recommended Next Step: Implement Read-Identity Clustering

Algorithm:
```
1. For each detected locus:
   a. Collect all reads mapped to it
   b. Bin the locus into windows (e.g., 10kb)
   c. For each bin, compute read set identity
   d. Find positions where Jaccard between adjacent bins < 0.3
   e. These are gene boundaries
```

## Missed Genes Analysis

### Category 1: Small genes (<1kb) absorbed into larger loci
- ID_14: AC005562.2 (112 bp) - 100% covered but IoU=0.0006
- ID_131: AC105272.1 (200 bp) - inside larger locus
- ID_213: AC239809.1 (218 bp) - inside larger locus
- ID_222: AC243829.6 (107 bp) - inside larger locus

**Status:** These ARE detected (inside larger loci) but have low IoU due to size difference.

### Category 2: Large genes partially detected (30-50% coverage)
- ID_163: GUSBP1 (346 kb) - only 130 kb detected (38%)
- ID_188: HERC2 (211 kb) - only 63 kb detected (34%)
- ID_215: PDE4DIP (239 kb) - only 137 kb detected (32%)

**Root cause:** Over-aggressive splitting of large genes.

### Category 3: Different chromosomes (transitive discovery working)
- ID_163, ID_175, ID_182, ID_24: Genes on multiple chromosomes - transitive discovery IS working
- Genes detected on all chromosomes, but some have partial coverage

### Category 4: Isolated genes (no shared reads)
- ID_213: AC245407.1 (2.4 kb) - 0% coverage, isolated gene with no shared reads

## Files Modified

### Rust Source Files (in /scratch/jxi21/segments/gene_family_detector/src/)

1. **main.rs** - Main entry point
2. **transitive_detector.rs** - Transitive detection with coverage splitting
   - Added `peak_count` and `confidence` fields to TransitiveLocus
   - Coverage valley splitting with dynamic thresholds
   - Junction-based splitting for large loci
3. **coverage_valley.rs** - Valley detection and splitting
   - `split_by_junctions()` function for splice-based splitting
   - `build_coverage_profile()` for coverage analysis
4. **smart_seed_selector.rs** - Connectivity-based seed selection
5. **boundary_refinement.rs** - Locus boundary refinement

### Commits Made
```
01391c5 - Add peak_count and confidence fields to TransitiveLocus struct
c52913d - Add coverage_validator module
3cd23aa - Improved junction splitting (seed reads + 20kb gaps)
7e8170a - Add merge_small_segments
888c98b - Make coverage splitting mandatory with dynamic thresholds
```

## Seed Files Created

83 seed files created in `/scratch/jxi21/segments/rust_multiseed_test/`:
- `ID_14_seeds.bed` through `ID_92_seeds.bed`

## Output Files

- `/scratch/jxi21/segments/rust_multiseed_test/results/` - Detection results (27 families)
- `/scratch/jxi21/segments/rust_multiseed_test/ID_14_test.bed` - Latest test output

## Commands to Resume

```bash
# Continue batch evaluation
cd /scratch/jxi21/segments/gene_family_detector
python3 run_all_families.py

# Summarize completed results
python3 /tmp/summarize_results.py

# Analyze specific family
./target/release/gene_family_detector \
  --bam /scratch/jxi21/segments/A119b.bam \
  --seeds-bed /scratch/jxi21/segments/rust_multiseed_test/ID_14_seeds.bed \
  --seed-chrom NC_060941.1 --seed-start 31562334 --seed-end 31562446 \
  --output test_output.bed --transitive --threads 4
```

## Outstanding Tasks

1. **Implement read-identity clustering for locus splitting**
   - This is the most important improvement
   - Will correctly separate adjacent genes with low read overlap
   - Will preserve genes with high read overlap

2. **Run full batch on all 83 families**
   - Currently 27/83 complete
   - Need to continue and complete

3. **Improve specificity**
   - Current specificity is very low (~11%)
   - Too many false positive segments
   - Read-identity splitting should help

4. **Handle very small genes (<500 bp)**
   - Currently absorbed into larger loci
   - May need special handling or post-processing

## Technical Details

### Coverage Splitting Parameters (current)
```rust
valley_frac: 0.03-0.12  // Dynamic based on max coverage
min_prominence_frac: 0.20-0.35  // Dynamic based on variability
min_gap_bp: 2000  // Minimum gap between segments
min_size_for_junction_split: 200_000  // Only split loci >200kb
min_junction_gap: 20_000  // 20kb gaps for junction clusters
```

### Read Overlap Analysis (from testing)
```
Adjacent gene pairs share 0.6-13% of reads (Jaccard 0.006-0.014)
Distant genes share 0% reads
Very close genes share 60-99% reads (same transcription unit)
```

This suggests a Jaccard threshold of 0.3 for splitting decisions.
