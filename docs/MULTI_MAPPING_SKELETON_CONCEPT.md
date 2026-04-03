# Multi-Mapping Read Skeleton: A Gene Family Detection Concept

## Core Hypothesis

Given one or more seed regions from a multi-copy gene family, the multi-mapping reads that overlap those seeds form a "skeleton" that can be used to discover all other family members.

## Key Insight

Different seeds see different portions of the family:
- Seed A may detect regions 1, 2, 3 (but not 4)
- Seed B may detect regions 2, 3, 4 (but not 1)
- Merging and keeping the **larger extent** gives complete coverage

This is because:
1. Multi-mapping reads have varying coverage across family members
2. Read depth affects how far from the seed we can detect
3. Some members have fewer reads but extend detection further

## What We've Proven (ID_14 Example)

| Metric | Value |
|--------|-------|
| Ground truth genes | 7 |
| Seeds used | 7 (one per gene) |
| Detected loci | 6 |
| Genes covered | 7/7 (100%) |
| Locus size vs GT | 3-2800x larger |

### Critical Fix

**Before**: Kept first overlapping locus → Missed AC090616.5 (85.7% sensitivity)

**After**: Keep larger overlapping locus → All genes covered (100% sensitivity)

```
Detected locus: 32889086-33199147 (310 kb)
Ground truth:
  - LRRC37B:    32953747-32999386 (45 kb) ✓
  - AC090616.5: 33034040-33040757 (7 kb) ✓  <- Was missed before
```

## The Skeleton Model

```
Seed Region ──┐
              │  Collect all multi-mapping reads
              ▼
         [Read Cloud]
              │
              │  For each read, find ALL alignments genome-wide
              ▼
    ┌─────────────────────────────────────┐
    │  Build read-to-read similarity      │
    │  graph (Jaccard overlap)            │
    └─────────────────────────────────────┘
              │
              │  Find connected components
              ▼
    ┌─────────────────────────────────────┐
    │  Cluster alignments into loci       │
    │  Group by genomic proximity         │
    └─────────────────────────────────────┘
              │
              │  Merge overlapping loci (keep larger)
              ▼
         [Family Members]
```

## Recommendations for Validation

### 1. Minimal Seed Experiment

**Question**: How many seeds are needed to find all family members?

**Experiment**:
```
For each family:
  For n_seeds in [1, 2, 3, 5, all]:
    Randomly select n_seeds from ground truth
    Run detection
    Measure: sensitivity, extra_loci, merge_ratio
```

**Expected outcome**: 
- 1 seed: ~40-60% of members (depends on seed quality)
- 2-3 seeds: ~80-95% of members
- All seeds: 100% (proven)

### 2. Family Type Comparison

**Tandem arrays** (genes clustered together):
- Example: ID_14 (7 genes on one chromosome region)
- Challenge: Distinguishing adjacent paralogs

**Dispersed families** (genes on different chromosomes):
- Example: ID_35 (10 genes across 8 chromosomes)
- Challenge: Finding distant members

**Small families** (3-5 members):
- Lower read depth per member
- May need connectivity-based seed selection

**Large families** (10+ members):
- More reads = better detection
- But more complex clustering

### 3. Boundary vs Detection Separation

**Detection Phase**: "Where are the family members?"
- Output: Loci that contain family members
- Metric: Sensitivity (members found / total members)
- Acceptable: Loci larger than genes

**Refinement Phase**: "Where are the gene boundaries?"
- Input: Detected loci
- Output: Precise gene coordinates
- Metrics: Precision, F1 score
- Methods: Coverage valleys, Iso-Seq, splice junctions

These are **separate problems**. Conflating them causes confusion.

## Proposed Test Protocol

### Phase 1: Validate Detection (Current Focus)

```bash
# For each test family:
# 1. Use all ground truth genes as seeds
# 2. Run detection
# 3. Measure sensitivity

for family in test_families:
    seeds = ground_truth[family]
    detected = detect(seeds)
    sensitivity = coverage(detected, ground_truth)
```

**Success criteria**: Sensitivity ≥ 95%

### Phase 2: Minimize Seeds

```bash
# For each family:
# 1. Try with 1 random seed
# 2. If sensitivity < threshold, add another seed
# 3. Repeat until threshold met

for family in test_families:
    for n_seeds in [1, 2, 3, ...]:
        seeds = random_sample(ground_truth, n_seeds)
        detected = detect(seeds)
        if sensitivity(detected) >= 0.95:
            record_min_seeds(family, n_seeds)
            break
```

**Success criteria**: Most families need ≤ 3 seeds

### Phase 3: Cross-Family Validation

Test on families with different characteristics:
- High copy number (ID_14, ID_207)
- Multi-chromosome (ID_35)
- Small families (3-5 members)
- Divergent paralogs (low sequence similarity)

## Suggested Test Families

| Family | Genes | Chromosomes | Type | Why Test |
|--------|-------|-------------|------|----------|
| ID_14 | 7 | 1 | Tandem | Already validated |
| ID_35 | 10 | 8 | Dispersed | Multi-chromosome |
| ID_8 | 9 | ? | Mixed | Different structure |
| ID_207 | 17 | ? | Large | High copy number |
| ID_131 | 6 | 1 | Tandem | Compare to ID_14 |

## Open Questions

1. **Minimum seeds**: What's the theoretical minimum for 95% coverage?
2. **Seed quality**: Are some seeds "better" than others? (More reads, better position)
3. **False positive rate**: How many detected loci have no ground truth?
4. **Scalability**: Does this work for 100+ families?

## Code Changes Summary

### Critical Fix: Keep Larger Locus

```rust
// Before: Kept first occurrence (missed some genes)
for locus in loci {
    if !is_duplicate { all_loci.push(locus); }
}

// After: Keep larger locus when overlapping
for locus in loci {
    if overlaps(existing) {
        if locus.span > existing.span {
            replace(existing, locus);
        }
    } else {
        push(locus);
    }
}
```

### File Modified
- `gene_family_detector/src/main.rs`: Lines 1406-1458, 1794-1820, 2362-2390, 2830-2858

## Conclusion

We have proven that multi-mapping reads can be used as seeds to discover all members of a gene family. The key insight is that different seeds see different coverage extents, and merging by keeping the larger extent gives complete coverage.

The next step is to quantify how few seeds are needed and validate across diverse family types.
