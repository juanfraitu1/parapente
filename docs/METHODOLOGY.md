# Methodology

## Overview

Parapente discovers paralogous gene families by analyzing multi-mapping RNA reads. This document describes the algorithm in detail.

## Core Concept

### The Multi-Mapping Signal

When a gene or lncRNA is duplicated, the copies (paralogs) share high sequence similarity. RNA-seq reads from these paralogs will:
1. Map uniquely to their source gene (primary alignment)
2. **Map to multiple paralogs** (secondary/supplementary alignments)

These "multi-mapping" reads are typically discarded in standard RNA-seq analysis. Parapente treats them as a **signature of paralogy**.

## Algorithm Details

### Step 1: Seed Read Collection

For each seed region:
- Query BAM for all alignments overlapping the seed
- Collect read names of all mapped reads
- Include both primary and supplementary alignments

```rust
seed_reads = query_bam(chrom, start, end)
```

### Step 2: Genome-Wide Alignment Discovery

Find **all** alignments of seed reads across the genome:
- Scan only chromosomes containing seed regions
- Build lightweight index: `read_name → [alignments]`
- Result: Set of all genomic locations touched by seed reads

### Step 3: Locus Clustering

Cluster alignments into candidate loci:
- Sort alignments by chromosome and position
- Merge alignments within `cluster_distance` (default: 50kb)
- Adaptive gap: Use 90th percentile of same-read inter-alignment gaps

```
Before clustering:  ████    ████    ████    ████    ████
After clustering:   ████████████████    ████████████████
```

### Step 4: Coverage Valley Splitting

Large merged regions often contain multiple genes. Split using coverage valleys:

1. **Build coverage profile**: Bin read coverage across locus
2. **Detect valleys**: Find local minima in coverage
3. **Split at valleys**: Separate regions at significant valleys
4. **Merge small segments**: Combine short segments (< 5kb)

```
Coverage:    /
            / \
           /   \    /\
          /     \  /  \
         /   ↓   \/    \
        /  (valley)     \
       /_________________\
       Split here ↑
```

### Step 5: Seed Preservation

Critical bug fix: Ensure all seed regions appear in output:
- Check if each seed overlaps any output locus
- If missing (lost during splitting), add back with its reads

### Step 6: Filtering

Apply filters to remove artifacts:

| Filter | Purpose | Default |
|--------|---------|---------|
| Coverage std | Remove flat/noisy coverage | 50 |
| Jaccard | Remove low-similarity hits | 0.3 |
| Multi-seed | Require support from N seeds | 1 |
| Core reads | Require shared reads across seeds | 1 |
| Coverage quality | Gene-like coverage pattern | 0.0 (off) |

## Multi-Seed Mode

For families with multiple known members:

```
Seed 1:  AAAAAA          (reads: R1, R2, R3)
Seed 2:      BBBBBB      (reads: R2, R3, R4)
Seed 3:          CCCCCC  (reads: R3, R4, R5)
                ↑
         R3 is "core read" (present in all seeds)
```

Loci sharing **core reads** are more likely to be true paralogs.

## Transitive Discovery

Find distant paralogs through **transitive connectivity**:

```
Seed reads → Locus A (jaccard=0.9)
                ↓
         Shared reads → Locus B (jaccard=0.1, but shares reads with A)
                ↓
         Shared reads → Locus C (jaccard=0.05, but shares reads with B)
```

Even with low direct Jaccard to seed, transitive loci can be discovered.

## Key Parameters

### Critical Parameters

| Parameter | Description | Recommended |
|-----------|-------------|-------------|
| `--cluster-distance` | Max gap for merging alignments | 50000 |
| `--min-std-coverage` | Min coverage variation | 50 |
| `--min-jaccard-locus` | Min similarity to seed | 0.3 |
| `--min-seeds` | Min seeds supporting locus | 2 |

### Coverage Valley Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--valley-frac` | Valley depth threshold | 0.08 |
| `--valley-min-segment-bp` | Min segment after split | 1500 |
| `--valley-merge-min-segment-bp` | Min segment to keep | 5000 |

## Complexity Analysis

- **Time**: O(R + A log A) where R = seed reads, A = alignments
- **Space**: O(A) for alignment index
- **Parallel**: Embarrassingly parallel per-locus processing

## Validation

Tested on 83 human gene families with segmental duplications:
- 100% sensitivity (all 362 genes found)
- Handles families with 99%+ sequence identity
- Robust to collapsed assemblies
