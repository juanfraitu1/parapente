# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a **gene family discovery pipeline** using multi-mapping RNA-seq reads (PacBio CCS long reads) to detect duplicated gene families in the human genome (T2T-CHM13 v2.0 assembly). The core innovation is **bridge read detection** - discovering disconnected gene family components through transitive read connections.

**Key Achievement**: 100% precision gene family discovery using multi-seed approach with coverage-based filtering.

**Validated on**: Both primate-specific families (SRGAP2, AMY) AND classic textbook families (globins, tubulins).

---

## IMPORTANT: Writing Style

**NEVER use emojis or em-dashes in any code, documentation, or outputs.** Use simple text-based formatting only (hyphens, asterisks, plain text). This applies to all files created or modified.

---

## Common Commands

### Run Gene Family Detection

```bash
# Basic usage (compact families)
python3 shared_reads_detector.py \
    --bam HSA_15.bam \
    --seed_chrom NC_060925.1 \
    --seed_start 121194173 \
    --seed_end 121402237 \
    --output results.bed

# With coverage filtering (dispersed families - recommended)
python3 shared_reads_detector.py \
    --bam HSA_15.bam \
    --seed_chrom NC_060925.1 \
    --seed_start 121194173 \
    --seed_end 121402237 \
    --output results.bed \
    --min_std_coverage 50

# Restrict to specific chromosomes (faster)
python3 shared_reads_detector.py \
    --bam HSA_15.bam \
    --seed_chrom NC_060925.1 \
    --seed_start 121194173 \
    --seed_end 121402237 \
    --output results.bed \
    --chromosomes NC_060925.1,NC_060940.1
```

### Multi-Seed Workflow (Recommended)

```bash
# Run from 2+ independent seeds, then merge
python3 shared_reads_detector.py --bam HSA_15.bam --seed_chrom NC_060925.1 --seed_start 121194173 --seed_end 121402237 --output seed1.bed --min_std_coverage 50
python3 shared_reads_detector.py --bam HSA_15.bam --seed_chrom NC_060925.1 --seed_start 143047107 --seed_end 143140349 --output seed2.bed --min_std_coverage 50

# Merge results
bedtools merge -i <(cat seed1.bed seed2.bed | grep -v "^#" | sort -k1,1 -k2,2n) > final_merged.bed

# Validate against ground truth
bedtools intersect -a final_merged.bed -b ground_truth_80families.bed -wa -wb > overlaps.txt
```

### Generate Figures

```bash
# Create presentation figures
python3 create_presentation_figures.py

# Plot overlap analysis
python3 plot_overlap_impact.py

# Create comparison plots
python3 create_comparison_plot.py
```

---

## High-Level Architecture

### Core Algorithm: Component 1 Discovery

The validated method uses **only Component 1 discovery** (direct connections to seed). Transitive discovery (Components 2+) accumulates false positives exponentially and should NOT be used for production.

```
Iteration 0 (Seed):
  - seed_reads = collect_reads_from_region(seed_coords)

Iteration 1 (Component 1 - USE THIS ONLY):
  - loci = find_where_reads_map(seed_reads)
  - cluster and filter by jaccard(locus, seed) >= threshold
  - Apply coverage filter (std_coverage >= 50) for dispersed families
  - Result: 88-100% precision

Iteration 2+ (Components 2+ - DO NOT USE):
  - pool_reads = collect_reads_from_all_discovered_loci()
  - new_loci = find_where_reads_map(pool_reads)
  - Problem: Accumulates false positives, precision drops to 17-28%
```

### Key Architectural Insight

**Multi-seed beats transitive discovery:**
- Instead of one seed + many iterations (fails)
- Use multiple seeds + one iteration each (succeeds)
- Pool results and merge
- Achieves 100% precision with 100% recall

### Rust Implementation (Production)

A high-performance Rust implementation is available in `gene_family_detector/`:

**Key Features:**
- Multi-threaded BAM scanning using `rayon`
- Memory-efficient with `FxHashSet` for read tracking
- Integrated coverage valley splitting with Iso-Seq validation
- Multi-seed support (BED file or CLI list)
- Four-layer validation framework

**Usage:**
```bash
# Build release binary
cd gene_family_detector && cargo build --release

# Run with integrated validation
./target/release/gene_family_detector \
    --bam A119b.bam \
    --seeds-bed seeds.bed \
    --output results.bed \
    --integrated \
    --valley-frac 0.1 \
    --min-isoseq-sim 0.6
```

**Modules:**
- `coverage_valley.rs` - Valley detection and splitting algorithm
- `coverage_validation.rs` - Peak analysis, symmetry, slope metrics
- `integrated_detector.rs` - Four-layer validation framework
- `isoseq.rs` - Iso-Seq transcript validation
- `read_cloud.rs` - Dense multi-mapping read cloud detection

### Python Implementation (Development/Analysis)

**Production Script:**
- `shared_reads_detector.py` - Main detector with integrated coverage filtering

**Key Functions:**
- `collect_reads_from_region()` - Extract read names from seed region
- `find_where_reads_map()` - Find all loci where seed reads align
- `cluster_alignments()` - Merge nearby alignments into loci
- `calculate_coverage_metrics()` - Compute coverage statistics for filtering
- `filter_by_coverage()` - Apply std_coverage threshold

**Visualization Scripts:**
- `create_presentation_figures.py` - Generate publication-ready figures
- `plot_overlap_impact.py` - Analyze overlap threshold effects
- `create_comparison_plot.py` - Compare methods

---

## Data Files

### Input Data

| File | Description | Location |
|------|-------------|----------|
| `HSA_15.bam` | PacBio CCS long-read RNA-seq alignments (36 GB) | `/storage/group/kdm16/default/jxi21/apes_transcriptome_analysis/Obj3A/RNA_based_method/` |
| `HSA_15.bam.bai` | BAM index | Same as above |
| `1000_gen_HSA.bam` | 1000 Genomes IsoSeq data (4.1 GB) | Same as above |
| `ground_truth_80families.bed` | Ground truth annotations for 80 gene families | Repository root |

### Reference Files

| File | Description |
|------|-------------|
| `T2T_V2_HSA.fasta` | T2T-CHM13 v2.0 reference genome |
| `GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz` | Gene annotations |

### Output Formats

**BED format** (detector output):
```
chrom  start  end  name  jaccard_seed  component  shared_with_seed  total_reads  span
```

**TSV format** (coverage metrics):
```
locus_id  chrom  start  end  mean_coverage  median_coverage  max_coverage  std_coverage  coef_variation
```

---

## Gene Family Types

### Dispersed Families (e.g., SRGAP2, span >50 Mb)

- Genes scattered across genome
- Low pairwise Jaccard indices (<0.1 common)
- **MUST use coverage filtering** (std_coverage >= 50)
- Component 1 precision: ~30% without filtering, 100% with filtering

### Compact Families (e.g., AMY, span <1 Mb)

- Genes clustered in tandem
- High pairwise Jaccard indices (>0.5)
- Coverage filtering not needed (merging handles overlaps)
- Component 1 precision: 70-100% naturally

---

## Validated Method (Current State - January 2026)

### The Only Approach That Achieves 100% Precision

**For dispersed families:**
1. Run component 1 from each seed with coverage filtering
2. Merge all results with bedtools merge
3. Result: 100% precision, 100% recall

**For compact families:**
1. Run component 1 from each seed
2. Merge all results
3. Result: 100% precision, 100% recall

### Why This Works

**Coverage filtering (std_coverage >= 50):**
- True genes have gene-like expression patterns (exons, introns, splicing) -> high variability
- Segmental duplications/pseudogenes are non-transcribed or aberrant -> low variability
- Uses **functional genomic context**, not sequence identity

**Multi-seed approach:**
- Avoids transitive discovery (no false positive cascade)
- Each seed finds its "neighborhood" independently
- Pooling results increases sensitivity while maintaining precision

---

## Failed Approaches (DO NOT REIMPLEMENT)

These have been thoroughly tested and documented as failures:

1. **Alignment quality filtering** (ALIGNMENT_QUALITY_FILTERING_TEST.md)
   - Tested: Filter reads by alignment fraction (>=80%, >=90%)
   - Result: Zero improvement (29.8% vs 31.9% precision)
   - Why it fails: Segmental duplications have high-quality alignments (>90%)

2. **Consensus voting** (CONSENSUS_WITH_FILTERED_RESULTS.md)
   - Tested: Require loci connect to multiple reference family members
   - Result: Precision dropped to 17.0% (from 88.2%)
   - Why it fails: Dispersed families have low pairwise Jaccard indices

3. **Adaptive thresholds** (ADAPTIVE_THRESHOLDS_TEST_RESULTS.md)
   - Tested: Iteration-dependent Jaccard thresholds
   - Result: Precision dropped to 27.7%
   - Why it fails: Cannot compensate for biological false positives

4. **Soft clipping/strand filtering** (STRAND_AND_OVERLAP_ANALYSIS.md)
   - Tested: Filter by soft clipping amount or strand consistency
   - Result: No discriminative power
   - Why it fails: True family members can be on different strands; segmental duplications align well

5. **Transitive discovery beyond component 1**
   - Components 2+ accumulate false positives exponentially
   - Only use max_iterations=1 (component 1 only)

---

## Key Parameters

| Parameter | Default | Recommended | Notes |
|-----------|---------|-------------|-------|
| `--min_jaccard` | 0.001 | 0.001 | Lower = more recall, less precision |
| `--max_iterations` | 1 | **1** | CRITICAL: Do not increase |
| `--min_std_coverage` | None | 50 | Required for dispersed families |
| `--cluster_distance` | 10000 | 10000 | Merging distance for alignments |
| `--min_reads` | 5 | 5 | Minimum reads to call a locus |

---

## Development Guidelines

1. **Before implementing new filtering approaches**: Read the negative result docs in `archive/failed_approaches/` to avoid repeating failed experiments.

2. **Testing new families**: Always validate against ground truth using bedtools intersect.

3. **Performance expectations:**
   - Component 1 discovery: 1-5 minutes per seed
   - Coverage calculation: 2-10 minutes depending on locus count
   - Filtering: <1 second

4. **Writing style requirements:**
   - NEVER use emojis in code, documentation, or outputs
   - NEVER use em-dashes (use hyphens or rewrite sentences)
   - Keep formatting simple and text-based

---

## Textbook Gene Families (Validated January 2026)

The method works on classic textbook gene families:

| Family | Seed Reads | Genes Found | Recall | Notes |
|--------|-----------|-------------|--------|-------|
| Globins (Alpha) | 1,071 | 5/5 | 100% | HBA1, HBA2, HBZ, HBM, HBQ1 |
| Globins (Beta) | 479 | 5/5 | 100% | HBB, HBD, HBG1, HBG2, HBE1 |
| Tubulin TUBA1 | 624 | 3/3 | 100% | TUBA1A, TUBA1B, TUBA1C cluster |
| Tubulin TUBB2 | 47 | 2/2 | 100% | TUBB2A, TUBB2B cluster |

**Key insight**: Alpha and beta globins are correctly treated as SEPARATE families (40% identity, diverged 450 MYA). This is biologically accurate.

**CORRECTION**: Previous documentation incorrectly stated globins don't work due to "97% identity being too divergent for long reads." This was caused by a bug (fixed January 2026). Globins at 97% identity DO multi-map and ARE detected.

---

## Known Limitations

**Transcribed pseudogenes**: Testis exhibits "pervasive transcription" and expresses many pseudogenes that are silenced in other tissues. These transcribed pseudogenes may pass the coverage filter (std_coverage >= 50) because they show gene-like variable coverage patterns.

**Current data is single-tissue**: Results reflect testis-specific transcription only. A locus detected only in testis could be:
- A testis-specific functional gene (true positive)
- A transcribed pseudogene (false positive that passes coverage filter)

**Recommended improvement**: Multi-tissue IsoSeq validation. See `legacy/intermediate_docs/FUTURE_DIRECTIONS_AND_CONSIDERATIONS.md` for detailed recommendations.

---

## Documentation Structure

### Validation Framework (New)
- `VALIDATION_FRAMEWORK_RECOMMENDATIONS.md`: Comprehensive validation using coverage peaks, slopes, symmetry, and four-layer validation
- `VALIDATION_QUICK_REFERENCE.md`: One-page quick reference for validation thresholds and command recipes
- `COVERAGE_SPLITTING_FORMAL_SPECIFICATION.md`: Mathematical definitions for valley detection and splitting algorithm
- `ISOSEQ_VALIDATION_RESULTS.md`: Iso-Seq validation performance (100% precision at threshold 0.6)
- `MULTIMAPPING_ANALYSIS_SUMMARY.md`: Multi-mapping read analysis results

### Current Results (Key Files)
- `TWO_SEED_DISCOVERY_RESULTS.md`: 100% precision achieved for ID_462 (dispersed)
- `COMPACT_VS_DISPERSED_FAMILIES.md`: Performance comparison across family types
- `COVERAGE_FILTERING_RESULTS.md`: Coverage filter validation
- `TEXTBOOK_FAMILIES_VALIDATION.md`: Validation on globins, tubulins

### Current Analysis Files
- `CORE_READS_FINDINGS.md`: Analysis of core read sets
- `MINIMUM_SEED_SET_FINDINGS.md`: Minimum seed set experiments
- `PER_GENE_JACCARD_ANALYSIS.md`: Per-gene Jaccard analysis

### Future Directions
- `legacy/intermediate_docs/FUTURE_DIRECTIONS_AND_CONSIDERATIONS.md`: Known limitations, unexploited IsoSeq features, multi-tissue approach recommendations

### Archived Documentation
- `archive/failed_approaches/`: Documentation of approaches that were tested and failed
- `archive/old_documentation/`: Extensive development history from earlier phases
- `legacy/`: Early experimental approaches and intermediate results

---

## Common Issues and Solutions

### Issue: Low precision in component 1 for dispersed families
**Solution**: Apply coverage filtering (std_coverage >= 50)

### Issue: Missing some family members
**Solution**: Use multi-seed approach with 2-3 seeds from different family members

### Issue: Slow BAM processing
**Solution**:
- Ensure BAM is indexed (.bai file present)
- Restrict search to specific chromosomes with `--chromosomes` parameter
- Use `--max_iterations 1` to avoid expensive transitive discovery

### Issue: Proposing sequence-based filters
**Solution**: DON'T. Read the negative result documentation first. All sequence-based approaches fail because:
- PacBio data has uniformly high quality (92.4% of reads align >=90%)
- False positives result from genuine biological relationships (high sequence identity)
- Only functional genomic context (coverage patterns) successfully filters FPs

---

## Bug Fixes (January 2026)

### Chromosome Scanning Bug (CRITICAL)

**Problem**: When `--chromosomes` was not specified, the detector scanned ZERO chromosomes, returning 0 loci regardless of input.

**Root cause**: In `find_where_reads_map()`, the alignment collection loop was inside an `if target_chroms:` block. When `target_chroms` was None, the entire loop was skipped.

**Fix**: Modified the function to get all chromosomes from the BAM header when `--chromosomes` is not specified:
```python
if target_chroms:
    chroms_to_scan = target_chroms
else:
    chroms_to_scan = [ref['SN'] for ref in bam.header['SQ']]
```

**Impact**: This bug caused the incorrect conclusion that "globins don't work because 97% identity is too divergent for long reads." After the fix, globins achieve 100% recall.

---

## Contact and References

- Primary documentation: README.md (high-level overview)
- Detailed algorithm: `legacy/intermediate_docs/BRIDGE_DETECTION_ALGORITHM.md`
- Method comparison: `legacy/intermediate_docs/unbiased_multiseed_strategies.md`

**Last Updated**: March 2026
