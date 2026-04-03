# Performance Evaluation

## Benchmark Dataset

**83 Human Gene Families** from T2T-CHM13 genome assembly

| Statistic | Value |
|-----------|-------|
| Total families | 83 |
| Total genes | 362 |
| Avg genes/family | 4.4 |
| Max genes/family | 17 (ID_207) |
| Chromosomes | 23 (chr1-chr22, chrX) |

### Dataset Characteristics

Families were selected to represent challenging cases:
- **Segmental duplications**: Genes in recent, high-identity duplications
- **Multi-chromosomal**: Families spanning multiple chromosomes
- **Tandem arrays**: Adjacent paralogs (e.g., TRIM genes)
- **Processed pseudogenes**: Retrocopied gene copies

## Sensitivity Results

### Overall Performance

| Metric | Value |
|--------|-------|
| Genes found | 362 / 362 |
| **Sensitivity** | **100%** |
| False negatives | 0 |

### By Family Size

| Family Size | Count | Sensitivity |
|-------------|-------|-------------|
| 1 gene | 8 | 100% |
| 2-3 genes | 25 | 100% |
| 4-6 genes | 31 | 100% |
| 7-10 genes | 15 | 100% |
| 11+ genes | 4 | 100% |

### Challenging Cases

| Family | Genes | Challenge | Status |
|--------|-------|-----------|--------|
| ID_14 (LRRC37) | 7 | 99%+ identity | ✅ Found all 7 |
| ID_163 (GUSB) | 6 | Segmental dup | ✅ Found all 6 |
| ID_215 (PDE4DIP) | 7 | Complex region | ✅ Found all 7 |
| ID_207 (SPDYE) | 17 | Large family | ✅ Found all 17 |

## Precision Analysis

Precision varies by family complexity:

| Family | True Positives | False Positives | Precision |
|--------|---------------|-----------------|-----------|
| ID_12 | 4 | 49 | 7.5% |
| ID_14 | 7 | 17 | 29.2% |
| ID_163 | 6 | 132 | 4.3% |
| ID_215 | 7 | 73 | 8.8% |

**Note**: Low precision is expected for discovery methods. False positives can be filtered using:
- Coverage quality thresholds
- Jaccard similarity cutoffs
- Multi-seed support requirements

## Runtime Performance

### Single Family (8 threads)

| Metric | Value |
|--------|-------|
| Median runtime | 35s |
| Mean runtime | 48s |
| Max runtime | 121s (ID_163) |

### Batch Processing (83 families)

| Configuration | Total Time |
|---------------|------------|
| Sequential | ~2 hours |
| Parallel (8 threads/family) | ~35 minutes |

### Scaling

| Threads | Speedup | Efficiency |
|---------|---------|------------|
| 1 | 1.0x | 100% |
| 4 | 3.2x | 80% |
| 8 | 5.5x | 69% |

*Efficiency limited by BAM I/O and coverage splitting synchronization.*

## Memory Usage

| Metric | Value |
|--------|-------|
| Peak memory (single family) | ~2 GB |
| Memory scaling | Sub-linear with family size |
| BAM index | Memory-mapped (~25 MB) |

## Comparison to Other Methods

| Method | Input | Sensitivity | Precision | Novel Discovery |
|--------|-------|-------------|-----------|-----------------|
| **Parapente** | RNA BAM | **100%** | ~10% | ✅ Yes |
| BLAST + Synteny | DNA | 85% | 90% | ❌ No (needs query) |
| Genome Assembly | DNA reads | 70% | 95% | ❌ No (collapses paralogs) |
| Iso-Seq + StringTie | Long reads | 90% | 85% | ⚠️ Partial |

*Note: Comparison methods estimated from published benchmarks on similar datasets.*

## Filter Impact

Applying filters improves precision while maintaining sensitivity:

| Filter | Loci Kept | Sensitivity | Precision |
|--------|-----------|-------------|-----------|
| None | 2687 | 100% | 13.5% |
| `--min-std-coverage 50` | 2410 | 100% | 15.0% |
| `--min-jaccard 0.3` | 1985 | 100% | 18.2% |
| `--min-seeds 2` | 1650 | 100% | 21.9% |
| Combined | 892 | 100% | 40.6% |

## Known Limitations

1. **Requires seed**: Needs at least one known family member
2. **RNA-only**: Cannot find non-transcribed pseudogenes
3. **Coverage-dependent**: Low-expression genes may be missed
4. **Large families**: Runtime increases with family size

## Future Improvements

- [ ] De novo seed finding (no prior knowledge needed)
- [ ] Support for single-cell RNA-seq
- [ ] Integration with gene annotation (GFF filtering)
- [ ] GPU acceleration for coverage analysis
