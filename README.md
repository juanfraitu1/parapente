# Parapente 🪂

**RNA-based paralog discovery using multi-mapping reads**

[![Rust](https://img.shields.io/badge/Rust-1.75%2B-orange)](https://www.rust-lang.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

## Overview

Parapente is a novel method for discovering paralogous gene families from RNA sequencing data. Unlike traditional approaches that rely on DNA sequence homology (which miss expressed pseudogenes) or genome assembly (which collapses highly similar sequences), Parapente leverages **multi-mapping RNA reads** as a fingerprint for paralog detection.

Parapente detects both **protein-coding genes** and **lncRNAs** (long non-coding RNAs), as both can occur in segmental duplications and produce multi-mapping reads.

### Key Innovation

Multi-mapping reads (reads that align to multiple genomic locations) are typically discarded in standard RNA-seq analysis. Parapente treats them as valuable signal: reads that map to multiple locations likely originate from recently duplicated genes or lncRNAs (segmental duplications) with high sequence similarity.

## How It Works

1. **Seed Selection**: Start with one or more known members of a gene family
2. **Read Collection**: Gather all multi-mapping reads from seed regions
3. **Transitive Discovery**: Find all genomic locations sharing these reads
4. **Coverage Splitting**: Separate merged regions into individual gene loci
5. **Validation**: Filter artifacts using coverage patterns and multi-seed support

### Algorithm Overview

```
Input: BAM file + seed regions
       ↓
Collect multi-mapping reads from seeds
       ↓
Find all alignments of seed reads across genome
       ↓
Cluster alignments into candidate loci
       ↓
Split by coverage valleys (separate merged genes)
       ↓
Filter by coverage quality, Jaccard similarity, seed support
       ↓
Output: Discovered paralog loci with metrics
```

## Installation

### From Source

```bash
git clone https://github.com/juanfraitu1/parapente.git
cd parapente
cargo build --release
```

The binary will be at `target/release/parapente`.

### Requirements

- Rust 1.75+
- Indexed BAM file (`.bam` + `.bai`)
- (Optional) Python 3.8+ for evaluation scripts

## Usage

### Basic Usage

```bash
# Single seed region
parapente \
  --bam sample.bam \
  --seeds-bed seeds.bed \
  --output results.bed \
  --threads 8
```

### Input Format

**Seeds BED file** (`seeds.bed`):
```
chr1    119682444    119734263    POM121
chr1    119839186    119871776    POM121C
```

### Output Format

**Results BED** with additional columns:
```
#chrom  start      end        jaccard  distance  component  reads
chr1    119682444  119734263  1.0000   0         1          2206
chr1    120137898  120187092  0.8913   0         1          2475
```

Columns:
- `jaccard`: Jaccard similarity to seed (0-1)
- `distance`: Graph distance from seed (0 = seed itself)
- `component`: Connected component ID
- `reads`: Number of supporting reads

## Performance

Validated on 83 challenging human gene families (T2T-CHM13 genome):

| Metric | Value |
|--------|-------|
| Families tested | 83 |
| Total genes | 362 |
| Sensitivity | **100%** |
| Avg runtime (8 threads) | ~40s/family |

### Challenging Cases

- **ID_14 (LRRC37)**: 7 genes on chromosome 1 with 99%+ identity ✅
- **ID_163 (GUSB)**: 6 genes across segmental duplications ✅
- **ID_215 (PDE4DIP)**: 7 genes in complex region ✅

## Key Features

### Multi-Seed Support
```bash
# Use multiple seeds for better discovery
parapente --seeds-bed family_seeds.bed --min-seeds 2 ...
```

### Filtering Options
```bash
# Coverage quality filter (gene-like patterns)
--min-coverage-quality 0.4

# Jaccard similarity threshold
--min-jaccard-locus 0.3

# Multi-seed support (loci must share reads with N seeds)
--min-seeds 2

# Coverage std filter (default: 50)
--min-std-coverage 50
```

### Parallel Processing
```bash
# Use 8 threads for faster processing
--threads 8
```

## Project Structure

```
parapente/
├── Cargo.toml          # Rust package config
├── src/
│   ├── main.rs         # CLI and main workflow
│   ├── transitive_detector.rs  # Core discovery algorithm
│   ├── coverage_valley.rs      # Coverage-based splitting
│   └── ...             # Additional modules
├── scripts/
│   ├── run_all_families.py     # Batch processing
│   └── generate_per_member_tsv.py  # Evaluation
├── data/
│   └── ground_truth_80families.bed  # Validation set
├── docs/
│   ├── METHODOLOGY.md
│   └── PERFORMANCE.md
└── results/
    └── all_families/   # Evaluation results
```

## Evaluation

Generate detailed per-family evaluation:

```bash
python scripts/generate_per_member_tsv.py
```

This creates `{family_id}_per_member.tsv` files with:
- Per-gene detection status
- IoU (Intersection over Union) scores
- Sensitivity and precision metrics
- Jaccard similarity and read counts

## Citation

If you use Parapente in your research, please cite:

```bibtex
@software{parapente2024,
  title={Parapente: RNA-based paralog discovery using multi-mapping reads},
  author={Iturralde Martinez, Francisco},
  year={2024},
  url={https://github.com/juanfraitu1/parapente}
}
```

## License

MIT License - see [LICENSE](LICENSE) file.

## Acknowledgments

- T2T-CHM13 human genome assembly
- PacBio Iso-Seq data from [source]
- Rust noodles crate for BAM parsing

## Contact

For questions or issues, please open a GitHub issue.
