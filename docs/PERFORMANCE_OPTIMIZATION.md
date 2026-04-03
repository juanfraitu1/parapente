# Performance Optimization for Low-Powered Computers

This document describes the performance optimizations implemented in parapente for running on low-powered computers (laptops, older hardware, systems with limited RAM).

## Overview

The main bottlenecks in parapente are:
1. **BAM I/O**: Creating multiple BAM readers and scanning large chromosome regions
2. **Memory allocations**: String allocations for read names, hash set resizing
3. **Peak memory usage**: Storing all alignments in memory at once

## Implemented Optimizations

### 1. Low-Memory Mode (`--low-memory`)

A new CLI flag enables memory-efficient code paths:

```bash
parapente \
  --bam A119b.bam \
  --seeds-bed seeds.bed \
  --output results.bed \
  --low-memory \
  --low-mem-batch-size 50000 \
  --low-mem-chunk-size 500000
```

### 2. Batch Read Collection

The `batch_collect_reads()` function in `optimized_detector.rs` uses a **single BAM reader** to collect reads from multiple seed regions, rather than creating a new reader for each seed.

**Before** (inefficient):
```rust
for seed in seeds {
    let mut reader = bam::io::indexed_reader::Builder::default()
        .build_from_path(bam_path)?;  // New reader per seed!
    let reads = collect_reads(&mut reader, seed)?;
}
```

**After** (optimized):
```rust
let mut reader = bam::io::indexed_reader::Builder::default()
    .build_from_path(bam_path)?;  // Single reader
for seed in seeds {
    let reads = collect_reads(&mut reader, seed)?;
}
```

### 3. Configuration Options

New CLI arguments for tuning memory usage:

| Argument | Default | Description |
|----------|---------|-------------|
| `--low-memory` | false | Enable low-memory optimizations |
| `--low-mem-batch-size` | 100000 | Batch size for BAM queries |
| `--low-mem-chunk-size` | 1000000 | Max chunk size for chromosome scanning (bp) |

Smaller values = less memory usage but potentially slower performance.

### 4. Data Structure Optimizations

- **FxHashSet**: Already used throughout the codebase for faster hashing
- **String interning**: Available via `ReadInterner` for reducing memory when read names are repeated

## Usage Recommendations

### For Laptops/Low-RAM Systems (< 8GB RAM)

```bash
parapente \
  --bam A119b.bam \
  --seeds-bed seeds.bed \
  --output results.bed \
  --low-memory \
  --low-mem-batch-size 25000 \
  --low-mem-chunk-size 250000 \
  --threads 2
```

### For Standard Desktops (8-16GB RAM)

```bash
parapente \
  --bam A119b.bam \
  --seeds-bed seeds.bed \
  --output results.bed \
  --low-memory
```

### For High-End Workstations (> 32GB RAM)

Use the default mode for maximum performance:

```bash
parapente \
  --bam A119b.bam \
  --seeds-bed seeds.bed \
  --output results.bed \
  --threads 8
```

## Future Optimizations

Potential future improvements:

1. **Streaming processing**: Process loci as they're found rather than collecting all in memory
2. **Async I/O**: Use async/await for overlapping I/O and computation
3. **Memory-mapped BAM**: For systems with sufficient address space
4. **Parallel chunked scanning**: Scan chromosome chunks in parallel
5. **Read name compression**: Use integer IDs instead of strings for read names

## Benchmarking

To measure the impact of optimizations:

```bash
# Without low-memory mode
time parapente --bam A119b.bam --seeds-bed test_seeds.bed --output standard.bed

# With low-memory mode
time parapente --bam A119b.bam --seeds-bed test_seeds.bed --output lowmem.bed --low-memory

# Compare memory usage
/usr/bin/time -v parapente ... 2>&1 | grep "Maximum resident"
```

## Implementation Details

The optimization module is in `src/optimized_detector.rs` and provides:

- `LowMemoryConfig`: Configuration struct for memory tuning
- `batch_collect_reads()`: Efficient multi-region read collection
- `detect_transitive_multi_seed_low_memory()`: Low-memory detection entry point

The `--low-memory` flag triggers the use of these optimized functions in `main.rs`.
