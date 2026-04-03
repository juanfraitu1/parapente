#!/usr/bin/env python3
"""
Batch process all 83 gene families with up to 5 iterations.
Uses existing results where available, runs detector where needed.
"""

import subprocess
import pandas as pd
from pathlib import Path
import re
from collections import defaultdict
import argparse


def parse_ground_truth(bed_file):
    """Parse ground truth and group by family."""
    families = defaultdict(list)
    
    with open(bed_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 4:
                continue
            
            chrom, start, end, name = parts[0], int(parts[1]), int(parts[2]), parts[3]
            match = re.search(r'ID_(\d+)', name)
            if match:
                family_id = f"ID_{match.group(1)}"
                families[family_id].append({
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'name': name
                })
    
    return families


def parse_seeds(bed_file):
    """Parse seed regions grouped by family."""
    seeds = defaultdict(list)
    
    with open(bed_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 4:
                continue
            
            chrom, start, end, name = parts[0], int(parts[1]), int(parts[2]), parts[3]
            match = re.search(r'ID_(\d+)', name)
            if match:
                family_id = f"ID_{match.group(1)}"
                seeds[family_id].append({
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'name': name
                })
    
    return seeds


def run_detector_iteration(family_id, seeds, bam_path, output_dir, iteration, binary_path, threads=4):
    """Run one iteration of the detector."""
    output_file = output_dir / f"{family_id}_iter{iteration}.bed"
    
    # Build seed arguments
    seed_args = []
    for seed in seeds:
        seed_args.extend([
            "--seed-chrom", seed['chrom'],
            "--seed-start", str(seed['start']),
            "--seed-end", str(seed['end'])
        ])
    
    cmd = [
        binary_path,
        "--bam", bam_path,
        *seed_args,
        "--output", str(output_file),
        "--min-jaccard", "0.001",
        "--min-reads", "10",
        "--transitive",
        "--refine-boundaries",
        "--max-trim-frac", "0.85",
        "--min-tx-density", "0.1",
        "--threads", str(threads)
    ]
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=600,  # 10 minute timeout
        )
        
        # Check for refined output
        refined_file = Path(str(output_file).replace('.bed', '_refined.bed'))
        if refined_file.exists():
            return refined_file
        elif output_file.exists():
            return output_file
        else:
            return None
            
    except subprocess.TimeoutExpired:
        print(f"    Timeout for {family_id} iteration {iteration}")
        return None
    except Exception as e:
        print(f"    Error: {e}")
        return None


def parse_detection_results(bed_file):
    """Parse detection results."""
    regions = []
    
    with open(bed_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue
            
            regions.append({
                'chrom': parts[0],
                'start': int(parts[1]),
                'end': int(parts[2])
            })
    
    return regions


def calculate_overlap(truth, detected):
    """Calculate overlap."""
    if truth['chrom'] != detected['chrom']:
        return 0
    
    overlap_start = max(truth['start'], detected['start'])
    overlap_end = min(truth['end'], detected['end'])
    
    if overlap_end <= overlap_start:
        return 0
    
    return overlap_end - overlap_start


def evaluate_family(family_id, truth_regions, detected_regions):
    """Evaluate detection."""
    found = 0
    total_overlap = 0
    
    for truth in truth_regions:
        best_overlap = 0
        truth_span = truth['end'] - truth['start']
        
        for det in detected_regions:
            overlap = calculate_overlap(truth, det)
            if overlap > best_overlap:
                best_overlap = overlap
        
        if best_overlap > 0:
            found += 1
            total_overlap += (best_overlap / truth_span * 100) if truth_span > 0 else 0
    
    sensitivity = (found / len(truth_regions) * 100) if truth_regions else 0
    avg_overlap = (total_overlap / found) if found > 0 else 0
    
    return found, len(truth_regions), sensitivity, avg_overlap


def main():
    parser = argparse.ArgumentParser(description='Batch process all gene families')
    parser.add_argument('--max-iterations', type=int, default=5,
                        help='Maximum iterations per family (default: 5)')
    parser.add_argument('--force-rerun', action='store_true',
                        help='Rerun even if results exist')
    parser.add_argument('--threads', type=int, default=4,
                        help='Threads per family')
    args = parser.parse_args()
    
    print("=" * 80)
    print(f"BATCH PROCESSING: All 83 families with up to {args.max_iterations} iterations")
    print("=" * 80)
    print()
    
    # Paths
    ground_truth = "/scratch/jxi21/segments/ground_truth_80families.bed"
    seeds_file = "/scratch/jxi21/segments/seeds_from_ground_truth.bed"
    bam_path = "/scratch/jxi21/segments/A119b.bam"
    binary_path = "/scratch/jxi21/segments/gene_family_detector/target/release/gene_family_detector"
    output_dir = Path("/scratch/jxi21/segments/detector_results/minimum_seed")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load data
    print("Loading ground truth and seeds...")
    families = parse_ground_truth(ground_truth)
    seeds = parse_seeds(seeds_file)
    print(f"Found {len(families)} families")
    print()
    
    # Process each family
    results = []
    
    for i, (family_id, truth_regions) in enumerate(sorted(families.items()), 1):
        print(f"[{i}/{len(families)}] Processing {family_id} ({len(truth_regions)} genes)...")
        
        # Check existing iterations
        existing = sorted(output_dir.glob(f"{family_id}_iter*.bed"))
        n_existing = len(existing)
        
        if n_existing >= args.max_iterations and not args.force_rerun:
            print(f"    Already has {n_existing} iterations, skipping")
            best_file = existing[-1]  # Use last iteration
        elif family_id not in seeds:
            print(f"    WARNING: No seeds found for {family_id}, skipping")
            continue
        else:
            # Run detector for needed iterations
            family_seeds = seeds[family_id]
            print(f"    Using {len(family_seeds)} seed(s)")
            
            start_iter = n_existing + 1 if not args.force_rerun else 1
            
            for iteration in range(start_iter, args.max_iterations + 1):
                print(f"    Running iteration {iteration}...")
                result_file = run_detector_iteration(
                    family_id, family_seeds, bam_path, output_dir, 
                    iteration, binary_path, args.threads
                )
                
                if result_file:
                    print(f"      Saved: {result_file.name}")
                else:
                    print(f"      Failed to produce output")
                    break
            
            # Get best result
            existing = sorted(output_dir.glob(f"{family_id}_iter*.bed"))
            best_file = existing[-1] if existing else None
        
        if best_file:
            detected = parse_detection_results(best_file)
            found, total, sensitivity, overlap = evaluate_family(family_id, truth_regions, detected)
            n_chroms = len(set(r['chrom'] for r in truth_regions))
            
            results.append({
                'family_id': family_id,
                'genes_found': found,
                'genes_total': total,
                'sensitivity': sensitivity,
                'avg_overlap': overlap,
                'n_chromosomes': n_chroms,
                'n_iterations': len(existing) if existing else 0,
                'n_detected_loci': len(detected)
            })
            
            print(f"    Result: {found}/{total} genes ({sensitivity:.1f}%), "
                  f"{len(detected)} loci detected")
        else:
            print(f"    ERROR: No output file")
        
        print()
    
    # Summary
    print("=" * 80)
    print("BATCH PROCESSING COMPLETE")
    print("=" * 80)
    print()
    
    df = pd.DataFrame(results)
    
    # Overall stats
    total_genes = df['genes_total'].sum()
    found_genes = df['genes_found'].sum()
    overall_sensitivity = found_genes / total_genes * 100
    
    print("OVERALL STATISTICS:")
    print(f"  Families processed: {len(df)}")
    print(f"  Total genes: {total_genes}")
    print(f"  Genes found: {found_genes}")
    print(f"  Overall sensitivity: {overall_sensitivity:.2f}%")
    print()
    
    # Distribution
    perfect = (df['sensitivity'] == 100).sum()
    good = ((df['sensitivity'] >= 75) & (df['sensitivity'] < 100)).sum()
    fair = ((df['sensitivity'] >= 50) & (df['sensitivity'] < 75)).sum()
    poor = (df['sensitivity'] < 50).sum()
    
    print("SENSITIVITY DISTRIBUTION:")
    print(f"  100%: {perfect} families ({perfect/len(df)*100:.1f}%)")
    print(f"  75-99%: {good} families ({good/len(df)*100:.1f}%)")
    print(f"  50-74%: {fair} families ({fair/len(df)*100:.1f}%)")
    print(f"  <50%: {poor} families ({poor/len(df)*100:.1f}%)")
    print()
    
    # Problematic
    if poor > 0:
        print("FAMILIES WITH <50% SENSITIVITY:")
        for _, row in df[df['sensitivity'] < 50].iterrows():
            print(f"  {row['family_id']}: {row['genes_found']}/{row['genes_total']} "
                  f"({row['sensitivity']:.1f}%) - {row['n_chromosomes']} chrom(s)")
        print()
    
    # Save results
    results_file = "/scratch/jxi21/segments/detector_evaluation_batch_5iter.csv"
    df.to_csv(results_file, index=False)
    print(f"Results saved to: {results_file}")
    
    # Print final summary table
    print()
    print("=" * 80)
    print("COMPLETE RESULTS TABLE")
    print("=" * 80)
    print(f"{'Family':<10} {'Found':>6} {'Total':>6} {'Sens%':>8} {'Chroms':>6} {'Iters':>6} {'Loci':>6}")
    print("-" * 80)
    for _, row in df.sort_values('family_id').iterrows():
        print(f"{row['family_id']:<10} {row['genes_found']:>6} {row['genes_total']:>6} "
              f"{row['sensitivity']:>8.1f} {row['n_chromosomes']:>6} {row['n_iterations']:>6} "
              f"{row['n_detected_loci']:>6}")
    print("=" * 80)


if __name__ == '__main__':
    main()
