#!/usr/bin/env python3
"""
Resume running gene family detector on remaining families.
"""

import subprocess
import os
import sys
from collections import defaultdict
import time

GT_FILE = "ground_truth_80families.bed"
BAM_FILE = "A119b.bam"
DETECTOR = "./gene_family_detector/target/release/gene_family_detector"
RESULTS_DIR = "results/all_families"
TIMEOUT = 600  # 10 minutes per family

def load_ground_truth():
    """Load ground truth and group by family."""
    families = defaultdict(list)
    with open(GT_FILE) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                chrom, start, end = parts[0], int(parts[1]), int(parts[2])
                name_id = parts[3]
                if '|' in name_id:
                    name, fam_id = name_id.split('|')
                else:
                    name = name_id
                    fam_id = name_id
                families[fam_id].append({
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'name': name
                })
    return families

def create_seeds_file(fam_id, genes):
    """Create seeds BED file."""
    fam_dir = os.path.join(RESULTS_DIR, fam_id)
    os.makedirs(fam_dir, exist_ok=True)
    seeds_file = os.path.join(fam_dir, f"{fam_id}_seeds.bed")
    with open(seeds_file, 'w') as f:
        for gene in genes:
            f.write(f"{gene['chrom']}\t{gene['start']}\t{gene['end']}\t{gene['name']}|{fam_id}\n")
    return seeds_file

def already_processed(fam_id):
    """Check if family already has predictions."""
    pred_file = os.path.join(RESULTS_DIR, fam_id, f"{fam_id}_predictions.bed")
    if os.path.exists(pred_file):
        with open(pred_file) as f:
            lines = sum(1 for line in f if not line.startswith('#') and line.strip())
        return lines > 0
    return False

def run_detection(fam_id, seeds_file):
    """Run gene family detection with new defaults."""
    fam_dir = os.path.join(RESULTS_DIR, fam_id)
    output_file = os.path.join(fam_dir, f"{fam_id}_predictions.bed")
    
    cmd = [
        DETECTOR,
        "--bam", BAM_FILE,
        "--seeds-bed", seeds_file,
        "--output", output_file,
        "--threads", "8"
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=TIMEOUT)
        if result.returncode != 0:
            print(f"  ERROR: {result.stderr[:200]}")
            return None, "error"
        
        # Count detected loci
        loci_count = 0
        if os.path.exists(output_file):
            with open(output_file) as f:
                for line in f:
                    if not line.startswith('#') and line.strip():
                        loci_count += 1
        
        return loci_count, "success"
    except subprocess.TimeoutExpired:
        print(f"  TIMEOUT after {TIMEOUT}s")
        return None, "timeout"
    except Exception as e:
        print(f"  ERROR: {e}")
        return None, "error"

def main():
    print("=" * 70)
    print("Resuming gene family detector run")
    print("=" * 70)
    
    # Load ground truth
    families = load_ground_truth()
    print(f"Loaded {len(families)} families from {GT_FILE}")
    
    # Statistics
    total = len(families)
    skipped = 0
    processed = 0
    success = 0
    failed = 0
    total_loci = 0
    
    start_time = time.time()
    
    for i, fam_id in enumerate(sorted(families.keys()), 1):
        genes = families[fam_id]
        
        # Skip if already processed
        if already_processed(fam_id):
            # Count existing loci
            pred_file = os.path.join(RESULTS_DIR, fam_id, f"{fam_id}_predictions.bed")
            with open(pred_file) as f:
                loci_count = sum(1 for line in f if not line.startswith('#') and line.strip())
            print(f"[{i}/{total}] {fam_id}: {len(genes)} genes... SKIPPED (already done, {loci_count} loci)")
            skipped += 1
            total_loci += loci_count
            continue
        
        processed += 1
        print(f"[{i}/{total}] {fam_id}: {len(genes)} genes... ", end='', flush=True)
        
        # Create seeds file
        seeds_file = create_seeds_file(fam_id, genes)
        
        # Run detection
        loci_count, status = run_detection(fam_id, seeds_file)
        
        if status == "success":
            success += 1
            total_loci += loci_count
            print(f"OK ({loci_count} loci)")
        else:
            failed += 1
            print(f"{status.upper()}")
    
    elapsed = time.time() - start_time
    
    print()
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Total families: {total}")
    print(f"Skipped (already done): {skipped}")
    print(f"Processed this run: {processed}")
    print(f"  Successful: {success}")
    print(f"  Failed: {failed}")
    print(f"Total loci detected: {total_loci}")
    print(f"Time: {elapsed:.1f}s ({elapsed/60:.1f} min)")
    print()
    print(f"Results updated in: {RESULTS_DIR}/")

if __name__ == "__main__":
    main()
