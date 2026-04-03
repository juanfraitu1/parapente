#!/usr/bin/env python3
"""
Batch evaluation of all 83 gene families using optimized multi-seed detection.
Produces a summary table with sensitivity, specificity, and IoU metrics.
"""

import subprocess
import os
import sys
from collections import defaultdict

GT_FILE = "ground_truth_80families.bed"
BAM_FILE = "A119b.bam"
DETECTOR = "./gene_family_detector/target/release/gene_family_detector"
OUTPUT_DIR = "batch_results"
TIMEOUT = 300  # 5 minutes per family

def load_ground_truth():
    """Load ground truth and group by family."""
    families = defaultdict(list)
    with open(GT_FILE) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                chrom, start, end = parts[0], int(parts[1]), int(parts[2])
                name_id = parts[3]
                # Extract family ID (e.g., "ID_14" from "LRRC37A|ID_14")
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
    """Create seeds BED file (using all genes as seeds)."""
    seeds_file = os.path.join(OUTPUT_DIR, f"{fam_id}_seeds.bed")
    with open(seeds_file, 'w') as f:
        for gene in genes:
            f.write(f"{gene['chrom']}\t{gene['start']}\t{gene['end']}\t{gene['name']}|{fam_id}\n")
    return seeds_file

def run_detection(fam_id, seeds_file):
    """Run gene family detection."""
    output_file = os.path.join(OUTPUT_DIR, f"{fam_id}_detected.bed")
    cmd = [
        DETECTOR,
        "--bam", BAM_FILE,
        "--seed-chrom", "NC_060925.1",  # Placeholder, seeds file takes precedence
        "--seed-start", "1",
        "--seed-end", "100",
        "--seeds-bed", seeds_file,
        "--output", output_file,
        "--transitive",
        "--threads", "4"
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=TIMEOUT)
        if result.returncode != 0:
            print(f"  ERROR: {result.stderr[:200]}")
            return None, "error"
        return output_file, "success"
    except subprocess.TimeoutExpired:
        print(f"  TIMEOUT after {TIMEOUT}s")
        return None, "timeout"
    except Exception as e:
        print(f"  ERROR: {e}")
        return None, "error"

def parse_detected(output_file):
    """Parse detected loci from output BED."""
    loci = []
    if not output_file or not os.path.exists(output_file):
        return loci
    with open(output_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                loci.append({
                    'chrom': parts[0],
                    'start': int(parts[1]),
                    'end': int(parts[2])
                })
    return loci

def compute_overlap(locus1, locus2):
    """Compute overlap between two loci."""
    if locus1['chrom'] != locus2['chrom']:
        return 0
    overlap_start = max(locus1['start'], locus2['start'])
    overlap_end = min(locus1['end'], locus2['end'])
    if overlap_start < overlap_end:
        return overlap_end - overlap_start
    return 0

def compute_iou(locus1, locus2):
    """Compute Intersection over Union."""
    overlap = compute_overlap(locus1, locus2)
    if overlap == 0:
        return 0
    union = (locus1['end'] - locus1['start']) + (locus2['end'] - locus2['start']) - overlap
    return overlap / union if union > 0 else 0

def bipartite_match_iou(gt_genes, detected_loci, threshold=0.1):
    """
    Compute bipartite matching IoU.
    Returns: matched_gt_count, matched_detected_count, total_iou
    """
    if not gt_genes or not detected_loci:
        return 0, 0, 0.0
    
    # Compute all IoU pairs
    pairs = []
    for i, gt in enumerate(gt_genes):
        for j, det in enumerate(detected_loci):
            iou = compute_iou(gt, det)
            if iou > 0:
                pairs.append((iou, i, j))
    
    # Sort by IoU descending
    pairs.sort(reverse=True)
    
    # Greedy matching
    matched_gt = set()
    matched_det = set()
    total_iou = 0.0
    
    for iou, i, j in pairs:
        if i not in matched_gt and j not in matched_det:
            matched_gt.add(i)
            matched_det.add(j)
            total_iou += iou
    
    return len(matched_gt), len(matched_det), total_iou

def evaluate_family(fam_id, gt_genes, detected_loci):
    """Evaluate detection results for a family."""
    # True positives: ground truth genes that overlap with detected loci
    tp = 0
    for gt in gt_genes:
        for det in detected_loci:
            if gt['chrom'] == det['chrom']:
                overlap_start = max(gt['start'], det['start'])
                overlap_end = min(gt['end'], det['end'])
                # Consider a gene "found" if at least 50% is covered
                gene_len = gt['end'] - gt['start']
                overlap_len = max(0, overlap_end - overlap_start)
                if overlap_len > gene_len * 0.5:
                    tp += 1
                    break
    
    # False negatives: ground truth genes not found
    fn = len(gt_genes) - tp
    
    # False positives: detected loci that don't overlap any ground truth
    fp = 0
    for det in detected_loci:
        found = False
        for gt in gt_genes:
            if gt['chrom'] == det['chrom']:
                overlap_start = max(gt['start'], det['start'])
                overlap_end = min(gt['end'], det['end'])
                overlap_len = max(0, overlap_end - overlap_start)
                if overlap_len > 0:
                    found = True
                    break
        if not found:
            fp += 1
    
    # Bipartite matching IoU
    matched_gt, matched_det, total_iou = bipartite_match_iou(gt_genes, detected_loci)
    avg_iou = total_iou / max(matched_gt, 1) if matched_gt > 0 else 0.0
    
    # Metrics
    sensitivity = tp / len(gt_genes) if gt_genes else 0.0
    specificity = 1 - (fp / len(detected_loci)) if detected_loci else 1.0
    
    return {
        'tp': tp,
        'fn': fn,
        'fp': fp,
        'sensitivity': sensitivity,
        'specificity': specificity,
        'iou': avg_iou,
        'matched_gt': matched_gt,
        'matched_det': matched_det
    }

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Load ground truth
    print("Loading ground truth...")
    families = load_ground_truth()
    print(f"Found {len(families)} families")
    
    # Process each family
    results = []
    fam_ids = sorted(families.keys())
    
    for i, fam_id in enumerate(fam_ids):
        genes = families[fam_id]
        print(f"\n[{i+1}/{len(fam_ids)}] {fam_id} ({len(genes)} genes)")
        
        # Create seeds file
        seeds_file = create_seeds_file(fam_id, genes)
        
        # Run detection
        output_file, status = run_detection(fam_id, seeds_file)
        
        if status != "success":
            results.append({
                'fam_id': fam_id,
                'n_seeds': len(genes),
                'n_gt': len(genes),
                'n_detected': 0,
                'status': status,
                'sensitivity': 0.0,
                'specificity': 1.0,
                'iou': 0.0
            })
            continue
        
        # Parse and evaluate
        detected_loci = parse_detected(output_file)
        metrics = evaluate_family(fam_id, genes, detected_loci)
        
        results.append({
            'fam_id': fam_id,
            'n_seeds': len(genes),
            'n_gt': len(genes),
            'n_detected': len(detected_loci),
            'status': status,
            'tp': metrics['tp'],
            'fn': metrics['fn'],
            'fp': metrics['fp'],
            'sensitivity': metrics['sensitivity'],
            'specificity': metrics['specificity'],
            'iou': metrics['iou']
        })
        
        print(f"  Detected: {len(detected_loci)} loci, Sens: {metrics['sensitivity']:.1%}, IoU: {metrics['iou']:.3f}")
    
    # Print summary table
    print("\n" + "="*100)
    print("SUMMARY TABLE")
    print("="*100)
    print(f"{'Family':<12} {'Seeds':>6} {'GT':>6} {'Found':>6} {'Extra':>6} {'Sens':>8} {'Spec':>8} {'IoU':>8} {'Status':<10}")
    print("-"*100)
    
    total_gt = 0
    total_found = 0
    total_extra = 0
    total_sens = 0
    total_spec = 0
    total_iou = 0
    n_success = 0
    
    for r in results:
        extra = r['n_detected'] - r['tp'] if 'tp' in r else r['n_detected']
        print(f"{r['fam_id']:<12} {r['n_seeds']:>6} {r['n_gt']:>6} {r['tp'] if 'tp' in r else 0:>6} {extra:>6} {r['sensitivity']:>7.1%} {r['specificity']:>7.1%} {r['iou']:>7.3f} {r['status']:<10}")
        
        if r['status'] == 'success':
            total_gt += r['n_gt']
            total_found += r['tp'] if 'tp' in r else 0
            total_extra += extra
            total_sens += r['sensitivity']
            total_spec += r['specificity']
            total_iou += r['iou']
            n_success += 1
    
    print("-"*100)
    if n_success > 0:
        avg_sens = total_sens / n_success
        avg_spec = total_spec / n_success
        avg_iou = total_iou / n_success
        print(f"{'TOTAL':<12} {'':<6} {total_gt:>6} {total_found:>6} {total_extra:>6}")
        print(f"{'AVERAGE':<12} {'':<6} {'':<6} {'':<6} {'':<6} {avg_sens:>7.1%} {avg_spec:>7.1%} {avg_iou:>7.3f}")
    
    # Save to file
    table_file = os.path.join(OUTPUT_DIR, "summary_table.tsv")
    with open(table_file, 'w') as f:
        f.write("Family\tSeeds\tGT\tFound\tExtra\tSensitivity\tSpecificity\tIoU\tStatus\n")
        for r in results:
            extra = r['n_detected'] - r['tp'] if 'tp' in r else r['n_detected']
            f.write(f"{r['fam_id']}\t{r['n_seeds']}\t{r['n_gt']}\t{r['tp'] if 'tp' in r else 0}\t{extra}\t{r['sensitivity']:.4f}\t{r['specificity']:.4f}\t{r['iou']:.4f}\t{r['status']}\n")
    
    print(f"\nResults saved to {table_file}")

if __name__ == "__main__":
    main()
