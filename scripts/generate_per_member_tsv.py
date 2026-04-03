#!/usr/bin/env python3
"""
Generate detailed per-member TSV files for all gene families.
Uses existing prediction files in results/all_families/.
"""

import os
import glob
import re
from collections import defaultdict

GT_FILE = "ground_truth_80families.bed"
RESULTS_DIR = "results/all_families"


def parse_ground_truth():
    """Parse ground truth BED file."""
    families = defaultdict(list)
    
    with open(GT_FILE) as f:
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
                gene_name = name.split('|')[0] if '|' in name else name
                families[family_id].append({
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'name': gene_name,
                    'full_name': name
                })
    
    return families


def parse_predictions(pred_file):
    """Parse prediction BED file."""
    predictions = []
    
    with open(pred_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue
            
            pred = {
                'chrom': parts[0],
                'start': int(parts[1]),
                'end': int(parts[2]),
                'jaccard': float(parts[3]) if len(parts) > 3 else 0.0,
                'distance': int(parts[4]) if len(parts) > 4 else 0,
                'component': int(parts[5]) if len(parts) > 5 else 0,
                'reads': int(parts[6]) if len(parts) > 6 else 0,
            }
            predictions.append(pred)
    
    return predictions


def compute_overlap(start1, end1, start2, end2):
    """Compute overlap between two intervals."""
    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)
    return max(0, overlap_end - overlap_start)


def compute_iou(truth, detected):
    """Compute IoU between truth and detected regions."""
    if truth['chrom'] != detected['chrom']:
        return 0.0
    
    overlap = compute_overlap(truth['start'], truth['end'], detected['start'], detected['end'])
    union = (truth['end'] - truth['start']) + (detected['end'] - detected['start']) - overlap
    
    return overlap / union if union > 0 else 0.0


def generate_per_member_tsv(family_id, genes, predictions):
    """Generate per-member TSV file for a family."""
    
    # Calculate family-level metrics
    copies_available = len(genes)
    
    # Find which genes were found
    found_genes = set()
    gene_best_iou = {}
    
    for gene in genes:
        gene_key = (gene['chrom'], gene['start'], gene['end'])
        best_iou = 0.0
        
        for pred in predictions:
            if pred['chrom'] == gene['chrom']:
                iou = compute_iou(gene, pred)
                if iou > best_iou:
                    best_iou = iou
                if iou > 0:  # Any overlap counts as found
                    found_genes.add(gene_key)
        
        gene_best_iou[gene_key] = best_iou
    
    copies_found = len(found_genes)
    copies_extra = len(predictions) - copies_found
    
    sensitivity = copies_found / copies_available if copies_available > 0 else 0.0
    precision = copies_found / len(predictions) if predictions else 0.0
    
    # Build rows
    rows = []
    
    # Add seed rows (genes)
    for i, gene in enumerate(genes):
        gene_key = (gene['chrom'], gene['start'], gene['end'])
        iou = gene_best_iou.get(gene_key, 0.0)
        
        # Find matching prediction for this gene (if any)
        matching_pred = None
        best_overlap = 0
        for pred in predictions:
            if pred['chrom'] == gene['chrom']:
                overlap = compute_overlap(gene['start'], gene['end'], pred['start'], pred['end'])
                if overlap > best_overlap:
                    best_overlap = overlap
                    matching_pred = pred
        
        row = {
            'location': f"{gene['chrom']}:{gene['start']}-{gene['end']}",
            'member_name': gene['name'],
            'family_name': family_id,
            'family_id': family_id,
            'copies_available': copies_available,
            'copies_found': copies_found,
            'copies_extra': copies_extra,
            'sensitivity': sensitivity,
            'precision': precision,
            'iou': iou,
            'pred_jaccard': matching_pred['jaccard'] if matching_pred else '',
            'pred_reads': matching_pred['reads'] if matching_pred else '',
            'pred_span_bp': (matching_pred['end'] - matching_pred['start']) if matching_pred else '',
            'n_peaks': '',  # Not available in current output
            'frac_high_bins': '',  # Would need coverage analysis
            'is_seed': 1
        }
        rows.append(row)
    
    # Add extra predictions (non-seeds)
    extra_count = 0
    for pred in predictions:
        # Check if this prediction overlaps any gene
        overlaps_gene = False
        for gene in genes:
            if pred['chrom'] == gene['chrom']:
                overlap = compute_overlap(pred['start'], pred['end'], gene['start'], gene['end'])
                if overlap > 0:
                    overlaps_gene = True
                    break
        
        if not overlaps_gene:
            extra_count += 1
            row = {
                'location': f"{pred['chrom']}:{pred['start']}-{pred['end']}",
                'member_name': f"extra_{extra_count}",
                'family_name': family_id,
                'family_id': family_id,
                'copies_available': copies_available,
                'copies_found': copies_found,
                'copies_extra': copies_extra,
                'sensitivity': sensitivity,
                'precision': precision,
                'iou': '',  # No truth to compare against
                'pred_jaccard': pred['jaccard'],
                'pred_reads': pred['reads'],
                'pred_span_bp': pred['end'] - pred['start'],
                'n_peaks': '',
                'frac_high_bins': '',
                'is_seed': 0
            }
            rows.append(row)
    
    return rows


def main():
    print("=" * 70)
    print("Generating per-member TSV files for all families")
    print("=" * 70)
    
    # Load ground truth
    families = parse_ground_truth()
    print(f"Loaded {len(families)} families from {GT_FILE}")
    print()
    
    # Process each family
    total_updated = 0
    
    for family_id in sorted(families.keys()):
        genes = families[family_id]
        
        # Load predictions
        pred_file = os.path.join(RESULTS_DIR, family_id, f"{family_id}_predictions.bed")
        if not os.path.exists(pred_file):
            print(f"[{family_id}] No predictions file found, skipping")
            continue
        
        predictions = parse_predictions(pred_file)
        
        # Generate TSV rows
        rows = generate_per_member_tsv(family_id, genes, predictions)
        
        # Write TSV file
        tsv_file = os.path.join(RESULTS_DIR, family_id, f"{family_id}_per_member.tsv")
        
        with open(tsv_file, 'w') as f:
            # Write header
            headers = ['location', 'member_name', 'family_name', 'family_id', 
                      'copies_available', 'copies_found', 'copies_extra',
                      'sensitivity', 'precision', 'iou', 'pred_jaccard', 
                      'pred_reads', 'pred_span_bp', 'n_peaks', 'frac_high_bins', 'is_seed']
            f.write('\t'.join(headers) + '\n')
            
            # Write rows
            for row in rows:
                values = [str(row.get(h, '')) for h in headers]
                f.write('\t'.join(values) + '\n')
        
        print(f"[{family_id}] {len(genes)} genes, {len(predictions)} predictions -> {tsv_file}")
        total_updated += 1
    
    print()
    print("=" * 70)
    print(f"Updated {total_updated} families")
    print("=" * 70)


if __name__ == "__main__":
    main()
