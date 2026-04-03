//! Fast gene family detector with valley splitting and artifact filtering
//! 
//! Usage:
//!   gene_family_detector_fast --bam input.bam --seed-chrom NC_060925.1 --seed-start 1000000 --seed-end 2000000 --output results.bed

use anyhow::{Context, Result};
use clap::Parser;
use hashbrown::{HashMap, HashSet};
use noodles::bam;
use noodles::sam::alignment::Record;
use std::fs::File;
use std::io::Write;
use std::path::Path;

const INDEX_VERSION: u32 = 3;  // Keep compatible with existing index
const SPATIAL_BIN_SIZE: i64 = 10000; // 10kb bins for spatial indexing

#[derive(Debug, Clone)]
struct Alignment {
    chrom_id: u32,
    start: i32,
    end: i32,
}

/// Spatial index for efficient region queries
/// Maps (chrom_id, bin_index) -> set of read names
struct SpatialIndex {
    bin_size: i64,
    bins: HashMap<(u32, i64), HashSet<String>>,
}

impl SpatialIndex {
    fn build(index: &HashMap<String, Vec<Alignment>>, bin_size: i64) -> Self {
        let mut bins: HashMap<(u32, i64), HashSet<String>> = HashMap::new();
        
        for (read_name, alignments) in index {
            for aln in alignments {
                let chrom_id = aln.chrom_id;
                let start_bin = aln.start as i64 / bin_size;
                let end_bin = aln.end as i64 / bin_size;
                
                // Add read to all bins it overlaps
                for bin_idx in start_bin..=end_bin {
                    bins.entry((chrom_id, bin_idx))
                        .or_default()
                        .insert(read_name.clone());
                }
            }
        }
        
        SpatialIndex { bin_size, bins }
    }
    
    /// Query reads overlapping a region
    fn query_region(&self, chrom_id: u32, start: i64, end: i64) -> HashSet<String> {
        let mut result = HashSet::new();
        let start_bin = start / self.bin_size;
        let end_bin = end / self.bin_size;
        
        for bin_idx in start_bin..=end_bin {
            if let Some(reads) = self.bins.get(&(chrom_id, bin_idx)) {
                result.extend(reads.iter().cloned());
            }
        }
        
        result
    }
    
    /// Query reads overlapping multiple regions (union)
    fn query_regions(&self, regions: &[(u32, i64, i64)]) -> HashSet<String> {
        let mut result = HashSet::new();
        for (chrom_id, start, end) in regions {
            result.extend(self.query_region(*chrom_id, *start, *end));
        }
        result
    }
}

struct BamIndex {
    index: HashMap<String, Vec<Alignment>>,
    chrom_names: Vec<String>,
    spatial: SpatialIndex,
}

impl BamIndex {
    fn build(bam_path: &str, include_secondary: bool, include_supplementary: bool) -> Result<Self> {
        println!("Building BAM index...");
        let start = std::time::Instant::now();

        let mut reader = bam::io::reader::Builder::default()
            .build_from_path(bam_path)
            .with_context(|| format!("Failed to open BAM: {}", bam_path))?;

        let header = reader.read_header()?;
        let ref_seqs = header.reference_sequences();

        let mut chrom_names: Vec<String> = Vec::with_capacity(ref_seqs.len());
        for i in 0..ref_seqs.len() {
            if let Some((name, _)) = ref_seqs.get_index(i) {
                chrom_names.push(std::str::from_utf8(name).unwrap_or("unknown").to_string());
            }
        }

        println!("  Indexing {} chromosomes...", chrom_names.len());

        let mut index: HashMap<String, Vec<Alignment>> = HashMap::with_capacity(20_000_000);
        let mut counter: u64 = 0;
        let mut last_report: u64 = 0;

        for result in reader.records() {
            let record = result?;
            let flags = record.flags();

            if flags.is_unmapped() { continue; }
            if !include_supplementary && flags.is_supplementary() { continue; }
            if !include_secondary && flags.is_secondary() { continue; }

            let name = match record.name() {
                Some(n) => match std::str::from_utf8(n) {
                    Ok(s) => s.to_string(),
                    Err(_) => continue,
                },
                None => continue,
            };

            let chrom_id = record.reference_sequence_id()
                .transpose()?
                .map(|id| usize::from(id) as u32)
                .unwrap_or(u32::MAX);

            if chrom_id == u32::MAX || chrom_id as usize >= chrom_names.len() { continue; }

            let aln_start = record.alignment_start()
                .transpose()?
                .map(|p| usize::from(p) as i32)
                .unwrap_or(0);
            let aln_end = record.alignment_end()
                .transpose()?
                .map(|p| usize::from(p) as i32)
                .unwrap_or(aln_start + 1);

            index.entry(name).or_default().push(Alignment {
                chrom_id, start: aln_start, end: aln_end,
            });

            counter += 1;
            if counter - last_report >= 50_000_000 {
                let elapsed = start.elapsed().as_secs_f64();
                println!("  Indexed {}M alignments ({:.1}M/sec)", 
                    counter / 1_000_000, counter as f64 / elapsed / 1_000_000.0);
                last_report = counter;
            }
        }

        println!("  Indexed {} alignments from {} reads in {:.1}s",
            counter, index.len(), start.elapsed().as_secs_f64());
        
        // Build spatial index (optimized inline during main build)
        println!("Building spatial index...");
        let spatial_start = std::time::Instant::now();
        // Build directly from the completed index
        let mut bins: HashMap<(u32, i64), HashSet<String>> = HashMap::with_capacity(100_000);
        for (read_name, alignments) in &index {
            for aln in alignments {
                let chrom_id = aln.chrom_id;
                let start_bin = aln.start as i64 / SPATIAL_BIN_SIZE;
                let end_bin = aln.end as i64 / SPATIAL_BIN_SIZE;
                for bin_idx in start_bin..=end_bin {
                    bins.entry((chrom_id, bin_idx))
                        .or_default()
                        .insert(read_name.clone());
                }
            }
        }
        let spatial = SpatialIndex { bin_size: SPATIAL_BIN_SIZE, bins };
        println!("  Built {} bins in {:.1}s", 
            spatial.bins.len(), spatial_start.elapsed().as_secs_f64());

        Ok(BamIndex { index, chrom_names, spatial })
    }

    fn save(&self, path: &str) -> Result<()> {
        println!("Saving index to {}...", path);
        let start = std::time::Instant::now();

        let mut file = File::create(path)?;
        file.write_all(&INDEX_VERSION.to_le_bytes())?;

        let chrom_count = self.chrom_names.len() as u32;
        file.write_all(&chrom_count.to_le_bytes())?;
        for name in &self.chrom_names {
            let bytes = name.as_bytes();
            file.write_all(&(bytes.len() as u32).to_le_bytes())?;
            file.write_all(bytes)?;
        }

        file.write_all(&(self.index.len() as u64).to_le_bytes())?;

        for (read_name, alignments) in &self.index {
            let bytes = read_name.as_bytes();
            file.write_all(&(bytes.len() as u32).to_le_bytes())?;
            file.write_all(bytes)?;
            file.write_all(&(alignments.len() as u32).to_le_bytes())?;

            for aln in alignments {
                file.write_all(&aln.chrom_id.to_le_bytes())?;
                file.write_all(&aln.start.to_le_bytes())?;
                file.write_all(&aln.end.to_le_bytes())?;
            }
        }

        println!("  Saved in {:.1}s ({} MB)", 
            start.elapsed().as_secs_f64(), std::fs::metadata(path)?.len() / 1_000_000);
        Ok(())
    }

    fn load(path: &str) -> Result<Self> {
        println!("Loading index from {}...", path);
        let start = std::time::Instant::now();

        let bytes = std::fs::read(path)?;
        let mut pos = 0;

        let version = u32::from_le_bytes(bytes[pos..pos+4].try_into()?);
        pos += 4;
        if version != INDEX_VERSION {
            anyhow::bail!("Index version mismatch: {} vs {}", version, INDEX_VERSION);
        }

        let chrom_count = u32::from_le_bytes(bytes[pos..pos+4].try_into()?) as usize;
        pos += 4;
        let mut chrom_names = Vec::with_capacity(chrom_count);
        for _ in 0..chrom_count {
            let len = u32::from_le_bytes(bytes[pos..pos+4].try_into()?) as usize;
            pos += 4;
            chrom_names.push(String::from_utf8(bytes[pos..pos+len].to_vec())?);
            pos += len;
        }

        let entry_count = u64::from_le_bytes(bytes[pos..pos+8].try_into()?) as usize;
        pos += 8;
        let mut index: HashMap<String, Vec<Alignment>> = HashMap::with_capacity(entry_count);

        for _ in 0..entry_count {
            let name_len = u32::from_le_bytes(bytes[pos..pos+4].try_into()?) as usize;
            pos += 4;
            let read_name = String::from_utf8(bytes[pos..pos+name_len].to_vec())?;
            pos += name_len;

            let aln_count = u32::from_le_bytes(bytes[pos..pos+4].try_into()?) as usize;
            pos += 4;

            let mut alignments = Vec::with_capacity(aln_count);
            for _ in 0..aln_count {
                let chrom_id = u32::from_le_bytes(bytes[pos..pos+4].try_into()?);
                pos += 4;
                let start = i32::from_le_bytes(bytes[pos..pos+4].try_into()?);
                pos += 4;
                let end = i32::from_le_bytes(bytes[pos..pos+4].try_into()?);
                pos += 4;
                alignments.push(Alignment { chrom_id, start, end });
            }
            index.insert(read_name, alignments);
        }

        println!("  Loaded {} reads in {:.1}s", index.len(), start.elapsed().as_secs_f64());
        
        // Build spatial index
        println!("Building spatial index...");
        let spatial_start = std::time::Instant::now();
        let spatial = SpatialIndex::build(&index, SPATIAL_BIN_SIZE);
        println!("  Built spatial index with {} bins in {:.1}s", 
            spatial.bins.len(), spatial_start.elapsed().as_secs_f64());
        
        Ok(BamIndex { index, chrom_names, spatial })
    }

    fn index_path(bam_path: &str) -> String {
        format!("{}.gfi", bam_path)
    }

    fn load_or_build(bam_path: &str, include_secondary: bool, include_supplementary: bool) -> Result<Self> {
        let index_path = Self::index_path(bam_path);

        if Path::new(&index_path).exists() {
            let bam_mtime = std::fs::metadata(bam_path)?.modified()?;
            let index_mtime = std::fs::metadata(&index_path)?.modified()?;
            if index_mtime >= bam_mtime {
                return Self::load(&index_path);
            }
            println!("Index outdated, rebuilding...");
        }

        let index = Self::build(bam_path, include_secondary, include_supplementary)?;
        index.save(&index_path)?;
        Ok(index)
    }

    fn get_alignments_by_chrom(&self, read_names: &HashSet<String>) -> HashMap<String, Vec<(i64, i64, String)>> {
        let mut by_chrom: HashMap<String, Vec<(i64, i64, String)>> = HashMap::new();
        for name in read_names {
            if let Some(alignments) = self.index.get(name) {
                for aln in alignments {
                    if aln.chrom_id as usize >= self.chrom_names.len() { continue; }
                    let chrom = &self.chrom_names[aln.chrom_id as usize];
                    by_chrom.entry(chrom.clone())
                        .or_default()
                        .push((aln.start as i64, aln.end as i64, name.clone()));
                }
            }
        }
        by_chrom
    }
    
    /// Query reads overlapping a genomic region using spatial index
    fn query_region(&self, chrom: &str, start: i64, end: i64) -> HashSet<String> {
        if let Some(chrom_id) = self.chrom_names.iter().position(|n| n == chrom) {
            self.spatial.query_region(chrom_id as u32, start, end)
        } else {
            HashSet::new()
        }
    }
    
    /// Query reads overlapping multiple regions (union)
    fn query_regions(&self, regions: &[(String, i64, i64)]) -> HashSet<String> {
        let chrom_regions: Vec<(u32, i64, i64)> = regions.iter()
            .filter_map(|(c, s, e)| {
                self.chrom_names.iter().position(|n| n == c)
                    .map(|id| (id as u32, *s, *e))
            })
            .collect();
        self.spatial.query_regions(&chrom_regions)
    }
    
    /// Find reads that overlap with any of the given loci regions using spatial index
    fn find_reads_overlapping_loci(&self, loci: &[(String, i64, i64, HashSet<String>)]) -> HashSet<String> {
        let mut result = HashSet::new();
        
        // Build locus regions by chromosome
        let mut regions_by_chrom: HashMap<u32, Vec<(i64, i64)>> = HashMap::new();
        for (chrom, start, end, _) in loci {
            if let Some(chrom_id) = self.chrom_names.iter().position(|n| n == chrom) {
                regions_by_chrom.entry(chrom_id as u32)
                    .or_default()
                    .push((*start, *end));
            }
        }
        
        // Scan index for reads overlapping any locus
        for (read_name, alignments) in &self.index {
            'aln_loop: for aln in alignments {
                let aln_start = aln.start as i64;
                let aln_end = aln.end as i64;
                
                if let Some(regions) = regions_by_chrom.get(&aln.chrom_id) {
                    for (r_start, r_end) in regions {
                        // Simple overlap check
                        if aln_start < *r_end && aln_end > *r_start {
                            result.insert(read_name.clone());
                            break 'aln_loop;
                        }
                    }
                }
            }
        }
        
        result
    }
}

/// Build coverage profile from index
fn build_coverage_profile(index: &BamIndex, chrom: &str, start: i64, end: i64, reads: &HashSet<String>, bin_size: i64) -> Vec<f64> {
    let chrom_id = match index.chrom_names.iter().position(|n| n == chrom) {
        Some(id) => id as u32,
        None => return vec![],
    };

    let num_bins = ((end - start) / bin_size + 1).max(1) as usize;
    let mut coverage = vec![0.0; num_bins];

    for read_name in reads {
        if let Some(alignments) = index.index.get(read_name) {
            for aln in alignments {
                if aln.chrom_id != chrom_id { continue; }
                let aln_start = aln.start as i64;
                let aln_end = aln.end as i64;
                if aln_end < start || aln_start > end { continue; }
                
                let s = ((aln_start.max(start) - start) / bin_size) as usize;
                let e = ((aln_end.min(end) - start) / bin_size) as usize;
                for b in s..=e.min(num_bins - 1) {
                    coverage[b] += 1.0;
                }
            }
        }
    }
    coverage
}

/// Find valleys in coverage profile
fn find_valleys(coverage: &[f64], valley_frac: f64) -> Vec<usize> {
    let max_cov = coverage.iter().copied().fold(0.0, f64::max);
    if max_cov < 2.0 { return vec![]; }
    
    let threshold = max_cov * valley_frac;
    let mut valleys = vec![];
    let mut in_valley = false;
    let mut valley_start = 0;
    
    for (i, &cov) in coverage.iter().enumerate() {
        if cov < threshold && !in_valley {
            in_valley = true;
            valley_start = i;
        } else if cov >= threshold && in_valley {
            in_valley = false;
            valleys.push((valley_start + i) / 2);
        }
    }
    valleys
}

/// Split locus by coverage valleys
fn split_locus(chrom: &str, start: i64, end: i64, reads: &HashSet<String>, index: &BamIndex, valley_frac: f64, min_seg_bp: i64, bin_size: i64) -> Vec<(String, i64, i64, HashSet<String>)> {
    let coverage = build_coverage_profile(index, chrom, start, end, reads, bin_size);
    if coverage.is_empty() {
        return vec![(chrom.to_string(), start, end, reads.clone())];
    }
    
    let valleys = find_valleys(&coverage, valley_frac);
    if valleys.is_empty() {
        return vec![(chrom.to_string(), start, end, reads.clone())];
    }
    
    let mut segments = vec![];
    let mut seg_start = start;
    
    for &valley_bin in &valleys {
        let valley_pos = start + valley_bin as i64 * bin_size;
        if valley_pos - seg_start >= min_seg_bp {
            let seg_reads: HashSet<String> = reads.iter()
                .filter(|r| index.index.get(*r).map_or(false, |alns| 
                    alns.iter().any(|a| a.chrom_id as usize == index.chrom_names.iter().position(|n| n == chrom).unwrap_or(999) 
                        && (a.start as i64) < valley_pos && (a.end as i64) > seg_start)))
                .cloned()
                .collect();
            if !seg_reads.is_empty() {
                segments.push((chrom.to_string(), seg_start, valley_pos, seg_reads));
            }
        }
        seg_start = valley_pos;
    }
    
    if end - seg_start >= min_seg_bp {
        let seg_reads: HashSet<String> = reads.iter()
            .filter(|r| index.index.get(*r).map_or(false, |alns| 
                alns.iter().any(|a| a.chrom_id as usize == index.chrom_names.iter().position(|n| n == chrom).unwrap_or(999) 
                    && (a.start as i64) < end && (a.end as i64) > seg_start)))
            .cloned()
            .collect();
        if !seg_reads.is_empty() {
            segments.push((chrom.to_string(), seg_start, end, seg_reads));
        }
    }
    
    if segments.is_empty() { vec![(chrom.to_string(), start, end, reads.clone())] } else { segments }
}

/// Cluster alignments into loci
fn cluster_loci(by_chrom: HashMap<String, Vec<(i64, i64, String)>>, cluster_distance: i64, min_reads: usize, min_size_bp: i64, max_size_bp: i64) -> Vec<(String, i64, i64, HashSet<String>)> {
    let mut loci = vec![];

    for (chrom, mut aligns) in by_chrom {
        if aligns.len() < min_reads { continue; }
        aligns.sort_by_key(|a| a.0);

        let mut current_reads: HashSet<String> = [aligns[0].2.clone()].iter().cloned().collect();
        let mut current_start = aligns[0].0;
        let mut current_end = aligns[0].1;

        for (start, end, name) in aligns.into_iter().skip(1) {
            if start - current_end <= cluster_distance {
                current_end = current_end.max(end);
                current_reads.insert(name);
            } else {
                let span = current_end - current_start;
                if current_reads.len() >= min_reads && span >= min_size_bp && span <= max_size_bp {
                    loci.push((chrom.clone(), current_start, current_end, current_reads.clone()));
                }
                current_reads = [name].iter().cloned().collect();
                current_start = start;
                current_end = end;
            }
        }

        let span = current_end - current_start;
        if current_reads.len() >= min_reads && span >= min_size_bp && span <= max_size_bp {
            loci.push((chrom, current_start, current_end, current_reads));
        }
    }
    loci
}

fn collect_seed_reads(index: &BamIndex, seed_chrom: &str, seed_start: i64, seed_end: i64) -> HashSet<String> {
    let mut seed_reads = HashSet::new();
    let chrom_id = index.chrom_names.iter().position(|n| n == seed_chrom).map(|i| i as u32);
    if chrom_id.is_none() { return seed_reads; }
    let chrom_id = chrom_id.unwrap();

    for (read_name, alignments) in &index.index {
        for aln in alignments {
            if aln.chrom_id == chrom_id && (aln.start as i64) < seed_end && (aln.end as i64) > seed_start {
                seed_reads.insert(read_name.clone());
                break;
            }
        }
    }
    seed_reads
}

/// Load multiple seeds from BED file
fn load_seeds_bed(path: &str) -> Result<Vec<(String, i64, i64, String)>> {
    let mut seeds = vec![];
    let content = std::fs::read_to_string(path)?;
    
    for line in content.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() >= 3 {
            let chrom = parts[0].to_string();
            let start: i64 = parts[1].parse()?;
            let end: i64 = parts[2].parse()?;
            let name = if parts.len() >= 4 { parts[3].to_string() } else { format!("seed_{}", seeds.len()) };
            seeds.push((chrom, start, end, name));
        }
    }
    Ok(seeds)
}

/// Collect seed reads from multiple seeds (pooled)
fn collect_multi_seed_reads(index: &BamIndex, seeds: &[(String, i64, i64, String)]) -> HashSet<String> {
    let mut all_reads = HashSet::new();
    for (chrom, start, end, _name) in seeds {
        let reads = collect_seed_reads(index, chrom, *start, *end);
        all_reads.extend(reads);
    }
    all_reads
}

/// Detect loci with transitive discovery
fn detect_with_transitive_discovery(
    index: &BamIndex,
    initial_seeds: &HashSet<String>,
    cluster_distance: i64,
    min_reads: usize,
    min_size: i64,
    max_size: i64,
    split_valleys: bool,
    valley_frac: f64,
    min_seg: i64,
    bin_size: i64,
    max_iterations: usize,
) -> Vec<(String, i64, i64, HashSet<String>)> {
    let mut all_loci: Vec<(String, i64, i64, HashSet<String>)> = vec![];
    let mut known_regions: Vec<(String, i64, i64)> = vec![];
    let mut current_reads = initial_seeds.clone();
    
    for iteration in 0..max_iterations {
        println!("\n  Iteration {}: {} seed reads", iteration + 1, current_reads.len());
        
        // Find alignments
        let by_chrom = index.get_alignments_by_chrom(&current_reads);
        
        // Filter out known regions
        let mut new_aligns: HashMap<String, Vec<(i64, i64, String)>> = HashMap::new();
        for (chrom, aligns) in by_chrom {
            let filtered: Vec<_> = aligns.into_iter()
                .filter(|(s, e, _)| {
                    !known_regions.iter().any(|(kc, ks, ke)| {
                        kc == &chrom && s < ke && e > ks
                    })
                })
                .collect();
            if !filtered.is_empty() {
                new_aligns.insert(chrom, filtered);
            }
        }
        
        if new_aligns.is_empty() {
            println!("    No new alignments found");
            break;
        }
        
        // Cluster new loci
        let new_loci = cluster_loci(new_aligns, cluster_distance, min_reads, min_size, max_size);
        println!("    Found {} new loci", new_loci.len());
        
        if new_loci.is_empty() {
            break;
        }
        
        // Add to known regions and all loci
        for (chrom, s, e, reads) in &new_loci {
            known_regions.push((chrom.clone(), *s, *e));
            all_loci.push((chrom.clone(), *s, *e, reads.clone()));
        }
        
        // TRANSITIVE DISCOVERY: Use spatial index to find reads near discovered loci
        // Query expansion regions around top loci to find connected reads
        let max_loci_to_expand = 5;
        let expansion_window = 20000i64; // 20kb expansion window
        
        // Collect top loci by read count
        let mut loci_sorted: Vec<_> = new_loci.iter()
            .map(|(c, s, e, r)| (r.len(), c.clone(), *s, *e))
            .collect();
        loci_sorted.sort_by_key(|(n, _, _, _)| std::cmp::Reverse(*n));
        
        // Build expansion regions from top loci
        let expansion_regions: Vec<(String, i64, i64)> = loci_sorted.iter().take(max_loci_to_expand)
            .map(|(_, c, s, e)| {
                let expand_start = (*s - expansion_window).max(0);
                let expand_end = *e + expansion_window;
                (c.clone(), expand_start, expand_end)
            })
            .collect();
        
        // Use spatial index to efficiently query reads in expansion regions
        let nearby_reads = index.query_regions(&expansion_regions);
        let old_count = current_reads.len();
        for read in nearby_reads {
            current_reads.insert(read);
        }
        let new_count = current_reads.len() - old_count;
        
        println!("    Added {} reads from spatial query ({}kb window, total pool: {})", 
                 new_count, expansion_window / 1000, current_reads.len());
        
        if new_count == 0 {
            println!("    No new reads discovered, stopping");
            break;
        }
    }
    
    // Apply valley splitting if requested
    if split_valleys {
        let mut split_loci = vec![];
        for (chrom, s, e, reads) in all_loci {
            let segs = split_locus(&chrom, s, e, &reads, index, valley_frac, min_seg, bin_size);
            split_loci.extend(segs);
        }
        all_loci = split_loci;
    }
    
    all_loci
}

#[derive(Parser, Debug)]
#[command(name = "gene_family_detector_fast")]
struct Args {
    #[arg(short, long)] bam: String,
    #[arg(long)] seed_chrom: Option<String>,
    #[arg(long)] seed_start: Option<i64>,
    #[arg(long)] seed_end: Option<i64>,
    #[arg(short, long)] output: String,
    #[arg(long, default_value_t = 2)] min_reads: usize,
    #[arg(long, default_value_t = 10000)] min_size_bp: i64,
    #[arg(long, default_value_t = 2000000)] max_size_bp: i64,
    #[arg(long, default_value_t = 10000)] cluster_distance: i64,
    #[arg(long)] include_secondary: bool,
    #[arg(long)] include_supplementary: bool,
    #[arg(long)] rebuild_index: bool,
    #[arg(long)] split_valleys: bool,
    #[arg(long, default_value_t = 0.15)] valley_frac: f64,
    #[arg(long, default_value_t = 5000)] min_segment_bp: i64,
    #[arg(long, default_value_t = 100)] bin_size: i64,
    /// BED file with multiple seeds (enables transitive discovery)
    #[arg(long)] seeds_bed: Option<String>,
    /// Disable transitive discovery (default is transitive on).
    #[arg(long = "no-transitive", default_value_t = false, action = clap::ArgAction::SetTrue)]
    no_transitive: bool,
    /// Ignored (transitive is default). Kept so older scripts that pass `--transitive` still run.
    #[allow(dead_code)]
    #[arg(long = "transitive", hide = true, action = clap::ArgAction::SetTrue)]
    _transitive_deprecated_noop: bool,
    /// Max iterations for transitive discovery
    #[arg(long, default_value_t = 2)] max_iterations: usize,
}

fn main() -> Result<()> {
    let args = Args::parse();

    // Validate arguments
    let has_single_seed = args.seed_chrom.is_some() && args.seed_start.is_some() && args.seed_end.is_some();
    let has_multi_seed = args.seeds_bed.is_some();
    
    if !has_single_seed && !has_multi_seed {
        anyhow::bail!("Either provide --seed-chrom/--seed-start/--seed-end OR --seeds-bed <file>");
    }

    println!("{}", "=".repeat(80));
    println!("FAST GENE FAMILY DETECTOR");
    println!("{}", "=".repeat(80));
    
    // Handle multi-seed mode
    let multi_seed = args.seeds_bed.is_some();
    let seeds: Vec<(String, i64, i64, String)> = if let Some(ref bed_path) = args.seeds_bed {
        println!("Multi-seed mode: loading from {}", bed_path);
        load_seeds_bed(bed_path)?
    } else {
        vec![(args.seed_chrom.clone().unwrap(), args.seed_start.unwrap(), args.seed_end.unwrap(), "seed".to_string())]
    };
    
    if multi_seed {
        println!("Seeds: {}", seeds.len());
        for (i, (chrom, start, end, name)) in seeds.iter().enumerate() {
            println!("  {}: {}:{}-{} ({})", i + 1, chrom, start, end, name);
        }
    } else {
        println!("Seed: {}:{}-{}", 
            args.seed_chrom.as_ref().unwrap(), 
            args.seed_start.unwrap(), 
            args.seed_end.unwrap());
    }
    
    println!("Min: {} reads, {} bp | Max: {} bp", args.min_reads, args.min_size_bp, args.max_size_bp);
    if !args.no_transitive {
        println!("Transitive discovery: ON ({} iterations)", args.max_iterations);
    }
    if args.split_valleys {
        println!("Valley splitting: ON (frac={}, min_seg={}bp)", args.valley_frac, args.min_segment_bp);
    }
    println!();

    let start = std::time::Instant::now();

    let index = if args.rebuild_index {
        let idx = BamIndex::build(&args.bam, args.include_secondary, args.include_supplementary)?;
        idx.save(&BamIndex::index_path(&args.bam))?;
        idx
    } else {
        BamIndex::load_or_build(&args.bam, args.include_secondary, args.include_supplementary)?
    };

    println!("\nCollecting seed reads...");
    let seed_reads = if multi_seed {
        collect_multi_seed_reads(&index, &seeds)
    } else {
        collect_seed_reads(&index, &args.seed_chrom.clone().unwrap(), args.seed_start.unwrap(), args.seed_end.unwrap())
    };
    println!("  Found {} seed reads", seed_reads.len());

    if seed_reads.is_empty() {
        std::fs::write(&args.output, "# No loci found\n")?;
        return Ok(());
    }

    // Detection with optional transitive discovery
    let loci = if !args.no_transitive {
        println!("\nRunning transitive discovery...");
        detect_with_transitive_discovery(
            &index, &seed_reads,
            args.cluster_distance, args.min_reads, args.min_size_bp, args.max_size_bp,
            args.split_valleys, args.valley_frac, args.min_segment_bp, args.bin_size,
            args.max_iterations
        )
    } else {
        println!("\nFinding shared reads...");
        let by_chrom = index.get_alignments_by_chrom(&seed_reads);
        println!("  Found alignments on {} chromosomes", by_chrom.len());

        println!("\nClustering loci...");
        let raw_loci = cluster_loci(by_chrom, args.cluster_distance, args.min_reads, args.min_size_bp, args.max_size_bp);
        println!("  Found {} raw loci", raw_loci.len());

        // Split by valleys
        if args.split_valleys {
            println!("\nSplitting by coverage valleys...");
            let mut split = vec![];
            let mut split_count = 0;
            for (chrom, s, e, reads) in raw_loci {
                let segs = split_locus(&chrom, s, e, &reads, &index, args.valley_frac, args.min_segment_bp, args.bin_size);
                if segs.len() > 1 { split_count += 1; }
                split.extend(segs);
            }
            println!("  Split {} loci into {} segments", split_count, split.len());
            split
        } else {
            raw_loci
        }
    };

    println!("\nWriting output...");
    let mut file = File::create(&args.output)?;
    writeln!(file, "# Fast Gene Family Detection Results")?;
    if multi_seed {
        writeln!(file, "# Seeds: {} (pooled)", seeds.len())?;
    } else {
        writeln!(file, "# Seed: {}:{}-{}", 
            args.seed_chrom.as_ref().unwrap(), 
            args.seed_start.unwrap(), 
            args.seed_end.unwrap())?;
    }
    writeln!(file, "# Seed reads: {}", seed_reads.len())?;
    writeln!(file, "# Loci found: {}", loci.len())?;
    if !args.no_transitive {
        writeln!(file, "# Transitive discovery: {} iterations", args.max_iterations)?;
    }
    writeln!(file, "#chrom\tstart\tend\tname\tshared_reads\tspan_bp")?;

    for (i, (chrom, s, e, reads)) in loci.iter().enumerate() {
        writeln!(file, "{}\t{}\t{}\tlocus_{}\t{}\t{}", chrom, s, e, i + 1, reads.len(), e - s)?;
    }

    println!("\n{}", "=".repeat(80));
    println!("COMPLETE: {} loci in {:.1}s", loci.len(), start.elapsed().as_secs_f64());
    println!("{}", "=".repeat(80));

    Ok(())
}
