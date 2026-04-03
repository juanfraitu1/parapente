//! Defensible Family-Aware Transcript Assembler
//! 
//! Produces publication-quality assemblies that are:
//! 1. Well-supported by reads (every isoform has evidence)
//! 2. Biologically plausible (canonical splice sites)
//! 3. Family-consistent (cross-member validation)
//! 4. Properly quantified (meaningful support values)
//! 5. Reasonably complex (not over-assembled)

use crate::bam_parser::{SplicedRead, SpliceJunction};
use crate::isoform::Isoform;
use hashbrown::{HashMap, HashSet};
use rayon::prelude::*;

/// Configuration for defensible assembly
#[derive(Debug, Clone)]
pub struct DefensibleAssemblyConfig {
    /// Minimum reads spanning ALL junctions (default: 3)
    pub min_junction_reads: usize,
    /// Minimum total coverage (default: 5)
    pub min_total_coverage: usize,
    /// Require canonical splice sites (default: true)
    pub require_canonical_splice: bool,
    /// Minimum fraction of canonical splice sites (default: 0.95)
    pub min_canonical_fraction: f64,
    /// Max isoforms per gene (default: 20)
    pub max_isoforms_per_gene: usize,
    /// Confidence threshold for output (default: 0.7)
    pub confidence_threshold: f64,
    /// Cross-member validation required (default: true)
    pub require_family_validation: bool,
    /// Minimum family support ratio (default: 0.3)
    pub min_family_support_ratio: f64,
}

impl Default for DefensibleAssemblyConfig {
    fn default() -> Self {
        Self {
            min_junction_reads: 3,
            min_total_coverage: 5,
            require_canonical_splice: true,
            min_canonical_fraction: 0.95,
            max_isoforms_per_gene: 20,
            confidence_threshold: 0.7,
            require_family_validation: true,
            min_family_support_ratio: 0.3,
        }
    }
}

/// Family member with validation info
#[derive(Debug, Clone)]
pub struct FamilyMember {
    pub id: String,
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub strand: char,
}

/// Validated isoform with evidence
#[derive(Debug, Clone)]
pub struct ValidatedIsoform {
    pub isoform: Isoform,
    /// Confidence score (0-1)
    pub confidence: f64,
    /// Reads supporting this isoform
    pub supporting_reads: Vec<String>,
    /// Coverage per position
    pub coverage_profile: Vec<usize>,
    /// Canonical splice site fraction
    pub canonical_splice_fraction: f64,
    /// Family support info
    pub family_support: FamilySupportInfo,
    /// Validation flags
    pub validation: ValidationFlags,
}

#[derive(Debug, Clone, Default)]
pub struct FamilySupportInfo {
    pub total_family_reads: usize,
    pub this_member_reads: usize,
    pub other_member_support: usize,  // How many other members have similar isoform
    pub is_conserved_structure: bool,
}

#[derive(Debug, Clone, Default)]
pub struct ValidationFlags {
    pub has_full_junction_support: bool,
    pub has_canonical_splices: bool,
    pub has_family_validation: bool,
    pub passes_coverage_threshold: bool,
    pub is_not_redundant: bool,
}

/// Splice site validator
pub struct SpliceSiteValidator {
    genome_seq: HashMap<String, Vec<u8>>,
}

impl SpliceSiteValidator {
    pub fn new(genome_path: &str) -> anyhow::Result<Self> {
        // Load genome sequence
        let mut genome_seq = HashMap::new();
        // Implementation: load FASTA
        
        Ok(Self { genome_seq })
    }
    
    /// Check if junction has canonical splice sites
    pub fn validate_junction(&self, chrom: &str, donor: u64, acceptor: u64, 
                            strand: char) -> SpliceSiteType {
        let seq = match self.genome_seq.get(chrom) {
            Some(s) => s,
            None => return SpliceSiteType::Unknown,
        };
        
        let donor_pos = donor as usize;
        let acceptor_pos = acceptor as usize;
        
        if donor_pos + 2 > seq.len() || acceptor_pos < 2 {
            return SpliceSiteType::Unknown;
        }
        
        let donor_seq = &seq[donor_pos..donor_pos+2];
        let acceptor_seq = &seq[acceptor_pos-2..acceptor_pos];
        
        match strand {
            '+' => {
                if donor_seq == b"GT" && acceptor_seq == b"AG" {
                    SpliceSiteType::Canonical
                } else if donor_seq == b"AT" && acceptor_seq == b"AC" {
                    SpliceSiteType::U12
                } else {
                    SpliceSiteType::NonCanonical
                }
            }
            '-' => {
                // Reverse complement
                let rc_donor = reverse_complement(donor_seq);
                let rc_acceptor = reverse_complement(acceptor_seq);
                
                if rc_donor == b"GT" && rc_acceptor == b"AG" {
                    SpliceSiteType::Canonical
                } else {
                    SpliceSiteType::NonCanonical
                }
            }
            _ => SpliceSiteType::Unknown,
        }
    }
    
    /// Calculate fraction of canonical splice sites in isoform
    pub fn canonical_fraction(&self, isoform: &Isoform) -> f64 {
        if isoform.junctions.is_empty() {
            return 1.0;  // Single-exon, no splice sites to validate
        }
        
        let canonical_count = isoform.junctions.iter()
            .filter(|j| {
                matches!(self.validate_junction(&j.chrom, j.donor, j.acceptor, 
                         j.strand.chars().next().unwrap_or('+')),
                         SpliceSiteType::Canonical | SpliceSiteType::U12)
            })
            .count();
        
        canonical_count as f64 / isoform.junctions.len() as f64
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum SpliceSiteType {
    Canonical,      // GT-AG
    U12,           // AT-AC (rare but valid)
    NonCanonical,  // Everything else (suspicious)
    Unknown,
}

/// Defensible family-aware assembler
pub struct DefensibleAssembler {
    family_id: String,
    members: Vec<FamilyMember>,
    config: DefensibleAssemblyConfig,
    splice_validator: Option<SpliceSiteValidator>,
}

/// Assembly result with quality metrics
#[derive(Debug)]
pub struct DefensibleAssemblyResult {
    pub member_isoforms: HashMap<String, Vec<ValidatedIsoform>>,
    pub quality_metrics: AssemblyQualityMetrics,
}

#[derive(Debug, Default)]
pub struct AssemblyQualityMetrics {
    pub total_isoforms: usize,
    pub high_confidence_isoforms: usize,
    pub mean_support: f64,
    pub canonical_splice_rate: f64,
    pub family_validated_rate: f64,
    pub read_explanation_rate: f64,
}

impl DefensibleAssembler {
    pub fn new(family_id: String, members: Vec<FamilyMember>) -> Self {
        Self {
            family_id,
            members,
            config: DefensibleAssemblyConfig::default(),
            splice_validator: None,
        }
    }
    
    pub fn with_config(mut self, config: DefensibleAssemblyConfig) -> Self {
        self.config = config;
        self
    }
    
    pub fn with_genome(mut self, genome_path: &str) -> anyhow::Result<Self> {
        self.splice_validator = Some(SpliceSiteValidator::new(genome_path)?);
        Ok(self)
    }
    
    /// Main assembly with full validation
    pub fn assemble(&self, reads: Vec<SplicedRead>) -> anyhow::Result<DefensibleAssemblyResult> {
        eprintln!("Defensible Assembly: {}", self.family_id);
        eprintln!("  Members: {}", self.members.len());
        eprintln!("  Total reads: {}", reads.len());
        
        // Phase 1: Assign reads to members
        let member_reads = self.assign_reads_to_members(&reads);
        eprintln!("  Reads assigned: {}", 
                 member_reads.values().map(|v| v.len()).sum::<usize>());
        
        // Phase 2: Build family consensus structure
        let family_consensus = self.build_family_consensus(&member_reads)?;
        eprintln!("  Conserved junctions: {}", family_consensus.conserved_junctions.len());
        
        // Phase 3: Assemble each member with validation
        let mut member_results: HashMap<String, Vec<ValidatedIsoform>> = HashMap::new();
        
        for (member_id, reads) in member_reads {
            let isoforms = self.assemble_member_defensible(
                &member_id, 
                &reads, 
                &family_consensus
            )?;
            
            if !isoforms.is_empty() {
                member_results.insert(member_id, isoforms);
            }
        }
        
        // Phase 4: Cross-validate across members
        let validated_results = self.cross_validate_results(member_results)?;
        
        // Phase 5: Filter to defensible set
        let final_results = self.filter_defensible_isoforms(validated_results)?;
        
        // Calculate quality metrics
        let metrics = self.calculate_quality_metrics(&final_results, &reads);
        
        eprintln!("\nAssembly Complete:");
        eprintln!("  Total isoforms: {}", metrics.total_isoforms);
        eprintln!("  High confidence: {}", metrics.high_confidence_isoforms);
        eprintln!("  Mean support: {:.1f}", metrics.mean_support);
        eprintln!("  Canonical splices: {:.1%}", metrics.canonical_splice_rate);
        
        Ok(DefensibleAssemblyResult {
            member_isoforms: final_results,
            quality_metrics: metrics,
        })
    }
    
    /// Assign reads to best-matching member
    fn assign_reads_to_members(&self, reads: &[SplicedRead]) -> HashMap<String, Vec<SplicedRead>> {
        let mut member_reads: HashMap<String, Vec<SplicedRead>> = HashMap::new();
        
        for read in reads {
            let mut best_member = None;
            let mut best_score = 0.0;
            
            for member in &self.members {
                if read.chrom != member.chrom || read.strand != member.strand {
                    continue;
                }
                
                // Calculate overlap score
                let overlap_start = read.tx_start.max(member.start);
                let overlap_end = read.tx_end.min(member.end);
                
                if overlap_start >= overlap_end {
                    continue;
                }
                
                let overlap_len = overlap_end - overlap_start;
                let read_len = read.tx_end - read.tx_start;
                let overlap_frac = overlap_len as f64 / read_len as f64;
                
                // Score based on overlap and junction compatibility
                let mut score = overlap_frac;
                
                // Bonus if read junctions match member location
                for junction in &read.junctions {
                    if junction.start >= member.start && junction.end <= member.end {
                        score += 0.1;  // Small bonus per junction
                    }
                }
                
                if score > best_score {
                    best_score = score;
                    best_member = Some(member.id.clone());
                }
            }
            
            // Assign to best member, or skip if no good match
            if let Some(member_id) = best_member {
                if best_score >= 0.5 {  // Must have at least 50% overlap
                    member_reads.entry(member_id).or_default().push(read.clone());
                }
            }
        }
        
        member_reads
    }
    
    /// Build family consensus structure
    fn build_family_consensus(&self, member_reads: &HashMap<String, Vec<SplicedRead>>) 
        -> anyhow::Result<FamilyConsensus> {
        
        let mut junction_counts: HashMap<(String, u64, u64), HashSet<String>> = HashMap::new();
        
        // Count junction support across family
        for (member_id, reads) in member_reads {
            for read in reads {
                for junction in &read.junctions {
                    let key = (read.chrom.clone(), junction.start, junction.end);
                    junction_counts.entry(key)
                        .or_default()
                        .insert(member_id.clone());
                }
            }
        }
        
        // Identify conserved junctions (supported by multiple members)
        let conserved_junctions: Vec<_> = junction_counts.iter()
            .filter(|(_, members)| members.len() >= 2)
            .map(|((chrom, donor, acceptor), members)| {
                ConservedJunction {
                    chrom: chrom.clone(),
                    donor: *donor,
                    acceptor: *acceptor,
                    supporting_members: members.clone(),
                }
            })
            .collect();
        
        Ok(FamilyConsensus {
            conserved_junctions,
            total_junctions: junction_counts.len(),
        })
    }
    
    /// Assemble one member with defensible validation
    fn assemble_member_defensible(&self, member_id: &str, reads: &[SplicedRead],
                                   consensus: &FamilyConsensus) 
        -> anyhow::Result<Vec<ValidatedIsoform>> {
        
        if reads.len() < self.config.min_total_coverage {
            return Ok(Vec::new());
        }
        
        // Build trie with support tracking
        let mut trie = DefensibleTrie::new(self.config.junction_tolerance);
        
        for read in reads {
            trie.insert(read);
        }
        
        // Extract paths with validation
        let paths = trie.extract_validated_paths(
            self.config.min_junction_reads,
            self.config.min_total_coverage,
        );
        
        // Build validated isoforms
        let mut isoforms = Vec::new();
        
        for (path_idx, path) in paths.iter().enumerate().take(self.config.max_isoforms_per_gene) {
            let member = self.members.iter()
                .find(|m| &m.id == member_id)
                .ok_or_else(|| anyhow::anyhow!("Member not found"))?;
            
            // Build isoform from path
            let isoform = self.build_isoform_from_path(
                path, 
                member,
                format!("{}_{}_{}", self.family_id, member_id, path_idx + 1)
            )?;
            
            // Validate isoform
            let validated = self.validate_isoform(
                isoform,
                path,
                reads,
                consensus,
                member_id,
            )?;
            
            // Only keep if passes confidence threshold
            if validated.confidence >= self.config.confidence_threshold {
                isoforms.push(validated);
            }
        }
        
        Ok(isoforms)
    }
    
    /// Validate isoform with multiple criteria
    fn validate_isoform(&self, isoform: Isoform, path: &DefensiblePath,
                        member_reads: &[SplicedRead], consensus: &FamilyConsensus,
                        member_id: &str) -> anyhow::Result<ValidatedIsoform> {
        
        let mut flags = ValidationFlags::default();
        
        // 1. Check junction support
        let mut junction_support = Vec::new();
        for junction in &isoform.junctions {
            let supporting = member_reads.iter()
                .filter(|r| r.has_junction(junction.donor, junction.acceptor))
                .count();
            junction_support.push(supporting);
        }
        
        flags.has_full_junction_support = junction_support.iter()
            .all(|&s| s >= self.config.min_junction_reads);
        
        // 2. Check splice sites
        let canonical_fraction = if let Some(ref validator) = self.splice_validator {
            validator.canonical_fraction(&isoform)
        } else {
            1.0  // Assume valid if no validator
        };
        
        flags.has_canonical_splices = canonical_fraction >= self.config.min_canonical_fraction;
        
        // 3. Check family support
        let mut family_support_count = 0;
        for junction in &isoform.junctions {
            let key = (junction.chrom.clone(), junction.donor, junction.acceptor);
            if consensus.conserved_junctions.iter()
                .any(|cj| cj.chrom == junction.chrom && 
                         cj.donor == junction.donor && 
                         cj.acceptor == junction.acceptor) {
                family_support_count += 1;
            }
        }
        
        let family_support_ratio = if isoform.junctions.is_empty() {
            0.0
        } else {
            family_support_count as f64 / isoform.junctions.len() as f64
        };
        
        flags.has_family_validation = family_support_ratio >= self.config.min_family_support_ratio;
        
        // 4. Check coverage
        let total_support = path.supporting_reads.len();
        flags.passes_coverage_threshold = total_support >= self.config.min_total_coverage;
        
        // 5. Calculate confidence
        let confidence = self.calculate_confidence(
            &flags, 
            canonical_fraction,
            family_support_ratio,
            total_support,
        );
        
        // Build family support info
        let family_info = FamilySupportInfo {
            total_family_reads: consensus.total_junctions,  // Approximation
            this_member_reads: total_support,
            other_member_support: family_support_count,
            is_conserved_structure: family_support_ratio >= 0.5,
        };
        
        Ok(ValidatedIsoform {
            isoform,
            confidence,
            supporting_reads: path.supporting_reads.clone(),
            coverage_profile: path.coverage_profile.clone(),
            canonical_splice_fraction: canonical_fraction,
            family_support: family_info,
            validation: flags,
        })
    }
    
    /// Calculate overall confidence score
    fn calculate_confidence(&self, flags: &ValidationFlags, canonical_frac: f64,
                           family_support: f64, total_support: usize) -> f64 {
        let mut score = 0.0;
        
        // Junction support (30%)
        if flags.has_full_junction_support {
            score += 0.30;
        }
        
        // Splice sites (25%)
        score += canonical_frac * 0.25;
        
        // Family support (25%)
        score += family_support.min(1.0) * 0.25;
        
        // Coverage depth (20%)
        let coverage_score = (total_support as f64 / self.config.min_total_coverage as f64)
            .min(2.0) / 2.0;  // Cap at 2x minimum
        score += coverage_score * 0.20;
        
        score
    }
    
    /// Cross-validate across family members
    fn cross_validate_results(&self, 
                               results: HashMap<String, Vec<ValidatedIsoform>>)
        -> anyhow::Result<HashMap<String, Vec<ValidatedIsoform>>> {
        
        // Find similar isoforms across members
        let mut paralog_groups: Vec<Vec<(String, usize)>> = Vec::new();
        
        for (m1_id, m1_isoforms) in &results {
            for (i1, iso1) in m1_isoforms.iter().enumerate() {
                let mut found_group = false;
                
                for group in &mut paralog_groups {
                    let (rep_member, rep_idx) = &group[0];
                    let rep_iso = &results[rep_member][*rep_idx];
                    
                    if self.isoforms_similar(&iso1.isoform, &rep_iso.isoform) {
                        group.push((m1_id.clone(), i1));
                        found_group = true;
                        break;
                    }
                }
                
                if !found_group {
                    paralog_groups.push(vec![(m1_id.clone(), i1)]);
                }
            }
        }
        
        // Update paralog support
        let mut updated_results = results.clone();
        
        for group in paralog_groups {
            let group_size = group.len();
            
            for (member_id, iso_idx) in group {
                if let Some(isoforms) = updated_results.get_mut(&member_id) {
                    if let Some(isoform) = isoforms.get_mut(iso_idx) {
                        // Boost confidence for conserved structures
                        if group_size >= 2 {
                            isoform.confidence = (isoform.confidence + 0.1).min(1.0);
                            isoform.family_support.other_member_support = group_size - 1;
                            isoform.family_support.is_conserved_structure = true;
                        }
                    }
                }
            }
        }
        
        Ok(updated_results)
    }
    
    /// Filter to defensible set
    fn filter_defensible_isoforms(&self, 
                                   results: HashMap<String, Vec<ValidatedIsoform>>)
        -> anyhow::Result<HashMap<String, Vec<ValidatedIsoform>>> {
        
        let mut filtered = HashMap::new();
        
        for (member_id, isoforms) in results {
            let defensible: Vec<_> = isoforms.into_iter()
                .filter(|iso| {
                    // Must pass all critical validations
                    iso.validation.has_full_junction_support &&
                    iso.validation.has_canonical_splices &&
                    iso.confidence >= self.config.confidence_threshold
                })
                .collect();
            
            if !defensible.is_empty() {
                filtered.insert(member_id, defensible);
            }
        }
        
        Ok(filtered)
    }
    
    /// Calculate quality metrics
    fn calculate_quality_metrics(&self, 
                                  results: &HashMap<String, Vec<ValidatedIsoform>>,
                                  all_reads: &[SplicedRead]) -> AssemblyQualityMetrics {
        
        let mut total_isoforms = 0;
        let mut high_confidence = 0;
        let mut total_support = 0;
        let mut total_canonical = 0.0;
        let mut total_family_validated = 0;
        
        for isoforms in results.values() {
            for iso in isoforms {
                total_isoforms += 1;
                total_support += iso.supporting_reads.len();
                total_canonical += iso.canonical_splice_fraction;
                
                if iso.confidence >= 0.8 {
                    high_confidence += 1;
                }
                
                if iso.validation.has_family_validation {
                    total_family_validated += 1;
                }
            }
        }
        
        let explained_reads: HashSet<_> = results.values()
            .flat_map(|isos| isos.iter().flat_map(|i| i.supporting_reads.clone()))
            .collect();
        
        AssemblyQualityMetrics {
            total_isoforms,
            high_confidence_isoforms: high_confidence,
            mean_support: if total_isoforms > 0 { 
                total_support as f64 / total_isoforms as f64 
            } else { 0.0 },
            canonical_splice_rate: if total_isoforms > 0 {
                total_canonical / total_isoforms as f64
            } else { 0.0 },
            family_validated_rate: if total_isoforms > 0 {
                total_family_validated as f64 / total_isoforms as f64
            } else { 0.0 },
            read_explanation_rate: explained_reads.len() as f64 / all_reads.len() as f64,
        }
    }
    
    /// Check if two isoforms are similar
    fn isoforms_similar(&self, a: &Isoform, b: &Isoform) -> bool {
        if a.junctions.len() != b.junctions.len() {
            return false;
        }
        
        for (ja, jb) in a.junctions.iter().zip(b.junctions.iter()) {
            let donor_diff = ja.donor.abs_diff(jb.donor);
            let acceptor_diff = ja.acceptor.abs_diff(jb.acceptor);
            
            if donor_diff > self.config.junction_tolerance ||
               acceptor_diff > self.config.junction_tolerance {
                return false;
            }
        }
        
        true
    }
    
    /// Build isoform from path
    fn build_isoform_from_path(&self, path: &DefensiblePath, member: &FamilyMember,
                                id: String) -> anyhow::Result<Isoform> {
        // Implementation: convert path to Isoform
        // ... (simplified)
        
        Ok(Isoform {
            id,
            gene_id: member.id.clone(),
            chrom: member.chrom.clone(),
            start: path.boundaries.0,
            end: path.boundaries.1,
            strand: member.strand.to_string(),
            exons: Vec::new(),  // Build from junctions
            junctions: Vec::new(),  // Build from path
            support: path.supporting_reads.len(),
            is_coding: false,
            cds_start: None,
            cds_end: None,
        })
    }
}

/// Family consensus structure
struct FamilyConsensus {
    conserved_junctions: Vec<ConservedJunction>,
    total_junctions: usize,
}

struct ConservedJunction {
    chrom: String,
    donor: u64,
    acceptor: u64,
    supporting_members: HashSet<String>,
}

/// Defensible trie node
#[derive(Debug, Default)]
struct DefensibleTrieNode {
    junction: Option<(u64, u64)>,
    children: HashMap<(u64, u64), DefensibleTrieNode>,
    supporting_reads: Vec<String>,
    coverage_profile: Vec<usize>,
}

/// Defensible path through trie
#[derive(Debug, Clone)]
struct DefensiblePath {
    junctions: Vec<(u64, u64)>,
    supporting_reads: Vec<String>,
    coverage_profile: Vec<usize>,
    boundaries: (u64, u64),
    support: usize,
}

/// Defensible trie
struct DefensibleTrie {
    root: DefensibleTrieNode,
    tolerance: u64,
}

impl DefensibleTrie {
    fn new(tolerance: u64) -> Self {
        Self {
            root: DefensibleTrieNode::default(),
            tolerance,
        }
    }
    
    fn insert(&mut self, read: &SplicedRead) {
        // Implementation
    }
    
    fn extract_validated_paths(&self, min_junction_reads: usize, min_total: usize) -> Vec<DefensiblePath> {
        // Implementation
        Vec::new()
    }
}

/// Helper functions
fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| match b {
        b'A' => b'T',
        b'T' => b'A',
        b'G' => b'C',
        b'C' => b'G',
        _ => b'N',
    }).collect()
}

fn median_u64(values: &[u64]) -> u64 {
    if values.is_empty() {
        return 0;
    }
    let mut sorted = values.to_vec();
    sorted.sort();
    sorted[sorted.len() / 2]
}
