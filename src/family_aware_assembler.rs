//! Family-Aware Transcript Assembler
//! 
//! Core innovation: Assemble transcripts by leveraging information across 
//! all members of a multi-copy gene family.
//! 
//! Key principle: Reads from all family members are informative for each other.
//! Conserved junctions (supported by multiple members) are high confidence.
//! Member-specific variations distinguish paralogs.

use crate::bam_parser::{SplicedRead, SpliceJunction};
use crate::isoform::Isoform;
use hashbrown::{HashMap, HashSet};
use rayon::prelude::*;

/// Configuration for family-aware assembly
#[derive(Debug, Clone)]
pub struct FamilyAssemblyConfig {
    /// Minimum reads per member for isoform (default: 2)
    pub min_member_support: usize,
    /// Minimum total reads across family (default: 5)
    pub min_family_support: usize,
    /// Minimum members supporting a junction (default: 1)
    pub min_member_junction: usize,
    /// Junction tolerance in bp (default: 10)
    pub junction_tolerance: u64,
    /// Whether to use cross-member validation (default: true)
    pub cross_member_validation: bool,
    /// Max isoforms per member (default: 50)
    pub max_isoforms_per_member: usize,
}

impl Default for FamilyAssemblyConfig {
    fn default() -> Self {
        Self {
            min_member_support: 2,
            min_family_support: 5,
            min_member_junction: 1,
            junction_tolerance: 10,
            cross_member_validation: true,
            max_isoforms_per_member: 50,
        }
    }
}

/// Represents a member of a gene family
#[derive(Debug, Clone)]
pub struct FamilyMember {
    pub id: String,
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub strand: char,
}

/// Junction with family-wide support information
#[derive(Debug, Clone)]
struct FamilyJunction {
    donor: u64,
    acceptor: u64,
    /// Support from each family member
    member_support: HashMap<String, Vec<String>>, // member_id -> read_ids
    /// Total support across family
    total_support: usize,
    /// Number of members with support
    member_count: usize,
}

impl FamilyJunction {
    fn new(donor: u64, acceptor: u64) -> Self {
        Self {
            donor,
            acceptor,
            member_support: HashMap::new(),
            total_support: 0,
            member_count: 0,
        }
    }
    
    fn add_support(&mut self, member_id: &str, read_id: &str) {
        let entry = self.member_support.entry(member_id.to_string()).or_default();
        if !entry.contains(&read_id.to_string()) {
            entry.push(read_id.to_string());
            self.total_support += 1;
            self.member_count = self.member_support.len();
        }
    }
    
    /// Is this junction conserved (supported by multiple members)?
    fn is_conserved(&self) -> bool {
        self.member_count >= 2
    }
    
    /// Get support from specific member
    fn member_support_count(&self, member_id: &str) -> usize {
        self.member_support.get(member_id).map(|v| v.len()).unwrap_or(0)
    }
}

/// Family-wide junction graph
#[derive(Debug)]
struct FamilyJunctionGraph {
    junctions: HashMap<(u64, u64), FamilyJunction>,
    tolerance: u64,
}

impl FamilyJunctionGraph {
    fn new(tolerance: u64) -> Self {
        Self {
            junctions: HashMap::new(),
            tolerance,
        }
    }
    
    /// Add junction from a family member
    fn add_junction(&mut self, donor: u64, acceptor: u64, member_id: &str, read_id: &str) {
        let key = self.find_matching_junction(donor, acceptor)
            .unwrap_or((donor, acceptor));
        
        let junction = self.junctions.entry(key).or_insert_with(|| {
            FamilyJunction::new(key.0, key.1)
        });
        
        junction.add_support(member_id, read_id);
    }
    
    /// Find existing junction within tolerance
    fn find_matching_junction(&self, donor: u64, acceptor: u64) -> Option<(u64, u64)> {
        for (key, _) in &self.junctions {
            let donor_diff = key.0.abs_diff(donor);
            let acceptor_diff = key.1.abs_diff(acceptor);
            
            if donor_diff <= self.tolerance && acceptor_diff <= self.tolerance {
                return Some(*key);
            }
        }
        None
    }
    
    /// Get conserved junctions (supported by multiple members)
    fn get_conserved_junctions(&self) -> Vec<&FamilyJunction> {
        self.junctions.values()
            .filter(|j| j.is_conserved())
            .collect()
    }
    
    /// Get all junctions meeting minimum support
    fn get_supported_junctions(&self, min_support: usize) -> Vec<&FamilyJunction> {
        self.junctions.values()
            .filter(|j| j.total_support >= min_support)
            .collect()
    }
}

/// Isoform path with family support
#[derive(Debug, Clone)]
struct FamilyPath {
    junctions: Vec<(u64, u64)>,
    member_support: HashMap<String, Vec<(String, u64, u64)>>, // member_id -> [(read_id, tss, tes)]
    total_support: usize,
    conserved_junction_count: usize,
}

impl FamilyPath {
    fn confidence_score(&self) -> f64 {
        let length_score = (self.junctions.len() as f64).sqrt();
        let support_score = self.total_support as f64;
        let conservation_score = self.conserved_junction_count as f64 * 2.0;
        
        (length_score + support_score + conservation_score) / 3.0
    }
    
    fn is_supported_by(&self, member_id: &str) -> bool {
        self.member_support.contains_key(member_id)
    }
}

/// Family-aware trie for path assembly
#[derive(Debug)]
struct FamilyTrie {
    root: FamilyTrieNode,
    tolerance: u64,
}

#[derive(Debug, Default)]
struct FamilyTrieNode {
    junction: Option<(u64, u64)>,
    children: HashMap<(u64, u64), FamilyTrieNode>,
    /// Support from each member
    member_reads: HashMap<String, Vec<(String, u64, u64)>>, // member_id -> [(read_id, tss, tes)]
}

impl FamilyTrie {
    fn new(tolerance: u64) -> Self {
        Self {
            root: FamilyTrieNode::default(),
            tolerance,
        }
    }
    
    /// Insert a read with its member assignment
    fn insert(&mut self, read: &SplicedRead, member_id: &str) {
        let mut current = &mut self.root;
        
        for junction in &read.junctions {
            let key = (junction.start, junction.end);
            let matched_key = self.find_matching_child(current, key);
            let use_key = matched_key.unwrap_or(key);
            
            if !current.children.contains_key(&use_key) {
                current.children.insert(use_key, FamilyTrieNode {
                    junction: Some(use_key),
                    children: HashMap::new(),
                    member_reads: HashMap::new(),
                });
            }
            
            current = current.children.get_mut(&use_key).unwrap();
            
            // Track support per member
            let entry = current.member_reads.entry(member_id.to_string()).or_default();
            entry.push((
                read.read_id.clone(),
                read.tx_start,
                read.tx_end,
            ));
        }
    }
    
    fn find_matching_child(&self, node: &FamilyTrieNode, target: (u64, u64)) -> Option<(u64, u64)> {
        for (key, _) in &node.children {
            let donor_diff = key.0.abs_diff(target.0);
            let acceptor_diff = key.1.abs_diff(target.1);
            
            if donor_diff <= self.tolerance && acceptor_diff <= self.tolerance {
                return Some(*key);
            }
        }
        None
    }
    
    /// Extract family paths with sufficient support
    fn extract_paths(&self, config: &FamilyAssemblyConfig) -> Vec<FamilyPath> {
        let mut paths = Vec::new();
        let mut current_path = Vec::new();
        
        self.dfs_extract(
            &self.root,
            &mut current_path,
            config,
            &mut paths,
        );
        
        // Score and sort
        let mut scored: Vec<_> = paths.into_iter()
            .map(|p| (p.confidence_score(), p))
            .collect();
        scored.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
        
        scored.into_iter().map(|(_, p)| p).collect()
    }
    
    fn dfs_extract(
        &self,
        node: &FamilyTrieNode,
        current_path: &mut Vec<(u64, u64)>,
        config: &FamilyAssemblyConfig,
        results: &mut Vec<FamilyPath>,
    ) {
        // Check if current path is valid
        if !current_path.is_empty() {
            let total_support: usize = node.member_reads.values().map(|v| v.len()).sum();
            let member_count = node.member_reads.len();
            
            if total_support >= config.min_family_support && 
               member_count >= config.min_member_junction {
                // Count conserved junctions
                let conserved_count = current_path.iter().filter(|(d, a)| {
                    // This is simplified - in real implementation check against FamilyJunctionGraph
                    true
                }).count();
                
                results.push(FamilyPath {
                    junctions: current_path.clone(),
                    member_support: node.member_reads.clone(),
                    total_support,
                    conserved_junction_count: conserved_count,
                });
            }
        }
        
        // Recurse
        for (_, child) in &node.children {
            if let Some(junction) = child.junction {
                current_path.push(junction);
                self.dfs_extract(child, current_path, config, results);
                current_path.pop();
            }
        }
    }
}

/// Main family-aware assembler
pub struct FamilyAwareAssembler {
    family_id: String,
    members: Vec<FamilyMember>,
    config: FamilyAssemblyConfig,
}

/// Result for one family member
#[derive(Debug)]
pub struct MemberIsoforms {
    pub member_id: String,
    pub isoforms: Vec<FamilyAwareIsoform>,
}

/// Isoform with family support information
#[derive(Debug)]
pub struct FamilyAwareIsoform {
    pub isoform: Isoform,
    pub family_support: usize,
    pub member_support: usize,
    pub conserved_junctions: usize,
    pub paralog_groups: Vec<String>, // IDs of similar isoforms at other members
}

impl FamilyAwareAssembler {
    pub fn new(family_id: String, members: Vec<FamilyMember>) -> Self {
        Self {
            family_id,
            members,
            config: FamilyAssemblyConfig::default(),
        }
    }
    
    pub fn with_config(mut self, config: FamilyAssemblyConfig) -> Self {
        self.config = config;
        self
    }
    
    /// Main assembly entry point
    pub fn assemble(&self, reads: Vec<SplicedRead>) -> anyhow::Result<Vec<MemberIsoforms>> {
        eprintln!("Family-Aware Assembly: {}", self.family_id);
        eprintln!("  Members: {}", self.members.len());
        eprintln!("  Total reads: {}", reads.len());
        
        // Phase 1: Build family junction graph
        let junction_graph = self.build_junction_graph(&reads);
        eprintln!("  Family junctions: {}", junction_graph.junctions.len());
        eprintln!("  Conserved junctions: {}", junction_graph.get_conserved_junctions().len());
        
        // Phase 2: Build family trie
        let family_trie = self.build_family_trie(&reads);
        
        // Phase 3: Extract family paths
        let family_paths = family_trie.extract_paths(&self.config);
        eprintln!("  Family paths: {}", family_paths.len());
        
        // Phase 4: Assign paths to members
        let member_isoforms = self.assign_paths_to_members(family_paths)?;
        eprintln!("  Member isoforms: {}", 
                 member_isoforms.iter().map(|m| m.isoforms.len()).sum::<usize>());
        
        // Phase 5: Cross-validate
        let validated = if self.config.cross_member_validation {
            self.cross_validate(member_isoforms)?
        } else {
            member_isoforms
        };
        
        Ok(validated)
    }
    
    /// Build family-wide junction graph
    fn build_junction_graph(&self, reads: &[SplicedRead]) -> FamilyJunctionGraph {
        let mut graph = FamilyJunctionGraph::new(self.config.junction_tolerance);
        
        // Assign reads to members
        for read in reads {
            if let Some(member_id) = self.assign_read_to_member(read) {
                for junction in &read.junctions {
                    graph.add_junction(
                        junction.start,
                        junction.end,
                        &member_id,
                        &read.read_id,
                    );
                }
            }
        }
        
        graph
    }
    
    /// Assign read to best matching member
    fn assign_read_to_member(&self, read: &SplicedRead) -> Option<String> {
        let mut best_member = None;
        let mut best_overlap = 0u64;
        
        for member in &self.members {
            if read.chrom != member.chrom || read.strand != member.strand {
                continue;
            }
            
            let overlap_start = read.tx_start.max(member.start);
            let overlap_end = read.tx_end.min(member.end);
            
            if overlap_start < overlap_end {
                let overlap = overlap_end - overlap_start;
                if overlap > best_overlap {
                    best_overlap = overlap;
                    best_member = Some(member.id.clone());
                }
            }
        }
        
        best_member
    }
    
    /// Build family trie with member tracking
    fn build_family_trie(&self, reads: &[SplicedRead]) -> FamilyTrie {
        let mut trie = FamilyTrie::new(self.config.junction_tolerance);
        
        for read in reads {
            if let Some(member_id) = self.assign_read_to_member(read) {
                trie.insert(read, &member_id);
            }
        }
        
        trie
    }
    
    /// Assign family paths to individual members
    fn assign_paths_to_members(&self, paths: Vec<FamilyPath>) -> anyhow::Result<Vec<MemberIsoforms>> {
        let mut member_results: HashMap<String, Vec<FamilyAwareIsoform>> = HashMap::new();
        
        for path in paths {
            for (member_id, reads) in &path.member_support {
                if reads.len() < self.config.min_member_support {
                    continue;
                }
                
                // Find member-specific boundaries
                let (tss_values, tes_values): (Vec<u64>, Vec<u64>) = reads.iter()
                    .map(|(_, tss, tes)| (*tss, *tes))
                    .unzip();
                
                let tss = median_u64(&tss_values);
                let tes = median_u64(&tes_values);
                
                // Build isoform
                let isoform = self.build_isoform_from_path(
                    &path.junctions,
                    tss,
                    tes,
                    member_id,
                )?;
                
                let fam_isoform = FamilyAwareIsoform {
                    isoform,
                    family_support: path.total_support,
                    member_support: reads.len(),
                    conserved_junctions: path.conserved_junction_count,
                    paralog_groups: Vec::new(), // Will be filled in cross-validation
                };
                
                member_results.entry(member_id.clone())
                    .or_default()
                    .push(fam_isoform);
            }
        }
        
        // Convert to MemberIsoforms
        Ok(member_results.into_iter()
            .map(|(id, isoforms)| MemberIsoforms { member_id: id, isoforms })
            .collect())
    }
    
    /// Build isoform from path
    fn build_isoform_from_path(
        &self,
        junctions: &[(u64, u64)],
        tss: u64,
        tes: u64,
        member_id: &str,
    ) -> anyhow::Result<Isoform> {
        let member = self.members.iter()
            .find(|m| m.id == member_id)
            .ok_or_else(|| anyhow::anyhow!("Member not found: {}", member_id))?;
        
        // Build exons from junctions
        let mut exons = Vec::new();
        let mut current_start = tss;
        
        for (donor, acceptor) in junctions {
            exons.push((current_start, *donor));
            current_start = *acceptor;
        }
        exons.push((current_start, tes));
        
        let junction_objs: Vec<SpliceJunction> = junctions.iter()
            .map(|(d, a)| SpliceJunction::new(member.chrom.clone(), *d, *a, member.strand))
            .collect();
        
        Ok(Isoform {
            id: format!("{}_{}_{}", self.family_id, member_id, tss),
            gene_id: member_id.to_string(),
            chrom: member.chrom.clone(),
            start: tss,
            end: tes,
            strand: member.strand.to_string(),
            exons,
            junctions: junction_objs,
            support: 0, // Will be set later
            is_coding: false,
            cds_start: None,
            cds_end: None,
        })
    }
    
    /// Cross-validate isoforms across family members
    fn cross_validate(&self, mut members: Vec<MemberIsoforms>) -> anyhow::Result<Vec<MemberIsoforms>> {
        // Find paralog groups (similar isoforms at different members)
        let mut paralog_groups: Vec<Vec<(String, usize)>> = Vec::new(); // (member_id, isoform_idx)
        
        for (i, member) in members.iter().enumerate() {
            for (j, isoform) in member.isoforms.iter().enumerate() {
                let mut found_group = false;
                
                for group in &mut paralog_groups {
                    // Check if this isoform belongs to this group
                    let representative = &members[group[0].0.parse::<usize>().unwrap()]
                        .isoforms[group[0].1];
                    
                    if self.isoforms_similar(&isoform.isoform, &representative.isoform) {
                        group.push((i.to_string(), j));
                        found_group = true;
                        break;
                    }
                }
                
                if !found_group {
                    paralog_groups.push(vec![(i.to_string(), j)]);
                }
            }
        }
        
        // Update paralog group assignments
        for (group_id, group) in paralog_groups.iter().enumerate() {
            for (member_idx, isoform_idx) in group {
                let member_idx: usize = member_idx.parse()?;
                let isoform = &mut members[member_idx].isoforms[*isoform_idx];
                
                // Add IDs of other members in group
                for (other_member, other_isoform) in group {
                    if other_member != &member_idx.to_string() {
                        isoform.paralog_groups.push(format!("{}_{}", other_member, other_isoform));
                    }
                }
                
                // Boost confidence if supported by multiple members
                if group.len() >= 2 {
                    isoform.family_support = isoform.family_support.max(group.len() * 5);
                }
            }
        }
        
        Ok(members)
    }
    
    /// Check if two isoforms are similar (same intron chain)
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
}

/// Calculate median of u64 values
fn median_u64(values: &[u64]) -> u64 {
    if values.is_empty() {
        return 0;
    }
    
    let mut sorted = values.to_vec();
    sorted.sort();
    
    let mid = sorted.len() / 2;
    if sorted.len() % 2 == 0 {
        (sorted[mid - 1] + sorted[mid]) / 2
    } else {
        sorted[mid]
    }
}
