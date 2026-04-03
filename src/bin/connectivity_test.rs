//! Standalone test for connectivity-based seed selection
#[path = "../smart_seed_selector.rs"]
mod smart_seed_selector;

use smart_seed_selector::{
    analyze_connectivity, ConnectivityAnalysis, ConnectivityParams, GeneRegion,
};
use smart_seed_selector::{get_selected_seeds, print_analysis_summary};

fn main() {
    let bam_path = "/scratch/jxi21/segments/A119b.bam";
    
    // Test families
    let families = vec![
        ("ID_14", vec![
            ("NC_060941.1", 31562334, 31562446, "AC005562.2|ID_14"),
            ("NC_060941.1", 31574334, 31581958, "LRRC37BP1|ID_14"),
            ("NC_060941.1", 32953747, 32999386, "LRRC37B|ID_14"),
            ("NC_060941.1", 33034040, 33040757, "AC090616.5|ID_14"),
            ("NC_060941.1", 47154742, 47199790, "LRRC37A|ID_14"),
            ("NC_060941.1", 47373158, 47417300, "LRRC37A2|ID_14"),
            ("NC_060941.1", 65724068, 65789437, "LRRC37A3|ID_14"),
        ]),
        ("ID_146", vec![
            ("NC_060927.1", 198342015, 198343130, "AC233280.19|ID_146"),
            ("NC_060927.1", 198403551, 198444967, "AC124944.2|ID_146"),
        ]),
        ("ID_127", vec![
            ("NC_060925.1", 143000000, 143100000, "Gene1|ID_127"),
            ("NC_060925.1", 143200000, 143300000, "Gene2|ID_127"),
            ("NC_060925.1", 143400000, 143500000, "Gene3|ID_127"),
        ]),
    ];
    
    println!("=== Connectivity-Based Seed Selection Test ===\n");
    
    for (fam_id, genes) in &families {
        println!("=== {} ({}) ===", fam_id, genes.len());
        
        let gene_regions: Vec<GeneRegion> = genes.iter().map(|(c, s, e, n)| {
            GeneRegion {
                chrom: c.to_string(),
                start: *s,
                end: *e,
                name: Some(n.to_string()),
            }
        }).collect();
        
        let gene_names: Vec<String> = genes.iter().map(|(_, _, _, n)| n.clone()).collect();
        
        let params = ConnectivityParams {
            min_jaccard_floor: 0.05,
            adaptive_quantile: Some(0.15),
            min_shared_reads: 1,
            hypergeom_max_pvalue: None,
        };
        match analyze_connectivity(bam_path, &gene_regions, &params) {
            Ok(analysis) => {
                print_analysis_summary(&analysis, &gene_names);
            },
            Err(e) => {
                eprintln!("  Error: {}", e);
            }
        }
        println!();
    }
}
