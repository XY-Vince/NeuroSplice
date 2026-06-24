import os
import csv
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
from scipy.stats import hypergeom

# Configuration
RESULTS_DIR = "results"
REGIONS = ["Hippocampus", "Prefrontal_Cortex", "Striatum"]
# Use consistent lowercase naming
REGION_DIRS = {
    "Hippocampus": "rmats_hippocampus_filtered",
    "Prefrontal_Cortex": "rmats_prefrontal_cortex_filtered",
    "Striatum": "rmats_striatum_filtered"
}

# Thresholds
FDR_CUTOFF = 0.05
DPSI_CUTOFF = 0.10

def get_sig_genes(region):
    rmats_dir = os.path.join(RESULTS_DIR, REGION_DIRS[region])
    file_path = os.path.join(rmats_dir, "SE.MATS.JCEC.txt")
    
    if not os.path.exists(file_path):
        print(f"Warning: Results not found for {region}")
        return set()
    
    genes = set()
    with open(file_path, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            try:
                fdr = float(row['FDR'])
                dpsi = abs(float(row['IncLevelDifference']))
                if fdr < FDR_CUTOFF and dpsi > DPSI_CUTOFF:
                    gene = row['geneSymbol'].strip('"')
                    genes.add(gene)
            except ValueError:
                continue
    
    print(f"{region}: {len(genes)} significant genes")
    # Save list
    with open(os.path.join(RESULTS_DIR, f"{region}_AS_genes.txt"), 'w') as out:
        out.write("\n".join(sorted(list(genes))))
        
    return genes

def main():
    print("--- Tier 2: Cross-Region Analysis ---")
    
    gene_sets = {}
    for region in REGIONS:
        gene_sets[region] = get_sig_genes(region)
        
    # Check if we have data for all 3 (for now, script might run partially)
    sets = [gene_sets[r] for r in REGIONS]
    labels = REGIONS
    
    if all(len(s) == 0 for s in sets):
        print("No significant genes found yet. Exiting.")
        return

    # Venn Diagram
    plt.figure(figsize=(10, 10))
    venn3(sets, set_labels=labels)
    plt.title("Overlap of AS Genes (FDR<0.05, dPSI>0.1)")
    plt.savefig(os.path.join(RESULTS_DIR, "venn_regions.png"))
    print("Venn diagram saved to venn_regions.png")
    
    # Core Signature (Intersection of all 3)
    core_genes = sets[0].intersection(sets[1]).intersection(sets[2])
    print(f"\nCore ADHD Signature ({len(core_genes)} genes):")
    print(", ".join(list(core_genes)))
    
    with open(os.path.join(RESULTS_DIR, "core_adhd_signature.txt"), 'w') as f:
        f.write("\n".join(sorted(list(core_genes))))

if __name__ == "__main__":
    main()
