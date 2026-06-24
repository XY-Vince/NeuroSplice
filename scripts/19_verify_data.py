import os
import csv

REGIONS = ["Hippocampus", "Prefrontal_Cortex", "Striatum"]
CORE_GENES = ["Cep63", "Pdlim7", "Carmil1", "Crem", "Spata5", "Dctn1", "Oprl1"]
FILES = {
    "Hippocampus": "results/rmats_hippocampus/SE.MATS.JCEC.txt",
    "Prefrontal_Cortex": "results/rmats_prefrontal_cortex/SE.MATS.JCEC.txt",
    "Striatum": "results/rmats_striatum/SE.MATS.JCEC.txt"
}

def check_directionality():
    print("--- Core Gene Directionality Check ---")
    data = {gene: {} for gene in CORE_GENES}
    
    for region, filepath in FILES.items():
        if not os.path.exists(filepath):
            print(f"File not found: {filepath}")
            continue
            
        with open(filepath, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                gene = row['geneSymbol'].replace('"', '')
                if gene in CORE_GENES:
                    dpsi = float(row['IncLevelDifference'])
                    fdr = float(row['FDR'])
                    if fdr < 0.05: # Only significant ones
                        data[gene][region] = dpsi

    # Print Report
    print(f"{'Gene':<10} {'Hippocampus':<15} {'PFC':<15} {'Striatum':<15} {'Consistent?'}")
    for gene in CORE_GENES:
        dirs = []
        row_str = f"{gene:<10} "
        for region in REGIONS:
            val = data[gene].get(region, "N/A")
            if val != "N/A":
                dirs.append(1 if val > 0 else -1)
                val_str = f"{val:.4f}"
            else:
                val_str = "-"
            row_str += f"{val_str:<15} "
        
        consistent = "YES" if (len(dirs) > 0 and abs(sum(dirs)) == len(dirs)) else "NO"
        # Logic: if all are 1, sum = len. if all -1, abs(sum) = len. Mixed will be less.
        
        print(row_str + consistent)

def count_total_genes():
    print("\n--- 983 Genes Verification ---")
    all_genes = set()
    rmats_dir = "results/rmats"
    for fname in os.listdir(rmats_dir):
        if fname.endswith("MATS.JCEC.txt"):
            path = os.path.join(rmats_dir, fname)
            with open(path, 'r') as f:
                reader = csv.DictReader(f, delimiter='\t')
                for row in reader:
                    try:
                        if float(row['FDR']) < 0.05:
                            all_genes.add(row['geneSymbol'].replace('"', ''))
                    except: pass
    
    print(f"Total Unique Significant Genes (FDR < 0.05) across all event types: {len(all_genes)}")

if __name__ == "__main__":
    check_directionality()
    count_total_genes()
