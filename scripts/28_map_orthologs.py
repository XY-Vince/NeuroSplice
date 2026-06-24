import os
import csv
import urllib.request
import pandas as pd

# Configuration
PROJECT_DIR = "/Volumes/Untitled/NeuroSplice"
MAGMA_DIR = os.path.join(PROJECT_DIR, "results", "magma")
INPUT_GENES_FILE = os.path.join(PROJECT_DIR, "results", "string", "significant_genes.txt")
OUTPUT_SET_FILE = os.path.join(MAGMA_DIR, "splicing_genes_human.txt")
MGI_URL = "http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt"

def ensure_dir(d):
    if not os.path.exists(d):
        os.makedirs(d)

def download_mgi_homology(dest_path):
    if not os.path.exists(dest_path):
        print("Downloading MGI Homology Report...")
        urllib.request.urlretrieve(MGI_URL, dest_path)
    return dest_path

def load_orthologs(rpt_path):
    # Columns: DB Class Key, Common Organism Name, Symbol, ...
    # We want to map Mouse Symbol -> Human Symbol
    # Algorithm: Group by DB Class Key. If cluster has Mouse and Human, map them.
    
    print("Parsing MGI Homology...")
    df = pd.read_csv(rpt_path, sep='\t')
    
    # Filter for Mouse and Human
    df = df[df['Common Organism Name'].isin(['mouse, laboratory', 'human'])]
    
    # Create dictionary: MouseSymbol -> HumanSymbol
    # Note: One-to-many or Many-to-one possible. We take all valid human symbols in the cluster.
    
    ortho_map = {} # Mouse -> [Human, Human...]
    
    groups = df.groupby('DB Class Key')
    
    count = 0
    for key, group in groups:
        mouse_genes = group[group['Common Organism Name'] == 'mouse, laboratory']['Symbol'].unique()
        human_genes = group[group['Common Organism Name'] == 'human']['Symbol'].unique()
        
        if len(mouse_genes) > 0 and len(human_genes) > 0:
            for m in mouse_genes:
                if m not in ortho_map:
                    ortho_map[m] = set()
                ortho_map[m].update(human_genes)
                count += 1
                
    print(f"Built map for {len(ortho_map)} mouse genes.")
    return ortho_map

def main():
    ensure_dir(MAGMA_DIR)
    
    # 1. Get Homology Data
    mgi_file = os.path.join(MAGMA_DIR, "HOM_MouseHumanSequence.rpt")
    download_mgi_homology(mgi_file)
    mapping = load_orthologs(mgi_file)
    
    # 2. Load Input Genes
    if not os.path.exists(INPUT_GENES_FILE):
        print(f"Error: Input file {INPUT_GENES_FILE} not found.")
        # Fallback: recover from rMATS directly if needed, but assuming file exists (generated in Step 6)
        return

    with open(INPUT_GENES_FILE, 'r') as f:
        mouse_genes = [line.strip() for line in f if line.strip()]
    
    print(f"Input Mouse Genes: {len(mouse_genes)}")
    
    # 3. Map
    human_genes = set()
    mapped_count = 0
    missing = []
    
    for m in mouse_genes:
        # Try exact match first
        if m in mapping:
            human_genes.update(mapping[m])
            mapped_count += 1
        else:
            # Try title case or upper case as fallback (sometimes symbols vary)
            m_cap = m.capitalize()
            if m_cap in mapping:
                human_genes.update(mapping[m_cap])
                mapped_count += 1
            else:
                missing.append(m)
                
    print(f"Mapped {mapped_count}/{len(mouse_genes)} genes to {len(human_genes)} human orthologs.")
    
    # 4. Write Output (MAGMA Set Format)
    # Format: SetName Gene1 Gene2 ...
    with open(OUTPUT_SET_FILE, 'w') as f:
        # Write single line
        line = "Splicing_Dominance_Genes\t" + "\t".join(sorted(list(human_genes)))
        f.write(line + "\n")
        
    print(f"Saved MAGMA gene set to {OUTPUT_SET_FILE}")

if __name__ == "__main__":
    main()
