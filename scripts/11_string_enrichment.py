import os
import csv
import urllib.request
import urllib.parse
import json
import time

# Configuration
PROJECT_DIR = "/Volumes/Untitled/NeuroSplice"
RMATS_DIR = os.path.join(PROJECT_DIR, "results", "rmats")
OUT_DIR = os.path.join(PROJECT_DIR, "results", "string")

FILES_TO_PROCESS = [
    "SE.MATS.JCEC.txt",
    "MXE.MATS.JCEC.txt",
    "A3SS.MATS.JCEC.txt",
    "A5SS.MATS.JCEC.txt",
    "RI.MATS.JCEC.txt"
]

FDR_THRESHOLD = 0.05
SPECIES_ID = 10090  # Mouse

def ensure_dir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

def get_significant_genes(rmats_dir, files, fdr_cutoff):
    sig_genes = set()
    print(f"Scanning for significant events (FDR < {fdr_cutoff})...")
    
    for filename in files:
        filepath = os.path.join(rmats_dir, filename)
        if not os.path.exists(filepath):
            print(f"Warning: File not found: {filepath}")
            continue
            
        print(f"Processing {filename}...")
        with open(filepath, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            # Check headers just in case
            if 'FDR' not in reader.fieldnames or 'geneSymbol' not in reader.fieldnames:
                 print(f"Error: Missing columns in {filename}. Available: {reader.fieldnames}")
                 continue
                 
            count = 0
            for row in reader:
                try:
                    fdr = float(row['FDR'])
                    if fdr < fdr_cutoff:
                        # geneSymbol is usually column 2 (0-indexed) or named 'geneSymbol'
                        gene = row['geneSymbol']
                        # Handle cases where gene might be quoted or have weird chars
                        gene = gene.strip('"')
                        if gene:
                            sig_genes.add(gene)
                            count += 1
                except ValueError:
                    continue # Skip header or bad parsing
            print(f"  Found {count} significant events.")
            
    return list(sig_genes)

def query_string_api(genes, species_id, out_dir):
    base_url = "https://string-db.org/api"
    gene_list_str = "\n".join(genes)
    
    print(f"Querying STRING for {len(genes)} genes...")
    
    # 1. Map Identifiers (Best Practice)
    map_url = f"{base_url}/json/get_string_ids"
    params = {
        "identifiers": gene_list_str,
        "species": species_id,
        "limit": 1,
        "echo_query": 1
    }
    
    # Use standard form data encoding
    data = urllib.parse.urlencode(params).encode('utf-8')
    try:
        req = urllib.request.Request(map_url, data=data)
        with urllib.request.urlopen(req) as response:
            mappings = json.loads(response.read().decode('utf-8'))
    except Exception as e:
        print(f"Error mapping IDs: {e}")
        return

    mapped_string_ids = list(set([item['stringId'] for item in mappings]))
    print(f"Successfully mapped {len(mapped_string_ids)} genes to STRING IDs.")
    
    if not mapped_string_ids:
        print("No genes mapped. Exiting.")
        return

    # 2. Get Network Image (SVG)
    network_url = f"{base_url}/svg/network"
    params = {
        "identifiers": "\n".join(mapped_string_ids),
        "species": species_id,
        "add_white_nodes": 0,
        "network_flavor": "confidence" # confidence links
    }
    data = urllib.parse.urlencode(params).encode('utf-8')
    
    print("Downloading network image...")
    try:
        req = urllib.request.Request(network_url, data=data)
        with urllib.request.urlopen(req) as response:
            svg_content = response.read()
            with open(os.path.join(out_dir, "network.svg"), 'wb') as f:
                f.write(svg_content)
        print("Saved network.svg")
    except Exception as e:
        print(f"Error downloading network image: {e}")

    # 3. Get Enrichment (JSON)
    enrichment_url = f"{base_url}/json/enrichment"
    params = {
        "identifiers": "\n".join(mapped_string_ids),
        "species": species_id
    }
    data = urllib.parse.urlencode(params).encode('utf-8')
    
    print("Retrieving enrichment data...")
    try:
        req = urllib.request.Request(enrichment_url, data=data)
        with urllib.request.urlopen(req) as response:
            enrichment_data = json.loads(response.read().decode('utf-8'))
            with open(os.path.join(out_dir, "enrichment.json"), 'w') as f:
                json.dump(enrichment_data, f, indent=2)
        print("Saved enrichment.json")
    except Exception as e:
        print(f"Error retrieving enrichment: {e}")

    # 4. Generate Interaction Table (TSV)
    interaction_url = f"{base_url}/tsv/network"
    params = {
        "identifiers": "\n".join(mapped_string_ids),
        "species": species_id,
        "required_score": 400 
    }
    data = urllib.parse.urlencode(params).encode('utf-8')
    
    print("Retrieving interaction table...")
    try:
        req = urllib.request.Request(interaction_url, data=data)
        with urllib.request.urlopen(req) as response:
            tsv_content = response.read()
            with open(os.path.join(out_dir, "interactions.tsv"), 'wb') as f:
                f.write(tsv_content)
        print("Saved interactions.tsv")
    except Exception as e:
        print(f"Error retrieving interactions: {e}")

def main():
    ensure_dir(OUT_DIR)
    
    significant_genes = get_significant_genes(RMATS_DIR, FILES_TO_PROCESS, FDR_THRESHOLD)
    print(f"Total unique significant genes: {len(significant_genes)}")
    
    if len(significant_genes) == 0:
        print("No significant genes found. Exiting.")
        return

    # Save gene list
    with open(os.path.join(OUT_DIR, "significant_genes.txt"), 'w') as f:
        f.write("\n".join(significant_genes))
        
    # Limit query size if too large (STRING has limits, usually 2000 is fine)
    if len(significant_genes) > 2000:
        print("Warning: Too many genes for one query. Taking top 2000.")
        significant_genes = significant_genes[:2000]

    query_string_api(significant_genes, SPECIES_ID, OUT_DIR)
    print("Analysis Complete.")

if __name__ == "__main__":
    main()
