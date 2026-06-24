import os
import json
import urllib.request
import urllib.parse
import sys

# Configuration
CORE_GENES = ["Dctn1", "Crem", "Pdlim7", "Spata5"]
SPECIES_ID = 10090  # Mouse
OUT_DIR = "results/string_core_network"

def ensure_dir(d):
    if not os.path.exists(d):
        os.makedirs(d)

def query_string_network(genes):
    base_url = "https://string-db.org/api"
    gene_str = "\n".join(genes)
    
    print(f"Querying STRING for Core 4: {genes}")
    
    # 1. Map IDs
    print("Mapping IDs...")
    params = {
        "identifiers": gene_str,
        "species": SPECIES_ID,
        "limit": 1,
        "echo_query": 1
    }
    data = urllib.parse.urlencode(params).encode('utf-8')
    try:
        req = urllib.request.Request(f"{base_url}/json/get_string_ids", data=data)
        with urllib.request.urlopen(req) as r:
            mappings = json.loads(r.read().decode('utf-8'))
    except Exception as e:
        print(f"Error mapping IDs: {e}")
        return

    string_ids = [m['stringId'] for m in mappings]
    print(f"Mapped to: {string_ids}")
    
    # 2. Get Interaction Partners (interaction_partners API or network with limit)
    # maximizing neighbors.
    print("Retrieving neighbors...")
    params = {
        "identifiers": "\n".join(string_ids),
        "species": SPECIES_ID,
        "add_nodes": 20, # Add top 20 neighbors
        "network_flavor": "confidence"
    }
    data = urllib.parse.urlencode(params).encode('utf-8')
    
    # Grab the network image
    try:
        req = urllib.request.Request(f"{base_url}/svg/network", data=data)
        with urllib.request.urlopen(req) as r:
            with open(os.path.join(OUT_DIR, "core4_network.svg"), 'wb') as f:
                f.write(r.read())
        print("Saved core4_network.svg")
    except Exception as e:
        print(f"Error getting network SVG: {e}")

    # 3. Get Enrichment for this cluster
    # Note: Enrichment API works on the input list. To include neighbors, we need their IDs.
    # We can get the text description of the network to find neighbors.
    
    # Let's get the interaction table to identify neighbors
    params['required_score'] = 400
    try:
        req = urllib.request.Request(f"{base_url}/tsv/network", data=data)
        with urllib.request.urlopen(req) as r:
            lines = r.read().decode('utf-8').split('\n')
            
        with open(os.path.join(OUT_DIR, "core4_interactions.tsv"), 'w') as f:
            f.write("\n".join(lines))
            
        # Extract all unique proteins from the interaction network
        all_proteins = set(string_ids)
        for line in lines[1:]:
            parts = line.split('\t')
            if len(parts) > 1:
                all_proteins.add(parts[0])
                all_proteins.add(parts[1])
        
        print(f"Expanded network has {len(all_proteins)} proteins.")
        
        # Now run enrichment on THIS expanded list
        print("Running enrichment on extended network...")
        enrich_params = {
            "identifiers": "\n".join(list(all_proteins)),
            "species": SPECIES_ID
        }
        data = urllib.parse.urlencode(enrich_params).encode('utf-8')
        
        req = urllib.request.Request(f"{base_url}/json/enrichment", data=data)
        with urllib.request.urlopen(req) as r:
            enrichment = json.loads(r.read().decode('utf-8'))
            with open(os.path.join(OUT_DIR, "core4_enrichment.json"), 'w') as f:
                json.dump(enrichment, f, indent=2)
            
            # Print top terms
            print("\nTop Enriched Pathways (FDR < 0.05):")
            count = 0
            for term in enrichment:
                if term['fdr'] < 0.05 and count < 10:
                    print(f"- {term['description']} (FDR: {term['fdr']:.2e})")
                    count += 1
                    
    except Exception as e:
        print(f"Error analyzing network: {e}")

if __name__ == "__main__":
    ensure_dir(OUT_DIR)
    query_string_network(CORE_GENES)
