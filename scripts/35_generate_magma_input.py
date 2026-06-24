import pandas as pd
import requests
import io
import os

# 1. Load Mouse Genes
df = pd.read_csv("results/final_gene_list.csv")
mouse_genes = df['GeneSymbol'].unique().tolist()
print(f"Loaded {len(mouse_genes)} mouse genes.")

# 2. Query MyGene.info for Orthologs
print("Querying MyGene.info for human orthologs...")
url = "https://mygene.info/v3/query"
headers = {'content-type': 'application/x-www-form-urlencoded'}
params = {
    'q': ",".join(mouse_genes),
    'scopes': 'symbol',
    'fields': 'homologene.genes',
    'species': 'mouse',
    'size': 1000
}

# MyGene via POST for large lists
res = requests.post(url, data={'q': ",".join(mouse_genes), 'scopes': 'symbol', 'fields': 'homologene.genes', 'species': 'mouse'}, headers=headers)
data = res.json()

human_entrez_ids = []

for item in data:
    if 'homologene' in item:
        # homologene can be a dict or list
        h = item['homologene']
        if isinstance(h, dict):
            h = [h]
        
        for entry in h:
            # Check for human (taxid 9606)
            if 'genes' in entry:
                genes = entry['genes']
                for gene in genes:
                    if gene[0] == 9606: # Human TaxID
                        human_entrez_ids.append(str(gene[1])) # Entrez ID

human_entrez_ids = sorted(list(set(human_entrez_ids)))
print(f"Mapped to {len(human_entrez_ids)} unique Human Entrez IDs.")

# 3. Write to File
out_file = "results/magma/magma_input_genes.txt"
with open(out_file, "w") as f:
    for gid in human_entrez_ids:
        f.write(f"{gid}\n")

print(f"Saved to {out_file}")
