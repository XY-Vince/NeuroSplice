import json
import os

try:
    with open('results/string/enrichment.json') as f:
        data = json.load(f)

    # STRING results are usually a list of dicts
    # keys: term, description, category, fdr, number_of_genes, etc.

    categories = ['Process', 'RCTM', 'WikiPathways']
    for cat in categories:
        print(f"\nTOP {cat}:")
        # Filter by category
        items = [x for x in data if x['category'] == cat]
        # Sort by FDR just in case (though usually sorted)
        items.sort(key=lambda x: x['fdr'])
        
        for item in items[:5]:
            print(f"- {item['description']} (FDR: {item['fdr']:.2e}, Genes: {item['number_of_genes']})")

except Exception as e:
    print(f"Error parsing: {e}")
