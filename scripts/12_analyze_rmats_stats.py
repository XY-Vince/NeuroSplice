import csv
import os
import math

# Configuration
RMATS_DIR = "results/rmats"
FILES = {
    "SE": "SE.MATS.JCEC.txt",
    "MXE": "MXE.MATS.JCEC.txt",
    "A3SS": "A3SS.MATS.JCEC.txt",
    "A5SS": "A5SS.MATS.JCEC.txt",
    "RI": "RI.MATS.JCEC.txt"
}

FDR_CUTOFF = 0.05

def analyze_file(event_type, filename):
    filepath = os.path.join(RMATS_DIR, filename)
    if not os.path.exists(filepath):
        print(f"Skipping {event_type} (File not found)")
        return None
    
    sig_events = []
    unique_genes = set()
    total_events = 0
    
    with open(filepath, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            total_events += 1
            try:
                fdr = float(row['FDR'])
                dpsi = abs(float(row['IncLevelDifference']))
                
                if fdr < FDR_CUTOFF:
                    gene = row['geneSymbol'].strip('"')
                    sig_events.append(dpsi)
                    unique_genes.add(gene)
            except ValueError:
                continue

    if not sig_events:
        return None
        
    # Calculate Stats
    count = len(sig_events)
    mean_dpsi = sum(sig_events) / count
    max_dpsi = max(sig_events)
    
    # Binning
    bins = {"<0.05": 0, "0.05-0.15": 0, "0.15-0.3": 0, ">0.3": 0}
    for d in sig_events:
        if d < 0.05: bins["<0.05"] += 1
        elif d < 0.15: bins["0.05-0.15"] += 1
        elif d < 0.3: bins["0.15-0.3"] += 1
        else: bins[">0.3"] += 1
        
    return {
        "Type": event_type,
        "Total_Events": total_events,
        "Sig_Events": count,
        "Unique_Genes": len(unique_genes),
        "dPSI_Mean": mean_dpsi,
        "dPSI_Max": max_dpsi,
        "dPSI_Dist": bins,
        "Genes_List": list(unique_genes)
    }

all_stats = []
all_sig_genes = []
total_sig_events = 0

print("--- Analysis Report ---")
print(f"FDR Threshold: {FDR_CUTOFF}")
print(f"Method: JCEC (Junction Count + Exon Count)\n")

for etype, fname in FILES.items():
    res = analyze_file(etype, fname)
    if res:
        all_stats.append(res)
        all_sig_genes.extend(res['Genes_List'])
        total_sig_events += res["Sig_Events"]
        
        print(f"[{etype}]")
        print(f"  Significant Events: {res['Sig_Events']}")
        print(f"  Unique Genes: {res['Unique_Genes']}")
        print(f"  Avg |dPSI|: {res['dPSI_Mean']:.3f}")
        print(f"  Max |dPSI|: {res['dPSI_Max']:.3f}")
        print(f"  dPSI Distribution: {res['dPSI_Dist']}")
        print("")

# Aggregate Stats
unique_genes_total = len(set(all_sig_genes))
avg_events_per_gene = total_sig_events / unique_genes_total if unique_genes_total > 0 else 0

print("--- Aggregated Stats ---")
print(f"Total Significant Events: {total_sig_events}")
print(f"Total Unique Genes: {unique_genes_total}")
print(f"Average Events per Gene: {avg_events_per_gene:.2f}")
