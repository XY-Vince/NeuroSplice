import os
import csv
import math
import glob
import numpy as np

# Configuration
PROJECT_DIR = "/Volumes/Untitled/NeuroSplice"
RMATS_DIR = os.path.join(PROJECT_DIR, "results/rmats_striatum_filtered") # Using Striatum as proxy/strongest region, or check all.
# Actually we should check all regions for consistency, but the User asked for "results"
# Let's check all 3 regions for the Core 4.

REGIONS = {
    "Hippocampus": os.path.join(PROJECT_DIR, "results/rmats_hippocampus_filtered"),
    "Prefrontal_Cortex": os.path.join(PROJECT_DIR, "results/rmats_prefrontal_cortex_filtered"),
    "Striatum": os.path.join(PROJECT_DIR, "results/rmats_striatum_filtered")
}

SALMON_DIR = os.path.join(PROJECT_DIR, "results/salmon_quant")
GTF_FILE = os.path.join(PROJECT_DIR, "reference/gencode.vM25.annotation.gtf")

CORE_GENES = ["Dctn1", "Crem", "Pdlim7", "Spata5"]
SPLICING_FACTORS = [
    "Rbofx1", "Rbfox2", "Rbfox3", 
    "Ptbp1", "Ptbp2",
    "Nova1", "Nova2",
    "Srsf1", "Srsf2", "Srsf3", "Srsf4", "Srsf5", "Srsf6", "Srsf7", "Srsf8", "Srsf9", "Srsf10", "Srsf11", "Srsf12"
]
# Note: Rbfox1 spelling -> Rbfox1

def load_t2g():
    print("Loading T2G...")
    t2g = {}
    g2name = {}
    with open(GTF_FILE, 'r') as f:
        for line in f:
            if line.startswith("#"): continue
            if "\ttranscript\t" not in line: continue
            parts = line.strip().split('\t')
            attr = parts[8]
            
            tid = ""
            gid = ""
            gname = ""
            if 'transcript_id "' in attr: tid = attr.split('transcript_id "')[1].split('"')[0]
            if 'gene_id "' in attr: gid = attr.split('gene_id "')[1].split('"')[0]
            if 'gene_name "' in attr: gname = attr.split('gene_name "')[1].split('"')[0]
            
            if tid and gid:
                t2g[tid] = gid
                if gname: g2name[gid] = gname
    return t2g, g2name

def apply_bh_correction(p_values):
    """Benjamini-Hochberg correction for FDR control."""
    n = len(p_values)
    # Sort p-values (keep original indices)
    sorted_indices = np.argsort(p_values)
    sorted_p = np.array(p_values)[sorted_indices]
    
    # Calculate BH critical values: (rank / n) * alpha? 
    # Or just adjust p-values: p_adj = p * n / rank
    # Monotonicity enforcement: p_adj_i = min(p_adj_i, p_adj_i+1)
    
    adj_p = np.zeros(n)
    for i, p in enumerate(sorted_p):
        rank = i + 1
        val = p * n / rank
        adj_p[i] = val
        
    # Enforce monotonicity (step-up) from end
    for i in range(n-2, -1, -1):
        adj_p[i] = min(adj_p[i], adj_p[i+1])
        
    # Map back to original order
    final_adj = np.zeros(n)
    final_adj[sorted_indices] = adj_p
    return final_adj

def check_core_splicing():
    print("\n--- Step 1: Core 4 Deep Dive (Splicing Effect Size) ---")
    print(f"{'Gene':<10} {'Region':<15} {'dPSI':<10} {'Category'}")
    
    for gene in CORE_GENES:
        for reg_name, reg_path in REGIONS.items():
            se_file = os.path.join(reg_path, "SE.MATS.JCEC.txt")
            if not os.path.exists(se_file): continue
            
            with open(se_file, 'r') as f:
                reader = csv.DictReader(f, delimiter='\t')
                for row in reader:
                    sym = row['geneSymbol'].replace('"', '')
                    if sym == gene:
                        dpsi = float(row['IncLevelDifference'])
                        fdr = float(row['FDR'])
                        
                        # Grading
                        abs_dpsi = abs(dpsi)
                        cat = "WEAK (<10%)"
                        if abs_dpsi >= 0.10: cat = "MODERATE (10-30%)"
                        if abs_dpsi >= 0.30: cat = "STRONG (>30%)"
                        
                        # Only print relevant hits (FDR < 0.05 or just the event)
                        # We print all for transparency
                        print(f"{gene:<10} {reg_name:<15} {dpsi:.4f}     {cat}")

def run_global_de(t2g, g2name):
    print("\n--- Step 3: Global Differential Expression (Stress Test) ---")
    print("Performing T-test on Log2(TPM+1) + Benjamini-Hochberg FDR Correction...")
    print("CRITERIA: FDR < 0.05 (No Fold-Change cutoff)")
    
    # Load metadata
    wt_samples = set()
    ko_samples = set()
    meta_dir = os.path.join(PROJECT_DIR, "data/metadata")
    if os.path.exists(meta_dir):
        for f in glob.glob(os.path.join(meta_dir, "wt_*.txt")):
            with open(f) as fh: wt_samples.update(fh.read().strip().split(','))
        for f in glob.glob(os.path.join(meta_dir, "ko_*.txt")):
            with open(f) as fh: ko_samples.update(fh.read().strip().split(','))
    
    # Clean IDs
    wt_samples = {os.path.basename(s.strip()).split('_')[0] for s in wt_samples if s.strip()}
    ko_samples = {os.path.basename(s.strip()).split('_')[0] for s in ko_samples if s.strip()}

    gene_data = {}
    quant_files = glob.glob(os.path.join(SALMON_DIR, "SRR*", "quant.sf"))
    
    for qf in quant_files:
        sname = os.path.basename(os.path.dirname(qf))
        group = "WT" if sname in wt_samples else ("KO" if sname in ko_samples else None)
        if group is None: continue
        
        sample_vals = {}
        with open(qf, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                tid = row['Name'].split('|')[0]
                tpm = float(row['TPM'])
                gid = t2g.get(tid)
                if not gid: continue
                sample_vals[gid] = sample_vals.get(gid, 0.0) + tpm
        
        for gid, val in sample_vals.items():
            if gid not in gene_data: gene_data[gid] = {'WT': [], 'KO': []}
            gene_data[gid][group].append(val)
            
    # Analyze
    from scipy import stats
    
    p_values = []
    gene_ids = []
    fc_values = []
    
    for gid, groups in gene_data.items():
        wts = groups['WT']
        kos = groups['KO']
        
        if len(wts) < 3 or len(kos) < 3: continue
        if np.mean(wts) < 1.0 and np.mean(kos) < 1.0: continue
        
        log_wts = np.log2(np.array(wts) + 1)
        log_kos = np.log2(np.array(kos) + 1)
        
        t_stat, p_val = stats.ttest_ind(log_wts, log_kos)
        
        mw = np.mean(wts)
        mk = np.mean(kos)
        fc = np.log2((mk + 0.01) / (mw + 0.01))
        
        if not np.isnan(p_val):
            p_values.append(p_val)
            gene_ids.append(gid)
            fc_values.append(fc)
            
    # Correction
    if not p_values:
        print("No genes passed filtering.")
        return

    adj_p_values = apply_bh_correction(p_values)
    
    sig_count = 0
    print(f"\n{'Gene':<15} {'Log2FC':<10} {'P-val':<10} {'FDR (BH)':<10}")
    
    for i, gid in enumerate(gene_ids):
        fdr = adj_p_values[i]
        fc = fc_values[i]
        
        # User Criterion: FDR < 0.05 ONLY (No FC cutoff)
        if fdr < 0.05:
            sig_count += 1
            if sig_count <= 20: # Top 20
                 name = g2name.get(gid, gid)
                 print(f"{name:<15} {fc:.3f}      {p_values[i]:.2e}      {fdr:.2e}")

    print(f"\n--- STRESS TEST RESULTS ---")
    print(f"Total Genes Analyzed: {len(gene_ids)}")
    print(f"Significant DE Genes (FDR < 0.05, NO FC cutoff): {sig_count}")
    
    ratio = "N/A"
    if sig_count > 0:
        ratio = f"{983/sig_count:.1f}:1" 
    
    print(f"Splicing (983) vs Expression ({sig_count}) Ratio: {ratio}")
    if sig_count > 400:
        print("RESULT: Splicing Dominance COLLAPSED (DE > 400).")
    else:
        print("RESULT: Splicing Dominance UPDHELD (DE < 400).")

if __name__ == "__main__":
    t2g, g2name = load_t2g()
    check_core_splicing()
    # check_sf_expression(t2g, g2name) # Already validated in Phase 4
    run_global_de(t2g, g2name)

