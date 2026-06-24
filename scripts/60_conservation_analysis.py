#!/usr/bin/env python3
"""
Phase 4: Virtual PCR & Conservation Analysis (Revised)
======================================================
For each Category A/B1 candidate gene, this script:
  1. Extracts the top rMATS splicing event (mm10 coordinates)
  2. Looks up human ortholog via Ensembl REST API (gene-level mapping)
  3. Verifies human gene has exon structure (transcript/exon count)
  4. Queries Open Targets for psychiatric GWAS associations
  5. Produces a unified conservation summary table

Note on LiftOver:
  The Ensembl REST API does not support cross-species coordinate mapping.
  We use orthology-based gene lookup instead. Exon-level coordinate mapping
  (mm10 exon → hg38 exon) requires UCSC liftOver binary + chain file,
  which is deferred to the HPC session.
"""

import os
import json
import time
import sys
import requests
import pandas as pd

PROJECT = os.path.expanduser("~/Desktop/WXY/NeuroSplice")
RMATS_DIR = os.path.join(PROJECT, "results/rmats_gse117357_hpc")
OUT_DIR = os.path.join(PROJECT, "results/conservation")
os.makedirs(OUT_DIR, exist_ok=True)

REGIONS_MAP = {"hippocampus": "Hippocampus", "prefrontal": "PFC", "striatum": "Striatum"}
EVENT_TYPES = ["SE", "A3SS", "A5SS", "MXE", "RI"]

TARGETS = {
    "Pts":    {"human": "PTS",    "category": "A"},
    "Lrp8":   {"human": "LRP8",   "category": "A"},
    "Myo9b":  {"human": "MYO9B",  "category": "A"},
    "Neil2":  {"human": "NEIL2",  "category": "B1"},
    "Unc13b": {"human": "UNC13B", "category": "B1"},
    "Tox3":   {"human": "TOX3",   "category": "B1"},
}

HEADERS = {"Content-Type": "application/json"}


def api_get(url, retries=2):
    """GET with retry and rate limiting."""
    for attempt in range(retries + 1):
        try:
            resp = requests.get(url, headers=HEADERS, timeout=20)
            time.sleep(0.5)
            if resp.status_code == 200:
                return resp.json()
            elif resp.status_code == 429:
                time.sleep(2)
                continue
            else:
                return None
        except requests.RequestException:
            if attempt < retries:
                time.sleep(1)
            else:
                return None
    return None


def api_post(url, payload, retries=2):
    """POST with retry and rate limiting (mirrors api_get for GraphQL endpoints)."""
    for attempt in range(retries + 1):
        try:
            resp = requests.post(url, json=payload, headers=HEADERS, timeout=20)
            time.sleep(0.5)
            if resp.status_code == 200:
                return resp.json()
            elif resp.status_code == 429:
                time.sleep(2)
                continue
            else:
                return None
        except requests.RequestException:
            if attempt < retries:
                time.sleep(1)
            else:
                return None
    return None


# ============================================================
# Step 4.1: Extract top rMATS event per gene
# ============================================================
def extract_top_events():
    print("=" * 70)
    print("Step 4.1: Extracting Top rMATS Events")
    print("=" * 70)

    all_events = []
    for rdir, rname in REGIONS_MAP.items():
        for ev in EVENT_TYPES:
            fpath = os.path.join(RMATS_DIR, rdir, f"{ev}.MATS.JCEC.txt")
            if not os.path.exists(fpath):
                continue
            df = pd.read_csv(fpath, sep="\t")
            df["geneSymbol"] = df["geneSymbol"].str.strip('"')
            for gene in TARGETS:
                hits = df[(df["geneSymbol"] == gene) & (df["FDR"] < 0.05)]
                for _, row in hits.iterrows():
                    entry = {
                        "gene": gene, "region": rname, "event_type": ev,
                        "chr": row["chr"], "strand": row["strand"],
                        "fdr": row["FDR"], "dpsi": row["IncLevelDifference"],
                        "abs_dpsi": abs(row["IncLevelDifference"]),
                    }
                    if ev == "SE":
                        entry["start"] = int(row["exonStart_0base"])
                        entry["end"] = int(row["exonEnd"])
                    elif ev in ["A3SS", "A5SS"]:
                        entry["start"] = int(row["longExonStart_0base"])
                        entry["end"] = int(row["longExonEnd"])
                    elif ev == "RI":
                        entry["start"] = int(row["riExonStart_0base"])
                        entry["end"] = int(row["riExonEnd"])
                    elif ev == "MXE":
                        entry["start"] = int(row["1stExonStart_0base"])
                        entry["end"] = int(row["1stExonEnd"])
                    all_events.append(entry)

    events_df = pd.DataFrame(all_events)
    if events_df.empty:
        import sys
        sys.exit("ERROR: No significant events found for any target gene. Check RMATS_DIR and gene symbols.")

    # Fix: use drop_duplicates instead of groupby().first() to safely preserve
    # the sort order established by sort_values — groupby does not guarantee it.
    top = (
        events_df
        .sort_values("abs_dpsi", ascending=False)
        .drop_duplicates(subset=["gene"], keep="first")
        .reset_index(drop=True)
        .drop(columns=["abs_dpsi"])
    )

    for _, row in top.iterrows():
        cat = TARGETS[row["gene"]]["category"]
        print(f"  [{cat}] {row['gene']:8s} | {row['region']:5s} {row['event_type']:4s} | "
              f"{row['chr']}:{row['start']}-{row['end']} ({row['strand']}) | "
              f"ΔΨ={row['dpsi']:+.3f}  FDR={row['fdr']:.2e}")

    return top


# ============================================================
# Step 4.2-4.3: Human Ortholog Lookup & Transcript Verification
# ============================================================
def lookup_human_orthologs():
    print(f"\n{'=' * 70}")
    print("Step 4.2-4.3: Human Ortholog Lookup & Transcript Verification")
    print("=" * 70)

    results = []
    for mouse_gene, info in TARGETS.items():
        human_gene = info["human"]

        # Lookup human gene
        data = api_get(f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{human_gene}?expand=1")
        if not data:
            print(f"  ✗ {mouse_gene:8s} → {human_gene:8s} — Ensembl lookup failed")
            results.append({
                "gene": mouse_gene, "human_gene": human_gene,
                "ensembl_id": None, "hg38_chr": None,
                "hg38_start": None, "hg38_end": None,
                "n_transcripts": 0, "n_exons": 0,
                "ortholog_status": "FAILED",
            })
            continue

        ens_id = data["id"]
        hg38_chr = f"chr{data['seq_region_name']}"
        hg38_start = data["start"]
        hg38_end = data["end"]
        strand = "+" if data["strand"] == 1 else "-"
        biotype = data.get("biotype", "unknown")

        # Get exon count
        exon_data = api_get(f"https://rest.ensembl.org/overlap/id/{ens_id}?feature=exon")
        n_exons = 0
        n_tx = 0
        if exon_data:
            n_exons = len(exon_data)
            # Fix: filter out empty strings before counting, to avoid inflating count
            # when a JSON entry is missing the "Parent" key.
            n_tx = len({ex.get("Parent") for ex in exon_data if ex.get("Parent")})

        marker = "✓" if biotype == "protein_coding" else "⚠"
        print(f"  {marker} {mouse_gene:8s} → {human_gene:8s} ({ens_id}) | "
              f"hg38: {hg38_chr}:{hg38_start}-{hg38_end} ({strand}) | "
              f"{n_exons} exons, {n_tx} transcripts | {biotype}")

        results.append({
            "gene": mouse_gene, "human_gene": human_gene,
            "ensembl_id": ens_id, "hg38_chr": hg38_chr,
            "hg38_start": hg38_start, "hg38_end": hg38_end,
            "hg38_strand": strand, "biotype": biotype,
            "n_transcripts": n_tx, "n_exons": n_exons,
            "ortholog_status": "CONFIRMED" if biotype == "protein_coding" else "NON_CODING",
        })

    return pd.DataFrame(results)


# ============================================================
# Step 4.5: Open Targets GWAS Check
# ============================================================
def check_gwas(ortholog_df):
    print(f"\n{'=' * 70}")
    print("Step 4.5: Psychiatric GWAS Association Check (Open Targets)")
    print("=" * 70)

    OT_URL = "https://api.platform.opentargets.org/api/v4/graphql"
    psych_keywords = [
        "adhd", "attention deficit", "hyperactiv", "schizophren", "bipolar",
        "depress", "anxiety", "autism", "psychiatric", "neurodevelop",
        "intellectual disabil", "cognitive",
    ]

    results = []
    for _, row in ortholog_df.iterrows():
        gene = row["gene"]
        ens_id = row.get("ensembl_id")
        if not ens_id:
            results.append({"gene": gene, "n_psych": 0, "psych_details": ""})
            continue

        query = """
        query($ensemblId: String!) {
          target(ensemblId: $ensemblId) {
            approvedSymbol
            associatedDiseases(page: {size: 500, index: 0}) {
              rows { disease { id name } score }
            }
          }
        }
        """
        # Fix: use api_post() retry wrapper instead of bare requests.post()
        # to handle transient 429/502 errors from the GraphQL API.
        try:
            data = api_post(OT_URL, {"query": query, "variables": {"ensemblId": ens_id}})

            if data is not None:
                target = data.get("data", {}).get("target")
                if target and target.get("associatedDiseases"):
                    diseases = target["associatedDiseases"]["rows"]
                    psych_hits = []
                    for d in diseases:
                        name = d["disease"]["name"].lower()
                        if any(kw in name for kw in psych_keywords):
                            psych_hits.append(f"{d['disease']['name']} (score={d['score']:.2f})")

                    if psych_hits:
                        print(f"  🟢 {row['human_gene']:8s} — {len(psych_hits)} psychiatric hit(s):")
                        for hit in psych_hits[:4]:
                            print(f"       • {hit}")
                    else:
                        print(f"  ⚪ {row['human_gene']:8s} — {len(diseases)} disease associations, 0 psychiatric")

                    results.append({
                        "gene": gene,
                        "total_disease_associations": len(diseases),
                        "n_psych": len(psych_hits),
                        "psych_details": "; ".join(psych_hits[:5]),
                    })
                else:
                    print(f"  ⚪ {row['human_gene']:8s} — no associations found")
                    results.append({"gene": gene, "total_disease_associations": 0,
                                    "n_psych": 0, "psych_details": ""})
            else:
                print(f"  ✗ {row['human_gene']:8s} — OT API failed after retries")
                results.append({"gene": gene, "n_psych": 0, "psych_details": ""})
        except requests.RequestException as e:
            print(f"  ✗ {row['human_gene']:8s} — Network/API error: {e}")
            results.append({"gene": gene, "n_psych": 0, "psych_details": ""})

    return pd.DataFrame(results)


# ============================================================
# Main
# ============================================================
def main():
    print("=" * 70)
    print("Phase 4: Virtual PCR & Conservation Analysis")
    print("=" * 70)

    if not os.path.exists(RMATS_DIR):
        sys.exit(f"ERROR: RMATS directory not found: {RMATS_DIR}")

    # Step 4.1
    top_events = extract_top_events()
    top_events.to_csv(os.path.join(OUT_DIR, "top_events_mm10.csv"), index=False)

    # Steps 4.2 & 4.3
    ortholog_df = lookup_human_orthologs()
    ortholog_df.to_csv(os.path.join(OUT_DIR, "human_orthologs.csv"), index=False)

    # Step 4.5
    gwas_df = check_gwas(ortholog_df)
    gwas_df.to_csv(os.path.join(OUT_DIR, "gwas_associations.csv"), index=False)

    # ---- Final Summary ----
    print(f"\n{'=' * 70}")
    print("CONSERVATION SUMMARY TABLE")
    print("=" * 70)

    merged = top_events[["gene", "region", "event_type", "dpsi", "fdr"]].copy()
    merged = merged.merge(ortholog_df[["gene", "human_gene", "ensembl_id", "hg38_chr",
                                        "hg38_start", "hg38_end", "n_transcripts",
                                        "n_exons", "ortholog_status"]], on="gene")
    merged = merged.merge(gwas_df[["gene", "n_psych", "psych_details"]], on="gene")

    print(f"\n{'Gene':8s} {'Cat':3s} {'Region':5s} {'Ev':4s} {'ΔΨ':>6s} "
          f"{'Human':8s} {'Ensembl':18s} {'Ortho':10s} {'Psych':5s}")
    print("-" * 80)
    for _, r in merged.iterrows():
        cat = TARGETS[r["gene"]]["category"]
        print(f"{r['gene']:8s} {cat:3s}  {r['region']:5s} {r['event_type']:4s} "
              f"{r['dpsi']:+.3f}  {r['human_gene']:8s} {r['ensembl_id']:18s} "
              f"{r['ortholog_status']:10s} {int(r['n_psych']):5d}")

    # Save merged
    merged.to_csv(os.path.join(OUT_DIR, "conservation_summary.csv"), index=False)
    print(f"\nAll outputs saved to: {OUT_DIR}")

    # Print key findings
    print(f"\n{'=' * 70}")
    print("KEY FINDINGS")
    print("=" * 70)
    confirmed = merged[merged["ortholog_status"] == "CONFIRMED"]
    print(f"  ✓ {len(confirmed)}/6 genes have confirmed human protein-coding orthologs")
    with_psych = merged[merged["n_psych"] > 0]
    print(f"  ✓ {len(with_psych)}/6 genes have psychiatric GWAS associations in humans")
    for _, r in with_psych.iterrows():
        details = r["psych_details"].split(";")[0].strip() if r["psych_details"] else ""
        print(f"       {r['gene']:8s} → {details}")
    no_psych = merged[merged["n_psych"] == 0]
    if len(no_psych) > 0:
        print(f"  ⚪ {len(no_psych)}/6 genes have NO psychiatric associations:")
        for _, r in no_psych.iterrows():
            print(f"       {r['gene']}")


if __name__ == "__main__":
    main()
