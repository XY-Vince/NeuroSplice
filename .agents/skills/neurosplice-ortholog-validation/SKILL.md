---
name: neurosplice-ortholog-validation
description: Validates human orthologs for the Adgrl3 Core 6 splicing candidates using VAST-DB and GTEx.
---

# NeuroSplice Ortholog Validation Skill

## Goal
Verify whether the specific exons identified in our "Core 6" mouse candidates are conserved in humans and whether they are known to be alternatively spliced in human neuronal contexts. This addresses Priority 6.

## The Core 6 Candidates
1. **Pts** (SE ID 15876, RI ID 1270) - C2_gene_multi_event_mixed
2. **Lrp8** (A3SS, SE) - A1b
3. **Myo9b** (A3SS, RI, SE) - A1b
4. **Bcl2l11** (RI) - C1
5. **Unc13b** (SE) - C2_gene_multi_event_mixed
6. **Neil2** - B (largest effect size, |ΔΨ| = 0.460)

## Workflow

1. **Ortholog Mapping**
   - For each mouse gene symbol above, identify the human ortholog (typically the same symbol, uppercase, e.g., *Pts* -> *PTS*).
   
2. **Coordinate & Event Translation**
   - Use Ensembl BioMart or UCSC Genome Browser (LiftOver) to map the specific mouse exon coordinates from the rMATS output (`results/rmats_gse117357/qc/rmats_candidate_events_ranked.csv`) to the human genome (GRCh38).

3. **VAST-DB Query**
   - Query VAST-DB (vastdb.crg.eu) using the human gene symbol.
   - Check if the orthologous exon is annotated as an alternative splicing event (e.g., cassette exon).
   - Record the `dPSI` or inclusion levels in human nervous system tissues (e.g., brain, cortex).

4. **GTEx Query**
   - Query the GTEx Portal (gtexportal.org) for the human gene.
   - Check the "Transcript" or "Exon" expression levels across brain tissues to confirm neuronal expression of the relevant isoforms.

5. **Reporting**
   - Output the findings into `results/conservation/core6_ortholog_validation.md`.
   - For each gene, explicitly state:
     - Is the exon conserved? (Yes/No - requires coordinate LiftOver or direct sequence mapping).
     - Is it alternatively spliced in humans? (Yes/No - requires VAST-DB or similar evidence).
     - Is there evidence of neuronal expression/regulation? (Yes/No - requires GTEx or VAST-DB brain tissue evidence with PSI > 5%).
     - Database source links.
