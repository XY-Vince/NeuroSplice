# AI Context Handoff: NeuroSplice Project
**Date:** 2026-06-24

## System Context
You are stepping into the "NeuroSplice" project, a multi-omic analysis of Alternative Splicing (AS) vs. Differential Expression (DEG) in an *Adgrl3* knockout mouse model (GSE117357). 
The core pipeline is hardened, and the primary analysis is complete. Your focus is on finalizing biological validation (human orthologs) and preparing the repository for publication.

**Experimental Design:** 10 WT vs 10 KO per brain region.
**Tissues Analyzed:** Hippocampus, Prefrontal Cortex (PFC), and Striatum (60 samples total).

## Core Scientific Thesis
**Adgrl3 loss drives a primary splicing phenotype that is mechanistically independent of transcriptional output.**
* **Burden Ratio:** There are 5- to 16-fold more AS events than DEGs across regions (e.g., Striatum has ~320 AS events vs 22 DEGs).
* **Independence:** Widespread splicing dysregulation occurs without significant changes to the background transcriptome (minimal DEG baseline, though *Adgrl3* itself is significantly downregulated as a positive control).
* **Example:** *Myo9b* harbors severe splicing alterations (Mean |ΔΨ| = 0.185 across its 3 events) despite showing near-zero differential expression (Striatum LFC = +0.052, padj = 0.804), strongly supporting the claim that AS is not just a secondary effect of expression shifts.

## Canonical Pipeline Spine
The local production path uses **Ensembl 102 GTF** (deprecated vM25 BAMs have been abandoned).
1. `scripts/58_deseq2_local_canonical.py`: DEG analysis (Wald Test). Design formula is strictly `~Genotype`. (Sex is unavailable in the metadata and cannot be modeled).
2. `scripts/59_rmats_qc_and_candidates.py`: AS quantification. Uses Likelihood Ratio Test (LRT). Gates candidates strictly on FDR < 0.05 and median junction coverage ≥ 10.
3. `scripts/51_final_categorization.py`: Splicing logic classification. Gates the `Strict_effect` tier at |ΔΨ| ≥ 0.10.

## Critical Technical & Math Details
* **Directionality:** $\Delta\Psi = \text{mean}(\Psi_{\text{WT}}) - \text{mean}(\Psi_{\text{KO}})$. Positive ΔΨ means higher inclusion in Wild Type (or lower inclusion in Knockout).
* **Metadata Quarantine:** The file `data/metadata/GSE117357_summary_DO_NOT_USE_wrong_dataset.csv` is quarantined due to contamination. The absolute source of truth is `data/metadata/GSE117357_complete_mapping.csv`.
* **Sashimi Plots:** Generated for 15/15 events across the "Core 6" candidates via a custom `pysam`+`matplotlib` script (bypassing `rmats2sashimiplot` due to sqlite3 issues). Plot files currently reside on an external drive. *(Eligibility required median junction depth ≥ 10 and ≥ 8/10 replicates with valid PSI per group).*
* **Namespace Join:** Ensembl `GeneID` suffixes (e.g., `.5`) must be stripped before joining with rMATS `geneSymbol`.

## Core 6 Candidates & 5-Tier Schema
Our candidates are classified by their cross-region behavior. Categories are ordered by cross-region scope, not effect size:
1.  **A1_pan_region_same_event**: No candidates. (Lost due to annotation stringency).
2.  **A1b_multi_region_same_event**: *Lrp8* (A3SS, SE), *Myo9b* (A3SS, RI, SE). (Supports systemic instability thesis).
3.  **B_region_restricted**: *Neil2* (|ΔΨ| = 0.460).
4.  **C1_region_opposite_same_event**: *Bcl2l11* (RI). *(Opposite-direction RI regulation across regions suggests tissue-specific regulatory biology, not noise).*
5.  **C2_gene_multi_event_mixed**: *Pts*, *Unc13b*. 
    * *Note on Pts*: *Pts* Striatum SE ID 15876 is the main visual anchor (clearest plot), but *Pts* Striatum RI ID 1270 is the supplementary mechanistic anchor due to its high support (ΔΨ = +0.230).

* **Note on Legacy Labels:** Script 51's CSV output uses legacy A/B/C/D groupings (e.g., *Pts* appears in `category_a` CSV). The 5-tier labels here were manually assigned based on event-level analysis in `results/final_categorization/category_c_event_resolution.csv`.
* **Note on Tox3:** Excluded due to insufficient median coverage depth, despite having the second-largest effect size in the dataset.

## Kill-Switches (Do Not Trust Results If:)
1. You use the quarantined summary CSV instead of `data/metadata/GSE117357_complete_mapping.csv`.
2. A sex covariate is included in DESeq2 (this metadata does not exist; its inclusion means you are using fake/imputed data).
3. WT and KO order is reversed in the $\Delta\Psi$ calculation.
4. DEG and AS joins fail to strip Ensembl `.version` suffixes.

## Immediate Next Steps (Your Tasks)
1. **Priority 6 (Human Ortholog Verification):** Use tools like VAST-DB or GTEx to confirm if the specific exons identified in the Core 6 are conserved in humans and known to be alternatively spliced in neuronal contexts.
2. **Priority 7 (Sashimi Rescue):** Ensure generated sashimi plots are copied from the external drive to the local `results/sashimi/` folder before the repository freeze.
3. **Priority 8 (GitHub Push):** Selectively stage and commit canonical scripts, summary CSVs, and documentation to the remote. (No BAMs/large .MATS files).
4. **Script 51 Refactor (Optional):** Update `51_final_categorization.py` to natively output the 5-tier schema (it currently outputs legacy A/B/C/D labels).

*Note: For the full audit trail and reviewer defense strategies, refer to `HANDOVER_SESSION_1.md` and `documents/Codebase_Ownership_Playbook.md` (both at v1.3).*
