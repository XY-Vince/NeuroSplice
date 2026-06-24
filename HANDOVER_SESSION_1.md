# HANDOVER: SESSION 1
**Project:** Adgrl3 Multi-omic Splicing Analysis (GSE117357)
**Version:** 1.3
**Date:** 2026-06-10
**Status:** Canonical processing complete. Moving to validation & dissemination.

*Document History: v1.3 resolves minor statistical representations, fixes the AS/DEG burden ratio, and documents Tox3 deprioritization. v1.2 patches factual errors regarding DEG stats, category ordering, and sashimi plot accounting. v1.1 integrates critical review corrections (metadata names, category label alignment, Core 6 definitions, and pipeline caveats).*

This document summarizes the current state of the NeuroSplice project, the key decisions and corrections made during this session, and the immediate next steps to bring the project to publication. Note: This document supersedes the candidate prioritization and category scheme in `NeuroSplice_Handoff_v7.0.md`, though pipeline specs and dataset metadata in that document remain canonical.

---

## 1. Current State of the Project

The computational pipeline has been completely refactored and hardened for local execution on the Mac. The pipeline now successfully relies on the canonical **Ensembl 102 GTF** reference, moving away from BAMs aligned against a deprecated annotation (GENCODE vM25).

### Completed Milestones:
*   **Canonical Pipeline Spine Locked:** We finalized the three core scripts:
    *   `58_deseq2_local_canonical.py` (Expression analysis)
    *   `59_rmats_qc_and_candidates.py` (Splicing quantification & filtering)
    *   `51_final_categorization.py` (Candidate classification)
*   **Genome-Wide Intersection (DEG vs AS):** We established a near-zero overlap between differential expression and alternative splicing. Only 2 genes show strict genome-wide overlap. None of our top prioritized AS candidates (*Pts*, *Lrp8*, *Myo9b*, *Bcl2l11*, *Unc13b*, *Neil2*) exhibit significant gene-level differential expression.
*   **Sashimi Plots Generated (15/15):** The standard `rmats2sashimiplot` tool failed due to MacOS `sqlite3` dependency issues. We successfully engineered a custom `pysam`+`matplotlib` script to bypass this and generated publication-quality plots for all 15 significant events across the Core 6 candidates. *(Note: My previous proofreading report had an arithmetic error—there are exactly 15 events: 5+3+3+1+2+1=15. All 15 generated successfully).*
*   **Documentation Hardened:** The `NeuroSplice_Handoff_v7.0.md` and `Codebase_Ownership_Playbook.md` have been meticulously updated to reflect the canonical pipeline and the factual corrections outlined below.

---

## 2. Key Decisions & Corrections Made

We made several critical, load-bearing scientific decisions that will defend the project against "Reviewer 2" level scrutiny:

### A. The "Sex Covariate" Correction
*   **Previous Assumption:** Older documentation stated the cohort was all-male and/or that `~sex + genotype` must be modeled.
*   **Decision:** We discovered that while older documentation and the original study claimed the cohort was "all male", sex is entirely unavailable in the verified GEO/SRA metadata for GSE117357. 
*   **Resolution:** The DESeq2 design formula has been locked to `~Genotype`. The absence of sex modeling is accurately documented as a known limitation of the public dataset, not a pipeline failure. Even if the cohort is indeed all-male, the absence of sex in the model would not be a confounding factor.

### B. Metadata Quarantine
*   **Finding:** The original summary file was found to be contaminated with samples from a different study (GSM5271282+), mixed genotypes, and a spurious/fake sex column.
*   **Resolution:** We renamed this file on disk to `GSE117357_summary_DO_NOT_USE_wrong_dataset.csv` and quarantined it. The absolute sole source of truth for the project's metadata is `GSE117357_complete_mapping.csv`.

### C. Identifier Harmonization Audit
*   **Finding:** An initial intersection of DEG and AS datasets yielded an incorrect "zero overlap" due to a namespace mismatch (joining Ensembl `GeneID` with rMATS `geneSymbol`).
*   **Resolution:** We implemented an audit trail and mapped Ensembl IDs (stripping `.version` suffixes like `.5` which break joins) to native Symbols, which confirmed 96.8–99.1% namespace coverage. This resulted in the accurate, defensible finding that the DEG/AS overlap is tiny (2 strict genes), but not absolutely zero.

### D. New Categorization Scheme
*   **Finding:** The old A1/A2/B/C/D scheme was inadequate for capturing complex splicing logic across tissues.
*   **Resolution:** We locked in a highly precise 5-tier scheme *(Note: Categories are ordered by cross-region scope, not by effect size or priority)*:
    1.  **A1_pan_region_same_event** (None currently assigned. *Annotation Update Impact:* *Lrp8* and *Myo9b* lost signal in one region after shifting to the Ensembl 102 GTF, dropping them to A1b. This reframes the instability as multi-region rather than strictly pan-tissue).
    2.  **A1b_multi_region_same_event** (*Lrp8*, *Myo9b*)
    3.  **B_region_restricted** (*Neil2*)
    4.  **C1_region_opposite_same_event** (*Bcl2l11* RI shows clear region-dependent opposite regulation)
    5.  **C2_gene_multi_event_mixed** (Different event loci in the same gene, e.g. *Pts*, *Unc13b*)
    *(Note: The current `51_final_categorization.py` script still outputs the legacy A1/A2/B/C/D labels. We manually mapped the top candidates to these fine-grained labels via `results/final_categorization/category_c_event_resolution.csv`. A future pipeline update must refactor script 51 to natively emit the 5-label schema).*

### E. Primary Visual Anchor Selected
*   **Decision:** While *Pts* Striatum RI (ID 1270) has the strongest statistical support (ΔΨ = +0.230) and directly supports the dopamine/BH4 mechanism, we selected **Pts Striatum SE ID 15876** as the main manuscript visual anchor because an exon-skipping event provides the clearest, most universally understood Sashimi plot for a general audience. This SE event belongs to the `Strict_effect` category (FDR < 0.05, |ΔΨ| ≥ 0.10), its JC/JCEC directions are completely consistent, and it satisfies all rigorous Sashimi plotting eligibility criteria. *Explicit action: SE 15876 will serve as the main figure visual anchor, while RI 1270 will be included as a supplementary figure to fully support the dopamine/BH4 narrative.*

### F. Tox3 Deprioritization
*   **Finding:** *Tox3* (PFC A3SS) had the second-largest effect size in the dataset (|ΔΨ| = 0.336) but lost support during the pipeline hardening.
*   **Resolution:** It was excluded from the Core 6 prioritization due to insufficient sashimi eligibility (median coverage depth was too low) after the annotation update to Ensembl 102.

---

## 3. The Core Scientific Thesis

**Adgrl3 loss drives a primary splicing phenotype that is mechanistically independent of transcriptional output.**
We see 5- to 16-fold more alternative splicing signals than differential expression signals (e.g., Striatum ~320 AS vs 22 DEG). The fact that genes like *Myo9b* show no genotype-dependent differential expression (e.g., in Striatum: LFC = +0.052, padj = 0.804) while harboring a strict splicing alteration (ΔΨ = +0.185 for its top event in the dataset) is a strong defense against claims that the splicing changes are just a secondary consequence of abundance shifts.

*(Note: Prior analyses included a Fisher's Exact Test against human digital phenotypes (see Handoff v7.0). This yielded 0/6 significant overlaps and has been deprioritized in favor of direct human ortholog transcript validation).*

---

## 4. Immediate Next Steps Checklist

To finish the project and formally wrap up the repository:

- [x] **Priorities 1–5 (Completed in prior sessions):** (1) Data recovery & backup, (2) Ensembl GTF refactoring, (3) Custom Sashimi plot generation, (4) DEG/AS genome-wide overlap analysis, (5) Handoff & Playbook documentation hardening.
- [ ] **Priority 6: Human Ortholog Verification:** Use tools like VAST-DB or GTEx to confirm whether the specific exons identified in our "Core 6" candidates (*Pts*, *Lrp8*, *Myo9b*, *Bcl2l11*, *Unc13b*, *Neil2*) are conserved in humans and known to be alternatively spliced in neuronal contexts. 
    *   *Neil2* is included due to having the highest dataset effect size (|ΔΨ| = 0.460).
    *   *Bcl2l11* (C1) is included because its clear region-dependent opposite regulation points to tissue-specific splicing biology rather than noise, making it a prime candidate for regulatory studies.
- [ ] **Priority 7: GitHub Push & Repository Freeze:** Selectively stage and commit the hardened canonical scripts, the summary CSVs, the Sashimi plots, and the updated documentation to the remote repository. (Avoid committing raw BAMs or massive `.MATS.JC.txt` files).
- [ ] **Priority 8: Copy Sashimi Plots:** Ensure that the generated sashimi plot files are copied from the external drive to the local `results/sashimi/` folder before the repository freeze.
- [ ] **Script 51 Refactor (Optional):** Update `51_final_categorization.py` to natively output the new 5-tier labeling schema.
- [ ] **Cross-Dataset Replication (Long-term):** Find independent *Adgrl3* or ADHD murine models to test if the Core 6 directionality holds true in other cohorts.
