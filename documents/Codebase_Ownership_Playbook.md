# Codebase Ownership Playbook
## From Orchestrator → Scientific Owner Without Reading the Whole Codebase
**Version:** 1.3

This playbook provides a systematic strategy to transition from a high-level pipeline manager to a confident, audit-ready scientific owner of the Adgrl3 splicing analysis code without getting lost in implementation boilerplate.

*Document History: v1.3 adds defensive Q&A for A1 events, corrects the LRT statistical nomenclature, and separates baseline vs control DE claims. v1.2 patched factual errors regarding statistical tests, tiering descriptions, and claim-to-file scope. v1.1 integrated critical review corrections (metadata names, category label alignment, threshold specifications).*

---

## Phase 0 — Verify Production Path

Before reading code, verify that these scripts are the current production path. Deeply understanding a stale script is false ownership.

* Confirm which scripts were actually run for current results.
* Confirm command lines, output folders, and timestamps.
* Confirm manuscript-facing outputs do not come from deprecated scripts.
* Confirm README / SCRIPT_INDEX / handoff documents point to the same canonical spine.

---

## Phase 1 — Data-Flow Map

Trace the canonical data flow to confirm the machinery of your pipeline. 

The canonical three-script spine is:
`58_deseq2_local_canonical.py`
→ `59_rmats_qc_and_candidates.py`
→ `51_final_categorization.py`

* **Script 58** = DEG contrast / expression baseline.
* **Script 59** = rMATS QC, candidate extraction, directionality, thresholds.
* **Script 51** = final category logic and manuscript-facing triage.

---

## Phase 2 — Load-Bearing Scientific Decisions

Only ~25–50 lines contain the main scientific decisions, but ownership also requires verifying that these lines are actually reached, executed, and consumed downstream.

### 1. PyDESeq2 Wald Test Contrast & Covariates
* **Script**: `scripts/58_deseq2_local_canonical.py`
* **Mathematical representation**: Contrast is defined as `Genotype: KO vs WT`. 
* **Critical Check**: The design formula is `~Genotype`. While older documentation claimed the cohort was "all male", sex is entirely unavailable in the verified GEO/SRA metadata and could not be modeled. This is a known limitation, not a pipeline defect.
* **Directionality**: Positive LFC = higher in KO; Negative LFC = higher in WT.

### 2. rMATS Group Directionality Audit
* **Script**: `scripts/59_rmats_qc_and_candidates.py` (`validate_directionality` logic)
* **Mathematical representation**: $\Delta\Psi = \text{mean}(\Psi_{\text{WT}}) - \text{mean}(\Psi_{\text{KO}})$
* **Scientific Defense**: "To guarantee directional accuracy, we implemented an automated audit parsing the input BAM file lists and cross-referencing their SRRs against metadata before assigning splicing phenotypes. Positive $\Delta\Psi$ = higher inclusion in WT."

### 3. Significance Tiering
* **Script**: `scripts/59_rmats_qc_and_candidates.py`
* **Scientific Defense**: "We tiered alternative splicing events to decouple biological effect size from statistical significance. Script 59 gates on FDR < 0.05; Script 51 further requires an absolute inclusion difference |ΔΨ| ≥ 0.10 for the Strict_effect categorization tier."

### 4. JC Primary, JCEC Sensitivity Standard
* **Script**: `scripts/59_rmats_qc_and_candidates.py`
* **Scientific Defense**: "JC is used as the primary evidence stream because splice-junction reads directly support exon-exon connectivity. JCEC is retained as a sensitivity stream incorporating exon-body evidence. JC/JCEC discordance should be treated as a caution flag, not ignored."

### 5. Sashimi Plot Eligibility Gating
* **Script**: `scripts/59_rmats_qc_and_candidates.py`
* **Scientific Defense**: "Sashimi plot eligibility was strictly gated. Events required a median junction depth of $\ge 10$ reads per group, with at least 80% (8/10) of individual group replicates meeting this depth and generating valid PSI estimates."

### 6. Core Recurrence Categories
* **Script**: `scripts/51_final_categorization.py`
* **Category Terminology Harmonization**:
  - **A1_pan_region_same_event** = same event, same direction across all three regions
  - **A1b_multi_region_same_event** = same event, same direction across two regions (e.g., *Lrp8* SE, *Myo9b* RI)
  - **B_region_restricted** = strong event limited to one region (e.g., *Neil2*)
  - **C1_region_opposite_same_event** = same event, opposite directions across ≥2 regions (e.g., *Bcl2l11* RI)
  - **C2_gene_multi_event_mixed** = different events within one gene show mixed directions
* *Action*: Ensure manuscript and defense strictly adhere to this resolved 5-label scheme. (Note: Script 51 requires refactoring to natively emit these labels instead of legacy A/B/C/D labels).

---

## Phase 2.5 — Manual Ownership Probes

For 2–3 representative events (e.g. *Pts* Striatum SE ID 15876, *Lrp8* Prefrontal SE (A1b), and *Bcl2l11* RI (C1)):
1. Open the raw `.MATS.JC.txt` row.
2. Parse `IncLevel1` and `IncLevel2`.
3. Compute mean WT PSI.
4. Compute mean KO PSI.
5. Confirm $\Delta\Psi$ = WT − KO.
6. Compute median junction depth per group, count replicates with reads $\ge 10$, and verify both criteria (median $\ge 10$, $n \ge 8$).
7. Confirm FDR / $\Delta\Psi$ tier.
8. Confirm whether the event enters candidate extraction.
9. Confirm how script 51 categorizes the gene.
10. **Confirm sashimi plot** for the anchor event visually matches the $\Delta\Psi$ direction and junction read counts.

*(Note: Environment verification is critical here. Ensure you are running Python with versions matching the canonical manifest: e.g., PyDESeq2, pandas, and numpy).*

---

## Phase 3 — Reviewer 2 Defense Drill

Prepare to defend your code's methodology with these standard reviewer challenges:

1. **Q: Why did you use JC instead of JCEC as the primary evidence stream?**
   * *A:* "Splice-junction reads (JC) directly support exon-exon connectivity, whereas exon-body reads (JCEC) can inflate counts for short exons. We utilize JC as the primary metric and JCEC purely as a sensitivity check."
2. **Q: How did you prevent sample swaps or genotype assignment errors?**
   * *A:* "We implemented a fail-loud directionality gate that parses the raw BAM lists and maps their SRR accessions directly to the metadata matrix before execution."
3. **Q: What is the denominator for your significant AS event fractions?**
   * *A:* "The denominators are the exact counts of testable events with sufficient coverage to perform a likelihood ratio test (LRT, as distinct from the Wald test used in DESeq2) in each region/event-type."
4. **Q: How does your categorization handle genes with multiple, distinct splicing events across tissues?**
   * *A:* "We explicitly distinguish A1/A1b (identical genomic coordinate splicing event across regions) from C2 (same gene, different coordinates/exons with mixed directions) to separate region-independent regulation from general splicing instability."
5. **Q: The Pts retained-intron event shows positive $\Delta\Psi$. Why does this justify a KO-induced phenotype?**
   * *A:* "A positive $\Delta\Psi$ indicates higher inclusion in WT (lower in KO). A loss of the retained-intron form disrupts the normal homeostatic ratio of transcript isoforms."
6. **Q: Why is the Pts SE (exon skipping) event chosen as the main figure, and how does it relate to the BH4 mechanism?**
   * *A:* "The SE event provides the clearest visual anchor for a general audience. While the RI event holds the strongest statistical support and direct link to the functional domain, both are presented—the SE in the main figure for clarity, and the RI in supplements to fully support the biological mechanism."
7. **Q: How do you interpret events where the same exon shows opposite inclusion direction in different tissues (e.g., Bcl2l11)?**
   * *A:* "We classify these as C1_region_opposite_same_event, interpreting them as tissue-specific regulatory context rather than a directional artifact. This is supported by our automated SRR-to-genotype audit confirming no sample swaps occurred."
8. **Q: How do you handle events where PSI cannot be computed for some samples?**
   * *A:* "We require that at least 8 out of 10 replicates per genotype have valid, non-NaN PSI values before an event is considered for visualization."
9. **Q: Could low-depth events bias your sashimi plot selections?**
   * *A:* "No, because we strictly gate sashimi eligibility on the *median* coverage ($\ge 10$ reads) across replicates, preventing single-replicate outliers from inflating eligibility."
10. **Q: If you have no true A1 (pan-tissue) events, how is the "widespread splicing dysregulation" thesis supported?**
    * *A:* "While strict A1 events dropped out due to annotation stringency, A1b events (recurrent in 2 of 3 regions, like Lrp8 and Myo9b) still demonstrate that splicing instability is a multi-region, systemic consequence of Adgrl3 loss, not a single-tissue artifact."

---

## Phase 4 — Claim-to-File Trace

For each major manuscript claim, maintain an unbroken trace to the file, script, and threshold:

| Claim | Source file | Script | Gene/event | Threshold | Figure/table | Caveat |
|---|---|---|---|---|---|---|
| Splicing vs DEG burden ratio (5- to 16-fold) | `rmats_thresholded_summary.csv` vs `[Region]_DESeq2_significant.csv` | 58 & 59 | N/A | FDR < 0.05 | Results Text | Ratio must be recalculated if DESeq2 padj thresholds change. |
| Gene-level non-overlap between DEG and AS signals | `genomewide_deg_as_overlap_union.csv` | 58+59 cross-analysis | N/A | Strict×Strict (2 genes overlap) | DEG vs AS comparison | DEG IDs have Ensembl version suffixes (e.g. `.5`) that must be stripped before join |
| *Pts* Striatal event altered | `rmats_candidate_events_ranked.csv` | 59 | *Pts* Striatum SE 15876 | FDR < 0.05 | Primary Candidate Sashimi | Explicit decision: SE 15876 is the main visual anchor; RI 1270 is the supplementary mechanistic anchor. |
| Multi-region splicing instability | `final_category_gene_summary.csv` | 51 | *Lrp8*, *Myo9b* | A1b criteria | Categorization Table | Do not mix older A1/A2 terminology in writing. |
| Region-dependent regulation | `results/final_categorization/category_c_event_resolution.csv` | 51 | *Bcl2l11* RI | C1 criteria | Categorization Table | Ensure directions match. |
| Baseline DEG is minimal | `[Region]_DESeq2_significant.csv` | 58 | N/A | padj < 0.05 | DEG vs AS comparison | Refers to the background transcriptome remaining largely stable. |
| *Adgrl3* KO Confirmation | `[Region]_DESeq2_significant.csv` | 58 | *Adgrl3* | padj < 0.05 | DEG vs AS comparison | Positive control confirming the knockout paradigm worked. |

*(Note: Files like `genomewide_deg_as_overlap_union.csv` were generated via ad-hoc cross-analysis and may not reside in the default `results/` output directory natively).*

---

## Phase 5 — What to Skip

Most of these lines are not part of the scientific argument. Verify once that they do not alter biological labels, signs, thresholds, or event identity; then skip implementation details:

- **Provenance manifests & SHA256 calculations** (Administrative bookkeeping)
- **CLI Argparse code & Logging** (Plumbing)
- **Column cleaning** (Stripping quotes, dropping duplicated IDs)

**Caution:** Do not fully ignore PSI/read parsing helpers. Verify once that comma-separated PSI/read fields and NaNs are handled correctly, then skip re-reading them.

---

## Phase 6 — Kill-Switches

Stop trusting the result if:
- Metadata file differs from canonical mapping.
- `GSE117357_summary_DO_NOT_USE_wrong_dataset.csv` is used (contaminated — sole metadata source is `GSE117357_complete_mapping.csv`).
- WT/KO order is reversed.
- Sex covariate is included without verified sex metadata available. If sex metadata becomes available in the future, re-evaluate the design formula.
- DEG and AS comparisons are joined on raw `GeneID` (Ensembl) without stripping `.version` suffixes.
- Regions are pooled unintentionally.
- Deprecated HISAT2 BAMs feed final figures.
- $\Delta\Psi$ sign cannot be reproduced manually.
- Candidate category cannot be traced to event evidence.
- Sashimi plot lacks read/junction support.
- Human translation claim lacks ortholog or supporting evidence.
