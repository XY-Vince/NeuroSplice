# NeuroSplice Script Index

This repository contains 74 scripts representing the complete history of the NeuroSplice project, from initial data download to final visualization.

To ensure reproducibility and prevent execution of deprecated logic, the scripts are strictly categorized below.

## 🟢 Canonical Spine (Production)
These are the ONLY scripts required to generate the final manuscript-facing results. They must be executed in this exact order:
1. `58_deseq2_local_canonical.py` - Performs canonical DEG analysis using PyDESeq2 (`~Genotype`).
2. `59_rmats_qc_and_candidates.py` - Performs rigorous rMATS QC, applies junction depth thresholds (median ≥ 10), and extracts candidate events.
3. `51_final_categorization.py` - Categorizes the candidate events into the 5-tier schema (A1, A1b, B, C1, C2). *(Note: Refactored to natively output the 5-tier labels).*

## 🟡 Upstream Preprocessing
These scripts are responsible for downloading data and building the raw alignment/quantification matrices. They are generally not re-run unless starting from scratch.
* `00_download_references.sh` through `06_process_existing.sh` - Reference and raw fastq download.
* `43_reprocess_gse117357_star.sh`, `44_star_docker.sh` - STAR alignment (Ensembl 102).
* `45_rmats_gse117357.sh`, `45_rmats_local_canonical.sh` - Raw rMATS execution.
* `56_featurecounts.sh`, `56_featurecounts_local_canonical.sh` - Gene-level read quantification for DESeq2.

## ⚪ Independent Datasets & Validations
Scripts used for cross-dataset validations and ortholog mapping.
* `37_fetch_gse173926_metadata.py` through `40_process_gse173926.sh` - Processing an independent dataset (GSE173926).
* `52_validate_gse173926.py`, `53_rmats_gse173926.sh`, `54_analyze_dataset_overlap.py` - GSE173926 comparisons.
* `60_conservation_analysis.py`, `55_human_conservation_analysis.py`, `47_human_translation_gwas.py` - Human ortholog and GWAS exploration.
* `27_setup_magma.sh` through `36_run_magma.sh` - MAGMA gene-set analysis.

## 🔴 Deprecated / Exploratory (Do Not Use for Final Results)
These scripts represent older pipeline versions (e.g. using GENCODE vM25 BAMs, legacy categorization logic, or alternative tools like salmon/stringtie).
* `10_rmats.sh` through `17_run_dapars.sh` - Legacy splicing and APA tools.
* `21_fix_stringtie.sh` through `26_disease_context.py` - Legacy network and context checks.
* `29_star_sparse_index.sh`, `30_run_star_pilot.sh`, `31_run_rmats_pilot.sh` - Pilot runs.
* `46_analyze_rmats_gse117357.py`, `48_verify_core_genes.py`, `49_find_strict_core_events.py`, `50_categorize_genes.py` - Precursors to the canonical 59->51 pipeline.
* `58_deseq2_analysis.py`, `59_fisher_test.py` - Precursors to the canonical scripts.
* `test_star_*.sh`, `debug_fifo.sh`, `transfer_to_bouchet.sh` - HPC and local debugging utilities.
