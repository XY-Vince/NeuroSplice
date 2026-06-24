# NeuroSplice Project Rules

These rules apply specifically to the NeuroSplice Adgrl3 multi-omic analysis and supersede any other global rules.

## Core Directives
1. **Metadata Quarantine**: The absolute source of truth is `data/metadata/GSE117357_complete_mapping.csv`. NEVER use the quarantined file `GSE117357_summary_DO_NOT_USE_wrong_dataset.csv`.
2. **Directionality**: ΔΨ = mean(WT) - mean(KO). Positive ΔΨ means higher inclusion in WT.
3. **Canonical Pipeline**: The core spine is scripts 58 -> 59 -> 51. Do not use legacy deprecated scripts.
4. **DESeq2 Design**: Formula must be `~Genotype`. DO NOT include a sex covariate as this metadata does not exist.
5. **Namespaces**: Ensembl `GeneID` suffixes (e.g., `.5`) must be stripped before joining with rMATS `geneSymbol`.

## Kill-Switches (Hard-Stop Conditions)
1. **Metadata Contamination**: Use of `GSE117357_summary_DO_NOT_USE_wrong_dataset.csv`.
2. **Fake Sex Covariates**: Including sex in the DESeq2 design.
3. **Reversed Directionality**: WT and KO order reversed in ΔΨ calculation.
4. **Namespace Collision**: Joining DEG and AS results without stripping `.version` Ensembl suffixes.

## Key Thresholds
- **Strict_effect**: Requires both FDR < 0.05 and an absolute inclusion difference |ΔΨ| ≥ threshold.
- *Note: Script 59 currently uses |ΔΨ| ≥ 0.05 for Strict_effect, but the handoff documentation cites |ΔΨ| ≥ 0.10. Clarification is pending.*

## Schema
When categorizing candidate events, strictly use this 5-tier classification scheme:
1. **A1_pan_region_same_event**: Same event recurrent at FDR < 0.05 in all 3 tissues.
2. **A1b_multi_region_same_event**: Same event recurrent in 2 of 3 tissues.
3. **B_region_restricted**: Event restricted to 1 tissue.
4. **C1_region_opposite_same_event**: Same event, opposite directions across tissues.
5. **C2_gene_multi_event_mixed**: Different events within one gene show mixed directions.

## Sashimi Plot Requirements
- Median junction depth ≥ 10 reads per group.
- At least 8/10 individual group replicates must meet this depth and generate valid PSI estimates.
