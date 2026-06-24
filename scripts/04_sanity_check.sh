#!/bin/bash
# =============================================================================
# GSE117357 Sanity Check: Validate Adgrl3 Knockout
# Critical validation before downstream analysis
# =============================================================================

set -euo pipefail

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate combio

# Configuration
PROJECT_DIR="/Users/guozhenghui/Desktop/WXY/ComBio/GSE117357_Analysis"
ALIGN_DIR="${PROJECT_DIR}/aligned/bam"
RESULTS_DIR="${PROJECT_DIR}/results"
META_DIR="${PROJECT_DIR}/data/metadata"
LOG_DIR="${PROJECT_DIR}/logs"

# Reference (UPDATE THIS)
GTF="/Users/guozhenghui/Desktop/WXY/ComBio/GSE117357_Analysis/reference/gencode.vM25.annotation.gtf"

# =============================================================================
# Step 1: Extract Adgrl3 Counts from STAR Output
# =============================================================================
echo "[$(date)] Extracting Adgrl3 counts from STAR ReadsPerGene files..."

mkdir -p "${RESULTS_DIR}/sanity_check"

# Find all ReadsPerGene files
echo "Sample,Adgrl3_unstranded,Adgrl3_sense,Adgrl3_antisense" > "${RESULTS_DIR}/sanity_check/adgrl3_star_counts.csv"

for gene_file in "${ALIGN_DIR}"/*_ReadsPerGene.out.tab; do
    if [[ -f "${gene_file}" ]]; then
        SAMPLE=$(basename "${gene_file}" _ReadsPerGene.out.tab)
        
        # Search for Adgrl3 or Lphn3 (old name)
        COUNTS=$(grep -E "Adgrl3|Lphn3" "${gene_file}" | head -1 || echo "N/A\tN/A\tN/A\tN/A")
        
        if [[ "${COUNTS}" != "N/A"* ]]; then
            echo "${SAMPLE},$(echo ${COUNTS} | awk '{print $2","$3","$4}')" >> \
                "${RESULTS_DIR}/sanity_check/adgrl3_star_counts.csv"
        else
            echo "${SAMPLE},NA,NA,NA" >> "${RESULTS_DIR}/sanity_check/adgrl3_star_counts.csv"
        fi
    fi
done

echo "[$(date)] Adgrl3 counts extracted"

# =============================================================================
# Step 2: Run featureCounts for Full Gene Counts (Optional)
# =============================================================================
run_featurecounts() {
    echo "[$(date)] Running featureCounts..."
    
    # Find all BAM files
    BAMS=$(ls "${ALIGN_DIR}"/*_Aligned.sortedByCoord.out.bam 2>/dev/null | tr '\n' ' ')
    
    if [[ -z "${BAMS}" ]]; then
        echo "[ERROR] No BAM files found in ${ALIGN_DIR}"
        exit 1
    fi
    
    featureCounts \
        -T 8 \
        -p \
        -a "${GTF}" \
        -o "${RESULTS_DIR}/sanity_check/gene_counts.txt" \
        ${BAMS} \
        2>&1 | tee "${LOG_DIR}/featurecounts.log"
    
    # Extract Adgrl3
    echo "[$(date)] Extracting Adgrl3 from featureCounts..."
    head -2 "${RESULTS_DIR}/sanity_check/gene_counts.txt" > \
        "${RESULTS_DIR}/sanity_check/adgrl3_featurecounts.txt"
    grep -E "Adgrl3|Lphn3" "${RESULTS_DIR}/sanity_check/gene_counts.txt" >> \
        "${RESULTS_DIR}/sanity_check/adgrl3_featurecounts.txt" || \
        echo "Adgrl3 not found in featureCounts output"
}

# =============================================================================
# Step 3: R Script for Validation Plot
# =============================================================================
cat > "${RESULTS_DIR}/sanity_check/validate_knockout.R" << 'RSCRIPT'
#!/usr/bin/env Rscript
# Adgrl3 Knockout Validation Script

library(ggplot2)
library(dplyr)

# ============================================================================
# Load Data
# ============================================================================
args <- commandArgs(trailingOnly = TRUE)
project_dir <- args[1]

counts_file <- file.path(project_dir, "results/sanity_check/gene_counts.txt")
meta_file <- file.path(project_dir, "data/metadata/sample_metadata.csv")
output_dir <- file.path(project_dir, "results/sanity_check")

# Read counts
counts <- read.table(counts_file, header=TRUE, row.names=1, skip=1, check.names=FALSE)

# Clean column names (remove path prefix)
colnames(counts) <- gsub(".*/(SRR[0-9]+)_.*", "\\1", colnames(counts))

# Find Adgrl3/Lphn3
adgrl3_row <- grep("Adgrl3|Lphn3", rownames(counts), value=TRUE)
if(length(adgrl3_row) == 0) {
  stop("Adgrl3/Lphn3 not found in count matrix!")
}

adgrl3_counts <- as.numeric(counts[adgrl3_row[1], -(1:5)])
sample_names <- colnames(counts)[-(1:5)]

# ============================================================================
# Load Metadata
# ============================================================================
if(file.exists(meta_file)) {
  meta <- read.csv(meta_file)
} else {
  # Create placeholder metadata
  meta <- data.frame(
    SRR = sample_names,
    condition = ifelse(grepl("KO|ko|knockout", sample_names, ignore.case=TRUE), "KO", "WT")
  )
  warning("Using inferred metadata - please verify condition assignments!")
}

# ============================================================================
# Calculate CPM
# ============================================================================
total_counts <- colSums(counts[, -(1:5)])
adgrl3_cpm <- (adgrl3_counts / total_counts) * 1e6

df <- data.frame(
  sample = sample_names,
  raw_counts = adgrl3_counts,
  cpm = adgrl3_cpm
)

# Merge with metadata
df <- merge(df, meta, by.x="sample", by.y="SRR", all.x=TRUE)

# ============================================================================
# Validation Check
# ============================================================================
if("condition" %in% colnames(df)) {
  ko_mean <- mean(df$cpm[df$condition == "KO"], na.rm=TRUE)
  wt_mean <- mean(df$cpm[df$condition == "WT"], na.rm=TRUE)
  
  ratio <- ko_mean / wt_mean
  
  cat("\n============ ADGRL3 KNOCKOUT VALIDATION ============\n")
  cat(sprintf("WT mean CPM: %.2f\n", wt_mean))
  cat(sprintf("KO mean CPM: %.2f\n", ko_mean))
  cat(sprintf("KO/WT ratio: %.4f (%.1f%% of WT)\n", ratio, ratio * 100))
  cat("\n")
  
  if(ratio < 0.05) {
    cat("✓ PASS: KO samples show <5% of WT expression\n")
    cat("  Knockout appears successful. Proceed with analysis.\n")
  } else if(ratio < 0.20) {
    cat("⚠ WARNING: KO samples show 5-20% of WT expression\n")
    cat("  This could indicate:\n")
    cat("  - Partial knockout (heterozygous?)\n")
    cat("  - Residual expression from truncated transcript\n")
    cat("  Review IGV tracks before proceeding.\n")
  } else {
    cat("✗ FAIL: KO samples show >20% of WT expression\n")
    cat("  This suggests:\n")
    cat("  - Metadata error (samples mislabeled?)\n")
    cat("  - Knockout failure\n")
    cat("  - Wrong gene targeted\n")
    cat("  DO NOT PROCEED. Investigate cause first.\n")
  }
  cat("=====================================================\n\n")
}

# ============================================================================
# Generate Plot
# ============================================================================
p <- ggplot(df, aes(x = condition, y = cpm, fill = condition)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 3, alpha = 0.8) +
  scale_fill_manual(values = c("KO" = "#E74C3C", "WT" = "#3498DB")) +
  labs(
    title = "Adgrl3 Expression: Knockout Validation",
    subtitle = sprintf("KO/WT ratio: %.2f%%", (ko_mean/wt_mean)*100),
    x = "",
    y = "CPM (Counts Per Million)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(output_dir, "adgrl3_knockout_validation.png"),
  p, width = 6, height = 5, dpi = 150
)

ggsave(
  file.path(output_dir, "adgrl3_knockout_validation.pdf"),
  p, width = 6, height = 5
)

# Save data
write.csv(df, file.path(output_dir, "adgrl3_expression_data.csv"), row.names=FALSE)

cat(sprintf("Plot saved: %s\n", file.path(output_dir, "adgrl3_knockout_validation.png")))
RSCRIPT

chmod +x "${RESULTS_DIR}/sanity_check/validate_knockout.R"

echo "[$(date)] R validation script created"

# =============================================================================
# Step 4: Run Validation
# =============================================================================
echo "[$(date)] Running validation..."

# Run featureCounts if not already done
if [[ ! -f "${RESULTS_DIR}/sanity_check/gene_counts.txt" ]]; then
    run_featurecounts
fi

# Run R script
Rscript "${RESULTS_DIR}/sanity_check/validate_knockout.R" "${PROJECT_DIR}"

echo "[$(date)] Sanity check complete!"
echo "Results: ${RESULTS_DIR}/sanity_check/"
echo ""
echo "NEXT STEPS:"
echo "1. Review adgrl3_knockout_validation.png"
echo "2. If PASS: Proceed to splicing analysis"
echo "3. If FAIL: Investigate data before continuing"
