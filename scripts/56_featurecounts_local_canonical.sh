#!/bin/bash
# ============================================================
# Canonical Local featureCounts Execution
# ============================================================
# Generates the canonical count matrix using:
#   GTF: Mus_musculus.GRCm38.102.chr.gtf (Ensembl 102 with chr prefix)
#   BAMs: Verified local STAR BAMs in data/aligned_star_reprocess/
# ============================================================
set -e

# Activate conda environment with Subread/featureCounts
source /opt/miniconda3/bin/activate combio

PROJECT_DIR="/Volumes/Untitled/NeuroSplice"
BAM_DIR="${PROJECT_DIR}/data/aligned_star_reprocess"
GTF="${PROJECT_DIR}/reference/Mus_musculus.GRCm38.102.chr.gtf"
OUT_DIR="${PROJECT_DIR}/results/featurecounts_canonical"
OUT_FILE="${OUT_DIR}/gene_counts.txt"
THREADS=8

echo "========================================="
echo "Running Local featureCounts (Canonical)"
echo "GTF: ${GTF}"
echo "BAMs: ${BAM_DIR}/*.bam"
echo "========================================="

mkdir -p "${OUT_DIR}"

# Find all 60 BAM files and pass them to featureCounts
BAM_FILES=$(ls ${BAM_DIR}/SRR*.bam)

featureCounts \
    -T ${THREADS} \
    -a "${GTF}" \
    -o "${OUT_FILE}" \
    -t exon \
    -g gene_id \
    -s 2 \
    ${BAM_FILES} \
    2>&1 | tee "${OUT_DIR}/featureCounts.log"

echo ""
echo "Done! Output matrix: ${OUT_FILE}"
