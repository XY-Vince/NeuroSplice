#!/bin/bash
# rMATS analysis for GSE173926 (MYT1L+/- Het vs WT)
# Cross-dataset validation for Category A genes

# Activate conda environment
source /opt/miniconda3/bin/activate combio-rmats

# Paths
PROJECT_DIR="/Volumes/Untitled/NeuroSplice"
ALIGNED_DIR="${PROJECT_DIR}/data/gse173926/aligned"
RMATS_INPUT="${PROJECT_DIR}/data/gse173926/rmats_input"
GTF="${PROJECT_DIR}/reference/gencode.vM25.annotation.gtf"
RMATS_OUT="${PROJECT_DIR}/results/gse173926_rmats_rerun"
THREADS=8
READ_LENGTH=76

mkdir -p "${RMATS_INPUT}"
mkdir -p "${RMATS_OUT}"

echo "Creating BAM lists..."

# Create temporary files with paths
TMPWT=$(mktemp)
TMPHET=$(mktemp)

# WT samples
echo "${ALIGNED_DIR}/SRR14494965.bam" >> "$TMPWT"
echo "${ALIGNED_DIR}/SRR14494966.bam" >> "$TMPWT"
echo "${ALIGNED_DIR}/SRR14494968.bam" >> "$TMPWT"
echo "${ALIGNED_DIR}/SRR14494970.bam" >> "$TMPWT"
echo "${ALIGNED_DIR}/SRR14494971.bam" >> "$TMPWT"
echo "${ALIGNED_DIR}/SRR14494974.bam" >> "$TMPWT"

# Het samples
echo "${ALIGNED_DIR}/SRR14494967.bam" >> "$TMPHET"
echo "${ALIGNED_DIR}/SRR14494969.bam" >> "$TMPHET"
echo "${ALIGNED_DIR}/SRR14494972.bam" >> "$TMPHET"
echo "${ALIGNED_DIR}/SRR14494973.bam" >> "$TMPHET"
echo "${ALIGNED_DIR}/SRR14494975.bam" >> "$TMPHET"
echo "${ALIGNED_DIR}/SRR14494976.bam" >> "$TMPHET"

# Convert to comma-separated on single line (rMATS format)
cat "$TMPWT" | tr '\n' ',' | sed 's/,$//' > "${RMATS_INPUT}/wt_bams.txt"
cat "$TMPHET" | tr '\n' ',' | sed 's/,$//' > "${RMATS_INPUT}/het_bams.txt"

rm -f "$TMPWT" "$TMPHET"

echo "BAM lists created (comma-separated format for rMATS):"
echo "  WT:  6 samples"
echo "  Het: 6 samples"

echo ""
echo "========================================="
echo "Running rMATS for GSE173926 (MYT1L Het vs WT)"
echo "========================================="
echo "WT BAMs:  ${RMATS_INPUT}/wt_bams.txt"
echo "Het BAMs: ${RMATS_INPUT}/het_bams.txt"
echo "Output:   ${RMATS_OUT}"
echo "========================================="

# Run rMATS
python /opt/miniconda3/envs/combio-rmats/bin/rmats.py \
    --b1 "${RMATS_INPUT}/wt_bams.txt" \
    --b2 "${RMATS_INPUT}/het_bams.txt" \
    --gtf "${GTF}" \
    --od "${RMATS_OUT}" \
    --tmp "${RMATS_OUT}/tmp" \
    -t single \
    --readLength ${READ_LENGTH} \
    --variable-read-length \
    --nthread ${THREADS} \
    --cstat 0.0001 \
    --tstat 6 \
    2>&1 | tee "${RMATS_OUT}/rmats.log"

STATUS=$?

echo ""
echo "========================================="
if [ $STATUS -eq 0 ]; then
    echo "✓ rMATS completed successfully for GSE173926"
    
    # Quick summary
    echo ""
    echo "Event Summary:"
    for event in SE A5SS A3SS MXE RI; do
        file="${RMATS_OUT}/${event}.MATS.JCEC.txt"
        if [ -f "$file" ]; then
            count=$(tail -n +2 "$file" | wc -l)
            echo "  ${event}: ${count} events"
        fi
    done
else
    echo "✗ rMATS failed with exit code ${STATUS}"
fi
echo "========================================="
echo "Output: ${RMATS_OUT}"
echo "Completed: $(date)"
