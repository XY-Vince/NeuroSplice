#!/bin/bash
# ============================================================
# DEPRECATED — DO NOT USE FOR CANONICAL ANALYSIS
# ============================================================
# The canonical rMATS results were generated on the Yale HPC using:
#   GTF: Mus_musculus.GRCm38.102.gtf
#   BAMs: /gpfs/gibbs/project/neurosplice_gse117357/alignments/
#   Output: /gpfs/gibbs/project/neurosplice_gse117357/rmats_out/
#
# KNOWN DELTA vs. HPC RUN: This local script is MISSING --allow-clipping,
# which was documented in NeuroSplice.md as a required parameter for
# 76bp single-end data to rescue soft-clipped junction reads.
#
# This script also references gencode.vM25.annotation.gtf (line 12),
# which is incompatible with the Ensembl 102 BAM coordinate space.
# ============================================================
# rMATS differential splicing analysis for GSE 117357
# WT vs KO comparison across three brain tissues
set -o pipefail

# Activate conda environment
source /opt/miniconda3/bin/activate combio-rmats

# Paths
PROJECT_DIR="/Volumes/Untitled/NeuroSplice"
BAM_LISTS="${PROJECT_DIR}/data/rmats_input"
GTF="${PROJECT_DIR}/reference/gencode.vM25.annotation.gtf"
RMATS_OUT="${PROJECT_DIR}/results/rmats_gse117357"
THREADS=8
READ_LENGTH=76  # From STAR logs: Average input read length = 76

mkdir -p "${RMATS_OUT}"

# Function to run rMATS for a tissue
run_rmats_tissue() {
    local tissue=$1
    local wt_bams="${BAM_LISTS}/wt_${tissue}_bams.txt"
    local ko_bams="${BAM_LISTS}/ko_${tissue}_bams.txt"
    local output="${RMATS_OUT}/${tissue}"
    
    echo "========================================="
    echo "Running rMATS for ${tissue}"
    echo "WT BAMs: ${wt_bams}"
    echo "KO BAMs: ${ko_bams}"
    echo "Output:  ${output}"
    echo "========================================="
    
    # Check files exist
    if [[ ! -f "$wt_bams" ]] || [[ ! -f "$ko_bams" ]]; then
        echo "ERROR: BAM list files not found for ${tissue}"
        return 1
    fi
    
    mkdir -p "${output}"
    
    # Note: --cstat 0.0001 is intentionally permissive; FDR < 0.05 applied downstream in 51_final_categorization.py
    python /opt/miniconda3/envs/combio-rmats/bin/rmats.py \
        --b1 "${wt_bams}" \
        --b2 "${ko_bams}" \
        --gtf "${GTF}" \
        --od "${output}" \
        --tmp "${output}/tmp" \
        -t single \
        --libType fr-firststrand \
        --novelSS 1 \
        --readLength ${READ_LENGTH} \
        --variable-read-length \
        --nthread ${THREADS} \
        --cstat 0.0001 \
        --tstat 6 \
        2>&1 | tee "${output}/rmats.log"
    
    local status=$?
    if [[ $status -eq 0 ]]; then
        echo "✓ rMATS completed successfully for ${tissue}"
    else
        echo "✗ rMATS failed for ${tissue} with exit code ${status}"
    fi
    
    return $status
}

# Run rMATS for each tissue
echo "Starting rMATS analysis for GSE117357 (WT vs KO)"
echo "Date: $(date)"

run_rmats_tissue "hippocampus"
HIPPO_STATUS=$?

run_rmats_tissue "prefrontal"
PFC_STATUS=$?

run_rmats_tissue "striatum"
STRIATUM_STATUS=$?

# Summary
echo ""
echo "========================================="
echo "rMATS Analysis Complete"
echo "========================================="
echo "Hippocampus:        $([ $HIPPO_STATUS -eq 0 ] && echo '✓ Success' || echo '✗ Failed')"
echo "Prefrontal Cortex:  $([ $PFC_STATUS -eq 0 ] && echo '✓ Success' || echo '✗ Failed')"
echo "Striatum:           $([ $STRIATUM_STATUS -eq 0 ] && echo '✓ Success' || echo '✗ Failed')"
echo ""

# Check for significant events
for tissue in hippocampus prefrontal striatum; do
    output="${RMATS_OUT}/${tissue}"
    if [[ -f "${output}/summary.txt" ]]; then
        echo "=== ${tissue} Summary ==="
        cat "${output}/summary.txt"
        echo ""
    fi
done

echo "Results directory: ${RMATS_OUT}"
echo "Completed: $(date)"
