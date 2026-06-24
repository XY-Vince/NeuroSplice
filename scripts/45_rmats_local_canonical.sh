#!/bin/bash
# ============================================================
# Canonical Local rMATS Execution
# ============================================================
# Generates canonical alternative splicing analysis using:
#   GTF:  Mus_musculus.GRCm38.102.chr.gtf (Ensembl 102 with chr prefix)
#   BAMs: Verified local STAR BAMs in data/aligned_star_reprocess/
#
# DIFFERENCES vs. the deprecated 45_rmats_gse117357.sh:
#   1. Uses Ensembl 102 GTF instead of gencode.vM25
#   2. Adds --allow-clipping to rescue ~11% of soft-clipped reads
#      (critical for 76bp single-end data)
#   3. Output goes to results/rmats_gse117357_canonical/ to avoid
#      overwriting the previous run
#
# Documented canonical parameters from NeuroSplice.md:
#   --readLength 76  --t single  --novelSS 1
#   --allow-clipping  --libType fr-firststrand
# ============================================================
# NOTE: Do NOT use 'set -e' here. Each tissue takes 6-12 hours.
# If one tissue fails, we must still attempt the remaining tissues.
# 'pipefail' alone ensures $? captures the rMATS exit code through | tee.
set -o pipefail

# Activate conda environment with rMATS
source /opt/miniconda3/bin/activate combio-rmats

# Paths
PROJECT_DIR="/Volumes/Untitled/NeuroSplice"
BAM_LISTS="${PROJECT_DIR}/data/rmats_input"
GTF="${PROJECT_DIR}/reference/Mus_musculus.GRCm38.102.chr.gtf"
RMATS_OUT="${PROJECT_DIR}/results/rmats_gse117357_canonical"
THREADS=8
READ_LENGTH=76  # From STAR logs: Average input read length = 76

# Pre-flight checks
if [[ ! -f "${GTF}" ]]; then
    echo "FATAL: GTF not found at ${GTF}"
    echo "Run Step 1 (download & transform Ensembl 102 GTF) first."
    exit 1
fi

mkdir -p "${RMATS_OUT}"

# Function to run rMATS for a single tissue
run_rmats_tissue() {
    local tissue=$1
    local wt_bams="${BAM_LISTS}/wt_${tissue}_bams.txt"
    local ko_bams="${BAM_LISTS}/ko_${tissue}_bams.txt"
    local output="${RMATS_OUT}/${tissue}"

    echo ""
    echo "========================================="
    echo "Running rMATS for ${tissue}"
    echo "WT BAMs: ${wt_bams}"
    echo "KO BAMs: ${ko_bams}"
    echo "GTF:     ${GTF}"
    echo "Output:  ${output}"
    echo "Started: $(date)"
    echo "========================================="

    # Check BAM list files exist
    if [[ ! -f "$wt_bams" ]] || [[ ! -f "$ko_bams" ]]; then
        echo "ERROR: BAM list files not found for ${tissue}"
        return 1
    fi

    # Clean previous tmp directory if it exists (rMATS will error otherwise)
    if [[ -d "${output}/tmp" ]]; then
        echo "Cleaning previous tmp directory for ${tissue}..."
        rm -rf "${output}/tmp"
    fi

    mkdir -p "${output}"

    # Note: --cstat 0.0001 is intentionally permissive;
    #       FDR < 0.05 is applied downstream in 51_final_categorization.py
    python /opt/miniconda3/envs/combio-rmats/bin/rmats.py \
        --b1 "${wt_bams}" \
        --b2 "${ko_bams}" \
        --gtf "${GTF}" \
        --od "${output}" \
        --tmp "${output}/tmp" \
        -t single \
        --libType fr-firststrand \
        --novelSS \
        --allow-clipping \
        --readLength ${READ_LENGTH} \
        --variable-read-length \
        --nthread ${THREADS} \
        --cstat 0.0001 \
        --tstat 6 \
        2>&1 | tee "${output}/rmats.log"

    local status=$?
    if [[ $status -eq 0 ]]; then
        echo "✓ rMATS completed successfully for ${tissue} at $(date)"
    else
        echo "✗ rMATS failed for ${tissue} with exit code ${status}"
    fi

    return $status
}

# Run rMATS for each tissue sequentially
echo "Starting Canonical rMATS analysis for GSE117357 (WT vs KO)"
echo "Date: $(date)"
echo "GTF: ${GTF}"
echo ""

run_rmats_tissue "hippocampus"
HIPPO_STATUS=$?

run_rmats_tissue "prefrontal"
PFC_STATUS=$?

run_rmats_tissue "striatum"
STRIATUM_STATUS=$?

# Summary
echo ""
echo "========================================="
echo "Canonical rMATS Analysis Complete"
echo "========================================="
echo "Hippocampus:        $([ $HIPPO_STATUS -eq 0 ] && echo '✓ Success' || echo '✗ Failed')"
echo "Prefrontal Cortex:  $([ $PFC_STATUS -eq 0 ] && echo '✓ Success' || echo '✗ Failed')"
echo "Striatum:           $([ $STRIATUM_STATUS -eq 0 ] && echo '✓ Success' || echo '✗ Failed')"
echo ""

# Show summary files if available
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
