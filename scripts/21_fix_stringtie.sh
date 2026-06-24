#!/bin/bash
set -euo pipefail

# Activate Environment
source /opt/miniconda3/etc/profile.d/conda.sh
conda activate stringtie_env
export PATH="/opt/miniconda3/envs/bedtools_env/bin:$PATH"

PROJECT_DIR="/Volumes/Untitled/NeuroSplice"
RESULTS_DIR="${PROJECT_DIR}/results/stringtie_tss"
GTF="${PROJECT_DIR}/reference/gencode.vM25.annotation.gtf"
CAGE_BED="${PROJECT_DIR}/data/cage/cage_peaks.bed"

echo "--- 1. Merging StringTie GTFs ---"
if [ ! -s "${RESULTS_DIR}/mergelist.txt" ]; then
    ls "${RESULTS_DIR}/SRR*.gtf" > "${RESULTS_DIR}/mergelist.txt"
fi

stringtie --merge \
    -G "${GTF}" \
    -o "${RESULTS_DIR}/stringtie_merged.gtf" \
    "${RESULTS_DIR}/mergelist.txt"

echo "Merge Complete. Size: $(du -h "${RESULTS_DIR}/stringtie_merged.gtf" | cut -f1)"

echo "--- 2. Extracting Predicted TSS ---"
awk '$3=="transcript" {if($7=="+") print $1"\t"$4-1"\t"$4"\t"$14"\t.\t+"; else print $1"\t"$5-1"\t"$5"\t"$14"\t.\t-"}' "${RESULTS_DIR}/stringtie_merged.gtf" | tr -d '";' > "${RESULTS_DIR}/predicted_tss.bed"

echo "--- 3. Intersecting with CAGE ---"
# Check if CAGE exists
if [ ! -f "${CAGE_BED}" ]; then
    echo "Unzipping CAGE..."
    zcat "${PROJECT_DIR}/data/cage/mm10_fair+new_CAGE_peaks_phase1and2.bed.gz" > "${CAGE_BED}"
fi

bedtools intersect \
    -a "${RESULTS_DIR}/predicted_tss.bed" \
    -b "${CAGE_BED}" \
    -wa -wb > "${RESULTS_DIR}/tss_cage_validated.txt"

echo "Validation Complete. Output lines:"
wc -l "${RESULTS_DIR}/tss_cage_validated.txt"
