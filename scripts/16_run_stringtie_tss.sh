#!/bin/bash
set -euo pipefail

# Activate Environment (StringTie should be in 'combio' or we install in specific env)
# If installed in base or combio, just use conda activate
source /opt/miniconda3/etc/profile.d/conda.sh
conda activate stringtie_env

# Directories
PROJECT_DIR="/Volumes/Untitled/NeuroSplice"
BAM_DIR="${PROJECT_DIR}/aligned/bam"
RESULTS_DIR="${PROJECT_DIR}/results/stringtie_tss"
CAGE_BED="${PROJECT_DIR}/data/cage/mm10_fair+new_CAGE_peaks_phase1and2.bed.gz" # Updated path
GTF="${PROJECT_DIR}/reference/gencode.vM25.annotation.gtf"

mkdir -p "${RESULTS_DIR}"

echo "--- Starting StringTie Assembly ---"

# Iterate all BAMs
# We use the metadata lists to get the BAMs
ALL_BAMS=$(cat ${PROJECT_DIR}/data/metadata/wt_*.txt ${PROJECT_DIR}/data/metadata/ko_*.txt | tr ',' '\n' | sort | uniq)

# Parallel Execution Config
MAX_JOBS=8
job_count=0

for BAM in ${ALL_BAMS}; do
    SAMPLE=$(basename "${BAM}" _Aligned.sortedByCoord.out.bam)
    OUT_GTF="${RESULTS_DIR}/${SAMPLE}.gtf"
    
    if [ -f "${OUT_GTF}" ]; then
        echo "Skipping ${SAMPLE} (Already exists)"
        continue
    fi
    
    echo "Assembling ${SAMPLE} [Job ${job_count}]..."
    (
        stringtie "${BAM}" \
            -G "${GTF}" \
            -o "${OUT_GTF}" \
            -p 2 \
            -v > "${RESULTS_DIR}/${SAMPLE}.log" 2>&1
    ) &
    
    ((job_count++))
    if (( job_count >= MAX_JOBS )); then
        wait
        job_count=0
    fi
done
wait # Wait for remaining jobs

# Since MacOS bash is old, let's use a simpler batch approach or check bash version.
# Assuming MacOS default bash 3.2.
# We will use the 'wait' at intervals pattern.


echo "--- StringTie Complete. Merging... ---"

# Merge all GTFs
ls "${RESULTS_DIR}"/*.gtf > "${RESULTS_DIR}/mergelist.txt"
stringtie --merge -G "${GTF}" -o "${RESULTS_DIR}/stringtie_merged.gtf" "${RESULTS_DIR}/mergelist.txt"

echo "--- Intersecting with CAGE ---"
# Extract TSS from Merged GTF (transcripts)
# Method: get 5' end of transcripts
# Strand +: Start
# Strand -: End
awk '$3=="transcript" {if($7=="+") print $1"\t"$4-1"\t"$4"\t"$14"\t.\t+"; else print $1"\t"$5-1"\t"$5"\t"$14"\t.\t-"}' "${RESULTS_DIR}/stringtie_merged.gtf" | tr -d '";' > "${RESULTS_DIR}/predicted_tss.bed"

# Intersect with CAGE
# CAGE is .bed.gz
zcat "${CAGE_BED}" > "${PROJECT_DIR}/data/cage/cage_peaks.bed"
bedtools intersect -a "${RESULTS_DIR}/predicted_tss.bed" -b "${PROJECT_DIR}/data/cage/cage_peaks.bed" -wa -wb > "${RESULTS_DIR}/tss_cage_validated.txt"

echo "Validation Complete. Results in ${RESULTS_DIR}/tss_cage_validated.txt"
