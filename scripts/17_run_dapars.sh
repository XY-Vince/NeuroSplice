#!/bin/bash
set -euo pipefail

# Activate Environment
source /opt/miniconda3/etc/profile.d/conda.sh
conda activate salmon_env
export PATH="/opt/miniconda3/envs/bedtools_env/bin:$PATH"

# Directories
PROJECT_DIR="/Volumes/Untitled/NeuroSplice"
BAM_DIR="${PROJECT_DIR}/aligned/bam"
RESULTS_DIR="${PROJECT_DIR}/results/dapars"
CONFIG_DIR="${RESULTS_DIR}/configs"
BEDGRAPH_DIR="${RESULTS_DIR}/bedgraph"
ANNOTATION="${PROJECT_DIR}/reference/dapars_annotation.bed"

mkdir -p "${RESULTS_DIR}" "${CONFIG_DIR}" "${BEDGRAPH_DIR}"

# Input Lists (Pairwise Comparison per Region)
REGIONS=("Hippocampus" "Prefrontal_Cortex" "Striatum")

# Loop Regions
for REGION in "${REGIONS[@]}"; do
    echo "Processing Region: ${REGION}"
    R_LOWER=$(echo "$REGION" | tr '[:upper:]' '[:lower:]')
    
    # We need to get the BAM files for this region
    # We can reuse metadata/wt_hippocampus.txt logic, but we need comma separated list for DaPars config
    # DaPars Config:
    # Group1_Tophat_BAMs = sample1.bedgraph,sample2.bedgraph
    # Group2_Tophat_BAMs = sample3.bedgraph,sample4.bedgraph
    
    # Original lists are BAM paths. We need to convert to BedGraph paths.
    WT_LIST_FILE="${PROJECT_DIR}/data/metadata/wt_${R_LOWER}.txt"
    KO_LIST_FILE="${PROJECT_DIR}/data/metadata/ko_${R_LOWER}.txt"
    
    if [ ! -f "${WT_LIST_FILE}" ]; then
        # Handle Prefrontal-Cortex vs Prefrontal_Cortex naming issue
        if [ "$REGION" == "Prefrontal_Cortex" ]; then
             WT_LIST_FILE="${PROJECT_DIR}/data/metadata/wt_prefrontal_cortex.txt"
             KO_LIST_FILE="${PROJECT_DIR}/data/metadata/ko_prefrontal_cortex.txt"
        fi
    fi

    # Convert BAMs to BedGraph (Parallelize this?)
    echo "  Generating BedGraph files..."
    
    process_bam() {
        BAM_PATH=$1
        SAMPLE=$(basename "$BAM_PATH" _Aligned.sortedByCoord.out.bam)
        BG_PATH="${BEDGRAPH_DIR}/${SAMPLE}.bedgraph"
        
        if [ ! -f "${BG_PATH}" ]; then
            echo "    Converting ${SAMPLE}..."
            bedtools genomecov -ibam "${BAM_PATH}" -bg -split > "${BG_PATH}"
        fi
        echo "${BG_PATH}"
    }
    
    # Extract BAM paths and convert
    # Reading comma separated list
    WT_BAMS=$(cat "${WT_LIST_FILE}" | tr ',' ' ')
    KO_BAMS=$(cat "${KO_LIST_FILE}" | tr ',' ' ')
    
    WT_BEDGRAPHS=""
    KO_BEDGRAPHS=""
    
    # We can iterate and run background jobs
    job_count=0
    for BAM in $WT_BAMS $KO_BAMS; do
        process_bam "$BAM" &
        ((job_count++))
        if (( job_count >= 8 )); then wait; job_count=0; fi
    done
    wait
    
    # Collect paths
    for BAM in $WT_BAMS; do
        SAMPLE=$(basename "$BAM" _Aligned.sortedByCoord.out.bam)
        WT_BEDGRAPHS="${WT_BEDGRAPHS}${BEDGRAPH_DIR}/${SAMPLE}.bedgraph,"
    done
    WT_BEDGRAPHS=${WT_BEDGRAPHS%,} # Remove trailing comma
    
    for BAM in $KO_BAMS; do
        SAMPLE=$(basename "$BAM" _Aligned.sortedByCoord.out.bam)
        KO_BEDGRAPHS="${KO_BEDGRAPHS}${BEDGRAPH_DIR}/${SAMPLE}.bedgraph,"
    done
    KO_BEDGRAPHS=${KO_BEDGRAPHS%,}

    # Generate Config
    CONFIG_FILE="${CONFIG_DIR}/dapars_config_${R_LOWER}.txt"
    echo "Annotated_3UTR=${ANNOTATION}" > "${CONFIG_FILE}"
    echo "Group1_Tophat_BAMs=${WT_BEDGRAPHS}" >> "${CONFIG_FILE}"
    echo "Group2_Tophat_BAMs=${KO_BEDGRAPHS}" >> "${CONFIG_FILE}"
    echo "Output_directory=${RESULTS_DIR}/${R_LOWER}" >> "${CONFIG_FILE}"
    echo "Output_result_file=dapars_result" >> "${CONFIG_FILE}"
    echo "Num_least_in_group1=5" >> "${CONFIG_FILE}" # Default coverage filters
    echo "Num_least_in_group2=5" >> "${CONFIG_FILE}"
    echo "Coverage_cutoff=30" >> "${CONFIG_FILE}"
    echo "FDR_cutoff=0.05" >> "${CONFIG_FILE}"
    echo "PDUI_cutoff=0.2" >> "${CONFIG_FILE}"
    echo "Fold_change_cutoff=0.59" >> "${CONFIG_FILE}"
    
    # Run DaPars
    echo "  Running DaPars for ${REGION}..."
    python tools/DaPars/src/DaPars_main.py "${CONFIG_FILE}" > "${RESULTS_DIR}/dapars_${R_LOWER}.log" 2>&1
    
    echo "  Finished ${REGION}."
done
