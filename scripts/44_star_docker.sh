#!/bin/bash
# Script: STAR Alignment via Docker (Native Linux speed on macOS)
# Uses biocontainers STAR image for faster alignment than Rosetta emulation

PROJECT_DIR="/Volumes/Untitled/NeuroSplice"
FASTQ_DIR="$PROJECT_DIR/data/trimmed"
ALIGNED_DIR="$PROJECT_DIR/data/aligned_star_reprocess"
STAR_INDEX="$PROJECT_DIR/reference/mm10_STAR_index"
GTF="$PROJECT_DIR/reference/gencode.vM25.annotation.gtf"
MAP_FILE="$PROJECT_DIR/data/metadata/sample_map.txt"
THREADS=8

# Docker image
STAR_IMAGE="quay.io/biocontainers/star:2.7.11b--h43eeafb_2"

# Activate conda for samtools
source /opt/miniconda3/bin/activate combio-rmats
export PATH=$PATH:/opt/miniconda3/envs/combio-rmats/bin:/opt/miniconda3/envs/combio/bin

mkdir -p $ALIGNED_DIR

echo "--- STAR Alignment via Docker ---"
echo "Using image: $STAR_IMAGE"

while read -r line; do
    SRR=$(echo "$line" | awk '{print $1}')
    SAMPLE_NAME=$(echo "$line" | awk '{print $2}')
    
    echo "Processing $SRR ($SAMPLE_NAME)..."
    
    R1="$FASTQ_DIR/${SRR}_trimmed.fastq.gz"
    OUT_PREFIX="$ALIGNED_DIR/${SRR}."
    
    if [[ ! -f "$R1" ]]; then
        echo "  [WARNING] Raw file $R1 not found. Skipping."
        continue
    fi
    
    # Check if BAM exists and is valid (size > 100MB)
    if [[ -f "${OUT_PREFIX}Aligned.sortedByCoord.out.bam" ]]; then
        SIZE=$(stat -f%z "${OUT_PREFIX}Aligned.sortedByCoord.out.bam")
        if (( SIZE > 100000000 )); then
            echo "  Already aligned and valid size ($SIZE bytes). Skipping."
            continue
        else
            echo "  BAM exists but too small ($SIZE bytes). Re-running."
        fi
    fi
    
    # Clean up any previous tmp directory
    rm -rf /tmp/${SRR}_STARtmp
    
    echo "  Starting STAR via Docker..."
    
    # Run STAR in Docker with volume mounts
    # Mount project directory to /data inside container
    docker run --rm \
        -v "$PROJECT_DIR:/data" \
        -v "/tmp:/tmp" \
        --memory=32g \
        "$STAR_IMAGE" \
        STAR --genomeDir /data/reference/mm10_STAR_index \
             --readFilesIn /data/data/trimmed/${SRR}_trimmed.fastq.gz \
             --readFilesCommand zcat \
             --outFileNamePrefix /data/data/aligned_star_reprocess/${SRR}. \
             --outSAMtype BAM SortedByCoordinate \
             --outTmpDir /tmp/${SRR}_STARtmp \
             --runThreadN $THREADS \
             --outSAMstrandField intronMotif \
             --outFilterMultimapNmax 10 \
             --alignSJoverhangMin 8 \
             --alignSJDBoverhangMin 1 \
             --outFilterMismatchNmax 999 \
             --outFilterMismatchNoverLmax 0.04 \
             --alignIntronMin 20 \
             --alignIntronMax 1000000 \
             --alignMatesGapMax 1000000
    
    STAR_EXIT=$?
    
    # Cleanup
    rm -rf /tmp/${SRR}_STARtmp
    
    if (( STAR_EXIT == 0 )); then
        echo "  Indexing BAM..."
        samtools index "${OUT_PREFIX}Aligned.sortedByCoord.out.bam"
        echo "  Done: $SRR"
    else
        echo "  STAR failed for $SRR. Exit code $STAR_EXIT."
    fi
done < "$MAP_FILE"

echo "--- Alignment Complete ---"
