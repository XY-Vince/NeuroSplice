#!/bin/bash
PROJECT_DIR="/Volumes/Untitled/NeuroSplice"
FASTQ_DIR="$PROJECT_DIR/data/trimmed"
ALIGNED_DIR="$PROJECT_DIR/data/aligned_star_file_test"
STAR_INDEX="$PROJECT_DIR/reference/mm10_STAR_index"
SRR="SRR7540620"
R1_PLAIN="$FASTQ_DIR/${SRR}_trimmed.fastq"

# Export PATH for STAR
export PATH=$PATH:/opt/miniconda3/envs/combio/bin:/opt/miniconda3/envs/star_env/bin

mkdir -p $ALIGNED_DIR

echo "Testing Plain File Input for $SRR..."
ls -lh "$R1_PLAIN"

# Run STAR on plain file
STAR --genomeDir $STAR_INDEX \
    --readFilesIn "$R1_PLAIN" \
    --outFileNamePrefix "$ALIGNED_DIR/${SRR}." \
    --outSAMtype BAM SortedByCoordinate \
    --runThreadN 8 \
    --outFilterMultimapNmax 10

echo "Exit Code: $?"
ls -l $ALIGNED_DIR
grep "Uniquely mapped reads" "$ALIGNED_DIR/${SRR}.Log.final.out"
