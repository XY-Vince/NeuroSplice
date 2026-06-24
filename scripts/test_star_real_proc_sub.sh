#!/bin/bash
PROJECT_DIR="/Volumes/Untitled/NeuroSplice"
FASTQ_DIR="$PROJECT_DIR/data/trimmed"
ALIGNED_DIR="$PROJECT_DIR/data/aligned_star_sub_test"
STAR_INDEX="$PROJECT_DIR/reference/mm10_STAR_index"
SRR="SRR7540620"
R1="$FASTQ_DIR/${SRR}_trimmed.fastq.gz"

# Export PATH for STAR
export PATH=$PATH:/opt/miniconda3/envs/combio/bin:/opt/miniconda3/envs/star_env/bin

mkdir -p $ALIGNED_DIR

echo "Testing Process Substitution <(...) for $SRR..."
echo "Input: $R1"

# Using Process Substitution
STAR --genomeDir $STAR_INDEX \
    --readFilesIn <(/usr/bin/gzip -dc "$R1") \
    --outFileNamePrefix "$ALIGNED_DIR/${SRR}." \
    --outSAMtype BAM SortedByCoordinate \
    --runThreadN 8 \
    --outFilterMultimapNmax 10

echo "Exit Code: $?"
ls -l $ALIGNED_DIR
grep "Uniquely mapped reads" "$ALIGNED_DIR/${SRR}.Log.final.out"
