#!/bin/bash
PROJECT_DIR="/Volumes/Untitled/NeuroSplice"
FASTQ_DIR="$PROJECT_DIR/data/trimmed"
ALIGNED_DIR="$PROJECT_DIR/data/aligned_star_test_proc"
STAR_INDEX="$PROJECT_DIR/reference/mm10_STAR_index"
SRR="SRR7540620"
R1="$FASTQ_DIR/${SRR}_trimmed.fastq.gz"

# Export PATH for STAR
export PATH=$PATH:/opt/miniconda3/envs/combio/bin:/opt/miniconda3/envs/star_env/bin

mkdir -p $ALIGNED_DIR
FILTER_OUT="$ALIGNED_DIR/${SRR}.Log.out"

echo "Testing Process Substitution for $SRR..."
echo "Input: $R1"

# Using Process Substitution <(...)
# Note: STAR needs to see a 'file', so we use --readFilesIn /dev/stdin and pipe it?
# OR STAR accepts standard input?
# STAR manual says: --readFilesIn Stdin
# So: gzip -dc file | STAR --readFilesIn Stdin ...

echo "Command: gzip -dc $R1 | STAR --readFilesIn /dev/stdin ..."

/usr/bin/gzip -dc "$R1" | STAR --genomeDir $STAR_INDEX \
    --readFilesIn /dev/stdin \
    --outFileNamePrefix "$ALIGNED_DIR/${SRR}." \
    --outSAMtype BAM SortedByCoordinate \
    --runThreadN 8 \
    --outFilterMultimapNmax 10

echo "Exit Code: $?"
ls -l $ALIGNED_DIR
grep "Uniquely mapped reads" "$ALIGNED_DIR/${SRR}.Log.final.out"
