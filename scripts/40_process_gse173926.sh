#!/bin/bash
# Using specific environment for rMATS if needed, but absolute path handles binary.


# Configuration
PROJECT_DIR="/Volumes/Untitled/NeuroSplice"
RAW_DIR="$PROJECT_DIR/data/gse173926/raw"
ALIGNED_DIR="$PROJECT_DIR/data/gse173926/aligned"
RMATS_OUT="$PROJECT_DIR/results/gse173926_rmats"
INDEX="$PROJECT_DIR/reference/mm10_hisat2_index/mm10"
GTF="$PROJECT_DIR/reference/gencode.vM25.annotation.gtf"
MAP_FILE="$PROJECT_DIR/data/metadata/GSE173926_srr_map.txt"
THREADS=4

# Activate Conda Env
source /opt/miniconda3/bin/activate combio-rmats
export PATH=$PATH:/opt/miniconda3/envs/combio/bin


mkdir -p $ALIGNED_DIR
# Clean previous results to avoid cache issues
rm -rf $RMATS_OUT
mkdir -p $RMATS_OUT
RMATS_TMP="$RMATS_OUT/tmp"
mkdir -p $RMATS_TMP

echo "--- Step 1: Alignment (HISAT2) ---"
# Read Map File skipping header
tail -n +2 $MAP_FILE | while read gsm srr group genotype; do
    echo "Processing $srr ($group)..."
    
    # Input Files
    fq1="$RAW_DIR/${srr}_1.fastq.gz"
    fq2="$RAW_DIR/${srr}_2.fastq.gz"
    bam="$ALIGNED_DIR/${srr}.bam"
    
    if [ ! -f $fq1 ]; then
        echo "Error: $fq1 not found. Run download script first."
        continue
    fi
    
    if [ -f $bam ]; then
        echo "BAM exists, skipping alignment."
    else
        hisat2 -p $THREADS -x $INDEX -1 $fq1 -2 $fq2 --dta | samtools sort -@ $THREADS -o $bam
        samtools index $bam
    fi
done

echo "--- Step 2: Prepare rMATS inputs ---"
# Create b1.txt (WT) and b2.txt (Het)
rm -f $RMATS_OUT/b1.txt $RMATS_OUT/b2.txt

tail -n +2 $MAP_FILE | while read gsm srr group genotype; do
    bam="$ALIGNED_DIR/${srr}.bam"
    # Mapping: WT -> b1 (Control), Het -> b2 (Case)
    # The user asked for "MYT1L+/-" which is Het.
    # rMATS is typically Case vs Control (b1 vs b2). 
    # Usually b1=Case, b2=Control OR b1=Test, b2=Control.
    # rMATS manual: --b1 <sample1> --b2 <sample2>. sample1 vs sample2.
    # Let's define b1 = Het (Case), b2 = WT (Control).
    
    if [ "$group" == "Het" ]; then
        echo -n "$bam," >> $RMATS_OUT/b1.txt # Case
    elif [ "$group" == "WT" ]; then
        echo -n "$bam," >> $RMATS_OUT/b2.txt # Control
    fi
done

# Remove trailing commas
sed -i '' 's/,$//' $RMATS_OUT/b1.txt
sed -i '' 's/,$//' $RMATS_OUT/b2.txt

echo "Case (Het): $(cat $RMATS_OUT/b1.txt)"
echo "Control (WT): $(cat $RMATS_OUT/b2.txt)"

echo "--- Step 3: Run rMATS ---"
rmats.py --b1 $RMATS_OUT/b1.txt --b2 $RMATS_OUT/b2.txt --gtf $GTF --od $RMATS_OUT --tmp $RMATS_TMP -t paired --readLength 151 --variable-read-length --nthread $THREADS --tstat 4 --libType fr-secondstrand

echo "--- Step 4: Check Core 4 Genes ---"
# Grep for Core 4 genes in the significant results
for file in $RMATS_OUT/*.JCEC.txt; do
    echo "Checking $file for Core 4..."
    head -n 1 $file > "${file}.core4.txt"
    grep -E "Crem|Dctn1|Pdlim7|Spata5" $file >> "${file}.core4.txt"
    cat "${file}.core4.txt"
done

echo "Processing Complete."
