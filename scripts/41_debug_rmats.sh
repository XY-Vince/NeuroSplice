#!/bin/bash
PROJECT_DIR="/Volumes/Untitled/NeuroSplice"
RMATS_OUT="$PROJECT_DIR/results/rmats_debug"
GTF="$PROJECT_DIR/reference/gencode.vM25.annotation.gtf"
BAM_WT="/Volumes/Untitled/NeuroSplice/data/gse173926/aligned/SRR14494965.bam"
BAM_HET="/Volumes/Untitled/NeuroSplice/data/gse173926/aligned/SRR14494966.bam"
RMATS_BIN="rmats.py" # Use from env

# Activate Conda Env
source /opt/miniconda3/bin/activate combio-rmats

rm -rf $RMATS_OUT
mkdir -p $RMATS_OUT/tmp

echo "$BAM_HET" > $RMATS_OUT/b1.txt
echo "$BAM_WT" > $RMATS_OUT/b2.txt

echo "Running rMATS Debug (fr-firststrand)..."
# Using fr-firststrand (common for dUTP methods)
rmats.py --b1 $RMATS_OUT/b1.txt --b2 $RMATS_OUT/b2.txt --gtf $GTF --od $RMATS_OUT --tmp $RMATS_OUT/tmp -t paired --readLength 151 --variable-read-length --nthread 2 --tstat 1 --libType fr-firststrand

echo "Counts:"
cat $RMATS_OUT/summary.txt
