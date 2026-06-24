#!/bin/bash
# Monitor GSE117357 Preprocessing Pipeline Progress
# Usage: bash scripts/monitor_preprocess.sh

PROJECT_DIR="/Volumes/Untitled/NeuroSplice"
TRIM_DIR="${PROJECT_DIR}/data/trimmed"
ALIGN_DIR="${PROJECT_DIR}/aligned/bam"
LOG_DIR="${PROJECT_DIR}/logs"
TOTAL_SAMPLES=25

echo "=== GSE117357 Preprocessing Pipeline Monitor ==="
echo "Time: $(date)"
echo ""

# Check if pipeline is running
FASTP_PID=$(pgrep -f "fastp.*GSE117357" | head -1)
STAR_PID=$(pgrep -f "STAR.*genomeDir" | head -1)
PIPELINE_PID=$(pgrep -f "03_preprocess.sh" | head -1)

if [ -n "$PIPELINE_PID" ]; then
    echo "✅ Pipeline running (PID: $PIPELINE_PID)"
else
    echo "❌ Pipeline NOT running"
fi

if [ -n "$FASTP_PID" ]; then
    echo "   └─ fastp active (PID: $FASTP_PID)"
fi
if [ -n "$STAR_PID" ]; then
    echo "   └─ STAR active (PID: $STAR_PID)"
fi

echo ""
echo "=== Progress ==="

# Count trimmed samples
TRIMMED=$(ls "${TRIM_DIR}"/*_1_trimmed.fastq.gz 2>/dev/null | wc -l | tr -d ' ')
echo "Trimmed: ${TRIMMED}/${TOTAL_SAMPLES}"

# Count aligned samples
ALIGNED=$(ls "${ALIGN_DIR}"/*_Aligned.sortedByCoord.out.bam 2>/dev/null | wc -l | tr -d ' ')
echo "Aligned: ${ALIGNED}/${TOTAL_SAMPLES}"

# Count indexed BAMs
INDEXED=$(ls "${ALIGN_DIR}"/*_Aligned.sortedByCoord.out.bam.bai 2>/dev/null | wc -l | tr -d ' ')
echo "Indexed: ${INDEXED}/${TOTAL_SAMPLES}"

echo ""
echo "=== Latest Log Activity ==="
LOG_FILE=$(ls -t "${LOG_DIR}"/preprocess_*.log 2>/dev/null | head -1)
if [ -n "$LOG_FILE" ]; then
    tail -5 "$LOG_FILE"
else
    echo "No log file found"
fi

echo ""
echo "=== Disk Usage ==="
du -sh "${TRIM_DIR}" "${ALIGN_DIR}" 2>/dev/null | awk '{print $2": "$1}'

echo ""
echo "💡 Monitor continuously: watch -n 30 bash scripts/monitor_preprocess.sh"
