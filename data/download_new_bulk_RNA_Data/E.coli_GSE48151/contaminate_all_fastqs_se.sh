#!/bin/bash
set -euo pipefail

############################################
# FIXED USER VARIABLES (DO NOT CHANGE)
############################################

HOST_DIR="/data/halimaakhter/download_new_bulk_RNA_Data/2_GSE114845_se_clean"
ECOLI_FASTQ="/data/halimaakhter/download_new_bulk_RNA_Data/7_E.coli_GSE48151/GSM1170025_LB_OD0p87_R1_001.fastq.gz"
CONTAM_PCT=70

OUT_DIR="${HOST_DIR}_CONTAM${CONTAM_PCT}"
mkdir -p "$OUT_DIR"

SEQKIT="$(dirname "$0")/seqkit"

echo "$SEQKIT"

############################################
# LOOP OVER ALL HOST FASTQS
############################################

for HOST_FASTQ in "$HOST_DIR"/*.fastq.gz; do
    BASENAME=$(basename "$HOST_FASTQ")
    OUT_FASTQ="$OUT_DIR/$BASENAME"

    echo "========================================"
    echo "Contaminating: $BASENAME"
    echo "========================================"

    TMPDIR=$(mktemp -d)

    HOST_READS=$(zcat "$HOST_FASTQ" | wc -l | awk '{print int($1/4)}')
    ECOLI_READS=$(zcat "$ECOLI_FASTQ" | wc -l | awk '{print int($1/4)}')

    ECOLI_N=$(( HOST_READS * CONTAM_PCT / 100 ))
    HOST_N=$(( HOST_READS - ECOLI_N ))

    echo "Host reads:   $HOST_READS"
    echo "E.coli reads: $ECOLI_READS"
    echo "Sampling host=$HOST_N  ecoli=$ECOLI_N"

    "$SEQKIT" sample -n "$HOST_N"  "$HOST_FASTQ"  > "$TMPDIR/host.fastq"
    "$SEQKIT" sample -n "$ECOLI_N" "$ECOLI_FASTQ" > "$TMPDIR/ecoli.fastq"

    "$SEQKIT" shuffle "$TMPDIR/host.fastq" "$TMPDIR/ecoli.fastq" \
        | gzip > "$OUT_FASTQ"

    rm -rf "$TMPDIR"

    echo "? Written: $OUT_FASTQ"
done

echo "========================================"
echo "ALL FASTQS CONTAMINATED SUCCESSFULLY"
echo "Output directory:"
echo "  $OUT_DIR"
echo "========================================"
