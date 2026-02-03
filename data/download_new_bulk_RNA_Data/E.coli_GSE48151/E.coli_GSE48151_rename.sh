#!/bin/bash
# Rename downloaded FASTQs for E. coli GSE48151 (contamination control) (single-end)
# Usage (run from anywhere in the repo):
#   bash data/download_new_bulk_RNA_Data/E.coli_GSE48151/E.coli_GSE48151_rename.sh

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
cd "$SCRIPT_DIR"

MAP="gsm_map_gse48151.tsv"

echo "Renaming FASTQs for E. coli GSE48151 (contamination control)..."
echo ""

while IFS=$'\t' read -r GSM LABEL; do
  [[ -z "$GSM" || -z "$LABEL" ]] && continue
  [[ "$GSM" =~ ^GSM ]] || continue

  IN="${GSM}.fastq.gz"
  OUT="${LABEL}_R1_001.fastq.gz"

  if [[ -f "$IN" ]]; then
    echo "$GSM -> $LABEL"
    mv "$IN" "$OUT"
  else
    echo "WARNING: $IN not found"
  fi
done < "$MAP"

echo ""
echo "Rename complete."
