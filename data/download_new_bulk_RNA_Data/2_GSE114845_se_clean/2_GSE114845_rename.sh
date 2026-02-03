#!/bin/bash
# Rename downloaded FASTQs for GSE114845 (Mouse sleep deprivation SE) (single-end)
# Usage (run from anywhere in the repo):
#   bash data/download_new_bulk_RNA_Data/2_GSE114845_se_clean/2_GSE114845_rename.sh

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
cd "$SCRIPT_DIR"

MAP="gsm_map_gse114845.tsv"

echo "Renaming FASTQs for GSE114845 (Mouse sleep deprivation SE)..."
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
