#!/bin/bash
# Rename downloaded FASTQs for GSE113754 (Mouse sleep deprivation PE) (paired-end)
# Usage (run from anywhere in the repo):
#   bash data/download_new_bulk_RNA_Data/3_GSE113754_pe_clean/3_GSE113754_rename.sh

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
cd "$SCRIPT_DIR"

MAP="gsm_map_gse113754.tsv"

echo "Renaming FASTQs for GSE113754 (Mouse sleep deprivation PE)..."
echo ""

while read -r GSM LABEL; do
  [[ -z "$GSM" || -z "$LABEL" ]] && continue
  [[ "$GSM" =~ ^GSM ]] || continue

  IN_R1="${GSM}_1.fastq.gz"
  IN_R2="${GSM}_2.fastq.gz"
  OUT_R1="${LABEL}_R1_001.fastq.gz"
  OUT_R2="${LABEL}_R2_001.fastq.gz"

  if [[ -f "$IN_R1" && -f "$IN_R2" ]]; then
    echo "$GSM -> $LABEL"
    mv "$IN_R1" "$OUT_R1"
    mv "$IN_R2" "$OUT_R2"
  else
    echo "WARNING: Missing $IN_R1 or $IN_R2 for $GSM"
  fi
done < "$MAP"

echo ""
echo "Rename complete."
