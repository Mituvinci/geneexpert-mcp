#!/bin/bash
# Download FASTQ files for GSE141496 (HeLa batch effect)
# Usage (run from anywhere in the repo):
#   bash data/download_new_bulk_RNA_Data/4_GSE141496_batch_effect/4_GSE141496_download.sh

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
cd "$SCRIPT_DIR"

REPO_ROOT="$SCRIPT_DIR/../../.."
DOWNLOAD_SCRIPT="$REPO_ROOT/bio_informatics/scripts/download_SRRs_from_GSM.sh"

if [[ ! -f "$DOWNLOAD_SCRIPT" ]]; then
  echo "Error: download_SRRs_from_GSM.sh not found at $DOWNLOAD_SCRIPT"
  exit 1
fi

TSV_FILE=$(ls gsm_map*.tsv 2>/dev/null | head -1)
if [[ -z "$TSV_FILE" ]]; then
  echo "Error: No gsm_map*.tsv found in $SCRIPT_DIR"
  exit 1
fi

echo "Dataset:  GSE141496 (HeLa batch effect)"
echo "Map file: $TSV_FILE"
echo ""

GSM_LIST=$(cut -f1 "$TSV_FILE" | grep "^GSM" | tr '\n' ' ')
if [[ -z "$GSM_LIST" ]]; then
  echo "Error: No GSM IDs found in $TSV_FILE"
  exit 1
fi

echo "GSM IDs: $GSM_LIST"
echo ""
echo "Starting download..."
echo ""

bash "$DOWNLOAD_SCRIPT" $GSM_LIST

echo ""
echo "Download complete. Next, run:"
echo "  bash data/download_new_bulk_RNA_Data/4_GSE141496_batch_effect/4_GSE141496_rename.sh"
