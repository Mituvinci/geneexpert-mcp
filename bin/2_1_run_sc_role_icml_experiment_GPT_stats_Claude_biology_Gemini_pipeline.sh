#!/bin/bash

#################################################
# ICML 2026 scRNA-seq Experiment Runner
#################################################
# Automates running all 5 systems on a scRNA dataset
#
# Usage:
#   bash bin/run_scrna_icml_experiment.sh pbmc_healthy_human human
#
# Arguments:
#   $1 = Dataset name (e.g., pbmc_healthy_human)
#   $2 = Organism (human, mouse)
#################################################

# NOTE: NO 'set -e' so that if one system fails/aborts, others still run

# Track which systems succeeded
SUCCESS_COUNT=0
FAIL_COUNT=0
declare -a FAILED_SYSTEMS

# Check arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <dataset> <organism>"
    echo "Example: $0 pbmc_healthy_human human"
    exit 1
fi

DATASET=$1
ORGANISM=$2

echo ""
echo "======================================================================"
echo "ICML 2026 scRNA-seq EXPERIMENT: $DATASET"
echo "======================================================================"
echo "Organism: $ORGANISM"
echo ""

# Find dataset directory
DATA_DIR="data/scRNA_data"
DATASET_PATH=$(find $DATA_DIR -type d -name "*${DATASET}*" | head -1)
ROLE="gp_st_cl_bl_gm_pl"

if [ -z "$DATASET_PATH" ]; then
    echo "❌ Error: Dataset $DATASET not found in $DATA_DIR"
    exit 1
fi

echo "✓ Found dataset: $DATASET_PATH"
echo ""

# Output base directory
RESULTS_DIR="experiments/scrna_results"
mkdir -p $RESULTS_DIR

echo "======================================================================"
echo "PHASE 1: RUN 5 EXPERIMENTAL SYSTEMS"
echo "======================================================================"
echo ""

# System 1: Multi-agent Parallel (default)
echo "--- System 1/5: Multi-agent Parallel ---"
if node bin/scrna_geneexpert.js analyze "$DATASET_PATH" \
  --output "$RESULTS_DIR/${DATASET}_parallel_${ROLE}" \
  --organism "$ORGANISM" \
  --auto-resolve auto \
  --gpt-role "stats" \
  --claude-role "biology" \
  --gemini-role "pipeline" \
  --verbose; then
  echo "✓ Multi-agent Parallel complete"
  ((SUCCESS_COUNT++))
else
  echo "✗ Multi-agent Parallel FAILED (aborted or error)"
  FAILED_SYSTEMS+=("Parallel")
  ((FAIL_COUNT++))
fi
echo ""


# System 2: Multi-agent Sequential
echo "--- System 2/5: Multi-agent Sequential ---"
if node bin/scrna_geneexpert.js analyze "$DATASET_PATH" \
  --output "$RESULTS_DIR/${DATASET}_sequential_${ROLE}" \
  --organism "$ORGANISM" \
  --sequential-chain \
  --auto-resolve auto \
  --gpt-role "stats" \
  --claude-role "biology" \
  --gemini-role "pipeline" \
  --verbose; then
  echo "✓ Multi-agent Sequential complete"
  ((SUCCESS_COUNT++))
else
  echo "✗ Multi-agent Sequential FAILED (aborted or error)"
  FAILED_SYSTEMS+=("Sequential")
  ((FAIL_COUNT++))
fi
echo ""


if [ $FAIL_COUNT -eq 0 ]; then
  echo "✓ All 5 systems completed successfully!"
  exit 0
else
  echo "⚠ $FAIL_COUNT system(s) failed or were aborted. Check logs above."
  exit 1
fi
