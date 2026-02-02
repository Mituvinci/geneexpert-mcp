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
  --output "$RESULTS_DIR/${DATASET}_parallel" \
  --organism "$ORGANISM" \
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
  --output "$RESULTS_DIR/${DATASET}_sequential" \
  --organism "$ORGANISM" \
  --sequential-chain \
  --verbose; then
  echo "✓ Multi-agent Sequential complete"
  ((SUCCESS_COUNT++))
else
  echo "✗ Multi-agent Sequential FAILED (aborted or error)"
  FAILED_SYSTEMS+=("Sequential")
  ((FAIL_COUNT++))
fi
echo ""


# System 3: Single-agent GPT-5.2
echo "--- System 3/5: Single-agent GPT-5.2 ---"
if node bin/scrna_geneexpert.js analyze "$DATASET_PATH" \
  --output "$RESULTS_DIR/${DATASET}_single_gpt" \
  --organism "$ORGANISM" \
  --single-agent gpt5.2 \
  --verbose; then
  echo "✓ Single-agent GPT-5.2 complete"
  ((SUCCESS_COUNT++))
else
  echo "✗ Single-agent GPT-5.2 FAILED (aborted or error)"
  FAILED_SYSTEMS+=("Single-GPT")
  ((FAIL_COUNT++))
fi
echo ""


# System 4: Single-agent Claude
echo "--- System 4/5: Single-agent Claude ---"
if node bin/scrna_geneexpert.js analyze "$DATASET_PATH" \
  --output "$RESULTS_DIR/${DATASET}_single_claude" \
  --organism "$ORGANISM" \
  --single-agent claude \
  --verbose; then
  echo "✓ Single-agent Claude complete"
  ((SUCCESS_COUNT++))
else
  echo "✗ Single-agent Claude FAILED (aborted or error)"
  FAILED_SYSTEMS+=("Single-Claude")
  ((FAIL_COUNT++))
fi
echo ""


# System 5: No-agent (force automation)
echo "--- System 5/5: No-agent (Force Automation) ---"
if node bin/scrna_geneexpert.js analyze "$DATASET_PATH" \
  --output "$RESULTS_DIR/${DATASET}_no_agent" \
  --organism "$ORGANISM" \
  --force-automation \
  --verbose; then
  echo "✓ No-agent complete"
  ((SUCCESS_COUNT++))
else
  echo "✗ No-agent FAILED (aborted or error)"
  FAILED_SYSTEMS+=("No-agent")
  ((FAIL_COUNT++))
fi
echo ""


echo "======================================================================"
echo "PHASE 2: CONVERT JSON → CSV (Easy Viewing)"
echo "======================================================================"
echo ""

if [ $SUCCESS_COUNT -gt 0 ]; then
  if node bin/json_to_csv.js convert --dir "$RESULTS_DIR" 2>/dev/null; then
    echo "✓ All JSON files converted to CSV"
  else
    echo "⚠ JSON→CSV conversion failed (non-critical)"
  fi
else
  echo "⚠ Skipping (no successful systems to convert)"
fi
echo ""

echo "======================================================================"
echo "PHASE 3: COST ANALYSIS"
echo "======================================================================"
echo ""

if [ $SUCCESS_COUNT -gt 0 ]; then
  if node bin/analyze_costs.js generate \
    --results "$RESULTS_DIR" \
    --output "$RESULTS_DIR/${DATASET}_cost_analysis.csv" 2>/dev/null; then
    echo "✓ Cost analysis complete"
  else
    echo "⚠ Cost analysis failed (non-critical)"
  fi
else
  echo "⚠ Skipping (no successful systems to analyze)"
fi
echo ""

echo ""
echo "======================================================================"
echo "EXPERIMENT COMPLETE: $DATASET"
echo "======================================================================"
echo ""
echo "Results saved in: $RESULTS_DIR"
echo ""
echo "System Completion Status:"
echo "  ✓ Successful: $SUCCESS_COUNT/5"
echo "  ✗ Failed:     $FAIL_COUNT/5"
if [ $FAIL_COUNT -gt 0 ]; then
  echo ""
  echo "Failed systems:"
  for system in "${FAILED_SYSTEMS[@]}"; do
    echo "  - $system"
  done
fi
echo ""
echo "Generated files (for successful systems):"
echo "  - ${DATASET}_parallel/        (Multi-agent parallel)"
echo "  - ${DATASET}_sequential/      (Multi-agent sequential)"
echo "  - ${DATASET}_single_gpt/      (Single GPT-5.2)"
echo "  - ${DATASET}_single_claude/   (Single Claude)"
echo "  - ${DATASET}_no_agent/        (No-agent baseline)"
echo "  - ${DATASET}_cost_analysis.csv (Cost breakdown)"
echo ""
echo "Each successful system output includes:"
echo "  - Seurat objects (stage1-5)"
echo "  - QC metrics, PCA results, cluster assignments"
echo "  - Marker genes (stage5)"
echo "  - Agent decision logs (JSON + CSV)"
echo ""
echo "======================================================================"
echo ""
if [ $FAIL_COUNT -eq 0 ]; then
  echo "✓ All 5 systems completed successfully!"
  exit 0
else
  echo "⚠ $FAIL_COUNT system(s) failed or were aborted. Check logs above."
  exit 1
fi
