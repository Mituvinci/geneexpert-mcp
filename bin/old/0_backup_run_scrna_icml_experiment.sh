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

set -e  # Exit on error

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



# System 2: Multi-agent Sequential
echo "--- System 2/5: Multi-agent Sequential ---"
node bin/scrna_geneexpert.js analyze "$DATASET_PATH" \
  --output "$RESULTS_DIR/${DATASET}_sequential" \
  --organism "$ORGANISM" \
  --sequential-chain \
  --verbose
echo "✓ Multi-agent Sequential complete"
echo ""

# rm -rf "$RESULTS_DIR/${DATASET}_sequential/stage"*"/"*.rds 2>/dev/null


# System 3: Single-agent GPT-5.2
echo "--- System 3/5: Single-agent GPT-5.2 ---"
node bin/scrna_geneexpert.js analyze "$DATASET_PATH" \
  --output "$RESULTS_DIR/${DATASET}_single_gpt" \
  --organism "$ORGANISM" \
  --single-agent gpt5.2 \
  --verbose
echo "✓ Single-agent GPT-5.2 complete"
echo ""

# rm -rf "$RESULTS_DIR/${DATASET}_single_gpt/stage"*"/"*.rds 2>/dev/null


# System 4: Single-agent Claude
echo "--- System 4/5: Single-agent Claude ---"
node bin/scrna_geneexpert.js analyze "$DATASET_PATH" \
  --output "$RESULTS_DIR/${DATASET}_single_claude" \
  --organism "$ORGANISM" \
  --single-agent claude \
  --verbose
echo "✓ Single-agent Claude complete"
echo ""

# rm -rf "$RESULTS_DIR/${DATASET}_single_claude/stage"*"/"*.rds 2>/dev/null


# System 5: No-agent (force automation)
echo "--- System 5/5: No-agent (Force Automation) ---"
node bin/scrna_geneexpert.js analyze "$DATASET_PATH" \
  --output "$RESULTS_DIR/${DATASET}_no_agent" \
  --organism "$ORGANISM" \
  --force-automation \
  --verbose
echo "✓ No-agent complete"
echo ""

# rm -rf "$RESULTS_DIR/${DATASET}_no_agent/stage"*"/"*.rds 2>/dev/null


echo "======================================================================"
echo "PHASE 2: CONVERT JSON → CSV (Easy Viewing)"
echo "======================================================================"
echo ""

node bin/json_to_csv.js convert --dir "$RESULTS_DIR"
echo "✓ All JSON files converted to CSV"
echo ""

echo "======================================================================"
echo "PHASE 3: COST ANALYSIS"
echo "======================================================================"
echo ""

node bin/analyze_costs.js generate \
  --results "$RESULTS_DIR" \
  --output "$RESULTS_DIR/${DATASET}_cost_analysis.csv"
echo "✓ Cost analysis complete"
echo ""

echo ""
echo "======================================================================"
echo "EXPERIMENT COMPLETE: $DATASET"
echo "======================================================================"
echo ""
echo "Results saved in: $RESULTS_DIR"
echo ""
echo "Generated files:"
echo "  - ${DATASET}_parallel/        (Multi-agent parallel)"
echo "  - ${DATASET}_sequential/      (Multi-agent sequential)"
echo "  - ${DATASET}_single_gpt/      (Single GPT-5.2)"
echo "  - ${DATASET}_single_claude/   (Single Claude)"
echo "  - ${DATASET}_no_agent/        (No-agent baseline)"
echo "  - ${DATASET}_cost_analysis.csv (Cost breakdown)"
echo ""
echo "Each system output includes:"
echo "  - Seurat objects (stage1-5)"
echo "  - QC metrics, PCA results, cluster assignments"
echo "  - Marker genes (stage5)"
echo "  - Agent decision logs (JSON + CSV)"
echo ""
echo "======================================================================"
