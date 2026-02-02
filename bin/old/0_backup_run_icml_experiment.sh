#!/bin/bash

#################################################
# ICML 2026 Experiment Runner
#################################################
# Automates running all 5 systems on a dataset
#
# Usage:
#   bash bin/run_icml_experiment.sh 1_GSE52778_pe_clean human untreated Dex
#
# Arguments:
#   $1 = Dataset name (e.g., GSE52778)
#   $2 = Organism (human, mouse, rat)
#   $3 = Control keyword (e.g., untreated, cont, ctrl)
#   $4 = Treatment keyword (e.g., Dex, treat, ips)
#################################################

set -e  # Exit on error

# Check arguments
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <dataset> <organism> <control_keyword> <treatment_keyword>"
    echo "Example: $0 GSE52778 human untreated Dex"
    exit 1
fi

DATASET=$1
ORGANISM=$2
CONTROL=$3
TREATMENT=$4

echo ""
echo "======================================================================"
echo "ICML 2026 EXPERIMENT: $DATASET"
echo "======================================================================"
echo "Organism: $ORGANISM"
echo "Control:  $CONTROL"
echo "Treatment: $TREATMENT"
echo ""

# Find dataset directory
DATA_DIR="data/download_new_bulk_RNA_Data"
DATASET_PATH=$(find $DATA_DIR -type d -name "*${DATASET}*" | head -1)

if [ -z "$DATASET_PATH" ]; then
    echo "❌ Error: Dataset $DATASET not found in $DATA_DIR"
    exit 1
fi

echo "✓ Found dataset: $DATASET_PATH"
echo ""

# Output base directory
RESULTS_DIR="experiments/results"
mkdir -p $RESULTS_DIR

echo "======================================================================"
echo "PHASE 1: RUN 5 EXPERIMENTAL SYSTEMS"
echo "======================================================================"
echo ""

# System 1: Multi-agent Parallel (default)
echo "--- System 1/5: Multi-agent Parallel ---"
node bin/geneexpert.js analyze "$DATASET_PATH" \
  --staged \
  --organism "$ORGANISM" \
  --comparison "$DATASET" \
  --control-keyword "$CONTROL" \
  --treatment-keyword "$TREATMENT" \
  --output "$RESULTS_DIR/${DATASET}_parallel"
echo "✓ Multi-agent Parallel complete"
echo ""

rm -f "$RESULTS_DIR/${DATASET}_parallel/stage2_alignment/bam_files/"*.{bam,vcf} 2>/dev/null



# System 2: Multi-agent Sequential
echo "--- System 2/5: Multi-agent Sequential ---"
node bin/geneexpert.js analyze "$DATASET_PATH" \
  --staged --sequential-chain \
  --organism "$ORGANISM" \
  --comparison "$DATASET" \
  --control-keyword "$CONTROL" \
  --treatment-keyword "$TREATMENT" \
  --output "$RESULTS_DIR/${DATASET}_sequential"
echo "✓ Multi-agent Sequential complete"
echo ""

# System 3: Single-agent GPT-5.2
echo "--- System 3/5: Single-agent GPT-5.2 ---"
node bin/geneexpert.js analyze "$DATASET_PATH" \
  --staged --single-agent gpt5.2 \
  --organism "$ORGANISM" \
  --comparison "$DATASET" \
  --control-keyword "$CONTROL" \
  --treatment-keyword "$TREATMENT" \
  --output "$RESULTS_DIR/${DATASET}_single_gpt"
echo "✓ Single-agent GPT-5.2 complete"
echo ""

# System 4: Single-agent Claude
echo "--- System 4/5: Single-agent Claude ---"
node bin/geneexpert.js analyze "$DATASET_PATH" \
  --staged --single-agent claude \
  --organism "$ORGANISM" \
  --comparison "$DATASET" \
  --control-keyword "$CONTROL" \
  --treatment-keyword "$TREATMENT" \
  --output "$RESULTS_DIR/${DATASET}_single_claude"
echo "✓ Single-agent Claude complete"
echo ""

# System 5: No-agent (force automation)
echo "--- System 5/5: No-agent (Force Automation) ---"
node bin/geneexpert.js analyze "$DATASET_PATH" \
  --staged --force-automation \
  --organism "$ORGANISM" \
  --comparison "$DATASET" \
  --control-keyword "$CONTROL" \
  --treatment-keyword "$TREATMENT" \
  --output "$RESULTS_DIR/${DATASET}_no_agent"
echo "✓ No-agent complete"
echo ""

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

echo "======================================================================"
echo "PHASE 4: EVALUATION (requires ground_truth.json)"
echo "======================================================================"
echo ""

if [ -f "experiments/ground_truth.json" ]; then
    node bin/evaluate.js metrics \
      --results "$RESULTS_DIR" \
      --truth experiments/ground_truth.json \
      --output "$RESULTS_DIR/${DATASET}_evaluation.csv"
    echo "✓ Evaluation complete"

    # Compare parallel vs sequential
    node bin/evaluate.js compare \
      --results "$RESULTS_DIR" \
      --truth experiments/ground_truth.json
    echo "✓ Parallel vs Sequential comparison complete"
else
    echo "⚠️  No ground_truth.json found - skipping evaluation"
    echo "   Create experiments/ground_truth.json to enable decision accuracy metrics"
fi

echo ""
echo "======================================================================"
echo "EXPERIMENT COMPLETE: $DATASET"
echo "======================================================================"
echo ""
echo "Results saved in: $RESULTS_DIR"
echo ""
echo "Generated files:"
echo "  - ${DATASET}_*_agent_decisions.json (raw logs)"
echo "  - ${DATASET}_*_metrics.csv (human-readable)"
echo "  - ${DATASET}_cost_analysis.csv (cost breakdown)"
echo "  - ${DATASET}_evaluation.csv (decision accuracy)"
echo ""
echo "======================================================================"
