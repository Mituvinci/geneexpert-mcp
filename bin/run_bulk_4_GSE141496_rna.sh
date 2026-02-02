#!/bin/bash

# NOTE: This script will continue running even if individual experiments fail
# Each experiment is independent and failures won't stop the batch

# Track success/failure
SUCCESS_COUNT=0
FAIL_COUNT=0

# Check arguments
if [ $# -lt 4 ]; then
  echo "Usage: bash bin/run_bulk_rna.sh <dataset> <organism> <control_kw> <treatment_kw>"
  echo ""
  echo "Example:"
  echo "  bash bin/run_bulk_rna.sh 1_GSE52778_pe_clean human untreated Dex"
  exit 1
fi


DATASET=$1
ORGANISM=$2
CONTROL=$3
TREATMENT=$4


echo "======================================================================"
echo "RUNNING ALL BULK RNA-SEQ ABLATION EXPERIMENTS"
echo "======================================================================"
echo "Dataset:   $DATASET"
echo "Organism:  $ORGANISM"
echo "Control:   $CONTROL"
echo "Treatment: $TREATMENT"
echo ""

echo "[1/10] No-Agent Baseline..."
if bash bin/0_1_run_bulk_no_agent_icml_experiment_optimized_gpt.sh "$DATASET" "$ORGANISM" "$CONTROL" "$TREATMENT"; then
  ((SUCCESS_COUNT++))
else
  echo "⚠️  Experiment 1 failed, continuing..."
  ((FAIL_COUNT++))
fi

echo "[2/10] Single-Agent GPT-5.2..."
if bash bin/0_2_run_bulk_single_agent_gpt_icml_experiment_optimized.sh "$DATASET" "$ORGANISM" "$CONTROL" "$TREATMENT"; then
  ((SUCCESS_COUNT++))
else
  echo "⚠️  Experiment 2 failed, continuing..."
  ((FAIL_COUNT++))
fi

echo "[3/10] Single-Agent Claude..."
if bash bin/0_3_run_bulk_single_agent_claude_icml_experiment_optimized.sh "$DATASET" "$ORGANISM" "$CONTROL" "$TREATMENT"; then
  ((SUCCESS_COUNT++))
else
  echo "⚠️  Experiment 3 failed, continuing..."
  ((FAIL_COUNT++))
fi

echo "[4/10] Single-Agent Gemini..."
if bash bin/0_4_run_bulk_single_agent_gemini_icml_experiment_optimized.sh "$DATASET" "$ORGANISM" "$CONTROL" "$TREATMENT"; then
  ((SUCCESS_COUNT++))
else
  echo "⚠️  Experiment 4 failed, continuing..."
  ((FAIL_COUNT++))
fi


echo "[10/10] Default Configuration..."
if bash bin/1_6_run_bulk_role_icml_experiment_optimized_default_parallel_sequential.sh "$DATASET" "$ORGANISM" "$CONTROL" "$TREATMENT"; then
  ((SUCCESS_COUNT++))
else
  echo "⚠️  Experiment 10 failed, continuing..."
  ((FAIL_COUNT++))
fi

echo ""
echo "======================================================================"
echo "ALL BULK RNA-SEQ EXPERIMENTS COMPLETE!"
echo "======================================================================"
echo "Results: experiments/results/${DATASET}_*"
echo ""
echo "Summary:"
echo "  ✓ Successful: $SUCCESS_COUNT/10"
if [ $FAIL_COUNT -gt 0 ]; then
  echo "  ✗ Failed:     $FAIL_COUNT/10"
  exit 1
else
  echo "  All experiments completed successfully!"
  exit 0
fi
