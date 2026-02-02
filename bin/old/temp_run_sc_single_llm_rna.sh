#!/bin/bash

# NOTE: This script will continue running even if individual experiments fail
# Each experiment is independent and failures won't stop the batch

# Track success/failure
SUCCESS_COUNT=0
FAIL_COUNT=0

# Check arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <dataset> <organism>"
    echo "Example: $0 pbmc_healthy_human human"
    exit 1
fi

DATASET=$1
ORGANISM=$2


echo "======================================================================"
echo "RUNNING ALL scRNA-SEQ ABLATION EXPERIMENTS"
echo "======================================================================"
echo "Dataset:  $DATASET"
echo "Organism: $ORGANISM"
echo ""


echo "[1/4] Single-Agent GPT-5.2..."
if bash bin/2_0_2_run_sc_sing_agent_gpt_icml_experiment.sh "$DATASET" "$ORGANISM"; then
  ((SUCCESS_COUNT++))
else
  echo "⚠️  Experiment 1 failed, continuing..."
  ((FAIL_COUNT++))
fi

echo "[2/4] Single-Agent Claude..."
if bash bin/2_0_3_run_sc_sing_agent_claude_icml_experiment.sh "$DATASET" "$ORGANISM"; then
  ((SUCCESS_COUNT++))
else
  echo "⚠️  Experiment 2 failed, continuing..."
  ((FAIL_COUNT++))
fi

echo "[3/4] Single-Agent Gemini..."
if bash bin/2_0_4_run_sc_sing_agent_gemini_icml_experiment.sh "$DATASET" "$ORGANISM"; then
  ((SUCCESS_COUNT++))
else
  echo "⚠️  Experiment 3 failed, continuing..."
  ((FAIL_COUNT++))
fi


echo "[4/4] Default Configuration..."
if bash bin/2_6_run_sc_default_role_icml_experiment_default_parallel_sequential.sh "$DATASET" "$ORGANISM"; then
  ((SUCCESS_COUNT++))
else
  echo "⚠️  Experiment 10 failed, continuing..."
  ((FAIL_COUNT++))
fi


echo ""
echo "======================================================================"
echo "ALL scRNA-SEQ EXPERIMENTS COMPLETE!"
echo "======================================================================"
echo "Results: experiments/scrna_results/${DATASET}_*"
echo ""
echo "Summary:"
echo "  ✓ Successful: $SUCCESS_COUNT/4"
if [ $FAIL_COUNT -gt 0 ]; then
  echo "  ✗ Failed:     $FAIL_COUNT/4"
  exit 1
else
  echo "  All experiments completed successfully!"
  exit 0
fi 
