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

echo "[5/10] Role Config 1: GPT=stats, Claude=biology, Gemini=pipeline..."
if bash bin/2_1_run_sc_role_icml_experiment_GPT_stats_Claude_biology_Gemini_pipeline.sh "$DATASET" "$ORGANISM"; then
  ((SUCCESS_COUNT++))
else
  echo "⚠️  Experiment 5 failed, continuing..."
  ((FAIL_COUNT++))
fi

echo "[6/10] Role Config 2: GPT=pipeline, Claude=stats, Gemini=biology..."
if bash bin/2_2_run_sc_role_icml_experiment_GPT_pipeline_Claude_stats_Gemini_biology.sh "$DATASET" "$ORGANISM"; then
  ((SUCCESS_COUNT++))
else
  echo "⚠️  Experiment 6 failed, continuing..."
  ((FAIL_COUNT++))
fi

echo "[7/10] Role Config 3: GPT=pipeline, Claude=biology, Gemini=stats..."
if bash bin/2_3_run_sc_role_icml_experiment_GPT_pipeline_Claude_biology_Gemini_stats.sh "$DATASET" "$ORGANISM"; then
  ((SUCCESS_COUNT++))
else
  echo "⚠️  Experiment 7 failed, continuing..."
  ((FAIL_COUNT++))
fi

echo "[8/10] Role Config 4: GPT=biology, Claude=stats, Gemini=pipeline..."
if bash bin/2_4_run_sc_role_icml_experiment_GPT_biology_Claude_stats_Gemini_pipeline.sh "$DATASET" "$ORGANISM"; then
  ((SUCCESS_COUNT++))
else
  echo "⚠️  Experiment 8 failed, continuing..."
  ((FAIL_COUNT++))
fi

echo "[9/10] Role Config 5: GPT=biology, Claude=pipeline, Gemini=stats..."
if bash bin/2_5_run_sc_role_icml_experiment_GPT_biology_Claude_pipeline_Gemini_stats.sh "$DATASET" "$ORGANISM"; then
  ((SUCCESS_COUNT++))
else
  echo "⚠️  Experiment 9 failed, continuing..."
  ((FAIL_COUNT++))
fi

echo "[10/10] Default Configuration..."
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
echo "  ✓ Successful: $SUCCESS_COUNT/6"
if [ $FAIL_COUNT -gt 0 ]; then
  echo "  ✗ Failed:     $FAIL_COUNT/6"
  exit 1
else
  echo "  All experiments completed successfully!"
  exit 0
fi 
