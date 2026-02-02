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


echo "[6/10] Role Config 2: GPT=pipeline, Claude=stats, Gemini=biology..."
if bash bin/1_2_run_bulk_role_icml_experiment_optimized_GPT_pipeline_Claude_stats_Gemini_biology.sh "$DATASET" "$ORGANISM" "$CONTROL" "$TREATMENT"; then
  ((SUCCESS_COUNT++))
else
  echo "⚠️  Experiment 6 failed, continuing..."
  ((FAIL_COUNT++))
fi

echo "[7/10] Role Config 3: GPT=pipeline, Claude=biology, Gemini=stats..."
if bash bin/1_3_run_bulk_role_icml_experiment_optimized_GPT_pipeline_Claude_biology_Gemini_stats.sh "$DATASET" "$ORGANISM" "$CONTROL" "$TREATMENT"; then
  ((SUCCESS_COUNT++))
else
  echo "⚠️  Experiment 7 failed, continuing..."
  ((FAIL_COUNT++))
fi

echo "[8/10] Role Config 4: GPT=biology, Claude=stats, Gemini=pipeline..."
if bash bin/1_4_run_bulk_role_icml_experiment_optimized_GPT_biology_Claude_stats_Gemini_pipeline.sh "$DATASET" "$ORGANISM" "$CONTROL" "$TREATMENT"; then
  ((SUCCESS_COUNT++))
else
  echo "⚠️  Experiment 8 failed, continuing..."
  ((FAIL_COUNT++))
fi

echo "[9/10] Role Config 5: GPT=biology, Claude=pipeline, Gemini=stats..."
if bash bin/1_5_run_bulk_role_icml_experiment_optimized_GPT_biology_Claude_pipeline_Gemini_stats.sh "$DATASET" "$ORGANISM" "$CONTROL" "$TREATMENT"; then
  ((SUCCESS_COUNT++))
else
  echo "⚠️  Experiment 9 failed, continuing..."
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
