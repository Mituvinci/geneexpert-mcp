#!/bin/bash
#
# Optimized ICML Experiment Runner
# Runs preprocessing ONCE, then 5 agent architectures on preprocessed data
# Saves ~70% time compared to redundant preprocessing
#
# Usage:
#   bash bin/run_icml_experiment_optimized.sh <dataset_name> <organism> <control_kw> <treatment_kw>
#
# Example:
#   bash bin/run_icml_experiment_optimized.sh 1_GSE52778_pe_clean human untreated Dex



# Check arguments
if [ $# -lt 4 ]; then
  echo "Usage: bash bin/run_icml_experiment_optimized.sh <dataset> <organism> <control_kw> <treatment_kw>"
  echo ""
  echo "Example:"
  echo "  bash bin/run_icml_experiment_optimized.sh 1_GSE52778_pe_clean human untreated Dex"
  exit 1
fi

DATASET=$1
ORGANISM=$2
CONTROL=$3
TREATMENT=$4

# Paths
DATA_DIR="data/download_new_bulk_RNA_Data/${DATASET}"
PREPROCESSING_DIR="experiments/preprocessing/${DATASET}"
RESULTS_DIR="experiments/results"
ROLE="gp_bl_cl_pl_gm_st"

echo "======================================================================"
echo "ICML 2026 EXPERIMENT (OPTIMIZED): ${DATASET}"
echo "======================================================================"
echo "Organism: $ORGANISM"
echo "Control:  $CONTROL"
echo "Treatment: $TREATMENT"
echo ""

# Validate dataset exists
if [ ! -d "$DATA_DIR" ]; then
  echo "❌ Error: Dataset not found: $DATA_DIR"
  exit 1
fi

echo "✓ Found dataset: $DATA_DIR"
echo ""

# ======================================================================
# PHASE 0: PREPROCESSING (Stage 1 + Stage 2 ONCE)
# ======================================================================

echo "======================================================================"
echo "PHASE 0: PREPROCESSING (Stage 1 FastQC + Stage 2 Alignment)"
echo "======================================================================"
echo ""

if [ -d "$PREPROCESSING_DIR/stage2_alignment/bam_files" ] && { [ "$(ls -A $PREPROCESSING_DIR/stage2_alignment/bam_files/*.bam 2>/dev/null)" ] || [ "$(ls -A $PREPROCESSING_DIR/stage2_alignment/bam_files/*.count.csv 2>/dev/null)" ]; }; then
  echo "✓ Preprocessing already exists, skipping..."
  echo "  FastQC: $PREPROCESSING_DIR/stage1_validation"
  echo "  BAM:    $PREPROCESSING_DIR/stage2_alignment/bam_files"
else
  echo "Running preprocessing (this takes ~30 min)..."
  echo ""

  node bin/preprocess_dataset.js \
    "$DATA_DIR" \
    "$PREPROCESSING_DIR" \
    --organism "$ORGANISM"

  echo ""
  echo "✓ Preprocessing complete!"
fi

echo ""



# System 4: Multi-Agent Parallel (default)
echo "--- System 4/5: Multi-Agent Parallel ---"
node bin/geneexpert.js analyze "$DATA_DIR" \
  --staged \
  --use-existing-fastqc "$PREPROCESSING_DIR/stage1_validation" \
  --use-existing-bam "$PREPROCESSING_DIR/stage2_alignment/bam_files" \
  --organism "$ORGANISM" \
  --comparison "$DATASET" \
  --control-keyword "$CONTROL" \
  --treatment-keyword "$TREATMENT" \
  --gpt-role "biology" \
  --claude-role "pipeline" \
  --gemini-role "stats" \
  --output "$RESULTS_DIR/${DATASET}_parallel_${ROLE}"
echo "✓ Multi-Agent Parallel complete"
echo ""

# Clean up BAM files
if [ -d "$RESULTS_DIR/${DATASET}_parallel_${ROLE}/stage2_alignment/bam_files" ]; then
  rm -f "$RESULTS_DIR/${DATASET}_parallel_${ROLE}/stage2_alignment/bam_files/"*.{bam,vcf} 2>/dev/null
fi
echo ""

# System 5: Multi-Agent Sequential Chain
echo "--- System 5/5: Multi-Agent Sequential Chain ---"
node bin/geneexpert.js analyze "$DATA_DIR" \
  --staged \
  --sequential-chain \
  --use-existing-fastqc "$PREPROCESSING_DIR/stage1_validation" \
  --use-existing-bam "$PREPROCESSING_DIR/stage2_alignment/bam_files" \
  --organism "$ORGANISM" \
  --comparison "$DATASET" \
  --control-keyword "$CONTROL" \
  --treatment-keyword "$TREATMENT" \
  --gpt-role "biology" \
  --claude-role "pipeline" \
  --gemini-role "stats" \
  --output "$RESULTS_DIR/${DATASET}_sequential_${ROLE}"
echo "✓ Multi-Agent Sequential Chain complete"
echo ""

# Clean up BAM files
if [ -d "$RESULTS_DIR/${DATASET}_sequential_${ROLE}/stage2_alignment/bam_files" ]; then
  rm -f "$RESULTS_DIR/${DATASET}_sequential_${ROLE}/stage2_alignment/bam_files/"*.{bam,vcf} 2>/dev/null
fi
echo ""

# ======================================================================
# SUMMARY
# ======================================================================

echo "======================================================================"
echo "EXPERIMENT COMPLETE: ${DATASET}"
echo "======================================================================"
echo ""
echo "Preprocessing:"
echo "  $PREPROCESSING_DIR/"
echo ""
echo "Results:"

echo "  4. Multi Parallel:  $RESULTS_DIR/${DATASET}_parallel_${ROLE}/"
echo "  5. Multi Sequential: $RESULTS_DIR/${DATASET}_sequential_${ROLE}/"
echo ""


