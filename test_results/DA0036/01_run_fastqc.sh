#!/bin/bash

# FastQC Quality Control
# Input: /destiny/fastqs/DA0036/ (READ-ONLY)
# Output: /data/halimaakhter/multi_llm_mcp/test_results/DA0036/fastqc/

echo "Running FastQC on DA0036 samples..."
echo "Input: /destiny/fastqs/DA0036/"
echo "Output: /data/halimaakhter/multi_llm_mcp/test_results/DA0036/fastqc/"
echo ""

# Output directory (in YOUR folder, not /destiny!)
OUTPUT_DIR="/data/halimaakhter/multi_llm_mcp/test_results/DA0036/fastqc"
INPUT_DIR="/destiny/fastqs/DA0036"

# Create output directory if it doesn't exist
mkdir -p ${OUTPUT_DIR}

# Run FastQC on all samples (READ-ONLY from /destiny, WRITE to your folder)
fastqc ${INPUT_DIR}/*.fastq.gz \
    -o ${OUTPUT_DIR} \
    -t 4

echo ""
echo "✅ FastQC complete!"
echo "Results saved to: ${OUTPUT_DIR}"
echo ""
echo "Check the HTML reports:"
ls -lh ${OUTPUT_DIR}/*.html
