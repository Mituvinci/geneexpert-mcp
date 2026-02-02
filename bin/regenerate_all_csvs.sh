#!/bin/bash
#
# Regenerate all CSV files with fixed Stage 3 extraction
#
# This script finds all *_agent_decisions.json files and regenerates
# their corresponding CSV files with the updated json_to_csv.js that
# properly extracts DE_Method and Outlier_Action from Stage 3 agent responses.

echo "======================================================================="
echo "REGENERATING ALL CSV FILES WITH FIXED STAGE 3 EXTRACTION"
echo "======================================================================="
echo ""

RESULTS_DIR="/users/ha00014/Halimas_projects/multi_llm_mcp/experiments/results"
COUNT=0
SUCCESS=0
FAIL=0

# Find all JSON decision files (bulk RNA-seq)
while IFS= read -r json_file; do
    COUNT=$((COUNT + 1))
    FOLDER=$(dirname "$json_file")
    BASENAME=$(basename "$json_file" .json)
    OUTPUT_CSV="${FOLDER}/${BASENAME}_metrics.csv"

    echo "[$COUNT] Processing: $(basename "$FOLDER")"
    echo "    Input:  $(basename "$json_file")"
    echo "    Output: $(basename "$OUTPUT_CSV")"

    # Run conversion
    if node bin/json_to_csv.js convert --input "$json_file" --output "$OUTPUT_CSV" > /dev/null 2>&1; then
        SUCCESS=$((SUCCESS + 1))
        echo "    ✓ Success"
    else
        FAIL=$((FAIL + 1))
        echo "    ✗ Failed"
    fi
    echo ""
done < <(find "$RESULTS_DIR" -name "*_agent_decisions.json" -type f | sort)

echo "======================================================================="
echo "SUMMARY"
echo "======================================================================="
echo "Total files processed: $COUNT"
echo "Successful conversions: $SUCCESS"
echo "Failed conversions: $FAIL"
echo ""
echo "Next steps:"
echo "1. Run: python bin/aggregate_experiments.py"
echo "2. Run: node bin/evaluate_bulk.js"
echo ""
