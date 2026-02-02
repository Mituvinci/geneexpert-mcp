#!/bin/bash
#
# Fix batch effect EdgeR Excel outputs
# Bug: merge_results.R had hardcoded threshold column position (24)
# Fix: Now dynamically calculated based on number of samples
# Date: 2026-01-27
#

set -e

SCRIPTS_PATH="/users/ha00014/Halimas_projects/multi_llm_mcp/bio_informatics/scripts"
RESULTS_BASE="/users/ha00014/Halimas_projects/multi_llm_mcp/experiments/results"

echo "========================================================================"
echo "Batch Effect EdgeR Excel Output Fixer"
echo "========================================================================"
echo ""
echo "Bug: Threshold columns were inserted in middle of sample columns"
echo "Fix: Threshold columns now appear AFTER all data columns"
echo ""

# Activate conda
echo "Activating conda environment..."
source $(conda info --base)/etc/profile.d/conda.sh
conda activate pytorch

# Find all batch effect analyses
echo ""
echo "Finding all batch effect EdgeR analyses..."
BATCH_DIRS=$(find "$RESULTS_BASE" -name "stage_4_*.sh" -type f | xargs grep -l "DE_METHOD=\"batch_effect_edger\"" | xargs dirname | sort -u)

COUNT=0
for DIR in $BATCH_DIRS; do
    COUNT=$((COUNT + 1))
done

echo "Found $COUNT directories with batch effect EdgeR analyses"
echo ""

if [ $COUNT -eq 0 ]; then
    echo "No batch effect analyses found. Exiting."
    exit 0
fi

# Ask for confirmation
read -p "Regenerate Excel files for all $COUNT analyses? (y/n): " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Cancelled."
    exit 0
fi

# Process each directory
FIXED=0
FAILED=0

for DIR in $BATCH_DIRS; do
    echo ""
    echo "========================================================================"
    echo "Processing: $DIR"
    echo "========================================================================"

    # Extract comparison name and keywords from stage 4 script
    STAGE4_SCRIPT=$(ls "$DIR"/stage_4_*.sh 2>/dev/null | head -1)

    if [ ! -f "$STAGE4_SCRIPT" ]; then
        echo "ERROR: No stage 4 script found"
        FAILED=$((FAILED + 1))
        continue
    fi

    # Extract variables from script
    COMPARISON=$(grep "^COMPARISON=" "$STAGE4_SCRIPT" | cut -d'"' -f2)
    CONTROL_KW=$(grep "^CONTROL_KEYWORD=" "$STAGE4_SCRIPT" | cut -d'"' -f2)
    TREATMENT_KW=$(grep "^TREATMENT_KEYWORD=" "$STAGE4_SCRIPT" | cut -d'"' -f2)
    STAGE3_DIR=$(grep "^STAGE3_DIR=" "$STAGE4_SCRIPT" | cut -d'"' -f2)

    if [ -z "$COMPARISON" ] || [ -z "$CONTROL_KW" ] || [ -z "$TREATMENT_KW" ]; then
        echo "ERROR: Could not extract parameters from stage 4 script"
        FAILED=$((FAILED + 1))
        continue
    fi

    echo "Comparison: $COMPARISON"
    echo "Control keyword: $CONTROL_KW"
    echo "Treatment keyword: $TREATMENT_KW"

    # Check if files exist
    RPKM_FILE="${STAGE3_DIR}/${COMPARISON}.rpkm.entrz.csv"
    DE_FILE="${DIR}/${COMPARISON}_DE.csv"
    FINAL_FILE="${DIR}/${COMPARISON}_final.xlsx"

    if [ ! -f "$RPKM_FILE" ]; then
        echo "ERROR: RPKM file not found: $RPKM_FILE"
        FAILED=$((FAILED + 1))
        continue
    fi

    if [ ! -f "$DE_FILE" ]; then
        echo "ERROR: DE file not found: $DE_FILE"
        FAILED=$((FAILED + 1))
        continue
    fi

    # Backup old Excel file
    if [ -f "$FINAL_FILE" ]; then
        BACKUP_FILE="${FINAL_FILE%.xlsx}_BROKEN_backup.xlsx"
        echo "Backing up: $BACKUP_FILE"
        cp "$FINAL_FILE" "$BACKUP_FILE"
    fi

    # Regenerate Excel file
    echo "Regenerating Excel file..."
    cd "$DIR"

    if Rscript "$SCRIPTS_PATH/merge_results.R" \
        "$RPKM_FILE" \
        "$DE_FILE" \
        "$COMPARISON" \
        "$CONTROL_KW" \
        "$TREATMENT_KW" > /dev/null 2>&1; then

        echo "✓ SUCCESS: $FINAL_FILE"
        FIXED=$((FIXED + 1))
    else
        echo "✗ FAILED: Could not regenerate Excel file"
        FAILED=$((FAILED + 1))
    fi
done

echo ""
echo "========================================================================"
echo "SUMMARY"
echo "========================================================================"
echo "Total analyses found: $COUNT"
echo "Successfully fixed: $FIXED"
echo "Failed: $FAILED"
echo ""

if [ $FIXED -gt 0 ]; then
    echo "✓ Excel files have been regenerated with correct column structure."
    echo "  Old files backed up with _BROKEN_backup.xlsx suffix."
fi

echo ""
echo "Done!"
