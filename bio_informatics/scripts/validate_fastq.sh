#!/usr/bin/env bash

# Usage:
# ./validate_fastq.sh <input_dir> <output_report.tsv>

INPUT_DIR="$1"
REPORT="$2"

echo -e "file\treads\tlines\tmod4_ok\tpaired\tstatus" > "$REPORT"

declare -A READ_COUNTS

# Count total files and estimate time
echo ""
echo "Calculating file sizes and estimating time..."
total_files=$(ls "$INPUT_DIR"/*.fastq.gz 2>/dev/null | wc -l)
total_size_bytes=0

for f in "$INPUT_DIR"/*.fastq.gz; do
    size_bytes=$(stat -c%s "$f" 2>/dev/null || stat -f%z "$f" 2>/dev/null)
    total_size_bytes=$((total_size_bytes + size_bytes))
done

# Convert to GB and estimate time using awk
total_size_gb=$(awk "BEGIN {printf \"%.1f\", $total_size_bytes / 1024 / 1024 / 1024}")
estimated_minutes=$(awk "BEGIN {printf \"%.0f\", $total_size_gb * 1.0}")
if [ "$estimated_minutes" -lt 1 ]; then
    estimated_minutes=1
fi

echo "Total files to validate: $total_files"
echo "Total data size: ${total_size_gb} GB"
echo "Estimated time: ~${estimated_minutes} minutes"
echo "Processing FASTQ integrity checks..."
echo ""

current_file=0
start_time=$(date +%s)

for f in "$INPUT_DIR"/*.fastq.gz; do
    fname=$(basename "$f")
    current_file=$((current_file + 1))

    # Get file size for this specific file
    size_bytes=$(stat -c%s "$f" 2>/dev/null || stat -f%z "$f" 2>/dev/null)
    size_gb=$(awk "BEGIN {printf \"%.1f\", $size_bytes / 1024 / 1024 / 1024}")
    file_est_min=$(awk "BEGIN {printf \"%.1f\", $size_gb * 1.0}")

    echo "[$current_file/$total_files] Processing $fname (${size_gb} GB, ~${file_est_min} min)..."

    # Count lines (streaming, safe)
    total_lines=$(zcat "$f" | wc -l)

    # FASTQ rule
    if (( total_lines % 4 != 0 )); then
        echo -e "$fname\tNA\t$total_lines\tFAIL\tNA\tINVALID_FASTQ" >> "$REPORT"
        echo "  ✗ INVALID: Line count not divisible by 4"
        continue
    fi

    reads=$(( total_lines / 4 ))

    # Store for pairing check
    base=${fname%%_R[12]_001.fastq.gz}
    READ_COUNTS["$base"]+="$fname:$reads "

    echo -e "$fname\t$reads\t$total_lines\tOK\tPENDING\tVALID" >> "$REPORT"
    echo "  ✓ Valid: $reads reads"
done

echo ""
end_time=$(date +%s)
elapsed=$((end_time - start_time))
elapsed_min=$((elapsed / 60))
elapsed_sec=$((elapsed % 60))
echo "Validation completed in ${elapsed_min}m ${elapsed_sec}s"
echo ""

# Paired-end consistency check
for base in "${!READ_COUNTS[@]}"; do
    entries=(${READ_COUNTS[$base]})
    if (( ${#entries[@]} == 2 )); then
        r1_reads=${entries[0]##*:}
        r2_reads=${entries[1]##*:}
        if (( r1_reads != r2_reads )); then
            sed -i "s/$base.*VALID/$base\tMISMATCH/g" "$REPORT"
        fi
    fi
done
