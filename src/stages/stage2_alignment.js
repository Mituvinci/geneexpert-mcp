/**
 * Stage 2: Alignment + Alignment QC
 *
 * This stage:
 * 1. Generates a script to run alignment (FASTQ -> BAM) using subread-align
 * 2. Runs alignment QC screening to check mapping rates
 * 3. Parses the output files
 * 4. Formats output for agent review
 *
 * Scripts used:
 * - subread-align: Align reads to reference genome
 * - alignment_qc_screening.R: Check mapping rates, detect possible contamination
 */

import fs from 'fs';
import path from 'path';
import dotenv from 'dotenv';

dotenv.config();

const SCRIPTS_PATH = process.env.SCRIPTS_PATH || './bio_informatics/scripts';
const CONDA_ENV = process.env.CONDA_ENV || 'pytorch';
const INDEX_PATH = process.env.INDEX_PATH || './bio_informatics/reference_data/index';

/**
 * Map organism name to genome build
 */
const GENOME_MAP = {
  'mouse': 'mm10',
  'human': 'hg38',
  'rat': 'rn6'
};

/**
 * Generate Stage 2 script for alignment and QC
 *
 * @param {Object} dataInfo - Dataset information (samples, pairedEnd, etc.)
 * @param {Object} config - Pipeline configuration
 * @param {string} outputPath - Path to save the script
 * @param {Object} stage1Decision - Decision from Stage 1 (proceed, warnings)
 * @returns {Object} - { scriptPath, scriptContent }
 */
export function generateStage2Script(dataInfo, config, outputPath, stage1Decision = {}) {
  console.log('[Stage 2] Generating alignment script...');

  const genomeBuild = GENOME_MAP[config.organism] || config.organism;
  const scriptLines = [];

  // Header
  scriptLines.push('#!/bin/bash');
  scriptLines.push('#');
  scriptLines.push('# GeneExpert Stage 2: Alignment + Alignment QC');
  scriptLines.push(`# Generated: ${new Date().toISOString()}`);
  scriptLines.push(`# Dataset: ${config.comparison}`);
  scriptLines.push(`# Organism: ${config.organism} (${genomeBuild})`);
  scriptLines.push(`# Sequencing: ${dataInfo.pairedEnd ? 'Paired-end' : 'Single-end'}`);
  scriptLines.push('#');
  scriptLines.push('set -e  # Exit on error');
  scriptLines.push('');

  // Activate conda environment
  scriptLines.push('# Activate conda environment');
  scriptLines.push('source $(conda info --base)/etc/profile.d/conda.sh');
  scriptLines.push(`conda activate ${CONDA_ENV}`);
  scriptLines.push('');

  // Variables
  scriptLines.push('# Paths');
  scriptLines.push(`INPUT_DIR="${config.input}"`);
  scriptLines.push(`OUTPUT_DIR="${config.output}"`);
  scriptLines.push(`SCRIPTS_DIR="${SCRIPTS_PATH}"`);
  scriptLines.push(`INDEX_PATH="${INDEX_PATH}"`);
  scriptLines.push(`GENOME="${genomeBuild}"`);
  scriptLines.push('');

  // Create output directories
  scriptLines.push('# Create output directories');
  scriptLines.push('mkdir -p "$OUTPUT_DIR/stage2_alignment"');
  scriptLines.push('mkdir -p "$OUTPUT_DIR/stage2_alignment/bam_files"');
  scriptLines.push('mkdir -p "$OUTPUT_DIR/stage2_alignment/alignment_qc"');
  scriptLines.push('');

  // Step 1: Alignment
  scriptLines.push('echo "=========================================="');
  scriptLines.push('echo "Stage 2: Alignment + Alignment QC"');
  scriptLines.push('echo "=========================================="');
  scriptLines.push('echo ""');
  scriptLines.push('');

  scriptLines.push('echo "Step 2.1: Aligning reads to reference genome..."');
  scriptLines.push(`echo "  Genome: $GENOME"`);
  scriptLines.push(`echo "  Index: $INDEX_PATH/$GENOME"`);
  scriptLines.push(`echo "  Input: $INPUT_DIR"`);
  scriptLines.push(`echo "  Output: $OUTPUT_DIR/stage2_alignment/bam_files/"`);
  scriptLines.push('echo ""');
  scriptLines.push('');

  if (dataInfo.pairedEnd) {
    // Paired-end alignment with progress tracking
    scriptLines.push('# Calculate total size and estimate time');
    scriptLines.push('echo "Calculating file sizes and estimating alignment time..."');
    scriptLines.push('total_samples=$(ls "$INPUT_DIR/"*_R1_001.fastq.gz 2>/dev/null | wc -l)');
    scriptLines.push('total_size_bytes=0');
    scriptLines.push('for i in "$INPUT_DIR/"*_R1_001.fastq.gz; do');
    scriptLines.push('  fname=$(basename "$i" _R1_001.fastq.gz)');
    scriptLines.push('  size_r1=$(stat -c%s "$INPUT_DIR/${fname}_R1_001.fastq.gz" 2>/dev/null || stat -f%z "$INPUT_DIR/${fname}_R1_001.fastq.gz" 2>/dev/null)');
    scriptLines.push('  size_r2=$(stat -c%s "$INPUT_DIR/${fname}_R2_001.fastq.gz" 2>/dev/null || stat -f%z "$INPUT_DIR/${fname}_R2_001.fastq.gz" 2>/dev/null)');
    scriptLines.push('  total_size_bytes=$((total_size_bytes + size_r1 + size_r2))');
    scriptLines.push('done');
    scriptLines.push('total_size_gb=$(awk "BEGIN {printf \\"%.1f\\", $total_size_bytes / 1024 / 1024 / 1024}")');
    scriptLines.push('estimated_minutes=$(awk "BEGIN {printf \\"%.0f\\", $total_size_gb * 2.5}")');
    scriptLines.push('if [ "$estimated_minutes" -lt 1 ]; then estimated_minutes=1; fi');
    scriptLines.push('echo "Total samples: $total_samples"');
    scriptLines.push('echo "Total data size: ${total_size_gb} GB"');
    scriptLines.push('echo "Estimated time: ~${estimated_minutes} minutes (running in parallel)"');
    scriptLines.push('echo ""');
    scriptLines.push('');
    scriptLines.push('start_time=$(date +%s)');
    scriptLines.push('current_sample=0');
    scriptLines.push('');
    scriptLines.push('# Run alignment for paired-end samples');
    scriptLines.push('for i in "$INPUT_DIR/"*_R1_001.fastq.gz; do');
    scriptLines.push('  fname=$(basename "$i" _R1_001.fastq.gz)');
    scriptLines.push('  current_sample=$((current_sample + 1))');
    scriptLines.push('  size_r1=$(stat -c%s "$INPUT_DIR/${fname}_R1_001.fastq.gz" 2>/dev/null || stat -f%z "$INPUT_DIR/${fname}_R1_001.fastq.gz" 2>/dev/null)');
    scriptLines.push('  size_r2=$(stat -c%s "$INPUT_DIR/${fname}_R2_001.fastq.gz" 2>/dev/null || stat -f%z "$INPUT_DIR/${fname}_R2_001.fastq.gz" 2>/dev/null)');
    scriptLines.push('  sample_size_gb=$(awk "BEGIN {printf \\"%.1f\\", ($size_r1 + $size_r2) / 1024 / 1024 / 1024}")');
    scriptLines.push('  sample_est_min=$(awk "BEGIN {printf \\"%.1f\\", $sample_size_gb * 2.5}")');
    scriptLines.push('  echo "[$current_sample/$total_samples] Aligning $fname (${sample_size_gb} GB, ~${sample_est_min} min)..."');
    scriptLines.push('  subread-align -t 0 \\');
    scriptLines.push('    -i "$INDEX_PATH/$GENOME" \\');
    scriptLines.push('    -r "$INPUT_DIR/${fname}_R1_001.fastq.gz" \\');
    scriptLines.push('    -R "$INPUT_DIR/${fname}_R2_001.fastq.gz" \\');
    scriptLines.push('    -T 8 \\');
    scriptLines.push('    -o "$OUTPUT_DIR/stage2_alignment/bam_files/${fname}.bam" &> "$OUTPUT_DIR/stage2_alignment/bam_files/${fname}.log" &');
    scriptLines.push('done');
    scriptLines.push('');
    scriptLines.push('# Wait for all alignment jobs to complete');
    scriptLines.push('echo ""');
    scriptLines.push('echo "Waiting for all $total_samples alignment jobs to complete..."');
    scriptLines.push('wait');
    scriptLines.push('end_time=$(date +%s)');
    scriptLines.push('elapsed=$((end_time - start_time))');
    scriptLines.push('elapsed_min=$((elapsed / 60))');
    scriptLines.push('elapsed_sec=$((elapsed % 60))');
    scriptLines.push('echo "All alignment jobs completed in ${elapsed_min}m ${elapsed_sec}s"');
    scriptLines.push('echo ""');
  } else {
    // Single-end alignment with progress tracking
    scriptLines.push('# Calculate total size and estimate time');
    scriptLines.push('echo "Calculating file sizes and estimating alignment time..."');
    scriptLines.push('total_samples=$(ls "$INPUT_DIR/"*.fastq.gz 2>/dev/null | grep -v "_R2_" | wc -l)');
    scriptLines.push('total_size_bytes=0');
    scriptLines.push('for i in "$INPUT_DIR/"*.fastq.gz; do');
    scriptLines.push('  fname=$(basename "$i" .fastq.gz)');
    scriptLines.push('  if [[ "$fname" == *"_R2_"* ]]; then continue; fi');
    scriptLines.push('  size=$(stat -c%s "$i" 2>/dev/null || stat -f%z "$i" 2>/dev/null)');
    scriptLines.push('  total_size_bytes=$((total_size_bytes + size))');
    scriptLines.push('done');
    scriptLines.push('total_size_gb=$(awk "BEGIN {printf \\"%.1f\\", $total_size_bytes / 1024 / 1024 / 1024}")');
    scriptLines.push('estimated_minutes=$(awk "BEGIN {printf \\"%.0f\\", $total_size_gb * 2.5}")');
    scriptLines.push('if [ "$estimated_minutes" -lt 1 ]; then estimated_minutes=1; fi');
    scriptLines.push('echo "Total samples: $total_samples"');
    scriptLines.push('echo "Total data size: ${total_size_gb} GB"');
    scriptLines.push('echo "Estimated time: ~${estimated_minutes} minutes (running in parallel)"');
    scriptLines.push('echo ""');
    scriptLines.push('');
    scriptLines.push('start_time=$(date +%s)');
    scriptLines.push('current_sample=0');
    scriptLines.push('');
    scriptLines.push('# Run alignment for single-end samples');
    scriptLines.push('for i in "$INPUT_DIR/"*.fastq.gz; do');
    scriptLines.push('  fname=$(basename "$i" .fastq.gz)');
    scriptLines.push('  if [[ "$fname" == *"_R2_"* ]]; then continue; fi');
    scriptLines.push('  current_sample=$((current_sample + 1))');
    scriptLines.push('  size=$(stat -c%s "$i" 2>/dev/null || stat -f%z "$i" 2>/dev/null)');
    scriptLines.push('  sample_size_gb=$(awk "BEGIN {printf \\"%.1f\\", $size / 1024 / 1024 / 1024}")');
    scriptLines.push('  sample_est_min=$(awk "BEGIN {printf \\"%.1f\\", $sample_size_gb * 2.5}")');
    scriptLines.push('  echo "[$current_sample/$total_samples] Aligning $fname (${sample_size_gb} GB, ~${sample_est_min} min)..."');
    scriptLines.push('  subread-align -t 0 \\');
    scriptLines.push('    -i "$INDEX_PATH/$GENOME" \\');
    scriptLines.push('    -r "$INPUT_DIR/${fname}.fastq.gz" \\');
    scriptLines.push('    -T 8 \\');
    scriptLines.push('    -o "$OUTPUT_DIR/stage2_alignment/bam_files/${fname}.bam" &> "$OUTPUT_DIR/stage2_alignment/bam_files/${fname}.log" &');
    scriptLines.push('done');
    scriptLines.push('');
    scriptLines.push('# Wait for all alignment jobs to complete');
    scriptLines.push('echo ""');
    scriptLines.push('echo "Waiting for all $total_samples alignment jobs to complete..."');
    scriptLines.push('wait');
    scriptLines.push('end_time=$(date +%s)');
    scriptLines.push('elapsed=$((end_time - start_time))');
    scriptLines.push('elapsed_min=$((elapsed / 60))');
    scriptLines.push('elapsed_sec=$((elapsed % 60))');
    scriptLines.push('echo "All alignment jobs completed in ${elapsed_min}m ${elapsed_sec}s"');
    scriptLines.push('echo ""');
  }

  scriptLines.push('');
  scriptLines.push('# Verify BAM files were created successfully');
  scriptLines.push('echo "Verifying BAM files..."');
  scriptLines.push('bam_count=$(ls "$OUTPUT_DIR/stage2_alignment/bam_files/"*.bam 2>/dev/null | wc -l)');
  scriptLines.push('if [ "$bam_count" -eq 0 ]; then');
  scriptLines.push('  echo "ERROR: No BAM files created by subread-align!"');
  scriptLines.push('  echo "Check alignment logs in: $OUTPUT_DIR/stage2_alignment/bam_files/"');
  scriptLines.push('  exit 1');
  scriptLines.push('fi');
  scriptLines.push('echo "Successfully created $bam_count BAM files in $OUTPUT_DIR/stage2_alignment/bam_files/"');
  scriptLines.push('');

  // Step 2: Alignment QC
  scriptLines.push('echo ""');
  scriptLines.push('echo "Step 2.2: Running alignment QC screening..."');
  scriptLines.push('');
  scriptLines.push('Rscript "$SCRIPTS_DIR/alignment_qc_screening.R" \\');
  scriptLines.push('  "$OUTPUT_DIR/stage2_alignment/bam_files" \\');
  scriptLines.push('  "$OUTPUT_DIR/stage2_alignment/alignment_qc"');
  scriptLines.push('');
  scriptLines.push('ALIGN_QC_EXIT=$?');
  scriptLines.push('');

  // Handle QC exit codes
  scriptLines.push('# Store exit code for agent review');
  scriptLines.push('echo "$ALIGN_QC_EXIT" > "$OUTPUT_DIR/stage2_alignment/alignment_qc/exit_code.txt"');
  scriptLines.push('');

  scriptLines.push('if [ $ALIGN_QC_EXIT -eq 0 ]; then');
  scriptLines.push('  echo "QC Status: All samples PASS"');
  scriptLines.push('elif [ $ALIGN_QC_EXIT -eq 1 ]; then');
  scriptLines.push('  echo "QC Status: Some samples have WARNINGS"');
  scriptLines.push('elif [ $ALIGN_QC_EXIT -eq 2 ]; then');
  scriptLines.push('  echo "QC Status: Some samples FAIL or SEVERE_FAIL"');
  scriptLines.push('  echo "Review: $OUTPUT_DIR/stage2_alignment/alignment_qc/alignment_qc_report.txt"');
  scriptLines.push('fi');
  scriptLines.push('');

  // Create JSON summary
  scriptLines.push('# Create JSON summary for agents');
  scriptLines.push('cat > "$OUTPUT_DIR/stage2_alignment/stage2_summary.json" << EOF');
  scriptLines.push('{');
  scriptLines.push('  "stage": 2,');
  scriptLines.push('  "stage_name": "Alignment + Alignment QC",');
  scriptLines.push('  "timestamp": "$(date -Iseconds)",');
  scriptLines.push('  "input_dir": "$INPUT_DIR",');
  scriptLines.push('  "output_dir": "$OUTPUT_DIR/stage2_alignment",');
  scriptLines.push('  "genome_build": "$GENOME",');
  scriptLines.push('  "bam_count": $bam_count,');
  scriptLines.push('  "qc_exit_code": $ALIGN_QC_EXIT,');
  scriptLines.push('  "alignment_qc_summary": "$OUTPUT_DIR/stage2_alignment/alignment_qc/alignment_qc_summary.csv",');
  scriptLines.push('  "alignment_qc_report": "$OUTPUT_DIR/stage2_alignment/alignment_qc/alignment_qc_report.txt"');
  scriptLines.push('}');
  scriptLines.push('EOF');
  scriptLines.push('');

  // Final status
  scriptLines.push('echo ""');
  scriptLines.push('echo "=========================================="');
  scriptLines.push('echo "Stage 2 Complete!"');
  scriptLines.push('echo "=========================================="');
  scriptLines.push('echo "Outputs:"');
  scriptLines.push('echo "  - BAM files: $OUTPUT_DIR/stage2_alignment/bam_files/"');
  scriptLines.push('echo "  - Alignment QC: $OUTPUT_DIR/stage2_alignment/alignment_qc/"');
  scriptLines.push('echo "  - Summary JSON: $OUTPUT_DIR/stage2_alignment/stage2_summary.json"');
  scriptLines.push('echo ""');
  scriptLines.push('echo "Ready for agent review."');
  scriptLines.push('');

  // Exit with QC code (so caller knows status)
  scriptLines.push('# Exit with QC status code');
  scriptLines.push('exit $ALIGN_QC_EXIT');

  const scriptContent = scriptLines.join('\n');

  // Write script to file
  fs.writeFileSync(outputPath, scriptContent, { mode: 0o755 });
  console.log(`[Stage 2] Script saved: ${outputPath}`);

  return {
    scriptPath: outputPath,
    scriptContent
  };
}

/**
 * Parse Stage 2 output files
 *
 * @param {string} outputDir - Output directory containing stage 2 results
 * @returns {Object} - Parsed stage 2 output
 */
export function parseStage2Output(outputDir) {
  console.log('[Stage 2] Parsing output files...');

  const stage2Dir = path.join(outputDir, 'stage2_alignment');
  const qcDir = path.join(stage2Dir, 'alignment_qc');

  const result = {
    overall_status: 'unknown',
    exit_code: null,
    samples_aligned: 0,
    per_sample_mapping: {},
    failed_samples: [],
    warned_samples: [],
    passed_samples: [],
    mean_mapping_rate: null,
    median_mapping_rate: null,
    min_mapping_rate: null,
    max_mapping_rate: null,
    warnings: [],
    errors: []
  };

  // Read exit code
  const exitCodePath = path.join(qcDir, 'exit_code.txt');
  if (fs.existsSync(exitCodePath)) {
    result.exit_code = parseInt(fs.readFileSync(exitCodePath, 'utf-8').trim());
  }

  // Parse alignment QC summary CSV
  const qcSummaryPath = path.join(qcDir, 'alignment_qc_summary.csv');
  if (fs.existsSync(qcSummaryPath)) {
    const csvContent = fs.readFileSync(qcSummaryPath, 'utf-8');
    const lines = csvContent.trim().split('\n');

    if (lines.length > 1) {
      // Parse header
      const header = lines[0].split(',').map(h => h.replace(/"/g, ''));

      // Find column indices
      const sampleIdx = header.indexOf('sample');
      const mappedPctIdx = header.indexOf('mapped_pct');
      const statusIdx = header.indexOf('status');
      const reasonsIdx = header.indexOf('reasons');
      const totalFragmentsIdx = header.indexOf('total_fragments');
      const properlyPairedPctIdx = header.indexOf('properly_paired_pct');

      const mappingRates = [];

      for (let i = 1; i < lines.length; i++) {
        const values = lines[i].split(',').map(v => v.replace(/"/g, ''));

        const sample = values[sampleIdx];
        const mappedPct = parseFloat(values[mappedPctIdx]);
        const status = values[statusIdx];
        const reasons = values[reasonsIdx];
        const totalFragments = parseInt(values[totalFragmentsIdx]) || 0;
        const properlyPairedPct = parseFloat(values[properlyPairedPctIdx]);

        result.per_sample_mapping[sample] = {
          mapping_rate: mappedPct / 100,  // Convert to decimal
          mapping_rate_pct: mappedPct,
          status,
          reasons,
          total_fragments: totalFragments,
          properly_paired_pct: properlyPairedPct
        };

        if (!isNaN(mappedPct)) {
          mappingRates.push(mappedPct / 100);
        }

        // Categorize by status
        if (status === 'PASS') {
          result.passed_samples.push(sample);
        } else if (status === 'WARN') {
          result.warned_samples.push(sample);
          result.warnings.push(`${sample}: ${reasons}`);
        } else if (status === 'FAIL' || status === 'SEVERE_FAIL') {
          result.failed_samples.push({ sample, mapping_rate: mappedPct / 100, status, reasons });
          result.errors.push(`${sample}: ${status} - ${reasons}`);
        }

        result.samples_aligned++;
      }

      // Calculate statistics
      if (mappingRates.length > 0) {
        result.mean_mapping_rate = mappingRates.reduce((a, b) => a + b, 0) / mappingRates.length;
        result.min_mapping_rate = Math.min(...mappingRates);
        result.max_mapping_rate = Math.max(...mappingRates);

        // Median
        const sorted = [...mappingRates].sort((a, b) => a - b);
        const mid = Math.floor(sorted.length / 2);
        result.median_mapping_rate = sorted.length % 2 !== 0
          ? sorted[mid]
          : (sorted[mid - 1] + sorted[mid]) / 2;
      }
    }
  } else {
    result.errors.push('Alignment QC summary not found');
  }

  // Parse alignment QC report (text)
  const qcReportPath = path.join(qcDir, 'alignment_qc_report.txt');
  if (fs.existsSync(qcReportPath)) {
    result.qc_report_text = fs.readFileSync(qcReportPath, 'utf-8');
  }

  // Count BAM files
  const bamDir = path.join(stage2Dir, 'bam_files');
  if (fs.existsSync(bamDir)) {
    const bamFiles = fs.readdirSync(bamDir).filter(f => f.endsWith('.bam'));
    result.bam_count = bamFiles.length;
  }

  // Determine overall status
  if (result.exit_code === 0) {
    result.overall_status = 'PASS';
    result.recommendation = 'All samples passed alignment QC. Proceed to Stage 3.';
  } else if (result.exit_code === 1) {
    result.overall_status = 'WARN';
    result.recommendation = 'Some samples have warnings. Review and decide whether to proceed.';
  } else if (result.exit_code === 2) {
    result.overall_status = 'FAIL';
    result.recommendation = 'Some samples failed alignment QC. Consider removing failed samples or investigating root cause.';
  } else {
    result.overall_status = 'UNKNOWN';
    result.recommendation = 'Could not determine alignment QC status.';
  }

  console.log(`[Stage 2] Parsing complete. Status: ${result.overall_status}`);
  return result;
}

/**
 * Format Stage 2 output for agent review
 *
 * @param {Object} parsedOutput - Parsed stage 2 output from parseStage2Output()
 * @param {Object} dataInfo - Dataset information
 * @returns {string} - Formatted string for agent consumption
 */
export function formatStage2ForAgents(parsedOutput, dataInfo) {
  const lines = [];

  lines.push('## Alignment QC Results - Stage 2 Output');
  lines.push('');

  // Dataset overview
  lines.push('### Dataset Overview');
  lines.push(`- Organism: ${dataInfo.organism || 'unknown'}`);
  lines.push(`- Genome Build: ${parsedOutput.genome_build || GENOME_MAP[dataInfo.organism] || 'unknown'}`);
  lines.push(`- Sequencing Type: ${dataInfo.pairedEnd ? 'Paired-end' : 'Single-end'}`);
  lines.push(`- Samples Aligned: ${parsedOutput.samples_aligned}`);
  lines.push(`- BAM Files Created: ${parsedOutput.bam_count || parsedOutput.samples_aligned}`);
  lines.push('');

  // Mapping rate statistics
  lines.push('### Mapping Rate Summary');
  lines.push(`- Mean Mapping Rate: ${parsedOutput.mean_mapping_rate ? (parsedOutput.mean_mapping_rate * 100).toFixed(1) + '%' : 'N/A'}`);
  lines.push(`- Median Mapping Rate: ${parsedOutput.median_mapping_rate ? (parsedOutput.median_mapping_rate * 100).toFixed(1) + '%' : 'N/A'}`);
  lines.push(`- Range: ${parsedOutput.min_mapping_rate ? (parsedOutput.min_mapping_rate * 100).toFixed(1) : '?'}% - ${parsedOutput.max_mapping_rate ? (parsedOutput.max_mapping_rate * 100).toFixed(1) : '?'}%`);
  lines.push('');

  // Thresholds reminder
  lines.push('### Mapping Rate Thresholds');
  lines.push('- >= 80%: PASS (good quality)');
  lines.push('- 70-80%: WARN (acceptable but borderline)');
  lines.push('- 60-70%: FAIL (poor quality)');
  lines.push('- < 60%: SEVERE_FAIL (possible contamination or wrong genome)');
  lines.push('');

  // Per-sample results
  lines.push('### Per-Sample Mapping Rates');
  lines.push('| Sample | Mapping Rate | Status | Notes |');
  lines.push('|--------|-------------|--------|-------|');

  for (const [sample, data] of Object.entries(parsedOutput.per_sample_mapping)) {
    const mappingRate = data.mapping_rate_pct ? data.mapping_rate_pct.toFixed(1) + '%' : 'N/A';
    const status = data.status || 'unknown';
    const notes = data.reasons || '';
    lines.push(`| ${sample} | ${mappingRate} | ${status} | ${notes.substring(0, 50)} |`);
  }
  lines.push('');

  // Status counts
  lines.push('### Status Summary');
  lines.push(`- PASS: ${parsedOutput.passed_samples.length} samples`);
  lines.push(`- WARN: ${parsedOutput.warned_samples.length} samples`);
  lines.push(`- FAIL/SEVERE_FAIL: ${parsedOutput.failed_samples.length} samples`);
  lines.push('');

  // Failed samples detail
  if (parsedOutput.failed_samples.length > 0) {
    lines.push('### Failed Samples (Require Decision)');
    for (const sample of parsedOutput.failed_samples) {
      lines.push(`- **${sample.sample}**: ${(sample.mapping_rate * 100).toFixed(1)}% mapping rate [${sample.status}]`);
      lines.push(`  - Reason: ${sample.reasons}`);
    }
    lines.push('');
  }

  // Warnings
  if (parsedOutput.warnings.length > 0) {
    lines.push('### Warnings');
    for (const warning of parsedOutput.warnings) {
      lines.push(`- ${warning}`);
    }
    lines.push('');
  }

  // Group balance check
  if (dataInfo.groups) {
    lines.push('### Group Balance After Potential Removal');
    const groupCounts = {};

    for (const [group, samples] of Object.entries(dataInfo.groups)) {
      const failedInGroup = samples.filter(s =>
        parsedOutput.failed_samples.some(f => f.sample.includes(s) || s.includes(f.sample))
      );
      const remaining = samples.length - failedInGroup.length;
      groupCounts[group] = { original: samples.length, remaining, failed: failedInGroup };
      lines.push(`- ${group}: ${remaining}/${samples.length} samples would remain`);
    }

    // Warning if removal would create imbalance
    const remainingCounts = Object.values(groupCounts).map(g => g.remaining);
    if (Math.min(...remainingCounts) < 2) {
      lines.push('');
      lines.push('**WARNING: Removing failed samples would leave <2 samples in a group!**');
    }
    lines.push('');
  }

  // Overall status
  lines.push('### Overall Status');
  lines.push(`- QC Exit Code: ${parsedOutput.exit_code} (0=PASS, 1=WARN, 2=FAIL)`);
  lines.push(`- Overall: **${parsedOutput.overall_status}**`);
  lines.push(`- Recommendation: ${parsedOutput.recommendation}`);
  lines.push('');

  // Decision prompt
  lines.push('### Decision Required');
  lines.push('Based on the above results, choose one:');
  lines.push('- **PASS_ALL**: Proceed with all samples');
  lines.push('- **REMOVE_SAMPLES**: Remove failed samples and proceed (list which ones)');
  lines.push('- **ABORT**: Stop analysis - too many failures or systematic issues');
  lines.push('');

  return lines.join('\n');
}

/**
 * Execute Stage 2 and return results
 *
 * @param {Object} dataInfo - Dataset information
 * @param {Object} config - Pipeline configuration
 * @param {Object} stage1Decision - Decision from Stage 1
 * @returns {Promise<Object>} - { success, scriptPath, output }
 */
export async function executeStage2(dataInfo, config, stage1Decision = {}) {
  const scriptPath = path.join(config.output, `stage_2_${config.comparison}.sh`);

  // Generate script
  generateStage2Script(dataInfo, config, scriptPath, stage1Decision);

  return {
    scriptPath,
    scriptName: `stage_2_${config.comparison}.sh`,
    stage: 2,
    stageName: 'Alignment + Alignment QC'
  };
}

export default {
  generateStage2Script,
  parseStage2Output,
  formatStage2ForAgents,
  executeStage2
};
