/**
 * Stage 1: FASTQ Validation
 *
 * This stage:
 * 1. Generates a script to run FASTQ validation and FastQC
 * 2. Parses the output files
 * 3. Formats output for agent review
 *
 * Scripts used:
 * - validate_fastq.sh: Check file integrity (lines % 4 == 0, paired consistency)
 * - fastqc: Generate quality reports
 */

import fs from 'fs';
import path from 'path';
import dotenv from 'dotenv';

dotenv.config();

const SCRIPTS_PATH = process.env.SCRIPTS_PATH || '/users/ha00014/Halimas_projects/multi_llm_mcp/bio_informatics/scripts';
const CONDA_ENV = process.env.CONDA_ENV || 'pytorch';

/**
 * Generate Stage 1 script for FASTQ validation and FastQC
 *
 * @param {Object} dataInfo - Dataset information (samples, pairedEnd, etc.)
 * @param {Object} config - Pipeline configuration (input, output, organism, etc.)
 * @param {string} outputPath - Path to save the script
 * @returns {Object} - { scriptPath, scriptContent }
 */
export function generateStage1Script(dataInfo, config, outputPath) {
  console.log('[Stage 1] Generating FASTQ validation script...');

  const scriptLines = [];

  // Header
  scriptLines.push('#!/bin/bash');
  scriptLines.push('#');
  scriptLines.push('# GeneExpert Stage 1: FASTQ Validation');
  scriptLines.push(`# Generated: ${new Date().toISOString()}`);
  scriptLines.push(`# Dataset: ${config.comparison}`);
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
  scriptLines.push('');

  // Create output directories
  scriptLines.push('# Create output directories');
  scriptLines.push('mkdir -p "$OUTPUT_DIR/stage1_validation"');
  scriptLines.push('mkdir -p "$OUTPUT_DIR/stage1_validation/fastqc_results"');
  scriptLines.push('');

  // Step 1: FASTQ Validation
  scriptLines.push('echo "=========================================="');
  scriptLines.push('echo "Stage 1: FASTQ Validation"');
  scriptLines.push('echo "=========================================="');
  scriptLines.push('echo ""');
  scriptLines.push('');

  scriptLines.push('echo "Step 1.1: Validating FASTQ file integrity..."');
  scriptLines.push('bash "$SCRIPTS_DIR/validate_fastq.sh" "$INPUT_DIR" "$OUTPUT_DIR/stage1_validation/validation_report.tsv"');
  scriptLines.push('');

  scriptLines.push('# Check if validation found any issues');
  scriptLines.push('if grep -q "INVALID_FASTQ\\|MISMATCH" "$OUTPUT_DIR/stage1_validation/validation_report.tsv"; then');
  scriptLines.push('  echo "WARNING: Some FASTQ files have validation issues!"');
  scriptLines.push('  echo "Check: $OUTPUT_DIR/stage1_validation/validation_report.tsv"');
  scriptLines.push('  # Don\'t exit - let agents decide what to do');
  scriptLines.push('fi');
  scriptLines.push('');

  // Step 2: FastQC
  scriptLines.push('echo ""');
  scriptLines.push('echo "Step 1.2: Running FastQC quality control..."');

  // Build FastQC command based on file type
  scriptLines.push('cd "$INPUT_DIR"');
  scriptLines.push('fastqc *.fastq.gz -o "$OUTPUT_DIR/stage1_validation/fastqc_results" -t 4');
  scriptLines.push('cd -');
  scriptLines.push('');

  // Step 3: Generate summary
  scriptLines.push('echo ""');
  scriptLines.push('echo "Step 1.3: Generating FastQC summary..."');
  scriptLines.push('');

  // Create a summary of FastQC results
  scriptLines.push('# Extract FastQC summary for each sample');
  scriptLines.push('echo "sample,basic_statistics,per_base_quality,per_sequence_quality,per_base_content,gc_content,per_base_n,sequence_length,duplication,overrepresented,adapter_content" > "$OUTPUT_DIR/stage1_validation/fastqc_summary.csv"');
  scriptLines.push('');
  scriptLines.push('for zipfile in "$OUTPUT_DIR/stage1_validation/fastqc_results"/*_fastqc.zip; do');
  scriptLines.push('  if [ -f "$zipfile" ]; then');
  scriptLines.push('    sample=$(basename "$zipfile" _fastqc.zip)');
  scriptLines.push('    # Extract summary.txt from zip');
  scriptLines.push('    unzip -p "$zipfile" "*/summary.txt" 2>/dev/null | \\');
  scriptLines.push('    awk -v sample="$sample" \'BEGIN{ORS=","} {print $1} END{print "\\n"}\' | \\');
  scriptLines.push('    sed "s/^/$sample,/" >> "$OUTPUT_DIR/stage1_validation/fastqc_summary.csv"');
  scriptLines.push('  fi');
  scriptLines.push('done');
  scriptLines.push('');

  // Step 4: Calculate read count statistics
  scriptLines.push('echo ""');
  scriptLines.push('echo "Step 1.4: Calculating read count statistics..."');
  scriptLines.push('');
  scriptLines.push('# Create read count summary');
  scriptLines.push('echo "sample,reads,status" > "$OUTPUT_DIR/stage1_validation/read_counts.csv"');
  scriptLines.push('tail -n +2 "$OUTPUT_DIR/stage1_validation/validation_report.tsv" | \\');
  scriptLines.push('  awk -F\'\\t\' \'{print $1","$2","$6}\' >> "$OUTPUT_DIR/stage1_validation/read_counts.csv"');
  scriptLines.push('');

  // Calculate summary statistics
  scriptLines.push('# Calculate summary statistics');
  scriptLines.push('total_reads=$(tail -n +2 "$OUTPUT_DIR/stage1_validation/validation_report.tsv" | awk -F\'\\t\' \'{sum+=$2} END{print sum}\')');
  scriptLines.push('num_samples=$(tail -n +2 "$OUTPUT_DIR/stage1_validation/validation_report.tsv" | wc -l)');
  scriptLines.push('avg_reads=$((total_reads / num_samples))');
  scriptLines.push('');
  scriptLines.push('echo "Total reads: $total_reads"');
  scriptLines.push('echo "Number of samples: $num_samples"');
  scriptLines.push('echo "Average reads per sample: $avg_reads"');
  scriptLines.push('');

  // Create JSON summary for agent consumption
  scriptLines.push('# Create JSON summary for agents');
  scriptLines.push('cat > "$OUTPUT_DIR/stage1_validation/stage1_summary.json" << EOF');
  scriptLines.push('{');
  scriptLines.push('  "stage": 1,');
  scriptLines.push('  "stage_name": "FASTQ Validation",');
  scriptLines.push('  "timestamp": "$(date -Iseconds)",');
  scriptLines.push('  "input_dir": "$INPUT_DIR",');
  scriptLines.push('  "output_dir": "$OUTPUT_DIR/stage1_validation",');
  scriptLines.push('  "total_samples": $num_samples,');
  scriptLines.push('  "total_reads": $total_reads,');
  scriptLines.push('  "avg_reads_per_sample": $avg_reads,');
  scriptLines.push('  "validation_report": "$OUTPUT_DIR/stage1_validation/validation_report.tsv",');
  scriptLines.push('  "fastqc_summary": "$OUTPUT_DIR/stage1_validation/fastqc_summary.csv",');
  scriptLines.push('  "read_counts": "$OUTPUT_DIR/stage1_validation/read_counts.csv"');
  scriptLines.push('}');
  scriptLines.push('EOF');
  scriptLines.push('');

  // Final status
  scriptLines.push('echo ""');
  scriptLines.push('echo "=========================================="');
  scriptLines.push('echo "Stage 1 Complete!"');
  scriptLines.push('echo "=========================================="');
  scriptLines.push('echo "Outputs:"');
  scriptLines.push('echo "  - Validation report: $OUTPUT_DIR/stage1_validation/validation_report.tsv"');
  scriptLines.push('echo "  - FastQC results: $OUTPUT_DIR/stage1_validation/fastqc_results/"');
  scriptLines.push('echo "  - FastQC summary: $OUTPUT_DIR/stage1_validation/fastqc_summary.csv"');
  scriptLines.push('echo "  - Read counts: $OUTPUT_DIR/stage1_validation/read_counts.csv"');
  scriptLines.push('echo "  - Summary JSON: $OUTPUT_DIR/stage1_validation/stage1_summary.json"');
  scriptLines.push('echo ""');
  scriptLines.push('echo "Ready for agent review."');
  scriptLines.push('');

  const scriptContent = scriptLines.join('\n');

  // Write script to file
  fs.writeFileSync(outputPath, scriptContent, { mode: 0o755 });
  console.log(`[Stage 1] Script saved: ${outputPath}`);

  return {
    scriptPath: outputPath,
    scriptContent
  };
}

/**
 * Parse Stage 1 output files
 *
 * @param {string} outputDir - Output directory containing stage 1 results
 * @returns {Object} - Parsed stage 1 output
 */
export function parseStage1Output(outputDir) {
  console.log('[Stage 1] Parsing output files...');

  const stage1Dir = path.join(outputDir, 'stage1_validation');
  const result = {
    validation_status: 'unknown',
    fastqc_status: 'unknown',
    per_sample: {},
    warnings: [],
    errors: []
  };

  // Parse validation report
  const validationReportPath = path.join(stage1Dir, 'validation_report.tsv');
  if (fs.existsSync(validationReportPath)) {
    const validationContent = fs.readFileSync(validationReportPath, 'utf-8');
    const lines = validationContent.trim().split('\n');
    const header = lines[0].split('\t');

    let allValid = true;
    let hasWarnings = false;

    for (let i = 1; i < lines.length; i++) {
      const values = lines[i].split('\t');
      const sample = values[0];
      const reads = parseInt(values[1]) || 0;
      const status = values[5];

      result.per_sample[sample] = {
        reads,
        lines: parseInt(values[2]) || 0,
        mod4_ok: values[3] === 'OK',
        paired_status: values[4],
        status
      };

      if (status === 'INVALID_FASTQ') {
        allValid = false;
        result.errors.push(`${sample}: Invalid FASTQ format (lines not divisible by 4)`);
      } else if (status === 'MISMATCH') {
        hasWarnings = true;
        result.warnings.push(`${sample}: Paired-end read count mismatch`);
      }
    }

    result.validation_status = allValid ? (hasWarnings ? 'PASS_WITH_WARNING' : 'PASS') : 'FAIL';
    result.validation_report = validationContent;
  } else {
    result.errors.push('Validation report not found');
  }

  // Parse FastQC summary
  const fastqcSummaryPath = path.join(stage1Dir, 'fastqc_summary.csv');
  if (fs.existsSync(fastqcSummaryPath)) {
    const fastqcContent = fs.readFileSync(fastqcSummaryPath, 'utf-8');
    const lines = fastqcContent.trim().split('\n');

    let hasFailures = false;
    let hasWarnings = false;
    const fastqcResults = {};

    for (let i = 1; i < lines.length; i++) {
      const values = lines[i].split(',');
      const sample = values[0];

      fastqcResults[sample] = {
        basic_statistics: values[1],
        per_base_quality: values[2],
        per_sequence_quality: values[3],
        per_base_content: values[4],
        gc_content: values[5],
        per_base_n: values[6],
        sequence_length: values[7],
        duplication: values[8],
        overrepresented: values[9],
        adapter_content: values[10]
      };

      // Check for failures/warnings
      for (let j = 1; j < values.length; j++) {
        if (values[j] === 'FAIL') hasFailures = true;
        if (values[j] === 'WARN') hasWarnings = true;
      }

      // Add to per_sample if exists
      if (result.per_sample[sample]) {
        result.per_sample[sample].fastqc = fastqcResults[sample];
      }
    }

    result.fastqc_status = hasFailures ? 'FAIL' : (hasWarnings ? 'WARN' : 'PASS');
    result.fastqc_summary = fastqcContent;
  } else {
    result.warnings.push('FastQC summary not found - may still be processing');
  }

  // Parse read counts
  const readCountsPath = path.join(stage1Dir, 'read_counts.csv');
  if (fs.existsSync(readCountsPath)) {
    const readCountsContent = fs.readFileSync(readCountsPath, 'utf-8');
    const lines = readCountsContent.trim().split('\n');

    let totalReads = 0;
    let minReads = Infinity;
    let maxReads = 0;
    let sampleCount = 0;

    for (let i = 1; i < lines.length; i++) {
      const values = lines[i].split(',');
      const reads = parseInt(values[1]) || 0;

      totalReads += reads;
      minReads = Math.min(minReads, reads);
      maxReads = Math.max(maxReads, reads);
      sampleCount++;
    }

    result.read_statistics = {
      total_reads: totalReads,
      sample_count: sampleCount,
      avg_reads: sampleCount > 0 ? Math.round(totalReads / sampleCount) : 0,
      min_reads: minReads === Infinity ? 0 : minReads,
      max_reads: maxReads
    };

    // Check for low read count warning
    if (result.read_statistics.min_reads < 5000000) {
      result.warnings.push(`Some samples have low read counts (<5M): minimum = ${result.read_statistics.min_reads.toLocaleString()}`);
    }
  }

  // Overall status
  if (result.errors.length > 0) {
    result.overall_status = 'FAIL';
  } else if (result.warnings.length > 0) {
    result.overall_status = 'PASS_WITH_WARNING';
  } else {
    result.overall_status = 'PASS';
  }

  console.log(`[Stage 1] Parsing complete. Status: ${result.overall_status}`);
  return result;
}

/**
 * Format Stage 1 output for agent review
 *
 * @param {Object} parsedOutput - Parsed stage 1 output from parseStage1Output()
 * @param {Object} dataInfo - Dataset information
 * @returns {string} - Formatted string for agent consumption
 */
export function formatStage1ForAgents(parsedOutput, dataInfo) {
  const lines = [];

  lines.push('## FASTQ Validation Results - Stage 1 Output');
  lines.push('');

  // Dataset overview
  lines.push('### Dataset Overview');
  lines.push(`- Organism: ${dataInfo.organism || 'unknown'}`);
  lines.push(`- Comparison: ${dataInfo.comparison || 'unknown'}`);
  lines.push(`- Total Samples: ${Object.keys(parsedOutput.per_sample).length}`);
  lines.push(`- Sequencing Type: ${dataInfo.pairedEnd ? 'Paired-end (R1 + R2)' : 'Single-end'}`);
  lines.push('');

  // Read count statistics
  if (parsedOutput.read_statistics) {
    lines.push('### Read Count Statistics');
    lines.push(`- Total Reads: ${parsedOutput.read_statistics.total_reads.toLocaleString()}`);
    lines.push(`- Average Reads/Sample: ${parsedOutput.read_statistics.avg_reads.toLocaleString()}`);
    lines.push(`- Min Reads: ${parsedOutput.read_statistics.min_reads.toLocaleString()}`);
    lines.push(`- Max Reads: ${parsedOutput.read_statistics.max_reads.toLocaleString()}`);
    lines.push(`- Read Depth Ratio (max/min): ${(parsedOutput.read_statistics.max_reads / parsedOutput.read_statistics.min_reads).toFixed(2)}x`);
    lines.push('');
  }

  // Validation status
  lines.push('### FASTQ Validation Status');
  lines.push(`- File Integrity: ${parsedOutput.validation_status}`);
  lines.push(`- FastQC Status: ${parsedOutput.fastqc_status}`);
  lines.push(`- Overall: ${parsedOutput.overall_status}`);
  lines.push('');

  // Per-sample results
  lines.push('### Per-Sample Results');
  lines.push('| Sample | Reads | Validation | FastQC Quality |');
  lines.push('|--------|-------|------------|----------------|');

  for (const [sample, data] of Object.entries(parsedOutput.per_sample)) {
    const reads = data.reads ? data.reads.toLocaleString() : 'N/A';
    const validation = data.status || 'unknown';
    const fastqc = data.fastqc?.per_base_quality || 'N/A';
    lines.push(`| ${sample} | ${reads} | ${validation} | ${fastqc} |`);
  }
  lines.push('');

  // Warnings
  if (parsedOutput.warnings.length > 0) {
    lines.push('### Warnings');
    for (const warning of parsedOutput.warnings) {
      lines.push(`- ${warning}`);
    }
    lines.push('');
  }

  // Errors
  if (parsedOutput.errors.length > 0) {
    lines.push('### Errors (CRITICAL)');
    for (const error of parsedOutput.errors) {
      lines.push(`- ${error}`);
    }
    lines.push('');
  }

  // Decision prompt
  lines.push('### Decision Required');
  lines.push('Based on the above results, should we proceed to Stage 2 (Alignment)?');
  lines.push('- PASS: All files valid, quality acceptable - proceed');
  lines.push('- PASS_WITH_WARNING: Minor issues but can proceed - note concerns');
  lines.push('- FAIL: Critical issues - stop and investigate');
  lines.push('');

  return lines.join('\n');
}

/**
 * Execute Stage 1 and return results
 *
 * @param {Object} dataInfo - Dataset information
 * @param {Object} config - Pipeline configuration
 * @returns {Promise<Object>} - { success, scriptPath, output }
 */
export async function executeStage1(dataInfo, config) {
  const scriptPath = path.join(config.output, `stage_1_${config.comparison}.sh`);

  // Generate script
  generateStage1Script(dataInfo, config, scriptPath);

  return {
    scriptPath,
    scriptName: `stage_1_${config.comparison}.sh`,
    stage: 1,
    stageName: 'FASTQ Validation'
  };
}

export default {
  generateStage1Script,
  parseStage1Output,
  formatStage1ForAgents,
  executeStage1
};
