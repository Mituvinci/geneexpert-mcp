/**
 * scRNA-seq Stage 1: Load + QC Metrics
 *
 * Loads 10x Genomics data and computes QC metrics.
 * No filtering, no agent decision required.
 *
 * Input: filtered_feature_bc_matrix.h5
 * Output: seurat_stage1_raw.rds, qc_metrics_stage1.csv, qc_summary_stage1.json
 */

import fs from 'fs';
import path from 'path';
import dotenv from 'dotenv';

dotenv.config();

const SCRIPTS_PATH = process.env.SCRIPTS_PATH || '/users/ha00014/Halimas_projects/multi_llm_mcp/bio_informatics/scripts';

/**
 * Generate Stage 1 script for loading 10x data and computing QC
 */
export function generateStage1Script(dataInfo, config, outputPath) {
  console.log('[scRNA Stage 1] Generating load + QC script...');

  const scriptLines = [];

  // Header
  scriptLines.push('#!/bin/bash');
  scriptLines.push('#');
  scriptLines.push('# scRNA-seq Stage 1: Load + QC Metrics');
  scriptLines.push(`# Generated: ${new Date().toISOString()}`);
  scriptLines.push('#');
  scriptLines.push('set -e  # Exit on error');
  scriptLines.push('');

  // Find input data (prefer directory > CSV > .h5 to avoid hdf5r dependency)
  const matrixDir = path.join(config.input, 'filtered_feature_bc_matrix');
  const h5File = path.join(config.input, 'filtered_feature_bc_matrix.h5');

  // Check if 10x files are directly in input directory
  let has10xFilesInRoot = false;
  if (fs.existsSync(config.input) && fs.statSync(config.input).isDirectory()) {
    const files = fs.readdirSync(config.input);
    const hasBarcodes = files.some(f => f === 'barcodes.tsv.gz' || f === 'barcodes.tsv');
    const hasFeatures = files.some(f => f === 'features.tsv.gz' || f === 'features.tsv' || f === 'genes.tsv.gz' || f === 'genes.tsv');
    const hasMatrix = files.some(f => f === 'matrix.mtx.gz' || f === 'matrix.mtx');
    has10xFilesInRoot = hasBarcodes && hasFeatures && hasMatrix;
  }

  // Check for CSV/TSV files in directory
  let csvFile = null;
  if (fs.existsSync(config.input) && fs.statSync(config.input).isDirectory()) {
    const files = fs.readdirSync(config.input);
    const dataFile = files.find(f => /\.(csv|txt|tsv)$/i.test(f));
    if (dataFile) {
      csvFile = path.join(config.input, dataFile);
    }
  } else if (fs.existsSync(config.input) && /\.(csv|txt|tsv)$/i.test(config.input)) {
    // Input path IS a CSV file
    csvFile = config.input;
  }

  let inputData;
  if (has10xFilesInRoot) {
    // 10x files directly in input directory
    inputData = config.input;
  } else if (fs.existsSync(matrixDir) && fs.statSync(matrixDir).isDirectory()) {
    // 10x files in filtered_feature_bc_matrix/ subdirectory
    inputData = matrixDir;
  } else if (csvFile) {
    inputData = csvFile;
  } else if (fs.existsSync(h5File)) {
    inputData = h5File;
  } else {
    throw new Error('No scRNA data found. Expected:\n  - barcodes.tsv.gz + features.tsv.gz + matrix.mtx.gz (directly in input dir or in filtered_feature_bc_matrix/ subdirectory)\n  - .csv/.txt/.tsv file (genes=rows, cells=columns)\n  - .h5 file');
  }

  // Variables
  scriptLines.push('# Paths');
  scriptLines.push(`INPUT_DATA="${inputData}"`);
  scriptLines.push(`OUTPUT_DIR="${config.output}/stage1_load_qc"`);
  scriptLines.push(`SCRIPTS_DIR="${SCRIPTS_PATH}"`);
  scriptLines.push('');

  // Create output directory
  scriptLines.push('mkdir -p "$OUTPUT_DIR"');
  scriptLines.push('');

  // Run R script
  scriptLines.push('echo "=========================================="');
  scriptLines.push('echo "scRNA-seq Stage 1: Load + QC Metrics"');
  scriptLines.push('echo "=========================================="');
  scriptLines.push('echo ""');
  scriptLines.push('echo "Loading 10x data and computing QC metrics..."');
  scriptLines.push('');
  scriptLines.push('Rscript "$SCRIPTS_DIR/stage1_load_qc.R" "$INPUT_DATA" "$OUTPUT_DIR"');
  scriptLines.push('');

  // Summary
  scriptLines.push('echo ""');
  scriptLines.push('echo "Stage 1 Complete!"');
  scriptLines.push('echo "Outputs:"');
  scriptLines.push('echo "  - Seurat object: $OUTPUT_DIR/seurat_stage1_raw.rds"');
  scriptLines.push('echo "  - QC metrics: $OUTPUT_DIR/qc_metrics_stage1.csv"');
  scriptLines.push('echo "  - QC summary: $OUTPUT_DIR/qc_summary_stage1.json"');
  scriptLines.push('');

  const scriptContent = scriptLines.join('\n');
  fs.writeFileSync(outputPath, scriptContent, { mode: 0o755 });
  console.log(`[scRNA Stage 1] Script saved: ${outputPath}`);

  return { scriptPath: outputPath, scriptContent };
}

/**
 * Parse Stage 1 output (QC summary)
 */
export function parseStage1Output(outputDir) {
  console.log('[scRNA Stage 1] Parsing output...');

  const stage1Dir = path.join(outputDir, 'stage1_load_qc');
  const summaryPath = path.join(stage1Dir, 'qc_summary_stage1.json');

  if (!fs.existsSync(summaryPath)) {
    throw new Error('Stage 1 QC summary not found');
  }

  const summary = JSON.parse(fs.readFileSync(summaryPath, 'utf-8'));

  console.log(`[scRNA Stage 1] Loaded ${summary.cells_total} cells, ${summary.genes_total} genes`);
  console.log(`[scRNA Stage 1] Median nFeature: ${summary.nFeature_median}, Median %MT: ${summary.percent_mt_median.toFixed(2)}`);

  return {
    ...summary,
    overall_status: 'PASS',  // Stage 1 always passes (no agent decision)
    stage1Dir
  };
}

/**
 * Format Stage 1 output for display (no agent needed for Stage 1)
 */
export function formatStage1ForDisplay(parsedOutput) {
  const lines = [];

  lines.push('## scRNA-seq Stage 1: Load + QC Metrics');
  lines.push('');
  lines.push('### Dataset Overview');
  lines.push(`- Total Cells: ${parsedOutput.cells_total.toLocaleString()}`);
  lines.push(`- Total Genes: ${parsedOutput.genes_total.toLocaleString()}`);
  lines.push(`- Median Features/Cell: ${parsedOutput.nFeature_median.toLocaleString()}`);
  lines.push(`- Median UMI Count/Cell: ${parsedOutput.nCount_median.toLocaleString()}`);
  lines.push(`- Median % Mitochondrial: ${parsedOutput.percent_mt_median.toFixed(2)}%`);
  lines.push('');
  lines.push('✓ Stage 1 complete. Proceeding to Stage 2 (QC filtering).');
  lines.push('');

  return lines.join('\n');
}

export default {
  generateStage1Script,
  parseStage1Output,
  formatStage1ForDisplay
};
