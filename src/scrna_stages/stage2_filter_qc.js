/**
 * scRNA-seq Stage 2: QC Filtering
 *
 * Applies QC thresholds to filter low-quality cells.
 * AGENT DECISION REQUIRED: Stats agent reviews filtering results.
 *
 * Input: seurat_stage1_raw.rds
 * Output: seurat_stage2_filtered.rds, qc_summary_stage2.json
 *
 * Agent Decisions: PROCEED / PROCEED_WITH_WARNING / STOP_AND_REVIEW
 */

import fs from 'fs';
import path from 'path';
import dotenv from 'dotenv';

dotenv.config();

const SCRIPTS_PATH = process.env.SCRIPTS_PATH || '/users/ha00014/Halimas_projects/multi_llm_mcp/bio_informatics/scripts';

/**
 * Generate Stage 2 script for QC filtering
 * @param {Object} stage1Output - Stage 1 output with QC metrics
 * @param {Object} thresholds - Agent-recommended thresholds
 * @param {Object} config - Configuration
 * @param {string} outputPath - Output script path
 */
export function generateStage2Script(stage1Output, thresholds, config, outputPath) {
  console.log('[scRNA Stage 2] Generating QC filtering script...');

  const scriptLines = [];

  // Header
  scriptLines.push('#!/bin/bash');
  scriptLines.push('#');
  scriptLines.push('# scRNA-seq Stage 2: QC Filtering');
  scriptLines.push(`# Generated: ${new Date().toISOString()}`);
  scriptLines.push('#');
  scriptLines.push('# Agent-Recommended Thresholds:');
  scriptLines.push(`#   nFeature: ${thresholds.nFeature_min} - ${thresholds.nFeature_max}`);
  scriptLines.push(`#   %MT: < ${thresholds.percent_mt_max}%`);
  scriptLines.push('#');
  scriptLines.push('set -e');
  scriptLines.push('');

  // Variables
  const stage1RDS = path.join(stage1Output.stage1Dir, 'seurat_stage1_raw.rds');

  scriptLines.push('# Paths');
  scriptLines.push(`RDS_IN="${stage1RDS}"`);
  scriptLines.push(`OUTPUT_DIR="${config.output}/stage2_filter_qc"`);
  scriptLines.push(`SCRIPTS_DIR="${SCRIPTS_PATH}"`);
  scriptLines.push('');

  // Threshold parameters
  scriptLines.push('# Agent-recommended QC thresholds');
  scriptLines.push(`NFEATURE_MIN=${thresholds.nFeature_min}`);
  scriptLines.push(`NFEATURE_MAX=${thresholds.nFeature_max}`);
  scriptLines.push(`PERCENT_MT_MAX=${thresholds.percent_mt_max}`);
  scriptLines.push('');

  // Create output directory
  scriptLines.push('mkdir -p "$OUTPUT_DIR"');
  scriptLines.push('');

  // Run R script with threshold parameters
  scriptLines.push('echo "=========================================="');
  scriptLines.push('echo "scRNA-seq Stage 2: QC Filtering"');
  scriptLines.push('echo "=========================================="');
  scriptLines.push('echo ""');
  scriptLines.push('echo "Agent-recommended thresholds:"');
  scriptLines.push('echo "  nFeature: $NFEATURE_MIN - $NFEATURE_MAX"');
  scriptLines.push('echo "  %MT: < $PERCENT_MT_MAX%"');
  scriptLines.push('echo ""');
  scriptLines.push('echo "Applying QC filters..."');
  scriptLines.push('');
  scriptLines.push('Rscript "$SCRIPTS_DIR/stage2_filter_qc.R" "$RDS_IN" "$OUTPUT_DIR" "$NFEATURE_MIN" "$NFEATURE_MAX" "$PERCENT_MT_MAX"');
  scriptLines.push('');

  // Summary
  scriptLines.push('echo ""');
  scriptLines.push('echo "Stage 2 Complete!"');
  scriptLines.push('echo "Outputs:"');
  scriptLines.push('echo "  - Filtered Seurat object: $OUTPUT_DIR/seurat_stage2_filtered.rds"');
  scriptLines.push('echo "  - QC summary: $OUTPUT_DIR/qc_summary_stage2.json"');
  scriptLines.push('echo ""');
  scriptLines.push('echo "Ready for agent review."');
  scriptLines.push('');

  const scriptContent = scriptLines.join('\n');
  fs.writeFileSync(outputPath, scriptContent, { mode: 0o755 });
  console.log(`[scRNA Stage 2] Script saved: ${outputPath}`);

  return { scriptPath: outputPath, scriptContent };
}

/**
 * Parse Stage 2 output (QC filtering summary)
 */
export function parseStage2Output(outputDir) {
  console.log('[scRNA Stage 2] Parsing output...');

  const stage2Dir = path.join(outputDir, 'stage2_filter_qc');
  const summaryPath = path.join(stage2Dir, 'qc_summary_stage2.json');

  if (!fs.existsSync(summaryPath)) {
    throw new Error('Stage 2 QC summary not found');
  }

  const summary = JSON.parse(fs.readFileSync(summaryPath, 'utf-8'));

  console.log(`[scRNA Stage 2] Filtered ${summary.cells_before} → ${summary.cells_after} cells (${summary.percent_removed}% removed)`);

  return {
    ...summary,
    stage2Dir,
    // Initial status - will be overridden by agent
    overall_status: summary.percent_removed > 50 ? 'WARN' : 'PASS'
  };
}

/**
 * Format Stage 2 output for agent review
 */
export function formatStage2ForAgents(parsedOutput) {
  const lines = [];

  lines.push('## scRNA-seq Stage 2: QC Filtering Results');
  lines.push('');
  lines.push('### Filtering Summary');
  lines.push(`- Cells Before: ${parsedOutput.cells_before.toLocaleString()}`);
  lines.push(`- Cells After: ${parsedOutput.cells_after.toLocaleString()}`);
  lines.push(`- Percent Removed: ${parsedOutput.percent_removed}%`);
  lines.push('');
  lines.push('### Applied Thresholds');
  lines.push(`- Min Features: ${parsedOutput.thresholds.min_features}`);
  lines.push(`- Max Features: ${parsedOutput.thresholds.max_features}`);
  lines.push(`- Max % Mitochondrial: ${parsedOutput.thresholds.max_percent_mt}%`);
  lines.push('');
  lines.push('### Decision Required');
  lines.push('Evaluate whether QC filtering is reasonable for typical 10x Genomics scRNA-seq data.');
  lines.push('');
  lines.push('**Allowed Decisions:**');
  lines.push('- `PROCEED`: Filtering is reasonable, continue to normalization');
  lines.push('- `PROCEED_WITH_WARNING`: Filtering is acceptable but note concerns');
  lines.push('- `STOP_AND_REVIEW`: Filtering appears problematic, halt for manual review');
  lines.push('');

  return lines.join('\n');
}

export default {
  generateStage2Script,
  parseStage2Output,
  formatStage2ForAgents
};
