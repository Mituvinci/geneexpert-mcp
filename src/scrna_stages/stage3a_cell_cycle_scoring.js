/**
 * scRNA-seq Stage 3A: Cell Cycle Scoring
 *
 * Normalizes data, identifies HVGs, scores cell cycle phases.
 * Generates cell_cycle_before.jpg for agent review.
 *
 * Agent decision required: REMOVE_CELL_CYCLE or SKIP_CELL_CYCLE
 *
 * Input: seurat_stage2_filtered.rds
 * Output: seurat_stage3a_scored.rds, cell_cycle_before.jpg, cell_cycle_summary_stage3a.json
 */

import fs from 'fs';
import path from 'path';
import dotenv from 'dotenv';

dotenv.config();

const SCRIPTS_PATH = process.env.SCRIPTS_PATH || './bio_informatics/scripts';

export function generateStage3AScript(stage2Output, config, outputPath) {
  console.log('[scRNA Stage 3A] Generating cell cycle scoring script...');

  const scriptLines = [];

  scriptLines.push('#!/bin/bash');
  scriptLines.push('# scRNA-seq Stage 3A: Cell Cycle Scoring');
  scriptLines.push('set -e');
  scriptLines.push('');

  const stage2RDS = path.join(stage2Output.stage2Dir, 'seurat_stage2_filtered.rds');

  scriptLines.push(`RDS_IN="${stage2RDS}"`);
  scriptLines.push(`OUTPUT_DIR="${config.output}/stage3_normalize_hvg"`);
  scriptLines.push(`SCRIPTS_DIR="${SCRIPTS_PATH}"`);
  scriptLines.push('');
  scriptLines.push('mkdir -p "$OUTPUT_DIR"');
  scriptLines.push('');
  scriptLines.push('echo "=========================================="');
  scriptLines.push('echo "scRNA-seq Stage 3A: Cell Cycle Scoring"');
  scriptLines.push('echo "=========================================="');
  scriptLines.push(`Rscript "$SCRIPTS_DIR/stage3a_cell_cycle_scoring.R" "$RDS_IN" "$OUTPUT_DIR" "${config.organism}"`);
  scriptLines.push('echo "Stage 3A Complete! (Cell cycle scored - ready for agent review)"');
  scriptLines.push('');

  fs.writeFileSync(outputPath, scriptLines.join('\n'), { mode: 0o755 });
  return { scriptPath: outputPath };
}

export function parseStage3AOutput(outputDir) {
  const stage3Dir = path.join(outputDir, 'stage3_normalize_hvg');
  const summaryPath = path.join(stage3Dir, 'cell_cycle_summary_stage3a.json');

  if (!fs.existsSync(summaryPath)) {
    throw new Error('Stage 3A cell cycle summary not found');
  }

  const summary = JSON.parse(fs.readFileSync(summaryPath, 'utf-8'));
  console.log(`[scRNA Stage 3A] Identified ${summary.n_hvgs} highly variable genes`);

  // Check for cell cycle plots (separate Phase and Scores plots)
  const cellCyclePlotPhaseJPG = path.join(stage3Dir, 'cell_cycle_phase_before.jpg');
  const cellCyclePlotScoresJPG = path.join(stage3Dir, 'cell_cycle_scores_before.jpg');
  const cellCyclePlotPDF = path.join(stage3Dir, 'cell_cycle_before.pdf');

  const cellCyclePlotPhaseExists = fs.existsSync(cellCyclePlotPhaseJPG);
  const cellCyclePlotScoresExists = fs.existsSync(cellCyclePlotScoresJPG);

  if (summary.cell_cycle_detected) {
    console.log('[scRNA Stage 3A] Cell cycle phases detected:');
    if (summary.phase_distribution) {
      console.log(`  G1: ${summary.phase_distribution.G1 || 0} cells`);
      console.log(`  S: ${summary.phase_distribution.S || 0} cells`);
      console.log(`  G2M: ${summary.phase_distribution.G2M || 0} cells`);
    }
    console.log(`  PC1 correlation with S.Score: ${summary.pc1_correlation_s_score?.toFixed(3)}`);
    console.log(`  PC1 correlation with G2M.Score: ${summary.pc1_correlation_g2m_score?.toFixed(3)}`);
    console.log(`  Recommendation: ${summary.recommendation}`);
  } else {
    console.log('[scRNA Stage 3A] Cell cycle not detected (insufficient markers)');
    console.log(`  Recommendation: ${summary.recommendation || 'SKIP_CELL_CYCLE'}`);
  }

  return {
    ...summary,
    stage3Dir,
    cell_cycle_plot_phase_jpg: cellCyclePlotPhaseExists ? cellCyclePlotPhaseJPG : null,
    cell_cycle_plot_scores_jpg: cellCyclePlotScoresExists ? cellCyclePlotScoresJPG : null,
    cell_cycle_plot_pdf: fs.existsSync(cellCyclePlotPDF) ? cellCyclePlotPDF : null,
    overall_status: 'PENDING_AGENT_REVIEW'
  };
}

export function formatStage3AForAgents(parsedOutput) {
  let message = `## scRNA-seq Stage 3A: Cell Cycle Scoring\n\n`;

  message += `**Dataset Summary:**\n`;
  message += `- Highly Variable Genes (HVGs): ${parsedOutput.n_hvgs}\n`;
  message += `- Organism: ${parsedOutput.organism || 'Not specified'}\n\n`;

  if (parsedOutput.cell_cycle_detected) {
    message += `**Cell Cycle Detection: ✓ DETECTED**\n\n`;

    message += `**Phase Distribution:**\n`;
    if (parsedOutput.phase_distribution) {
      const total = Object.values(parsedOutput.phase_distribution).reduce((a, b) => a + b, 0);
      message += `- G1 phase: ${parsedOutput.phase_distribution.G1 || 0} cells (${((parsedOutput.phase_distribution.G1 || 0) / total * 100).toFixed(1)}%)\n`;
      message += `- S phase: ${parsedOutput.phase_distribution.S || 0} cells (${((parsedOutput.phase_distribution.S || 0) / total * 100).toFixed(1)}%)\n`;
      message += `- G2M phase: ${parsedOutput.phase_distribution.G2M || 0} cells (${((parsedOutput.phase_distribution.G2M || 0) / total * 100).toFixed(1)}%)\n\n`;
    }

    message += `**Cell Cycle Markers:**\n`;
    message += `- S phase markers present: ${parsedOutput.s_markers_present}/${parsedOutput.s_markers_present || 43}\n`;
    message += `- G2M phase markers present: ${parsedOutput.g2m_markers_present}/${parsedOutput.g2m_markers_present || 54}\n\n`;

    message += `**PCA Correlation Analysis:**\n`;
    message += `- PC1 correlation with S.Score: ${parsedOutput.pc1_correlation_s_score?.toFixed(3)}\n`;
    message += `- PC1 correlation with G2M.Score: ${parsedOutput.pc1_correlation_g2m_score?.toFixed(3)}\n\n`;

    message += `**Visual Evidence:**\n`;
    message += `- See attached: cell_cycle_before.jpg\n`;
    message += `  - PCA plot colored by cell cycle phase (G1, S, G2M)\n`;
    message += `  - Feature plots showing S.Score and G2M.Score distribution\n\n`;

    message += `**Interpretation Guide:**\n`;
    message += `- If PCA shows strong clustering by cell cycle phase (G1/S/G2M clearly separated):\n`;
    message += `  → Cell cycle is a confounding factor\n`;
    message += `  → Recommend: REMOVE_CELL_CYCLE\n`;
    message += `- If PCA shows mixed phases (no clear separation):\n`;
    message += `  → Cell cycle is NOT driving variance\n`;
    message += `  → Recommend: SKIP_CELL_CYCLE\n`;
    message += `- Correlation threshold: |r| > 0.3 suggests cell cycle effect\n\n`;

    message += `**Decision Required:**\n`;
    message += `Based on the cell_cycle_before.jpg plot and correlation analysis, should we:\n`;
    message += `1. REMOVE_CELL_CYCLE - Regress out S.Score and G2M.Score (proceed to Stage 3B)\n`;
    message += `2. SKIP_CELL_CYCLE - Skip regression, just scale data (skip Stage 3B)\n`;
  } else {
    message += `**Cell Cycle Detection: ✗ NOT DETECTED**\n\n`;
    message += `**Reason:** ${parsedOutput.reason}\n\n`;
    message += `**Recommendation:** SKIP_CELL_CYCLE (insufficient markers for regression)\n`;
  }

  return message;
}

export function formatStage3AForDisplay(parsedOutput) {
  let message = `## scRNA-seq Stage 3A: Cell Cycle Scoring\n\n`;
  message += `- HVGs Identified: ${parsedOutput.n_hvgs}\n`;

  if (parsedOutput.cell_cycle_detected) {
    message += `- Cell Cycle Phases Detected: ✓ YES\n`;
    if (parsedOutput.phase_distribution) {
      message += `  - G1: ${parsedOutput.phase_distribution.G1 || 0} cells\n`;
      message += `  - S: ${parsedOutput.phase_distribution.S || 0} cells\n`;
      message += `  - G2M: ${parsedOutput.phase_distribution.G2M || 0} cells\n`;
    }
    message += `- PC1 Correlations:\n`;
    message += `  - S.Score: ${parsedOutput.pc1_correlation_s_score?.toFixed(3)}\n`;
    message += `  - G2M.Score: ${parsedOutput.pc1_correlation_g2m_score?.toFixed(3)}\n`;
    message += `- Plots saved: cell_cycle_before.pdf/.jpg\n`;
  } else {
    message += `- Cell Cycle Detection: ⚠ SKIPPED (${parsedOutput.reason})\n`;
  }

  message += `\n→ Stage 3A complete. Agents will review cell cycle plot and decide next step.\n`;
  return message;
}

export default {
  generateStage3AScript,
  parseStage3AOutput,
  formatStage3AForAgents,
  formatStage3AForDisplay
};
