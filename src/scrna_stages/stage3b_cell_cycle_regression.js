/**
 * scRNA-seq Stage 3B: Cell Cycle Regression (Conditional)
 *
 * Either:
 * 1. Regress out cell cycle (REMOVE_CELL_CYCLE) - generates cell_cycle_after.jpg
 * 2. Skip regression (SKIP_CELL_CYCLE) - just scales data
 *
 * Agent decision from Stage 3A determines which path.
 *
 * Input: seurat_stage3a_scored.rds, hvg_genes.csv
 * Output: seurat_stage3_norm.rds (ready for Stage 4 PCA)
 */

import fs from 'fs';
import path from 'path';
import dotenv from 'dotenv';

dotenv.config();

const SCRIPTS_PATH = process.env.SCRIPTS_PATH || './bio_informatics/scripts';

export function generateStage3BScript(stage3AOutput, agentDecision, config, outputPath) {
  console.log('[scRNA Stage 3B] Generating cell cycle regression/skip script...');
  console.log(`[scRNA Stage 3B] Agent decision: ${agentDecision}`);

  const scriptLines = [];

  scriptLines.push('#!/bin/bash');
  scriptLines.push('# scRNA-seq Stage 3B: Cell Cycle Regression or Skip');
  scriptLines.push('set -e');
  scriptLines.push('');

  const stage3ARDS = path.join(stage3AOutput.stage3Dir, 'seurat_stage3a_scored.rds');
  const outputDir = stage3AOutput.stage3Dir;

  scriptLines.push(`RDS_IN="${stage3ARDS}"`);
  scriptLines.push(`OUTPUT_DIR="${outputDir}"`);
  scriptLines.push(`SCRIPTS_DIR="${SCRIPTS_PATH}"`);
  scriptLines.push('');

  if (agentDecision === 'REMOVE_CELL_CYCLE') {
    scriptLines.push('echo "=========================================="');
    scriptLines.push('echo "scRNA-seq Stage 3B: CELL CYCLE REGRESSION"');
    scriptLines.push('echo "=========================================="');
    scriptLines.push('echo "Agents determined: REMOVE cell cycle effects"');
    scriptLines.push('');
    scriptLines.push(`Rscript "$SCRIPTS_DIR/stage3b_cell_cycle_regression.R" "$RDS_IN" "$OUTPUT_DIR"`);
    scriptLines.push('echo "Stage 3B Complete! (Cell cycle regressed)"');
  } else {
    // SKIP_CELL_CYCLE
    scriptLines.push('echo "=========================================="');
    scriptLines.push('echo "scRNA-seq Stage 3B: SKIP CELL CYCLE (Scale Only)"');
    scriptLines.push('echo "=========================================="');
    scriptLines.push('echo "Agents determined: NO significant cell cycle effect"');
    scriptLines.push('');
    scriptLines.push(`Rscript "$SCRIPTS_DIR/stage3b_skip_cell_cycle.R" "$RDS_IN" "$OUTPUT_DIR"`);
    scriptLines.push('echo "Stage 3B Complete! (Scaled without regression)"');
  }

  scriptLines.push('');

  fs.writeFileSync(outputPath, scriptLines.join('\n'), { mode: 0o755 });
  return { scriptPath: outputPath, agentDecision };
}

export function parseStage3BOutput(outputDir, agentDecision) {
  const stage3Dir = path.join(outputDir, 'stage3_normalize_hvg');

  let summary = {};
  let summaryPath;

  if (agentDecision === 'REMOVE_CELL_CYCLE') {
    summaryPath = path.join(stage3Dir, 'regression_summary_stage3b.json');
    if (fs.existsSync(summaryPath)) {
      summary = JSON.parse(fs.readFileSync(summaryPath, 'utf-8'));
      console.log('[scRNA Stage 3B] Cell cycle regression completed');
      console.log(`  PC1 correlation with S.Score (after): ${summary.pc1_correlation_s_score_after?.toFixed(3)}`);
      console.log(`  PC1 correlation with G2M.Score (after): ${summary.pc1_correlation_g2m_score_after?.toFixed(3)}`);
    }
  } else {
    summaryPath = path.join(stage3Dir, 'skip_summary_stage3b.json');
    if (fs.existsSync(summaryPath)) {
      summary = JSON.parse(fs.readFileSync(summaryPath, 'utf-8'));
      console.log('[scRNA Stage 3B] Cell cycle regression skipped (no effect detected)');
    }
  }

  // Check for final RDS file
  const finalRDS = path.join(stage3Dir, 'seurat_stage3_norm.rds');
  if (!fs.existsSync(finalRDS)) {
    throw new Error('Stage 3B final RDS file not found');
  }

  // Check for cell_cycle_after plot (only if regression was done)
  const cellCycleAfterJPG = path.join(stage3Dir, 'cell_cycle_after.jpg');
  const cellCycleAfterPDF = path.join(stage3Dir, 'cell_cycle_after.pdf');
  const cellCycleAfterExists = fs.existsSync(cellCycleAfterJPG);

  return {
    ...summary,
    stage3Dir,
    agentDecision,
    cell_cycle_after_jpg: cellCycleAfterExists ? cellCycleAfterJPG : null,
    cell_cycle_after_pdf: cellCycleAfterExists ? cellCycleAfterPDF : null,
    overall_status: 'PASS'
  };
}

export function formatStage3BForDisplay(parsedOutput) {
  let message = `## scRNA-seq Stage 3B: Cell Cycle Regression\n\n`;

  if (parsedOutput.agentDecision === 'REMOVE_CELL_CYCLE') {
    message += `- Agent Decision: ✓ REMOVE CELL CYCLE\n`;
    message += `- Action Taken: Regressed out S.Score and G2M.Score\n`;
    if (parsedOutput.pc1_correlation_s_score_after !== undefined) {
      message += `- PC1 Correlations (AFTER regression):\n`;
      message += `  - S.Score: ${parsedOutput.pc1_correlation_s_score_after.toFixed(3)}\n`;
      message += `  - G2M.Score: ${parsedOutput.pc1_correlation_g2m_score_after.toFixed(3)}\n`;
    }
    message += `- Plots saved: cell_cycle_after.pdf/.jpg\n`;
    message += `\n→ Cell cycle effects removed. Data ready for PCA (Stage 4).\n`;
  } else {
    message += `- Agent Decision: ✓ SKIP CELL CYCLE\n`;
    message += `- Action Taken: Scaled data without regression\n`;
    message += `- Reason: ${parsedOutput.reason || 'No significant cell cycle effect detected'}\n`;
    message += `\n→ Data scaled (no regression). Ready for PCA (Stage 4).\n`;
  }

  return message;
}

export default {
  generateStage3BScript,
  parseStage3BOutput,
  formatStage3BForDisplay
};
