/**
 * scRNA-seq Stage 4: PCA + PC Selection
 *
 * Runs PCA and generates elbow plot.
 * AGENT DECISION REQUIRED: Stats agent selects number of PCs.
 *
 * Input: seurat_stage3_normalized.rds
 * Output: seurat_stage4_pca.rds, pca_summary_stage4.json
 *
 * Agent Decisions: USE_DEFAULT / SELECT_PC_RANGE / STOP_AND_REVIEW
 */

import fs from 'fs';
import path from 'path';
import dotenv from 'dotenv';

dotenv.config();

const SCRIPTS_PATH = process.env.SCRIPTS_PATH || '/users/ha00014/Halimas_projects/multi_llm_mcp/bio_informatics/scripts';

export function generateStage4Script(stage3Output, config, outputPath) {
  console.log('[scRNA Stage 4] Generating PCA script...');

  const scriptLines = [];
  scriptLines.push('#!/bin/bash');
  scriptLines.push('# scRNA-seq Stage 4: PCA');
  scriptLines.push('set -e');
  scriptLines.push('');

  const stage3RDS = path.join(stage3Output.stage3Dir, 'seurat_stage3_norm.rds');

  scriptLines.push(`RDS_IN="${stage3RDS}"`);
  scriptLines.push(`OUTPUT_DIR="${config.output}/stage4_pca"`);
  scriptLines.push(`SCRIPTS_DIR="${SCRIPTS_PATH}"`);
  scriptLines.push('');
  scriptLines.push('mkdir -p "$OUTPUT_DIR"');
  scriptLines.push('');
  scriptLines.push('echo "=========================================="');
  scriptLines.push('echo "scRNA-seq Stage 4: PCA"');
  scriptLines.push('echo "=========================================="');
  scriptLines.push('Rscript "$SCRIPTS_DIR/stage4_pca.R" "$RDS_IN" "$OUTPUT_DIR"');
  scriptLines.push('echo "Stage 4 Complete!"');
  scriptLines.push('echo "Ready for agent review (PC selection)."');
  scriptLines.push('');

  fs.writeFileSync(outputPath, scriptLines.join('\n'), { mode: 0o755 });
  return { scriptPath: outputPath };
}

export function parseStage4Output(outputDir) {
  const stage4Dir = path.join(outputDir, 'stage4_pca');
  const summaryPath = path.join(stage4Dir, 'pca_variance.json');

  if (!fs.existsSync(summaryPath)) {
    throw new Error('Stage 4 PCA summary not found');
  }

  const summary = JSON.parse(fs.readFileSync(summaryPath, 'utf-8'));
  console.log(`[scRNA Stage 4] PCA complete. Variance explained by PC 20: ${summary.variance_explained_pc20}%`);

  return { ...summary, stage4Dir };
}

export function formatStage4ForAgents(parsedOutput) {
  const lines = [];

  lines.push('## scRNA-seq Stage 4: PCA Results');
  lines.push('');
  lines.push('### PCA Variance Summary');
  lines.push(`- Total PCs Computed: ${parsedOutput.total_pcs}`);
  lines.push(`- Variance Explained (PC 1-10): ${parsedOutput.variance_explained_pc10}%`);
  lines.push(`- Variance Explained (PC 1-20): ${parsedOutput.variance_explained_pc20}%`);
  lines.push(`- Variance Explained (PC 1-30): ${parsedOutput.variance_explained_pc30}%`);
  lines.push('');
  lines.push('### Decision Required');
  lines.push('Based on the variance explained, select an appropriate PC range for downstream analysis.');
  lines.push('');
  lines.push('**Allowed Decisions:**');
  lines.push('- `USE_DEFAULT`: Use PCs 1-20 (standard for 10x data)');
  lines.push('- `SELECT_PC_RANGE`: Specify min_pc and max_pc (e.g., 1-15 or 1-30)');
  lines.push('- `STOP_AND_REVIEW`: PCA results appear problematic');
  lines.push('');

  return lines.join('\n');
}

export default { generateStage4Script, parseStage4Output, formatStage4ForAgents };
