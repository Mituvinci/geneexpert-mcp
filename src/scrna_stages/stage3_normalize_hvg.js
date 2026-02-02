/**
 * scRNA-seq Stage 3: Normalization + HVGs
 *
 * Normalizes data, identifies highly variable genes, and scales.
 * No agent decision required (deterministic).
 *
 * Input: seurat_stage2_filtered.rds
 * Output: seurat_stage3_normalized.rds, hvg_summary_stage3.json
 */

import fs from 'fs';
import path from 'path';
import dotenv from 'dotenv';

dotenv.config();

const SCRIPTS_PATH = process.env.SCRIPTS_PATH || './bio_informatics/scripts';

export function generateStage3Script(stage2Output, config, outputPath) {
  console.log('[scRNA Stage 3] Generating normalization + HVG script...');

  const scriptLines = [];

  scriptLines.push('#!/bin/bash');
  scriptLines.push('# scRNA-seq Stage 3: Normalization + HVGs');
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
  scriptLines.push('echo "scRNA-seq Stage 3: Normalization + HVGs + CELL CYCLE ANALYSIS"');
  scriptLines.push('echo "=========================================="');
  scriptLines.push(`Rscript "$SCRIPTS_DIR/stage3_normalize_hvg.R" "$RDS_IN" "$OUTPUT_DIR" "${config.organism}"`);
  scriptLines.push('echo "Stage 3 Complete!"');
  scriptLines.push('');

  fs.writeFileSync(outputPath, scriptLines.join('\n'), { mode: 0o755 });
  return { scriptPath: outputPath };
}

export function parseStage3Output(outputDir) {
  const stage3Dir = path.join(outputDir, 'stage3_normalize_hvg');
  const summaryPath = path.join(stage3Dir, 'hvg_summary_stage3.json');

  if (!fs.existsSync(summaryPath)) {
    throw new Error('Stage 3 HVG summary not found');
  }

  const summary = JSON.parse(fs.readFileSync(summaryPath, 'utf-8'));
  console.log(`[scRNA Stage 3] Identified ${summary.n_hvgs} highly variable genes`);

  // Check for cell cycle analysis outputs
  const cellCycleSummaryPath = path.join(stage3Dir, 'cell_cycle_summary.txt');
  const cellCycleBeforePlot = path.join(stage3Dir, 'cell_cycle_before.pdf');
  const cellCycleAfterPlot = path.join(stage3Dir, 'cell_cycle_after.pdf');

  const cellCycleAnalyzed = fs.existsSync(cellCycleSummaryPath) &&
                             fs.existsSync(cellCycleBeforePlot) &&
                             fs.existsSync(cellCycleAfterPlot);

  if (cellCycleAnalyzed) {
    console.log('[scRNA Stage 3] Cell cycle analysis completed (plots saved)');
  } else {
    console.log('[scRNA Stage 3] Cell cycle analysis skipped (markers not detected or data too small)');
  }

  return {
    ...summary,
    stage3Dir,
    overall_status: 'PASS',
    cell_cycle_analyzed: cellCycleAnalyzed
  };
}

export function formatStage3ForDisplay(parsedOutput) {
  let message = `## scRNA-seq Stage 3: Normalization + HVGs + CELL CYCLE ANALYSIS\n\n`;
  message += `- HVGs Identified: ${parsedOutput.n_hvgs}\n`;

  if (parsedOutput.cell_cycle_analyzed) {
    message += `- Cell Cycle Analysis: ✓ COMPLETED\n`;
    message += `  - Phases detected (S, G2M)\n`;
    message += `  - Plots saved: cell_cycle_before.pdf, cell_cycle_after.pdf\n`;
    message += `  - Cell cycle effects removed from data\n`;
  } else {
    message += `- Cell Cycle Analysis: ⚠ SKIPPED (markers not detected or data too small)\n`;
  }

  message += `\n✓ Stage 3 complete. Proceeding to PCA.\n`;
  return message;
}

export default { generateStage3Script, parseStage3Output, formatStage3ForDisplay };
