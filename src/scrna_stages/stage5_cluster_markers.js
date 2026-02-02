/**
 * scRNA-seq Stage 5: Clustering + Marker Discovery
 *
 * Performs graph-based clustering and identifies marker genes.
 * AGENT DECISION REQUIRED: Pipeline agent validates clustering structure.
 *
 * Input: seurat_stage4_pca.rds, selected PCs
 * Output: seurat_stage5_clustered.rds, cluster_summary_stage5.json, markers_stage5.csv
 *
 * Agent Decisions: ACCEPT_CLUSTERING / ADJUST_RESOLUTION / FLAG_SUSPICIOUS
 */

import fs from 'fs';
import path from 'path';
import dotenv from 'dotenv';

dotenv.config();

const SCRIPTS_PATH = process.env.SCRIPTS_PATH || './bio_informatics/scripts';

export function generateStage5Script(stage4Output, pcSelection, config, outputPath) {
  console.log('[scRNA Stage 5] Generating clustering script...');

  const scriptLines = [];
  scriptLines.push('#!/bin/bash');
  scriptLines.push('# scRNA-seq Stage 5: Clustering + Markers');
  scriptLines.push('set -e');
  scriptLines.push('');

  const stage4RDS = path.join(stage4Output.stage4Dir, 'seurat_stage4_pca.rds');
  const minPC = pcSelection.min_pc || 1;
  const maxPC = pcSelection.max_pc || 20;
  const resolution = config.resolution || 0.5;  // Use custom resolution if provided

  scriptLines.push(`RDS_IN="${stage4RDS}"`);
  scriptLines.push(`OUTPUT_DIR="${config.output}/stage5_cluster_markers"`);
  scriptLines.push(`SCRIPTS_DIR="${SCRIPTS_PATH}"`);
  scriptLines.push(`MIN_PC=${minPC}`);
  scriptLines.push(`MAX_PC=${maxPC}`);
  scriptLines.push(`RESOLUTION=${resolution}`);
  scriptLines.push('');
  scriptLines.push('mkdir -p "$OUTPUT_DIR"');
  scriptLines.push('');
  scriptLines.push('echo "=========================================="');
  scriptLines.push('echo "scRNA-seq Stage 5: Clustering + Markers"');
  scriptLines.push('echo "=========================================="');
  scriptLines.push('echo "Using PCs: $MIN_PC-$MAX_PC"');
  scriptLines.push('echo "Resolution: $RESOLUTION"');
  scriptLines.push('Rscript "$SCRIPTS_DIR/stage5_cluster_markers.R" "$RDS_IN" "$OUTPUT_DIR" "$MIN_PC" "$MAX_PC" "$RESOLUTION"');
  scriptLines.push('echo "Stage 5 Complete!"');
  scriptLines.push('echo "Ready for agent review (clustering validation)."');
  scriptLines.push('');

  fs.writeFileSync(outputPath, scriptLines.join('\n'), { mode: 0o755 });
  return { scriptPath: outputPath };
}

export function parseStage5Output(outputDir) {
  const stage5Dir = path.join(outputDir, 'stage5_cluster_markers');
  const summaryPath = path.join(stage5Dir, 'cluster_summary_stage5.json');

  if (!fs.existsSync(summaryPath)) {
    throw new Error('Stage 5 cluster summary not found');
  }

  const summary = JSON.parse(fs.readFileSync(summaryPath, 'utf-8'));
  console.log(`[scRNA Stage 5] Identified ${summary.n_clusters} clusters`);
  console.log(`[scRNA Stage 5] Total markers detected: ${summary.total_markers}`);

  // Check for UMAP plot
  const umapPlotJPG = path.join(stage5Dir, 'umap_plot.jpg');
  const umapPlotPDF = path.join(stage5Dir, 'umap_plot.pdf');

  return {
    ...summary,
    stage5Dir,
    umap_plot_jpg: fs.existsSync(umapPlotJPG) ? umapPlotJPG : null,
    umap_plot_pdf: fs.existsSync(umapPlotPDF) ? umapPlotPDF : null
  };
}

export function formatStage5ForAgents(parsedOutput) {
  const lines = [];

  lines.push('## scRNA-seq Stage 5: Clustering + Marker Discovery');
  lines.push('');
  lines.push('### Clustering Summary');
  lines.push(`- Number of Clusters: ${parsedOutput.n_clusters}`);
  lines.push(`- Total Markers Detected: ${parsedOutput.total_markers}`);
  lines.push(`- Median Markers per Cluster: ${parsedOutput.median_markers_per_cluster}`);
  lines.push(`- Resolution Used: ${parsedOutput.resolution}`);
  lines.push('');
  lines.push('### Cluster Size Distribution');
  for (const [cluster, size] of Object.entries(parsedOutput.cluster_sizes || {})) {
    lines.push(`- Cluster ${cluster}: ${size} cells`);
  }
  lines.push('');
  lines.push('### Decision Required');
  lines.push('Evaluate whether clustering appears technically reasonable.');
  lines.push('');
  lines.push('**Allowed Decisions:**');
  lines.push('- `ACCEPT_CLUSTERING`: Cluster sizes balanced, markers reasonable');
  lines.push('- `ADJUST_RESOLUTION`: Rerun Stage 5 with different resolution');
  lines.push('- `FLAG_SUSPICIOUS`: Clustering appears problematic, halt for review');
  lines.push('');

  return lines.join('\n');
}

export default { generateStage5Script, parseStage5Output, formatStage5ForAgents };
