/**
 * Stage 3c: Clustering + UMAP Visualization (scRNA-seq)
 *
 * This stage:
 * 1. Loads Seurat object from Stage 3b (with regressed data)
 * 2. Finds neighbors and clusters cells
 * 3. Runs UMAP dimensionality reduction
 * 4. Generates UMAP plots (colored by cluster and by Phase)
 * 5. Finds marker genes for each cluster
 *
 * Outputs:
 * - Seurat object with clusters + UMAP (.rds)
 * - UMAP plots (clusters, cell cycle phase) - PDF + JPG
 * - Cluster markers CSV
 * - Clustering summary (cluster sizes, top markers)
 *
 * Agent decision: APPROVE_CLUSTERS / ADJUST_RESOLUTION / REQUEST_REANALYSIS
 *
 * Generates script → Parses output → Formats for agents
 */

import fs from 'fs';
import path from 'path';

/**
 * Generate Stage 3c script: Clustering + UMAP
 *
 * @param {Object} dataInfo - Dataset information (samples, organism, etc.)
 * @param {Object} config - Configuration (paths, comparison name, etc.)
 * @param {string} outputPath - Path to write the generated script
 * @param {Object} stage3bResult - Results from Stage 3b (for RDS file location)
 * @returns {Object} - { scriptPath, scriptContent }
 */
export function generateStage3cScript(dataInfo, config, outputPath, stage3bResult = {}) {
  const scriptsPath = process.env.SCRIPTS_PATH || '/users/ha00014/bin';
  const condaEnv = process.env.CONDA_ENV || 'singlecell';

  // Paths
  const outputDir = path.dirname(outputPath);
  const stage3cDir = path.join(outputDir, 'stage3c_clustering_umap');
  const stage3bDir = path.join(outputDir, 'stage3b_cell_cycle_regression');

  // Comparison name
  const comparison = config.comparison || 'comparison';

  // Clustering resolution (default 0.8, can be adjusted)
  const resolution = config.clusteringResolution || 0.8;

  // Seurat RDS path from Stage 3b
  const seuratRdsPath = stage3bResult.seurat_rds_path || path.join(stage3bDir, `${comparison}_seurat_regressed.rds`);

  // Build script
  let script = `#!/bin/bash
#
# Stage 3c: Clustering + UMAP Visualization
# Generated: ${new Date().toISOString()}
#

set -e  # Exit on error

echo "========================================"
echo "Stage 3c: Clustering + UMAP"
echo "========================================"
echo ""

# Environment
SCRIPTS_PATH="${scriptsPath}"
CONDA_ENV="${condaEnv}"
OUTPUT_DIR="${stage3cDir}"
COMPARISON="${comparison}"
RESOLUTION="${resolution}"
SEURAT_RDS="${seuratRdsPath}"

# Activate conda environment
echo "[Stage 3c] Activating conda environment: \${CONDA_ENV}"
set +u
source $(conda info --base)/etc/profile.d/conda.sh
conda activate \${CONDA_ENV}
set -u

# Create output directory
mkdir -p "\${OUTPUT_DIR}"
cd "\${OUTPUT_DIR}"

echo ""
echo "Running Clustering + UMAP..."
echo "========================================"
echo ""

# Run seurat_clustering_umap.R
Rscript \${SCRIPTS_PATH}/seurat_clustering_umap.R \\
  "\${SEURAT_RDS}" \\
  "\${RESOLUTION}" \\
  "\${COMPARISON}"

if [ $? -ne 0 ]; then
  echo "ERROR: Clustering + UMAP failed"
  exit 1
fi

echo ""
echo "========================================"
echo "Stage 3c Complete"
echo "========================================"
echo ""
echo "Output files:"
echo "  - \${COMPARISON}_seurat_clustered.rds (Seurat object with clusters + UMAP)"
echo "  - \${COMPARISON}_umap_clusters.pdf (UMAP by cluster - user)"
echo "  - \${COMPARISON}_umap_clusters.jpg (UMAP by cluster - agents)"
echo "  - \${COMPARISON}_umap_phase.pdf (UMAP by phase - user)"
echo "  - \${COMPARISON}_umap_phase.jpg (UMAP by phase - agents)"
echo "  - \${COMPARISON}_cluster_markers.csv (marker genes)"
echo "  - \${COMPARISON}_clustering_summary.txt (summary statistics)"
echo ""
echo "Next: Stage 4 (Differential Expression Analysis)"
echo ""

exit 0
`;

  // Write script to file
  fs.writeFileSync(outputPath, script, { mode: 0o755 });

  return {
    scriptPath: outputPath,
    scriptContent: script
  };
}

/**
 * Parse Stage 3c output: Read clustering + UMAP results
 *
 * @param {string} outputDir - Output directory containing stage3c_clustering_umap/
 * @param {string} comparison - Comparison name
 * @returns {Object} - Parsed clustering results
 */
export function parseStage3cOutput(outputDir, comparison) {
  const stage3cDir = path.join(outputDir, 'stage3c_clustering_umap');

  // Initialize result structure
  const result = {
    overall_status: 'UNKNOWN',
    exit_code: 0,
    cells_processed: 0,
    n_clusters: 0,
    cluster_sizes: {},
    cluster_percentages: {},
    phase_by_cluster: {},
    top_markers: {},
    umap_clusters_plot_path: null,
    umap_phase_plot_path: null,
    seurat_rds_path: null,
    markers_csv_path: null,
    warnings: [],
    errors: []
  };

  try {
    // 1. Check if output directory exists
    if (!fs.existsSync(stage3cDir)) {
      result.errors.push(`Stage 3c output directory not found: ${stage3cDir}`);
      result.overall_status = 'ERROR';
      return result;
    }

    // 2. Parse summary file
    const summaryPath = path.join(stage3cDir, `${comparison}_clustering_summary.txt`);
    if (!fs.existsSync(summaryPath)) {
      result.errors.push(`Summary file not found: ${summaryPath}`);
      result.overall_status = 'ERROR';
      return result;
    }

    const summaryContent = fs.readFileSync(summaryPath, 'utf-8');

    // Extract cells processed
    const cellsMatch = summaryContent.match(/Cells processed:\s*(\d+)/);
    if (cellsMatch) {
      result.cells_processed = parseInt(cellsMatch[1]);
    }

    // Extract number of clusters
    const nClustersMatch = summaryContent.match(/Number of clusters:\s*(\d+)/);
    if (nClustersMatch) {
      result.n_clusters = parseInt(nClustersMatch[1]);
    }

    // Extract cluster sizes
    const clusterSizesSection = summaryContent.match(/Cluster Sizes:([\s\S]*?)Cell Cycle Phase Distribution/);
    if (clusterSizesSection) {
      const clusterLines = clusterSizesSection[1].match(/Cluster\s+(\d+):\s+(\d+)\s+cells\s+\(([\d.]+)%\)/g);
      if (clusterLines) {
        clusterLines.forEach(line => {
          const match = line.match(/Cluster\s+(\d+):\s+(\d+)\s+cells\s+\(([\d.]+)%\)/);
          if (match) {
            const clusterNum = match[1];
            result.cluster_sizes[clusterNum] = parseInt(match[2]);
            result.cluster_percentages[clusterNum] = parseFloat(match[3]);
          }
        });
      }
    }

    // Extract phase distribution per cluster
    const phaseSection = summaryContent.match(/Cell Cycle Phase Distribution per Cluster:([\s\S]*?)Top Marker Genes/);
    if (phaseSection) {
      const phaseLines = phaseSection[1].match(/Cluster\s+(\d+):\s+G1=(\d+),\s+S=(\d+),\s+G2M=(\d+)/g);
      if (phaseLines) {
        phaseLines.forEach(line => {
          const match = line.match(/Cluster\s+(\d+):\s+G1=(\d+),\s+S=(\d+),\s+G2M=(\d+)/);
          if (match) {
            const clusterNum = match[1];
            result.phase_by_cluster[clusterNum] = {
              G1: parseInt(match[2]),
              S: parseInt(match[3]),
              G2M: parseInt(match[4])
            };
          }
        });
      }
    }

    // Extract top marker genes per cluster
    const markersSection = summaryContent.match(/Top Marker Genes per Cluster \(Top 5\):([\s\S]*?)Output Files:/);
    if (markersSection) {
      const markerLines = markersSection[1].match(/Cluster\s+(\d+):\s+(.*)/g);
      if (markerLines) {
        markerLines.forEach(line => {
          const match = line.match(/Cluster\s+(\d+):\s+(.*)/);
          if (match) {
            const clusterNum = match[1];
            const genes = match[2].trim().split(',').map(g => g.trim());
            result.top_markers[clusterNum] = genes;
          }
        });
      }
    }

    // 3. Check for output files
    const umapClustersJpg = path.join(stage3cDir, `${comparison}_umap_clusters.jpg`);
    if (fs.existsSync(umapClustersJpg)) {
      result.umap_clusters_plot_path = umapClustersJpg;
    } else {
      result.warnings.push('UMAP clusters plot (JPG) not found');
    }

    const umapPhaseJpg = path.join(stage3cDir, `${comparison}_umap_phase.jpg`);
    if (fs.existsSync(umapPhaseJpg)) {
      result.umap_phase_plot_path = umapPhaseJpg;
    } else {
      result.warnings.push('UMAP phase plot (JPG) not found');
    }

    const rdsPath = path.join(stage3cDir, `${comparison}_seurat_clustered.rds`);
    if (fs.existsSync(rdsPath)) {
      result.seurat_rds_path = rdsPath;
    } else {
      result.errors.push('Seurat RDS file not found');
    }

    const markersCsv = path.join(stage3cDir, `${comparison}_cluster_markers.csv`);
    if (fs.existsSync(markersCsv)) {
      result.markers_csv_path = markersCsv;
    }

    // 4. Determine overall status
    if (result.errors.length === 0) {
      result.overall_status = 'SUCCESS';
      result.exit_code = 0;
    } else {
      result.overall_status = 'ERROR';
      result.exit_code = 1;
    }

  } catch (error) {
    result.errors.push(`Error parsing Stage 3c output: ${error.message}`);
    result.overall_status = 'ERROR';
    result.exit_code = 1;
  }

  return result;
}

/**
 * Format Stage 3c output for agent consumption
 *
 * @param {Object} parsedOutput - Parsed clustering results
 * @param {Object} dataInfo - Dataset information
 * @returns {string} - Markdown formatted output for agents
 */
export function formatStage3cForAgents(parsedOutput, dataInfo) {
  let formatted = `# Stage 3c: Clustering + UMAP Results (scRNA-seq)\n\n`;

  formatted += `**IMPORTANT:** Two UMAP plots will be provided as image attachments:\n`;
  formatted += `1. **UMAP colored by cluster** - Shows how cells are grouped\n`;
  formatted += `2. **UMAP colored by cell cycle phase** - Verifies cell cycle effect removal\n\n`;
  formatted += `You must VISUALLY EXAMINE both plots to assess clustering quality.\n\n`;

  // Dataset overview
  formatted += `## Dataset Overview\n`;
  formatted += `- Organism: ${dataInfo.organism || 'unknown'}\n`;
  formatted += `- Cells Processed: ${parsedOutput.cells_processed.toLocaleString()}\n`;
  formatted += `- Number of Clusters: ${parsedOutput.n_clusters}\n`;
  formatted += `- Data Type: Single-cell RNA-seq\n`;
  formatted += `\n`;

  // Cluster sizes
  formatted += `## Cluster Sizes\n`;
  for (const [clusterNum, size] of Object.entries(parsedOutput.cluster_sizes)) {
    const pct = parsedOutput.cluster_percentages[clusterNum];
    formatted += `- **Cluster ${clusterNum}**: ${size.toLocaleString()} cells (${pct.toFixed(1)}%)\n`;
  }
  formatted += `\n`;

  // Cell cycle phase distribution per cluster
  if (Object.keys(parsedOutput.phase_by_cluster).length > 0) {
    formatted += `## Cell Cycle Phase Distribution per Cluster\n`;
    formatted += `This shows whether cell cycle effect was successfully removed:\n\n`;
    for (const [clusterNum, phases] of Object.entries(parsedOutput.phase_by_cluster)) {
      formatted += `- **Cluster ${clusterNum}**: G1=${phases.G1}, S=${phases.S}, G2M=${phases.G2M}\n`;
    }
    formatted += `\n`;
    formatted += `**Note:** If regression was successful, clusters should NOT be dominated by a single phase.\n`;
    formatted += `Mixed phase distribution within clusters indicates cell cycle effect has been removed.\n\n`;
  }

  // Top marker genes
  if (Object.keys(parsedOutput.top_markers).length > 0) {
    formatted += `## Top Marker Genes per Cluster (Top 5)\n`;
    for (const [clusterNum, genes] of Object.entries(parsedOutput.top_markers)) {
      formatted += `- **Cluster ${clusterNum}**: ${genes.join(', ')}\n`;
    }
    formatted += `\n`;
  }

  // UMAP plots
  formatted += `## UMAP Plots\n`;
  if (parsedOutput.umap_clusters_plot_path) {
    formatted += `- Clusters: \`${parsedOutput.umap_clusters_plot_path}\`\n`;
  }
  if (parsedOutput.umap_phase_plot_path) {
    formatted += `- Cell Cycle Phase: \`${parsedOutput.umap_phase_plot_path}\`\n`;
  }
  formatted += `\n`;

  // Warnings and Errors
  if (parsedOutput.warnings.length > 0) {
    formatted += `## Warnings\n`;
    for (const warning of parsedOutput.warnings) {
      formatted += `- ⚠️ ${warning}\n`;
    }
    formatted += `\n`;
  }

  if (parsedOutput.errors.length > 0) {
    formatted += `## Errors\n`;
    for (const error of parsedOutput.errors) {
      formatted += `- ❌ ${error}\n`;
    }
    formatted += `\n`;
  }

  formatted += `---\n\n`;
  formatted += `## Your Task\n\n`;
  formatted += `Based on your VISUAL INSPECTION of the UMAP plots, determine:\n\n`;
  formatted += `1. **Are clusters biologically meaningful?**\n`;
  formatted += `   - Look at the UMAP (clusters): Are clusters well-separated and compact?\n`;
  formatted += `   - Check marker genes: Do they suggest distinct cell types or states?\n\n`;
  formatted += `2. **Was cell cycle effect successfully removed?**\n`;
  formatted += `   - Look at the UMAP (phase): Are phases mixed within clusters?\n`;
  formatted += `   - If phases form distinct clusters, cell cycle still dominates (regression may have failed)\n\n`;
  formatted += `3. **Is clustering resolution appropriate?**\n`;
  formatted += `   - Too low: Large heterogeneous clusters (under-clustering)\n`;
  formatted += `   - Too high: Many small similar clusters (over-clustering)\n`;
  formatted += `   - Just right: Distinct, homogeneous clusters with biological interpretation\n\n`;
  formatted += `4. **Decision Options:**\n`;
  formatted += `   - **APPROVE_CLUSTERS**: Clustering looks good, proceed to Stage 4 (DE analysis)\n`;
  formatted += `   - **ADJUST_RESOLUTION**: Clustering resolution needs adjustment (specify higher/lower)\n`;
  formatted += `   - **REQUEST_REANALYSIS**: Major issues detected, need to revisit earlier stages\n\n`;
  formatted += `5. **Quality Assessment Criteria:**\n`;
  formatted += `   - ✅ GOOD: Tight, well-separated clusters with mixed phases and distinct markers\n`;
  formatted += `   - ⚠️ ACCEPTABLE: Clusters somewhat diffuse but interpretable, phases mostly mixed\n`;
  formatted += `   - ❌ POOR: Poorly separated clusters, phase-driven clustering, or no distinct markers\n\n`;

  return formatted;
}

export default {
  generateStage3cScript,
  parseStage3cOutput,
  formatStage3cForAgents
};
