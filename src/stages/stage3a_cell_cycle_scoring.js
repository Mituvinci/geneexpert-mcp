/**
 * Stage 3a: Cell Cycle Scoring + Visualization (scRNA-seq)
 *
 * This stage:
 * 1. Loads Seurat object from count matrix
 * 2. Performs basic QC filtering
 * 3. Scores cells for cell cycle phase (S, G2M, G1)
 * 4. Runs PCA and generates DimPlot colored by Phase
 * 5. Calculates correlation between PCs and cell cycle scores
 *
 * Outputs:
 * - Seurat object with cell cycle scores (.rds)
 * - PCA DimPlot (before regression) - PDF + JPG
 * - Summary statistics (cell cycle distribution, PCA variance, correlations)
 *
 * Agent decision: INFORMATIONAL ONLY - always proceeds to Stage 3b (regression)
 *
 * Generates script → Parses output → Formats for agents
 */

import fs from 'fs';
import path from 'path';

/**
 * Generate Stage 3a script: Cell Cycle Scoring
 *
 * @param {Object} dataInfo - Dataset information (samples, organism, etc.)
 * @param {Object} config - Configuration (paths, comparison name, etc.)
 * @param {string} outputPath - Path to write the generated script
 * @returns {Object} - { scriptPath, scriptContent }
 */
export function generateStage3aScript(dataInfo, config, outputPath) {
  const scriptsPath = process.env.SCRIPTS_PATH || '/users/ha00014/bin';
  const condaEnv = process.env.CONDA_ENV || 'singlecell';

  // Paths
  const outputDir = path.dirname(outputPath);
  const stage3aDir = path.join(outputDir, 'stage3a_cell_cycle_scoring');

  // Comparison name
  const comparison = config.comparison || 'comparison';
  const organism = dataInfo.organism || 'human';

  // Count matrix path (from previous stage or provided)
  const countMatrixPath = config.countMatrixPath || path.join(outputDir, 'counts.csv');

  // Build script
  let script = `#!/bin/bash
#
# Stage 3a: Cell Cycle Scoring + Visualization
# Generated: ${new Date().toISOString()}
#

set -e  # Exit on error

echo "========================================"
echo "Stage 3a: Cell Cycle Scoring"
echo "========================================"
echo ""

# Environment
SCRIPTS_PATH="${scriptsPath}"
CONDA_ENV="${condaEnv}"
OUTPUT_DIR="${stage3aDir}"
COMPARISON="${comparison}"
ORGANISM="${organism}"
COUNT_MATRIX="${countMatrixPath}"

# Activate conda environment
echo "[Stage 3a] Activating conda environment: \${CONDA_ENV}"
set +u
source $(conda info --base)/etc/profile.d/conda.sh
conda activate \${CONDA_ENV}
set -u

# Create output directory
mkdir -p "\${OUTPUT_DIR}"
cd "\${OUTPUT_DIR}"

echo ""
echo "Running Cell Cycle Scoring..."
echo "========================================"
echo ""

# Run seurat_cell_cycle_scoring.R
Rscript \${SCRIPTS_PATH}/seurat_cell_cycle_scoring.R \\
  "\${COUNT_MATRIX}" \\
  "\${ORGANISM}" \\
  "\${COMPARISON}"

if [ $? -ne 0 ]; then
  echo "ERROR: Cell cycle scoring failed"
  exit 1
fi

echo ""
echo "========================================"
echo "Stage 3a Complete"
echo "========================================"
echo ""
echo "Output files:"
echo "  - \${COMPARISON}_seurat_scored.rds (Seurat object with scores)"
echo "  - \${COMPARISON}_cell_cycle_before.pdf (PCA plot - user)"
echo "  - \${COMPARISON}_cell_cycle_before.jpg (PCA plot - agents)"
echo "  - \${COMPARISON}_cell_cycle_summary.txt (summary statistics)"
echo ""
echo "Next: Stage 3b (Cell Cycle Regression)"
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
 * Parse Stage 3a output: Read cell cycle scoring results
 *
 * @param {string} outputDir - Output directory containing stage3a_cell_cycle_scoring/
 * @param {string} comparison - Comparison name
 * @returns {Object} - Parsed cell cycle scoring results
 */
export function parseStage3aOutput(outputDir, comparison) {
  const stage3aDir = path.join(outputDir, 'stage3a_cell_cycle_scoring');

  // Initialize result structure
  const result = {
    overall_status: 'UNKNOWN',
    exit_code: 0,
    cells_analyzed: 0,
    phase_distribution: {
      G1: 0,
      S: 0,
      G2M: 0
    },
    phase_percentages: {
      G1: 0,
      S: 0,
      G2M: 0
    },
    pca_variance: {
      PC1: null,
      PC2: null,
      PC3: null,
      PC4: null,
      PC5: null
    },
    correlations: {
      PC1_S: null,
      PC1_G2M: null,
      PC2_S: null,
      PC2_G2M: null
    },
    pca_plot_path: null,  // Path to JPEG for agents
    seurat_rds_path: null,
    warnings: [],
    errors: []
  };

  try {
    // 1. Check if output directory exists
    if (!fs.existsSync(stage3aDir)) {
      result.errors.push(`Stage 3a output directory not found: ${stage3aDir}`);
      result.overall_status = 'ERROR';
      return result;
    }

    // 2. Parse summary file
    const summaryPath = path.join(stage3aDir, `${comparison}_cell_cycle_summary.txt`);
    if (!fs.existsSync(summaryPath)) {
      result.errors.push(`Summary file not found: ${summaryPath}`);
      result.overall_status = 'ERROR';
      return result;
    }

    const summaryContent = fs.readFileSync(summaryPath, 'utf-8');

    // Extract cell cycle distribution
    const g1Match = summaryContent.match(/G1:\s*(\d+)\s*cells\s*\(([\d.]+)%\)/);
    const sMatch = summaryContent.match(/S:\s*(\d+)\s*cells\s*\(([\d.]+)%\)/);
    const g2mMatch = summaryContent.match(/G2M:\s*(\d+)\s*cells\s*\(([\d.]+)%\)/);

    if (g1Match) {
      result.phase_distribution.G1 = parseInt(g1Match[1]);
      result.phase_percentages.G1 = parseFloat(g1Match[2]);
    }
    if (sMatch) {
      result.phase_distribution.S = parseInt(sMatch[1]);
      result.phase_percentages.S = parseFloat(sMatch[2]);
    }
    if (g2mMatch) {
      result.phase_distribution.G2M = parseInt(g2mMatch[1]);
      result.phase_percentages.G2M = parseFloat(g2mMatch[2]);
    }

    result.cells_analyzed = result.phase_distribution.G1 + result.phase_distribution.S + result.phase_distribution.G2M;

    // Extract PCA variance
    for (let i = 1; i <= 5; i++) {
      const pcMatch = summaryContent.match(new RegExp(`PC${i}:\\s*([\\d.]+)%`));
      if (pcMatch) {
        result.pca_variance[`PC${i}`] = parseFloat(pcMatch[1]);
      }
    }

    // Extract correlations
    const corS_PC1 = summaryContent.match(/PC1 vs S\.Score:\s*r=([-\d.]+)/);
    const corG2M_PC1 = summaryContent.match(/PC1 vs G2M\.Score:\s*r=([-\d.]+)/);
    const corS_PC2 = summaryContent.match(/PC2 vs S\.Score:\s*r=([-\d.]+)/);
    const corG2M_PC2 = summaryContent.match(/PC2 vs G2M\.Score:\s*r=([-\d.]+)/);

    if (corS_PC1) result.correlations.PC1_S = parseFloat(corS_PC1[1]);
    if (corG2M_PC1) result.correlations.PC1_G2M = parseFloat(corG2M_PC1[1]);
    if (corS_PC2) result.correlations.PC2_S = parseFloat(corS_PC2[1]);
    if (corG2M_PC2) result.correlations.PC2_G2M = parseFloat(corG2M_PC2[1]);

    // 3. Check for output files
    const jpgPath = path.join(stage3aDir, `${comparison}_cell_cycle_before.jpg`);
    if (fs.existsSync(jpgPath)) {
      result.pca_plot_path = jpgPath;
    } else {
      result.warnings.push('PCA plot (JPG) not found');
    }

    const rdsPath = path.join(stage3aDir, `${comparison}_seurat_scored.rds`);
    if (fs.existsSync(rdsPath)) {
      result.seurat_rds_path = rdsPath;
    } else {
      result.errors.push('Seurat RDS file not found');
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
    result.errors.push(`Error parsing Stage 3a output: ${error.message}`);
    result.overall_status = 'ERROR';
    result.exit_code = 1;
  }

  return result;
}

/**
 * Format Stage 3a output for agent consumption
 *
 * @param {Object} parsedOutput - Parsed cell cycle scoring results
 * @param {Object} dataInfo - Dataset information
 * @returns {string} - Markdown formatted output for agents
 */
export function formatStage3aForAgents(parsedOutput, dataInfo) {
  let formatted = `# Stage 3a: Cell Cycle Scoring Results (scRNA-seq)\n\n`;

  formatted += `**IMPORTANT:** A PCA plot (PC1 vs PC2) will be provided as an image attachment.\n`;
  formatted += `The plot shows cells colored by cell cycle phase (G1, S, G2M).\n`;
  formatted += `This visualization shows cell cycle effect BEFORE regression.\n\n`;

  // Dataset overview
  formatted += `## Dataset Overview\n`;
  formatted += `- Organism: ${dataInfo.organism || 'unknown'}\n`;
  formatted += `- Cells Analyzed: ${parsedOutput.cells_analyzed.toLocaleString()}\n`;
  formatted += `- Data Type: Single-cell RNA-seq\n`;
  formatted += `\n`;

  // Cell cycle distribution
  formatted += `## Cell Cycle Phase Distribution\n`;
  formatted += `- **G1**: ${parsedOutput.phase_distribution.G1.toLocaleString()} cells (${parsedOutput.phase_percentages.G1.toFixed(1)}%)\n`;
  formatted += `- **S**: ${parsedOutput.phase_distribution.S.toLocaleString()} cells (${parsedOutput.phase_percentages.S.toFixed(1)}%)\n`;
  formatted += `- **G2M**: ${parsedOutput.phase_distribution.G2M.toLocaleString()} cells (${parsedOutput.phase_percentages.G2M.toFixed(1)}%)\n`;
  formatted += `\n`;

  // PCA variance
  formatted += `## PCA Variance Explained\n`;
  for (let i = 1; i <= 5; i++) {
    const variance = parsedOutput.pca_variance[`PC${i}`];
    if (variance !== null) {
      formatted += `- PC${i}: ${variance.toFixed(2)}%\n`;
    }
  }
  formatted += `\n`;

  // Correlations with cell cycle scores
  formatted += `## Correlation: PCs vs Cell Cycle Scores\n`;
  formatted += `These correlations indicate how much cell cycle dominates the principal components:\n\n`;
  formatted += `- **PC1 vs S.Score**: r=${parsedOutput.correlations.PC1_S?.toFixed(3) || 'N/A'}\n`;
  formatted += `- **PC1 vs G2M.Score**: r=${parsedOutput.correlations.PC1_G2M?.toFixed(3) || 'N/A'}\n`;
  formatted += `- **PC2 vs S.Score**: r=${parsedOutput.correlations.PC2_S?.toFixed(3) || 'N/A'}\n`;
  formatted += `- **PC2 vs G2M.Score**: r=${parsedOutput.correlations.PC2_G2M?.toFixed(3) || 'N/A'}\n`;
  formatted += `\n`;
  formatted += `**Interpretation:** |r| > 0.5 indicates strong correlation (cell cycle dominates that PC).\n\n`;

  // PCA plot
  if (parsedOutput.pca_plot_path) {
    formatted += `## PCA Plot (Before Regression)\n`;
    formatted += `\`${parsedOutput.pca_plot_path}\`\n`;
    formatted += `\n`;
  }

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
  formatted += `## Your Task (INFORMATIONAL ONLY)\n\n`;
  formatted += `This stage shows cell cycle effect BEFORE regression. Review the results:\n\n`;
  formatted += `1. **Cell Cycle Distribution**: Are cells balanced across phases or heavily skewed?\n`;
  formatted += `2. **PCA Correlations**: Do cell cycle scores strongly correlate with PC1/PC2?\n`;
  formatted += `3. **Visual Assessment**: In the PCA plot, do you see distinct clustering by cell cycle phase?\n\n`;
  formatted += `**NOTE:** The system will AUTOMATICALLY proceed to Stage 3b (Cell Cycle Regression).\n`;
  formatted += `No decision is required at this stage - this is for your information only.\n\n`;

  return formatted;
}

export default {
  generateStage3aScript,
  parseStage3aOutput,
  formatStage3aForAgents
};
