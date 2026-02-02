/**
 * Stage 3b: Cell Cycle Regression + Verification (scRNA-seq)
 *
 * This stage:
 * 1. Loads Seurat object from Stage 3a (with cell cycle scores)
 * 2. Performs cell cycle regression (ScaleData with vars.to.regress)
 * 3. Re-runs PCA on regressed data
 * 4. Generates DimPlot (after regression) and comparison plots
 * 5. Calculates regression success metrics (correlation reduction)
 *
 * Outputs:
 * - Seurat object with regressed data (.rds)
 * - PCA DimPlot (after regression) - PDF + JPG
 * - Before/After comparison plot - PDF + JPG
 * - Regression summary (correlation reduction, success status)
 *
 * Agent decision: REGRESSION_SUCCESS / REGRESSION_PARTIAL / REGRESSION_FAILED
 *
 * Generates script → Parses output → Formats for agents
 */

import fs from 'fs';
import path from 'path';

/**
 * Generate Stage 3b script: Cell Cycle Regression
 *
 * @param {Object} dataInfo - Dataset information (samples, organism, etc.)
 * @param {Object} config - Configuration (paths, comparison name, etc.)
 * @param {string} outputPath - Path to write the generated script
 * @param {Object} stage3aResult - Results from Stage 3a (for RDS file location)
 * @returns {Object} - { scriptPath, scriptContent }
 */
export function generateStage3bScript(dataInfo, config, outputPath, stage3aResult = {}) {
  const scriptsPath = process.env.SCRIPTS_PATH || './bio_informatics/scripts';
  const condaEnv = process.env.CONDA_ENV || 'singlecell';

  // Paths
  const outputDir = path.dirname(outputPath);
  const stage3bDir = path.join(outputDir, 'stage3b_cell_cycle_regression');
  const stage3aDir = path.join(outputDir, 'stage3a_cell_cycle_scoring');

  // Comparison name
  const comparison = config.comparison || 'comparison';

  // Seurat RDS path from Stage 3a
  const seuratRdsPath = stage3aResult.seurat_rds_path || path.join(stage3aDir, `${comparison}_seurat_scored.rds`);

  // Build script
  let script = `#!/bin/bash
#
# Stage 3b: Cell Cycle Regression + Verification
# Generated: ${new Date().toISOString()}
#

set -e  # Exit on error

echo "========================================"
echo "Stage 3b: Cell Cycle Regression"
echo "========================================"
echo ""

# Environment
SCRIPTS_PATH="${scriptsPath}"
CONDA_ENV="${condaEnv}"
OUTPUT_DIR="${stage3bDir}"
COMPARISON="${comparison}"
SEURAT_RDS="${seuratRdsPath}"

# Activate conda environment
echo "[Stage 3b] Activating conda environment: \${CONDA_ENV}"
set +u
source $(conda info --base)/etc/profile.d/conda.sh
conda activate \${CONDA_ENV}
set -u

# Create output directory
mkdir -p "\${OUTPUT_DIR}"
cd "\${OUTPUT_DIR}"

echo ""
echo "Running Cell Cycle Regression..."
echo "========================================"
echo ""

# Run seurat_cell_cycle_regression.R
Rscript \${SCRIPTS_PATH}/seurat_cell_cycle_regression.R \\
  "\${SEURAT_RDS}" \\
  "\${COMPARISON}"

if [ $? -ne 0 ]; then
  echo "ERROR: Cell cycle regression failed"
  exit 1
fi

echo ""
echo "========================================"
echo "Stage 3b Complete"
echo "========================================"
echo ""
echo "Output files:"
echo "  - \${COMPARISON}_seurat_regressed.rds (Seurat object with regressed data)"
echo "  - \${COMPARISON}_cell_cycle_after.pdf (PCA plot after - user)"
echo "  - \${COMPARISON}_cell_cycle_after.jpg (PCA plot after - agents)"
echo "  - \${COMPARISON}_cell_cycle_comparison.pdf (before/after comparison - user)"
echo "  - \${COMPARISON}_cell_cycle_comparison.jpg (before/after comparison - agents)"
echo "  - \${COMPARISON}_regression_summary.txt (regression metrics)"
echo ""
echo "Next: Stage 3c (Clustering + UMAP)"
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
 * Parse Stage 3b output: Read cell cycle regression results
 *
 * @param {string} outputDir - Output directory containing stage3b_cell_cycle_regression/
 * @param {string} comparison - Comparison name
 * @returns {Object} - Parsed regression results
 */
export function parseStage3bOutput(outputDir, comparison) {
  const stage3bDir = path.join(outputDir, 'stage3b_cell_cycle_regression');

  // Initialize result structure
  const result = {
    overall_status: 'UNKNOWN',
    exit_code: 0,
    cells_processed: 0,
    correlations_before: {
      PC1_S: null,
      PC1_G2M: null,
      PC2_S: null,
      PC2_G2M: null
    },
    correlations_after: {
      PC1_S: null,
      PC1_G2M: null,
      PC2_S: null,
      PC2_G2M: null
    },
    correlation_reduction: {
      PC1_S: null,
      PC1_G2M: null,
      PC2_S: null,
      PC2_G2M: null
    },
    regression_status: 'UNKNOWN',  // SUCCESS, PARTIAL, or FAILED
    max_correlation_after: null,
    comparison_plot_path: null,  // Path to comparison JPG for agents
    seurat_rds_path: null,
    warnings: [],
    errors: []
  };

  try {
    // 1. Check if output directory exists
    if (!fs.existsSync(stage3bDir)) {
      result.errors.push(`Stage 3b output directory not found: ${stage3bDir}`);
      result.overall_status = 'ERROR';
      return result;
    }

    // 2. Parse summary file
    const summaryPath = path.join(stage3bDir, `${comparison}_regression_summary.txt`);
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

    // Extract BEFORE correlations
    const beforeSection = summaryContent.match(/BEFORE Regression:([\s\S]*?)AFTER Regression:/);
    if (beforeSection) {
      const beforeText = beforeSection[1];
      const corS_PC1_before = beforeText.match(/PC1 vs S\.Score:\s*r=([-\d.]+)/);
      const corG2M_PC1_before = beforeText.match(/PC1 vs G2M\.Score:\s*r=([-\d.]+)/);
      const corS_PC2_before = beforeText.match(/PC2 vs S\.Score:\s*r=([-\d.]+)/);
      const corG2M_PC2_before = beforeText.match(/PC2 vs G2M\.Score:\s*r=([-\d.]+)/);

      if (corS_PC1_before) result.correlations_before.PC1_S = parseFloat(corS_PC1_before[1]);
      if (corG2M_PC1_before) result.correlations_before.PC1_G2M = parseFloat(corG2M_PC1_before[1]);
      if (corS_PC2_before) result.correlations_before.PC2_S = parseFloat(corS_PC2_before[1]);
      if (corG2M_PC2_before) result.correlations_before.PC2_G2M = parseFloat(corG2M_PC2_before[1]);
    }

    // Extract AFTER correlations and reductions
    const afterSection = summaryContent.match(/AFTER Regression:([\s\S]*?)Regression Status:/);
    if (afterSection) {
      const afterText = afterSection[1];

      // Parse "r=X.XXX (reduced by Y.YYY)" format
      const corS_PC1_after = afterText.match(/PC1 vs S\.Score:\s*r=([-\d.]+)\s*\(reduced by ([-\d.]+)\)/);
      const corG2M_PC1_after = afterText.match(/PC1 vs G2M\.Score:\s*r=([-\d.]+)\s*\(reduced by ([-\d.]+)\)/);
      const corS_PC2_after = afterText.match(/PC2 vs S\.Score:\s*r=([-\d.]+)\s*\(reduced by ([-\d.]+)\)/);
      const corG2M_PC2_after = afterText.match(/PC2 vs G2M\.Score:\s*r=([-\d.]+)\s*\(reduced by ([-\d.]+)\)/);

      if (corS_PC1_after) {
        result.correlations_after.PC1_S = parseFloat(corS_PC1_after[1]);
        result.correlation_reduction.PC1_S = parseFloat(corS_PC1_after[2]);
      }
      if (corG2M_PC1_after) {
        result.correlations_after.PC1_G2M = parseFloat(corG2M_PC1_after[1]);
        result.correlation_reduction.PC1_G2M = parseFloat(corG2M_PC1_after[2]);
      }
      if (corS_PC2_after) {
        result.correlations_after.PC2_S = parseFloat(corS_PC2_after[1]);
        result.correlation_reduction.PC2_S = parseFloat(corS_PC2_after[2]);
      }
      if (corG2M_PC2_after) {
        result.correlations_after.PC2_G2M = parseFloat(corG2M_PC2_after[1]);
        result.correlation_reduction.PC2_G2M = parseFloat(corG2M_PC2_after[2]);
      }
    }

    // Extract regression status
    const statusMatch = summaryContent.match(/Regression Status:\s*(\w+)/);
    if (statusMatch) {
      result.regression_status = statusMatch[1];  // SUCCESS or PARTIAL
    }

    // Extract max correlation after
    const maxCorMatch = summaryContent.match(/Max correlation after regression:\s*([\d.]+)/);
    if (maxCorMatch) {
      result.max_correlation_after = parseFloat(maxCorMatch[1]);
    }

    // 3. Check for output files
    const comparisonJpgPath = path.join(stage3bDir, `${comparison}_cell_cycle_comparison.jpg`);
    if (fs.existsSync(comparisonJpgPath)) {
      result.comparison_plot_path = comparisonJpgPath;
    } else {
      result.warnings.push('Comparison plot (JPG) not found');
    }

    const rdsPath = path.join(stage3bDir, `${comparison}_seurat_regressed.rds`);
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
    result.errors.push(`Error parsing Stage 3b output: ${error.message}`);
    result.overall_status = 'ERROR';
    result.exit_code = 1;
  }

  return result;
}

/**
 * Format Stage 3b output for agent consumption
 *
 * @param {Object} parsedOutput - Parsed regression results
 * @param {Object} dataInfo - Dataset information
 * @returns {string} - Markdown formatted output for agents
 */
export function formatStage3bForAgents(parsedOutput, dataInfo) {
  let formatted = `# Stage 3b: Cell Cycle Regression Results (scRNA-seq)\n\n`;

  formatted += `**IMPORTANT:** A comparison plot (Before vs After) will be provided as an image attachment.\n`;
  formatted += `The plot shows PCA colored by cell cycle phase BEFORE and AFTER regression.\n`;
  formatted += `You must VISUALLY EXAMINE both plots to assess regression effectiveness.\n\n`;

  // Dataset overview
  formatted += `## Dataset Overview\n`;
  formatted += `- Organism: ${dataInfo.organism || 'unknown'}\n`;
  formatted += `- Cells Processed: ${parsedOutput.cells_processed.toLocaleString()}\n`;
  formatted += `- Data Type: Single-cell RNA-seq\n`;
  formatted += `\n`;

  // Regression results
  formatted += `## Regression Results\n\n`;
  formatted += `### BEFORE Regression:\n`;
  formatted += `- PC1 vs S.Score: r=${parsedOutput.correlations_before.PC1_S?.toFixed(3) || 'N/A'}\n`;
  formatted += `- PC1 vs G2M.Score: r=${parsedOutput.correlations_before.PC1_G2M?.toFixed(3) || 'N/A'}\n`;
  formatted += `- PC2 vs S.Score: r=${parsedOutput.correlations_before.PC2_S?.toFixed(3) || 'N/A'}\n`;
  formatted += `- PC2 vs G2M.Score: r=${parsedOutput.correlations_before.PC2_G2M?.toFixed(3) || 'N/A'}\n`;
  formatted += `\n`;

  formatted += `### AFTER Regression:\n`;
  formatted += `- PC1 vs S.Score: r=${parsedOutput.correlations_after.PC1_S?.toFixed(3) || 'N/A'} (reduced by ${parsedOutput.correlation_reduction.PC1_S?.toFixed(3) || 'N/A'})\n`;
  formatted += `- PC1 vs G2M.Score: r=${parsedOutput.correlations_after.PC1_G2M?.toFixed(3) || 'N/A'} (reduced by ${parsedOutput.correlation_reduction.PC1_G2M?.toFixed(3) || 'N/A'})\n`;
  formatted += `- PC2 vs S.Score: r=${parsedOutput.correlations_after.PC2_S?.toFixed(3) || 'N/A'} (reduced by ${parsedOutput.correlation_reduction.PC2_S?.toFixed(3) || 'N/A'})\n`;
  formatted += `- PC2 vs G2M.Score: r=${parsedOutput.correlations_after.PC2_G2M?.toFixed(3) || 'N/A'} (reduced by ${parsedOutput.correlation_reduction.PC2_G2M?.toFixed(3) || 'N/A'})\n`;
  formatted += `\n`;

  // Regression status
  formatted += `### Regression Status: **${parsedOutput.regression_status}**\n`;
  if (parsedOutput.max_correlation_after !== null) {
    formatted += `- Max correlation after regression: ${parsedOutput.max_correlation_after.toFixed(3)}\n`;
    formatted += `- Threshold: 0.3 (correlations should be < 0.3 for successful regression)\n`;
  }
  formatted += `\n`;

  // Interpretation
  formatted += `### Interpretation:\n`;
  if (parsedOutput.regression_status === 'SUCCESS') {
    formatted += `✅ **SUCCESS**: Cell cycle effect has been successfully removed from the data.\n`;
    formatted += `The maximum correlation is below 0.3, indicating minimal residual cell cycle influence.\n`;
  } else if (parsedOutput.regression_status === 'PARTIAL') {
    formatted += `⚠️ **PARTIAL**: Cell cycle effect has been reduced but not completely eliminated.\n`;
    formatted += `The maximum correlation is still ≥ 0.3, indicating some residual cell cycle influence.\n`;
  } else {
    formatted += `❌ **FAILED**: Regression did not effectively reduce cell cycle correlations.\n`;
  }
  formatted += `\n`;

  // Comparison plot
  if (parsedOutput.comparison_plot_path) {
    formatted += `## Comparison Plot (Before vs After)\n`;
    formatted += `\`${parsedOutput.comparison_plot_path}\`\n`;
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
  formatted += `## Your Task\n\n`;
  formatted += `Based on your VISUAL INSPECTION of the comparison plot and the correlation metrics, determine:\n\n`;
  formatted += `1. **Did regression successfully remove cell cycle effect?**\n`;
  formatted += `   - Look at the before/after plots: Are phase-based clusters less distinct after regression?\n`;
  formatted += `   - Check correlations: Were they reduced to < 0.3 (absolute value)?\n\n`;
  formatted += `2. **Decision Options:**\n`;
  formatted += `   - **REGRESSION_SUCCESS**: Cell cycle effect removed, proceed to Stage 3c (Clustering)\n`;
  formatted += `   - **REGRESSION_PARTIAL**: Effect reduced but not eliminated - may still proceed cautiously\n`;
  formatted += `   - **REGRESSION_FAILED**: Regression ineffective, may need to adjust parameters or reconsider approach\n\n`;
  formatted += `3. **Recommendation:**\n`;
  formatted += `   - If max correlation after < 0.3: Recommend REGRESSION_SUCCESS\n`;
  formatted += `   - If 0.3 ≤ max correlation < 0.5: Recommend REGRESSION_PARTIAL (proceed with caution)\n`;
  formatted += `   - If max correlation ≥ 0.5: Recommend REGRESSION_FAILED (re-analysis needed)\n\n`;

  return formatted;
}

export default {
  generateStage3bScript,
  parseStage3bOutput,
  formatStage3bForAgents
};
