/**
 * Stage 4: Differential Expression Analysis
 *
 * This stage:
 * 1. Runs DE analysis using either:
 *    - simpleEdger3.R (no batch correction)
 *    - batch_effect_edgeR_v3.R (with batch correction)
 * 2. Merges RPKM + DE results into final output file
 *
 * Generates script → Parses output → Formats for agents
 */

import fs from 'fs';
import path from 'path';

/**
 * Generate Stage 4 script: Differential Expression Analysis
 *
 * @param {Object} dataInfo - Dataset information (samples, organism, etc.)
 * @param {Object} config - Configuration (paths, comparison name, etc.)
 * @param {string} outputPath - Path to write the generated script
 * @param {Object} stage3Result - Results from Stage 3 (DE method decision)
 * @returns {Object} - { scriptPath, scriptContent }
 */
export function generateStage4Script(dataInfo, config, outputPath, stage3Result = {}) {
  const scriptsPath = process.env.SCRIPTS_PATH || '/users/ha00014/bin';
  const condaEnv = process.env.CONDA_ENV || 'rnaseq';

  // Get DE method decision from Stage 3
  const deMethod = stage3Result.deMethod || 'simpleEdger';
  const batchSpecification = stage3Result.batchSpecification || 'auto';

  // Paths
  const outputDir = path.dirname(outputPath);
  const stage4Dir = path.join(outputDir, 'stage4_de_analysis');
  const stage3Dir = path.join(outputDir, 'stage3_quantification');

  // Comparison name
  const comparison = config.comparison || 'comparison';
  const controlKeyword = config.controlKeyword || 'cont';
  const treatmentKeyword = config.treatmentKeyword || 'ips';

  // Choose DE script based on Stage 3 decision
  let deScriptName;
  let deScriptArgs;

  if (deMethod === 'batch_effect_edger') {
    deScriptName = 'batch_effect_edgeR_v3.R';
    deScriptArgs = `"\${STAGE3_DIR}/\${COMPARISON}.count.filtered.csv" "\${CONTROL_KEYWORD}" "\${TREATMENT_KEYWORD}" "\${BATCH_SPEC}"`;
  } else {
    deScriptName = 'simpleEdger3.R';
    deScriptArgs = `"\${STAGE3_DIR}/\${COMPARISON}.count.filtered.csv" "\${CONTROL_KEYWORD}" "\${TREATMENT_KEYWORD}"`;
  }

  // Build script
  let script = `#!/bin/bash
#
# Stage 4: Differential Expression Analysis
# Generated: ${new Date().toISOString()}
#

set -e  # Exit on error

echo "========================================"
echo "Stage 4: Differential Expression Analysis"
echo "========================================"
echo ""

# Environment
SCRIPTS_PATH="${scriptsPath}"
CONDA_ENV="${condaEnv}"
OUTPUT_DIR="${stage4Dir}"
STAGE3_DIR="${stage3Dir}"
COMPARISON="${comparison}"
CONTROL_KEYWORD="${controlKeyword}"
TREATMENT_KEYWORD="${treatmentKeyword}"
DE_METHOD="${deMethod}"
BATCH_SPEC="${batchSpecification}"

# Activate conda environment (disable -u temporarily to avoid conda activation errors)
echo "[Stage 4] Activating conda environment: \${CONDA_ENV}"
set +u
source $(conda info --base)/etc/profile.d/conda.sh
conda activate \${CONDA_ENV}
set -u

# Create output directory
mkdir -p "\${OUTPUT_DIR}"
cd "\${OUTPUT_DIR}"

echo ""
echo "DE Method: \${DE_METHOD}"
echo "Control Keyword: \${CONTROL_KEYWORD}"
echo "Treatment Keyword: \${TREATMENT_KEYWORD}"
`;

  // Add batch specification info if using batch correction
  if (deMethod === 'batch_effect_edger') {
    script += `echo "Batch Specification: \${BATCH_SPEC}"
`;
  }

  script += `echo ""
echo "Step 1: Run Differential Expression Analysis"
echo "========================================"
echo ""

# Run DE analysis
Rscript \${SCRIPTS_PATH}/${deScriptName} \\
  ${deScriptArgs}

if [ $? -ne 0 ]; then
  echo "ERROR: Differential expression analysis failed"
  exit 2
fi

echo ""
echo "Step 2: Merge RPKM + DE Results"
echo "========================================"
echo ""

# Merge RPKM and DE results
Rscript \${SCRIPTS_PATH}/merge_results.R \\
  "\${STAGE3_DIR}/\${COMPARISON}.rpkm.entrz.csv" \\
  "\${COMPARISON}_DE.csv" \\
  "\${COMPARISON}"

if [ $? -ne 0 ]; then
  echo "ERROR: merge_results.R failed"
  exit 2
fi

echo ""
echo "========================================"
echo "Stage 4 Complete"
echo "========================================"
echo ""
echo "Output files:"
echo "  - \${COMPARISON}_DE.csv (DE results)"
echo "  - \${COMPARISON}_final.xlsx (merged RPKM + DE)"
echo ""
echo "DE Analysis: SUCCESS"

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
 * Parse Stage 4 output: Read DE analysis results
 *
 * @param {string} outputDir - Output directory containing stage4_de_analysis/
 * @returns {Object} - Parsed DE results with gene counts, top genes, etc.
 */
export function parseStage4Output(outputDir) {
  const stage4Dir = path.join(outputDir, 'stage4_de_analysis');

  // Initialize result structure
  const result = {
    overall_status: 'UNKNOWN',
    de_results_file: null,
    final_results_file: null,
    num_degs_up: 0,
    num_degs_down: 0,
    total_degs: 0,
    fc_range: [null, null],
    fdr_threshold: 0.05,
    top_up_genes: [],
    top_down_genes: [],
    total_genes_tested: 0,
    warnings: [],
    errors: []
  };

  try {
    // 1. Check if output directory exists
    if (!fs.existsSync(stage4Dir)) {
      result.errors.push(`Stage 4 output directory not found: ${stage4Dir}`);
      result.overall_status = 'ERROR';
      return result;
    }

    // 2. Find DE results file
    const deFiles = fs.readdirSync(stage4Dir).filter(f => f.endsWith('_DE.csv'));
    if (deFiles.length === 0) {
      result.errors.push('No DE results file found (*_DE.csv)');
      result.overall_status = 'ERROR';
      return result;
    }

    const deFile = deFiles[0];
    const comparison = deFile.replace('_DE.csv', '');
    result.de_results_file = path.join(stage4Dir, deFile);

    // 3. Find final results file (Excel)
    const finalFiles = fs.readdirSync(stage4Dir).filter(f => f.endsWith('_final.xlsx'));
    if (finalFiles.length > 0) {
      result.final_results_file = path.join(stage4Dir, finalFiles[0]);
    }

    // 4. Parse DE results CSV
    const deData = fs.readFileSync(result.de_results_file, 'utf-8');
    const lines = deData.split('\n').filter(l => l.trim().length > 0);

    if (lines.length < 2) {
      result.warnings.push('DE results file is empty or has no data rows');
      result.overall_status = 'WARNING';
      return result;
    }

    // Parse header to find column indices
    const headers = lines[0].split(',').map(h => h.trim().replace(/^"|"$/g, ''));
    const logFCIndex = headers.findIndex(h => h.toLowerCase().includes('logfc'));
    const fdrIndex = headers.findIndex(h => h.toLowerCase().includes('fdr') || h.toLowerCase().includes('adj.p'));
    const geneIndex = 0; // Usually first column

    if (logFCIndex === -1 || fdrIndex === -1) {
      result.warnings.push('Could not find logFC or FDR columns in DE results');
    }

    // Parse data rows
    const genes = [];
    for (let i = 1; i < lines.length; i++) {
      const values = lines[i].split(',').map(v => v.trim().replace(/^"|"$/g, ''));

      if (values.length < headers.length) continue;

      const geneName = values[geneIndex];
      const logFC = logFCIndex !== -1 ? parseFloat(values[logFCIndex]) : null;
      const fdr = fdrIndex !== -1 ? parseFloat(values[fdrIndex]) : null;

      if (!isNaN(logFC) && !isNaN(fdr)) {
        genes.push({ name: geneName, logFC, fdr });
      }
    }

    result.total_genes_tested = genes.length;

    // 5. Filter DEGs (FDR < 0.05 by default)
    const degs = genes.filter(g => g.fdr < result.fdr_threshold);
    const upGenes = degs.filter(g => g.logFC > 0).sort((a, b) => b.logFC - a.logFC);
    const downGenes = degs.filter(g => g.logFC < 0).sort((a, b) => a.logFC - b.logFC);

    result.num_degs_up = upGenes.length;
    result.num_degs_down = downGenes.length;
    result.total_degs = degs.length;

    // 6. Get top genes (top 10 each)
    result.top_up_genes = upGenes.slice(0, 10).map(g => ({
      name: g.name,
      logFC: g.logFC.toFixed(2),
      fdr: g.fdr.toExponential(2)
    }));

    result.top_down_genes = downGenes.slice(0, 10).map(g => ({
      name: g.name,
      logFC: g.logFC.toFixed(2),
      fdr: g.fdr.toExponential(2)
    }));

    // 7. Calculate FC range for all DEGs
    if (degs.length > 0) {
      const allLogFCs = degs.map(g => g.logFC);
      result.fc_range = [Math.min(...allLogFCs), Math.max(...allLogFCs)];
    }

    // 8. Set overall status
    if (result.total_degs === 0) {
      result.warnings.push('No differentially expressed genes found (FDR < 0.05)');
      result.overall_status = 'WARNING';
    } else {
      result.overall_status = 'SUCCESS';
    }

  } catch (error) {
    result.errors.push(`Error parsing Stage 4 output: ${error.message}`);
    result.overall_status = 'ERROR';
  }

  return result;
}

/**
 * Format Stage 4 output for agent consumption
 *
 * @param {Object} parsedOutput - Parsed DE results
 * @param {Object} dataInfo - Dataset information
 * @returns {string} - Markdown formatted output for agents
 */
export function formatStage4ForAgents(parsedOutput, dataInfo) {
  let formatted = `# Stage 4: Differential Expression Analysis Results\n\n`;

  // Dataset overview
  formatted += `## Dataset Overview\n`;
  formatted += `- Organism: ${dataInfo.organism || 'unknown'}\n`;
  formatted += `- Comparison: Control vs Treatment\n`;
  formatted += `- Total Genes Tested: ${parsedOutput.total_genes_tested.toLocaleString()}\n`;
  formatted += `\n`;

  // DE Results Summary
  formatted += `## Differential Expression Summary\n`;
  formatted += `- **Total DEGs (FDR < ${parsedOutput.fdr_threshold}):** ${parsedOutput.total_degs}\n`;
  formatted += `- **Up-regulated:** ${parsedOutput.num_degs_up}\n`;
  formatted += `- **Down-regulated:** ${parsedOutput.num_degs_down}\n`;

  if (parsedOutput.fc_range[0] !== null) {
    formatted += `- **LogFC Range:** ${parsedOutput.fc_range[0].toFixed(2)} to ${parsedOutput.fc_range[1].toFixed(2)}\n`;
  }
  formatted += `\n`;

  // Top Up-regulated Genes
  if (parsedOutput.top_up_genes.length > 0) {
    formatted += `## Top Up-Regulated Genes\n`;
    formatted += `| Gene | LogFC | FDR |\n`;
    formatted += `|------|-------|-----|\n`;

    for (const gene of parsedOutput.top_up_genes) {
      formatted += `| ${gene.name} | ${gene.logFC} | ${gene.fdr} |\n`;
    }
    formatted += `\n`;
  }

  // Top Down-regulated Genes
  if (parsedOutput.top_down_genes.length > 0) {
    formatted += `## Top Down-Regulated Genes\n`;
    formatted += `| Gene | LogFC | FDR |\n`;
    formatted += `|------|-------|-----|\n`;

    for (const gene of parsedOutput.top_down_genes) {
      formatted += `| ${gene.name} | ${gene.logFC} | ${gene.fdr} |\n`;
    }
    formatted += `\n`;
  }

  // Output Files
  formatted += `## Output Files\n`;
  formatted += `- DE Results: \`${path.basename(parsedOutput.de_results_file || 'N/A')}\`\n`;
  formatted += `- Final Results: \`${path.basename(parsedOutput.final_results_file || 'N/A')}\`\n`;
  formatted += `\n`;

  // Status
  formatted += `## Analysis Status\n`;
  formatted += `- Overall Status: **${parsedOutput.overall_status}**\n`;

  if (parsedOutput.warnings.length > 0) {
    formatted += `\n### Warnings\n`;
    for (const warning of parsedOutput.warnings) {
      formatted += `- ⚠️ ${warning}\n`;
    }
  }

  if (parsedOutput.errors.length > 0) {
    formatted += `\n### Errors\n`;
    for (const error of parsedOutput.errors) {
      formatted += `- ❌ ${error}\n`;
    }
  }

  formatted += `\n`;

  // Decision Prompt
  formatted += `---\n\n`;
  formatted += `## Your Decision (REQUIRED FORMAT)\n\n`;
  formatted += `Please review the DE analysis results above and provide your final assessment:\n\n`;

  formatted += `**Final_Decision:** [APPROVE / REQUEST_REANALYSIS]\n`;
  formatted += `- APPROVE: Results look good, analysis complete\n`;
  formatted += `- REQUEST_REANALYSIS: Issues detected, need to re-run with different parameters\n\n`;

  formatted += `**Confidence:** [HIGH / MEDIUM / LOW]\n\n`;
  formatted += `**Reasoning:** [2-3 sentences explaining your assessment]\n\n`;

  formatted += `**Biological_Assessment:** [Brief biological interpretation of the results]\n`;

  return formatted;
}

export default {
  generateStage4Script,
  parseStage4Output,
  formatStage4ForAgents
};
