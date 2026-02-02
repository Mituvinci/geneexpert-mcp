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
  const scriptsPath = process.env.SCRIPTS_PATH || './bio_informatics/scripts';
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

  // Custom thresholds for re-analysis (optional)
  const fdrThreshold = config.fdr_threshold || 0.05;
  const logfcThreshold = config.logfc_threshold || 0.585;

  // Choose DE script based on Stage 3 decision
  let deScriptName;
  let deScriptArgs;

  if (deMethod === 'batch_effect_edger') {
    deScriptName = 'batch_effect_edgeR_v3.R';
    deScriptArgs = `"\${STAGE3_DIR}/\${COMPARISON}.count.filtered.csv" "\${CONTROL_KEYWORD}" "\${TREATMENT_KEYWORD}" "\${BATCH_SPEC}" "\${FDR_THRESHOLD}" "\${LOGFC_THRESHOLD}"`;
  } else {
    deScriptName = 'simpleEdger3.R';
    deScriptArgs = `"\${STAGE3_DIR}/\${COMPARISON}.count.filtered.csv" "\${CONTROL_KEYWORD}" "\${TREATMENT_KEYWORD}" "\${FDR_THRESHOLD}" "\${LOGFC_THRESHOLD}"`;
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
FDR_THRESHOLD="${fdrThreshold}"
LOGFC_THRESHOLD="${logfcThreshold}"

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

  script += `echo "FDR Threshold: \${FDR_THRESHOLD}"
echo "LogFC Threshold: \${LOGFC_THRESHOLD}"
echo ""
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

# Rename DE results file to standard name
# simpleEdger3.R creates files like "Group1vsGroup2.csv"
# We need "\${COMPARISON}_DE.csv" for the merge step
echo ""
echo "Renaming DE results file to standard name..."
DE_FILE=\$(ls *vs*.csv 2>/dev/null | head -1)
if [ -n "\$DE_FILE" ]; then
  mv "\$DE_FILE" "\${COMPARISON}_DE.csv"
  echo "Renamed: \$DE_FILE -> \${COMPARISON}_DE.csv"
else
  echo "ERROR: Could not find DE results file"
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
  "\${COMPARISON}" \\
  "\${CONTROL_KEYWORD}" \\
  "\${TREATMENT_KEYWORD}"

if [ $? -ne 0 ]; then
  echo "ERROR: merge_results.R failed"
  exit 2
fi

echo ""
echo "Step 3: Generate MA Plot for Visualization"
echo "========================================"
echo ""

# Generate MA plot with lab thresholds: logCPM=0, logFC=0.585, FDR=0.05
Rscript \${SCRIPTS_PATH}/maplot.R \\
  "\${COMPARISON}_DE.csv" \\
  0 \\
  0.585 \\
  0.05

if [ $? -ne 0 ]; then
  echo "WARNING: MA plot generation failed (non-critical)"
else
  echo "MA plot generated: \${COMPARISON}_DE.csvMaPlot.pdf"

  # Convert PDF to JPEG for agent viewing (GPT-5.2 and Claude can't read PDF)
  echo "Converting MA plot PDF to JPEG for agent review..."
  if command -v pdftoppm &> /dev/null; then
    pdftoppm -jpeg -r 150 -singlefile "\${COMPARISON}_DE.csvMaPlot.pdf" "\${COMPARISON}_maplot"
    if [ -f "\${COMPARISON}_maplot.jpg" ]; then
      echo "JPEG created: \${COMPARISON}_maplot.jpg"
    else
      echo "WARNING: JPEG conversion failed, using ImageMagick fallback..."
      convert -density 150 "\${COMPARISON}_DE.csvMaPlot.pdf" "\${COMPARISON}_maplot.jpg" 2>/dev/null || echo "WARNING: Could not convert PDF to JPEG"
    fi
  else
    echo "WARNING: pdftoppm not found, trying ImageMagick..."
    convert -density 150 "\${COMPARISON}_DE.csvMaPlot.pdf" "\${COMPARISON}_maplot.jpg" 2>/dev/null || echo "WARNING: Could not convert PDF to JPEG"
  fi
fi

echo ""
echo "========================================"
echo "Stage 4 Complete"
echo "========================================"
echo ""
echo "Output files:"
echo "  - \${COMPARISON}_DE.csv (DE results)"
echo "  - \${COMPARISON}_final.xlsx (merged RPKM + DE with formulas)"
echo "  - \${COMPARISON}_DE.csvMaPlot.pdf (MA plot visualization)"
echo "  - \${COMPARISON}_maplot.jpg (MA plot for agents)"
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
    maplot_pdf: null,
    maplot_jpg: null,
    num_degs_up: 0,
    num_degs_down: 0,
    total_degs: 0,
    fc_range: [null, null],
    fdr_threshold: 0.05,
    top_up_genes: [],
    top_down_genes: [],
    top_genes_by_pvalue: [],
    total_genes_tested: 0,
    classifications: {
      failed2DownRegulate: 0,
      failed2UpRegulate: 0,
      nchg: 0
    },
    fdr_distribution: {},
    logfc_distribution: {},
    logfc_range: [null, null],
    logfc_mean_abs: null,
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

    // 3b. Find MA plot files
    const maplotPdfFiles = fs.readdirSync(stage4Dir).filter(f => f.includes('MaPlot.pdf'));
    if (maplotPdfFiles.length > 0) {
      result.maplot_pdf = path.join(stage4Dir, maplotPdfFiles[0]);
    }

    const maplotJpgFiles = fs.readdirSync(stage4Dir).filter(f => f.includes('maplot.jpg'));
    if (maplotJpgFiles.length > 0) {
      result.maplot_jpg = path.join(stage4Dir, maplotJpgFiles[0]);
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
    let totalGenesInFile = 0;

    for (let i = 1; i < lines.length; i++) {
      const values = lines[i].split(',').map(v => v.trim().replace(/^"|"$/g, ''));

      if (values.length < headers.length) continue;

      totalGenesInFile++; // Count ALL genes in file

      const geneName = values[geneIndex];
      const logFC = logFCIndex !== -1 ? parseFloat(values[logFCIndex]) : null;
      const fdr = fdrIndex !== -1 ? parseFloat(values[fdrIndex]) : null;

      // Only add to genes array if we have valid numeric values for analysis
      if (!isNaN(logFC) && !isNaN(fdr)) {
        genes.push({ name: geneName, logFC, fdr });
      }
    }

    // Total genes tested = all genes in the file (including those with NA values)
    result.total_genes_tested = totalGenesInFile;

    // 5. Parse additional columns for classification (logCPM, PValue)
    const logCPMIndex = headers.findIndex(h => h.toLowerCase().includes('logcpm'));
    const pvalIndex = headers.findIndex(h => h.toLowerCase() === 'pvalue');

    // Re-parse with all columns for classification
    const fullGenes = [];
    for (let i = 1; i < lines.length; i++) {
      const values = lines[i].split(',').map(v => v.trim().replace(/^"|"$/g, ''));
      if (values.length < headers.length) continue;

      const geneName = values[geneIndex];
      const logFC = logFCIndex !== -1 ? parseFloat(values[logFCIndex]) : null;
      const fdr = fdrIndex !== -1 ? parseFloat(values[fdrIndex]) : null;
      const logCPM = logCPMIndex !== -1 ? parseFloat(values[logCPMIndex]) : null;
      const pvalue = pvalIndex !== -1 ? parseFloat(values[pvalIndex]) : null;

      if (!isNaN(logFC) && !isNaN(fdr)) {
        fullGenes.push({ name: geneName, logFC, fdr, logCPM, pvalue });
      }
    }

    // 6. Apply lab classification thresholds (matching Excel formulas)
    const LAB_FDR_THRESHOLD = 0.05;
    const LAB_LOGFC_THRESHOLD = 0.585;  // log2(1.5) for 1.5-fold change
    const LAB_LOGCPM_THRESHOLD = 0;

    // Note: We're simplifying Expr check - in Excel it checks max(avg(control), avg(treatment)) > 2
    // Here we just check if logCPM > 0 (gene is expressed)
    const classifiedGenes = fullGenes.map(g => {
      const isExpressed = (g.logCPM !== null && g.logCPM > LAB_LOGCPM_THRESHOLD);
      let classification = 'nchg';  // no change

      if (isExpressed && g.fdr < LAB_FDR_THRESHOLD) {
        if (g.logFC > LAB_LOGFC_THRESHOLD) {
          classification = 'failed2DownRegulate';  // Note: Keeping your lab's naming convention
        } else if (g.logFC < -LAB_LOGFC_THRESHOLD) {
          classification = 'failed2UpRegulate';
        }
      }

      return { ...g, classification, isExpressed };
    });

    // Count classifications
    const failed2Down = classifiedGenes.filter(g => g.classification === 'failed2DownRegulate').length;
    const failed2Up = classifiedGenes.filter(g => g.classification === 'failed2UpRegulate').length;
    const noChange = classifiedGenes.filter(g => g.classification === 'nchg').length;

    result.classifications = {
      failed2DownRegulate: failed2Down,
      failed2UpRegulate: failed2Up,
      nchg: noChange
    };

    // 7. Calculate distribution stats (for agent decision-making)
    result.fdr_distribution = {
      fdr_0_05: fullGenes.filter(g => g.fdr < 0.05).length,
      fdr_0_10: fullGenes.filter(g => g.fdr < 0.10).length,
      fdr_0_20: fullGenes.filter(g => g.fdr < 0.20).length,
      fdr_0_50: fullGenes.filter(g => g.fdr < 0.50).length
    };

    result.logfc_distribution = {
      abs_logfc_0_5: fullGenes.filter(g => Math.abs(g.logFC) > 0.5).length,
      abs_logfc_1_0: fullGenes.filter(g => Math.abs(g.logFC) > 1.0).length,
      abs_logfc_2_0: fullGenes.filter(g => Math.abs(g.logFC) > 2.0).length
    };

    if (fullGenes.length > 0) {
      const allLogFCs = fullGenes.map(g => g.logFC);
      result.logfc_range = [Math.min(...allLogFCs).toFixed(2), Math.max(...allLogFCs).toFixed(2)];
      result.logfc_mean_abs = (fullGenes.reduce((sum, g) => sum + Math.abs(g.logFC), 0) / fullGenes.length).toFixed(2);
    }

    // 8. Get top genes by p-value (not FDR) - shows most promising genes
    const genesWithPval = fullGenes.filter(g => g.pvalue !== null && !isNaN(g.pvalue));
    const topByPvalue = genesWithPval.sort((a, b) => a.pvalue - b.pvalue).slice(0, 10);
    result.top_genes_by_pvalue = topByPvalue.map(g => ({
      name: g.name,
      logFC: g.logFC.toFixed(2),
      pvalue: g.pvalue.toExponential(2),
      fdr: g.fdr.toExponential(2)
    }));

    // 9. Filter DEGs (FDR < 0.05 by default - for backward compatibility)
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

    // 9. Parse threshold distribution (for re-analysis agent review)
    const thresholdDistFile = path.join(stage4Dir, 'threshold_distribution.json');
    if (fs.existsSync(thresholdDistFile)) {
      try {
        const thresholdDistData = fs.readFileSync(thresholdDistFile, 'utf8');
        result.threshold_distribution = JSON.parse(thresholdDistData);
      } catch (error) {
        result.warnings.push('Could not parse threshold_distribution.json');
      }
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
