/**
 * Stage 3: Quantification + QC Assessment
 *
 * This stage:
 * 1. Runs featureCounts to quantify gene expression
 * 2. Filters problematic gene IDs
 * 3. Normalizes with RPKM
 * 4. Adds gene symbols (Entrez)
 * 5. Runs PCA/QC assessment to detect outliers and batch effects
 *
 * Generates script → Parses output → Formats for agents
 */

import fs from 'fs';
import path from 'path';

/**
 * Generate Stage 3 script: Quantification + QC Assessment
 *
 * @param {Object} dataInfo - Dataset information (samples, organism, etc.)
 * @param {Object} config - Configuration (paths, comparison name, etc.)
 * @param {string} outputPath - Path to write the generated script
 * @param {Object} stage2Result - Results from Stage 2 (for BAM file locations)
 * @returns {Object} - { scriptPath, scriptContent }
 */
export function generateStage3Script(dataInfo, config, outputPath, stage2Result = {}) {
  const scriptsPath = process.env.SCRIPTS_PATH || '/users/ha00014/bin';
  const condaEnv = process.env.CONDA_ENV || 'rnaseq';

  // Determine organism annotation file
  const annotationMap = {
    mouse: 'mm10',
    human: 'hg38',
    rat: 'rn6'
  };
  const genomeBuild = annotationMap[dataInfo.organism?.toLowerCase()] || 'mm10';

  // Paths
  const outputDir = path.dirname(outputPath);
  const stage3Dir = path.join(outputDir, 'stage3_quantification');
  const bamDir = path.join(outputDir, 'stage2_alignment', 'bam_files');

  // Choose featureCounts script based on paired-end
  const featureCountsScript = dataInfo.pairedEnd ? 'featurecounts.R' : 'featurecounts_unpaired.R';

  // Comparison name
  const comparison = config.comparison || 'comparison';
  const controlKeyword = config.controlKeyword;
  const treatmentKeyword = config.treatmentKeyword;

  // Build script
  let script = `#!/bin/bash
#
# Stage 3: Quantification + QC Assessment
# Generated: ${new Date().toISOString()}
#

set -e  # Exit on error

echo "========================================"
echo "Stage 3: Quantification + QC Assessment"
echo "========================================"
echo ""

# Environment
SCRIPTS_PATH="${scriptsPath}"
CONDA_ENV="${condaEnv}"
OUTPUT_DIR="${stage3Dir}"
BAM_DIR="${bamDir}"
COMPARISON="${comparison}"
CONTROL_KEYWORD="${controlKeyword}"
TREATMENT_KEYWORD="${treatmentKeyword}"
ORGANISM="${dataInfo.organism}"
GENOME_BUILD="${genomeBuild}"

# Activate conda environment (disable -u temporarily to avoid conda activation errors)
echo "[Stage 3] Activating conda environment: \${CONDA_ENV}"
set +u
source $(conda info --base)/etc/profile.d/conda.sh
conda activate \${CONDA_ENV}
set -u

# Create output directory
mkdir -p "\${OUTPUT_DIR}"
cd "\${OUTPUT_DIR}"

echo ""
echo "Step 1: Feature Counts (Quantification)"
echo "========================================"
echo ""

# Run featureCounts
Rscript \${SCRIPTS_PATH}/${featureCountsScript} \\
  "\${GENOME_BUILD}" \\
  "\${COMPARISON}" \\
  "\${CONTROL_KEYWORD}" \\
  "\${TREATMENT_KEYWORD}" \\
  "\${BAM_DIR}"/*.bam

if [ $? -ne 0 ]; then
  echo "ERROR: featureCounts failed"
  exit 2
fi

echo ""
echo "Step 2: Filter Bad Gene IDs"
echo "========================================"
echo ""

# Run filterIDS.R
Rscript \${SCRIPTS_PATH}/filterIDS.R \\
  "\${COMPARISON}.count.csv"

if [ $? -ne 0 ]; then
  echo "ERROR: filterIDS.R failed"
  exit 2
fi

echo ""
echo "Step 3: RPKM Normalization"
echo "========================================"
echo ""

# Run RPKM normalization
Rscript \${SCRIPTS_PATH}/RPKM.R \\
  "\${COMPARISON}.count.filtered.csv"

if [ $? -ne 0 ]; then
  echo "ERROR: RPKM normalization failed"
  exit 2
fi

echo ""
echo "Step 4: Add Gene Symbols (Entrez)"
echo "========================================"
echo ""

# Add Entrez gene symbols
Rscript \${SCRIPTS_PATH}/entrz.R \\
  "\${COMPARISON}.rpkm.csv" \\
  "\${GENOME_BUILD}"

if [ $? -ne 0 ]; then
  echo "ERROR: entrz.R failed"
  exit 2
fi

echo ""
echo "Step 5: QC Assessment (PCA Analysis)"
echo "========================================"
echo ""

# Run PCA/QC assessment
# This detects:
# - Outliers (robust distance-based detection)
# - Batch effects (median distance to centroid)
# Exit codes:
#   0 = All clean (no outliers, no batch effects)
#   2 = Issues detected (outliers or batch effects found)

Rscript \${SCRIPTS_PATH}/qc_assessment_pca.R \\
  "\${COMPARISON}.count.filtered.csv" \\
  "\${COMPARISON}"

QC_EXIT_CODE=$?

echo ""
echo "========================================"
echo "Stage 3 Complete"
echo "========================================"
echo ""
echo "Output files:"
echo "  - \${COMPARISON}.count.csv (raw counts)"
echo "  - \${COMPARISON}.count.filtered.csv (filtered counts)"
echo "  - \${COMPARISON}.rpkm.csv (RPKM normalized)"
echo "  - \${COMPARISON}.rpkm.entrz.csv (with gene symbols)"
echo "  - \${COMPARISON}_qc_summary.txt (QC report)"
echo "  - \${COMPARISON}_pca_plot.pdf (PCA visualization)"
echo ""

if [ \$QC_EXIT_CODE -eq 0 ]; then
  echo "QC Status: PASS (no issues detected)"
  exit 0
elif [ \$QC_EXIT_CODE -eq 2 ]; then
  echo "QC Status: ISSUES DETECTED (outliers or batch effects found)"
  echo "Review QC summary and decide: keep all samples or remove outliers"
  exit 2
else
  echo "QC Status: ERROR (unexpected exit code: \$QC_EXIT_CODE)"
  exit 2
fi
`;

  // Write script to file
  fs.writeFileSync(outputPath, script, { mode: 0o755 });

  return {
    scriptPath: outputPath,
    scriptContent: script
  };
}

/**
 * Parse Stage 3 output: Read QC assessment results
 *
 * @param {string} outputDir - Output directory containing stage3_quantification/
 * @returns {Object} - Parsed QC results with outliers, batch effects, etc.
 */
export function parseStage3Output(outputDir) {
  const stage3Dir = path.join(outputDir, 'stage3_quantification');

  // Initialize result structure
  const result = {
    overall_status: 'UNKNOWN',
    exit_code: 0,
    outliers_detected: [],
    batch_effect_detected: false,
    batch_effect_severity: 'none',
    batch_effect_evidence: '',
    per_sample_qc: {},
    pc1_variance: null,
    pc2_variance: null,
    samples_analyzed: 0,
    genes_detected: 0,
    recommendation: '',
    warnings: [],
    errors: []
  };

  try {
    // 1. Check if output directory exists
    if (!fs.existsSync(stage3Dir)) {
      result.errors.push(`Stage 3 output directory not found: ${stage3Dir}`);
      result.overall_status = 'ERROR';
      return result;
    }

    // 2. Find comparison name from count file
    const countFiles = fs.readdirSync(stage3Dir).filter(f => f.endsWith('.count.filtered.csv'));
    if (countFiles.length === 0) {
      result.errors.push('No count file found (*.count.filtered.csv)');
      result.overall_status = 'ERROR';
      return result;
    }

    const comparison = countFiles[0].replace('.count.filtered.csv', '');

    // 3. Parse QC summary file
    const qcSummaryPath = path.join(stage3Dir, `${comparison}_qc_summary.txt`);
    if (fs.existsSync(qcSummaryPath)) {
      const qcSummary = fs.readFileSync(qcSummaryPath, 'utf-8');

      // Extract outliers
      const outlierMatch = qcSummary.match(/Outliers detected:\s*([^\n]+)/i);
      if (outlierMatch && outlierMatch[1].trim().toLowerCase() !== 'none') {
        result.outliers_detected = outlierMatch[1].split(',').map(s => s.trim());
      }

      // Extract batch effect detection
      const batchEffectMatch = qcSummary.match(/Batch effects detected:\s*(yes|no)/i);
      if (batchEffectMatch) {
        result.batch_effect_detected = batchEffectMatch[1].toLowerCase() === 'yes';
      }

      // Extract batch effect severity
      const severityMatch = qcSummary.match(/Batch effect severity:\s*([^\n]+)/i);
      if (severityMatch) {
        result.batch_effect_severity = severityMatch[1].trim().toLowerCase();
      }

      // Extract evidence
      const evidenceMatch = qcSummary.match(/Evidence:\s*([^\n]+)/i);
      if (evidenceMatch) {
        result.batch_effect_evidence = evidenceMatch[1].trim();
      }

      // Extract PCA variance
      const pc1Match = qcSummary.match(/PC1 variance:\s*([\d.]+)%/i);
      if (pc1Match) {
        result.pc1_variance = parseFloat(pc1Match[1]) / 100;
      }

      const pc2Match = qcSummary.match(/PC2 variance:\s*([\d.]+)%/i);
      if (pc2Match) {
        result.pc2_variance = parseFloat(pc2Match[1]) / 100;
      }

      // Extract recommendation
      const recMatch = qcSummary.match(/Recommendation:\s*([^\n]+)/i);
      if (recMatch) {
        result.recommendation = recMatch[1].trim();
      }

      // Determine exit code from summary
      if (result.outliers_detected.length > 0 || result.batch_effect_detected) {
        result.exit_code = 2;
        result.overall_status = 'ISSUES_DETECTED';
      } else {
        result.exit_code = 0;
        result.overall_status = 'PASS';
      }
    } else {
      result.warnings.push('QC summary file not found - using defaults');
    }

    // 4. Parse QC metrics CSV (if exists)
    const qcMetricsPath = path.join(stage3Dir, `${comparison}_qc_metrics.csv`);
    if (fs.existsSync(qcMetricsPath)) {
      result.per_sample_qc = parseQcMetricsCsv(qcMetricsPath);
      result.samples_analyzed = Object.keys(result.per_sample_qc).length;
    }

    // 5. Count genes from count file
    const countFilePath = path.join(stage3Dir, `${comparison}.count.filtered.csv`);
    if (fs.existsSync(countFilePath)) {
      const countData = fs.readFileSync(countFilePath, 'utf-8');
      const lines = countData.split('\n').filter(l => l.trim().length > 0);
      result.genes_detected = Math.max(0, lines.length - 1); // Subtract header
    }

    // 6. Set default recommendation if not found
    if (!result.recommendation) {
      if (result.overall_status === 'PASS') {
        result.recommendation = 'No issues detected. Proceed with simpleEdger (standard DE analysis).';
      } else if (result.batch_effect_detected) {
        result.recommendation = 'Batch effects detected. Use batch_effect_edger with batch correction.';
      } else if (result.outliers_detected.length > 0) {
        result.recommendation = 'Outliers detected. Consider removing outliers or proceed with caution.';
      }
    }

  } catch (error) {
    result.errors.push(`Error parsing Stage 3 output: ${error.message}`);
    result.overall_status = 'ERROR';
  }

  return result;
}

/**
 * Parse QC metrics CSV file
 */
function parseQcMetricsCsv(csvPath) {
  const metrics = {};

  try {
    const csvData = fs.readFileSync(csvPath, 'utf-8');
    const lines = csvData.split('\n').filter(l => l.trim().length > 0);

    if (lines.length < 2) return metrics;

    const headers = lines[0].split(',').map(h => h.trim());

    for (let i = 1; i < lines.length; i++) {
      const values = lines[i].split(',').map(v => v.trim());
      const sampleName = values[0];

      if (!sampleName) continue;

      metrics[sampleName] = {};
      for (let j = 1; j < headers.length; j++) {
        const key = headers[j];
        const value = values[j];

        // Try to parse as number
        const numValue = parseFloat(value);
        metrics[sampleName][key] = isNaN(numValue) ? value : numValue;
      }
    }
  } catch (error) {
    console.error(`Error parsing QC metrics CSV: ${error.message}`);
  }

  return metrics;
}

/**
 * Format Stage 3 output for agent consumption
 *
 * @param {Object} parsedOutput - Parsed QC assessment results
 * @param {Object} dataInfo - Dataset information
 * @returns {string} - Markdown formatted output for agents
 */
export function formatStage3ForAgents(parsedOutput, dataInfo) {
  const samples = dataInfo.samples || [];
  const groups = dataInfo.groups || {};

  let formatted = `# Stage 3: Quantification + QC Assessment Results\n\n`;

  // Dataset overview
  formatted += `## Dataset Overview\n`;
  formatted += `- Organism: ${dataInfo.organism || 'unknown'}\n`;
  formatted += `- Samples Analyzed: ${parsedOutput.samples_analyzed || samples.length}\n`;
  formatted += `- Genes Detected: ${parsedOutput.genes_detected?.toLocaleString() || 'unknown'}\n`;
  formatted += `- Sequencing Type: ${dataInfo.pairedEnd ? 'Paired-end' : 'Single-end'}\n`;
  formatted += `\n`;

  // PCA Results
  formatted += `## PCA Results\n`;
  if (parsedOutput.pc1_variance !== null) {
    formatted += `- PC1 Variance Explained: ${(parsedOutput.pc1_variance * 100).toFixed(1)}%\n`;
  }
  if (parsedOutput.pc2_variance !== null) {
    formatted += `- PC2 Variance Explained: ${(parsedOutput.pc2_variance * 100).toFixed(1)}%\n`;
  }
  formatted += `\n`;

  // Batch Effect Detection
  formatted += `## Batch Effect Detection\n`;
  formatted += `- Batch Effect Detected: **${parsedOutput.batch_effect_detected ? 'YES' : 'NO'}**\n`;
  if (parsedOutput.batch_effect_detected) {
    formatted += `- Severity: ${parsedOutput.batch_effect_severity}\n`;
    formatted += `- Evidence: ${parsedOutput.batch_effect_evidence}\n`;
  }
  formatted += `\n`;

  // Outlier Detection
  formatted += `## Outlier Detection\n`;
  if (parsedOutput.outliers_detected.length > 0) {
    formatted += `- Outliers Detected: **${parsedOutput.outliers_detected.length}**\n`;
    formatted += `- Outlier Samples:\n`;
    for (const sample of parsedOutput.outliers_detected) {
      // Find which group this sample belongs to
      let sampleGroup = 'unknown';
      for (const [group, sampleNames] of Object.entries(groups)) {
        if (sampleNames.includes(sample)) {
          sampleGroup = group;
          break;
        }
      }
      formatted += `  - ${sample} (${sampleGroup})\n`;
    }
  } else {
    formatted += `- Outliers Detected: **None**\n`;
  }
  formatted += `\n`;

  // Per-sample QC metrics (if available)
  if (Object.keys(parsedOutput.per_sample_qc).length > 0) {
    formatted += `## Per-Sample QC Metrics\n`;
    formatted += `| Sample | PC1 | PC2 | Distance to Centroid | Outlier |\n`;
    formatted += `|--------|-----|-----|---------------------|----------|\n`;

    for (const [sample, metrics] of Object.entries(parsedOutput.per_sample_qc)) {
      const isOutlier = parsedOutput.outliers_detected.includes(sample);
      formatted += `| ${sample} | `;
      formatted += `${metrics.PC1?.toFixed(2) || 'N/A'} | `;
      formatted += `${metrics.PC2?.toFixed(2) || 'N/A'} | `;
      formatted += `${metrics.distance?.toFixed(2) || 'N/A'} | `;
      formatted += `${isOutlier ? '⚠️ YES' : 'No'} |\n`;
    }
    formatted += `\n`;
  }

  // Group Balance Check
  formatted += `## Group Balance After Potential Outlier Removal\n`;
  const currentCounts = {};
  for (const [group, sampleNames] of Object.entries(groups)) {
    const remaining = sampleNames.filter(s => !parsedOutput.outliers_detected.includes(s));
    currentCounts[group] = remaining.length;
  }

  for (const [group, count] of Object.entries(currentCounts)) {
    formatted += `- ${group}: ${count} samples`;
    if (count < 2) {
      formatted += ` ⚠️ **WARNING: Less than 2 samples!**`;
    }
    formatted += `\n`;
  }
  formatted += `\n`;

  // QC Status and Recommendation
  formatted += `## QC Status\n`;
  formatted += `- Overall Status: **${parsedOutput.overall_status}**\n`;
  formatted += `- Exit Code: ${parsedOutput.exit_code} (0=PASS, 2=Issues)\n`;
  formatted += `- Recommendation: ${parsedOutput.recommendation}\n`;
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

  // Decision Prompt
  formatted += `---\n\n`;
  formatted += `## Your Decision (REQUIRED FORMAT)\n\n`;
  formatted += `Please review the QC assessment results above and provide your decision:\n\n`;

  formatted += `**DE_Method:** [simpleEdger / batch_effect_edger]\n`;
  formatted += `- Choose simpleEdger if NO batch effects detected\n`;
  formatted += `- Choose batch_effect_edger if batch effects detected\n\n`;

  if (parsedOutput.batch_effect_detected) {
    formatted += `**Batch_Specification:** [auto / paired / explicit like "1,1,2,2"]\n`;
    formatted += `- Use "auto" to let the pipeline infer batch structure\n`;
    formatted += `- Use "paired" if samples were processed in pairs (ctrl1+trt1, ctrl2+trt2)\n`;
    formatted += `- Use explicit labels for custom batch structure\n\n`;
  } else {
    formatted += `**Batch_Specification:** [N/A]\n\n`;
  }

  formatted += `**Outlier_Action:** [KEEP_ALL / REMOVE_OUTLIERS]\n`;
  formatted += `- KEEP_ALL: Keep all samples (even outliers)\n`;
  formatted += `- REMOVE_OUTLIERS: Remove detected outliers\n\n`;

  if (parsedOutput.outliers_detected.length > 0) {
    formatted += `**Outliers_to_Remove:** [${parsedOutput.outliers_detected.join(', ')} / None]\n`;
    formatted += `- List specific outliers to remove, or "None" to keep all\n\n`;
  } else {
    formatted += `**Outliers_to_Remove:** [None]\n\n`;
  }

  formatted += `**Confidence:** [HIGH / MEDIUM / LOW]\n\n`;
  formatted += `**Reasoning:** [2-3 sentences explaining your assessment]\n`;

  return formatted;
}

export default {
  generateStage3Script,
  parseStage3Output,
  formatStage3ForAgents
};
