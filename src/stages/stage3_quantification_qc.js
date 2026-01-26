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

# OPTIMIZATION: Check if count matrix already exists in BAM directory (from previous run)
# featureCounts is deterministic, so we can reuse existing count files to save time
EXISTING_COUNT_FILE="\${BAM_DIR}/\${COMPARISON}.count.csv"
if [ -f "\${EXISTING_COUNT_FILE}" ]; then
  echo "[Stage 3] Found existing count matrix: \${EXISTING_COUNT_FILE}"
  echo "[Stage 3] Copying existing count file (deterministic, no need to regenerate)"
  cp "\${EXISTING_COUNT_FILE}" "\${OUTPUT_DIR}/\${COMPARISON}.count.csv"

  if [ $? -ne 0 ]; then
    echo "WARNING: Failed to copy existing count file, will regenerate"
    # Fall through to run featureCounts
  else
    echo "[Stage 3] ✓ Count matrix copied successfully, skipping featureCounts"
    # Skip featureCounts since we copied the file
    echo "" > /dev/null  # Placeholder to maintain script structure
  fi
fi

# Run featureCounts only if count file doesn't exist yet
if [ ! -f "\${OUTPUT_DIR}/\${COMPARISON}.count.csv" ]; then
  echo "[Stage 3] Running featureCounts to generate count matrix..."
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
else
  echo "[Stage 3] Count matrix already exists, skipping featureCounts"
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
  "\${COMPARISON}" \\
  "\${CONTROL_KEYWORD}" \\
  "\${TREATMENT_KEYWORD}"

QC_EXIT_CODE=$?

# Note: R script now saves both PDF and JPEG directly
# No conversion needed - both formats generated by qc_assessment_pca.R

# Cleanup: Delete BAM files to save disk space
# BAM files are no longer needed after quantification (featureCounts) is complete
echo ""
echo "========================================"
echo "Cleanup: Deleting BAM files to save space"
echo "========================================"
BAM_DIR="${outputDir}/stage2_alignment/bam_files"
echo "[DEBUG] BAM_DIR set to: \$BAM_DIR"

if [ -d "\$BAM_DIR" ]; then
  echo "[DEBUG] BAM directory exists"
  BAM_COUNT=\$(find "\$BAM_DIR" -name "*.bam" | wc -l)
  BAM_SIZE=\$(du -sh "\$BAM_DIR" 2>/dev/null | cut -f1)

  echo "Found \$BAM_COUNT BAM files (\$BAM_SIZE total)"
  echo "Deleting BAM files and VCF files..."

  echo "[DEBUG] Running: find \"\$BAM_DIR\" -name \"*.bam\" -delete"
  find "\$BAM_DIR" -name "*.bam" -delete
  echo "[DEBUG] BAM files deleted"

  echo "[DEBUG] Running: find \"\$BAM_DIR\" -name \"*.bam.bai\" -delete"
  find "\$BAM_DIR" -name "*.bam.bai" -delete 2>/dev/null
  echo "[DEBUG] BAI files deleted"

  echo "[DEBUG] Running: find \"\$BAM_DIR\" -name \"*.indel.vcf\" -delete"
  find "\$BAM_DIR" -name "*.indel.vcf" -delete 2>/dev/null
  echo "[DEBUG] VCF files deleted"

  # Verify deletion
  REMAINING_BAM=\$(find "\$BAM_DIR" -name "*.bam" | wc -l)
  echo "[DEBUG] Remaining BAM files after deletion: \$REMAINING_BAM"

  echo "✓ BAM and VCF files deleted successfully"
else
  echo "[ERROR] BAM directory not found: \$BAM_DIR"
  echo "[DEBUG] Checking if parent directory exists..."
  ls -ld "\$(dirname \"\$BAM_DIR\")" 2>&1 || echo "[DEBUG] Parent directory also doesn't exist"
fi

echo "[DEBUG] Cleanup section completed"

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
echo "  - \${COMPARISON}_pca_plot.pdf (PCA visualization - for user)"
echo "  - \${COMPARISON}_pca_plot_qc.jpg (PCA visualization - for agents)"
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
    pca_plot_path: null,  // Path to PCA plot JPEG for agent review (PDF kept for user)
    sample_dictionary: null,  // Sample name dictionary (A1->full_name, B1->full_name, etc.)
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

      // Extract outliers (new format from qc_assessment_pca.R)
      // Looks for lines like: "- GSM1275870_N080611_untreated (Control)"
      // Or any sample name format: "- Sample1_control_rep1 (Control)"
      const outlierLines = qcSummary.match(/^\s*-\s+(\S+)\s+\(/gm);
      if (outlierLines && outlierLines.length > 0) {
        result.outliers_detected = outlierLines.map(line => {
          const match = line.match(/^\s*-\s+(\S+)\s+\(/);
          return match ? match[1].trim() : null;
        }).filter(Boolean);
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

    // 6. Find PCA plot path (JPEG for agents, PDF kept for user)
    const pcaPlotPathJPEG = path.join(stage3Dir, `${comparison}_pca_plot_qc.jpg`);
    const pcaPlotPathPDF = path.join(stage3Dir, `${comparison}_pca_plot.pdf`);

    if (fs.existsSync(pcaPlotPathJPEG)) {
      result.pca_plot_path = pcaPlotPathJPEG;
    } else if (fs.existsSync(pcaPlotPathPDF)) {
      // Fallback to PDF if JPEG conversion failed
      result.pca_plot_path = pcaPlotPathPDF;
      result.warnings.push('JPEG conversion failed - using PDF (may not work for all agents)');
    }

    // 7. Read sample name dictionary (A1->full_name, B1->full_name)
    const dictionaryPath = path.join(stage3Dir, `${comparison}_sample_dictionary.txt`);
    if (fs.existsSync(dictionaryPath)) {
      result.sample_dictionary = fs.readFileSync(dictionaryPath, 'utf-8');
    } else {
      result.warnings.push('Sample dictionary not found - agents will see full sample names in plot');
    }

    // 8. Set default recommendation if not found
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

  formatted += `**IMPORTANT:** A PCA plot (PC1 vs PC2) will be provided as an image attachment.\n`;
  formatted += `You must VISUALLY EXAMINE the PCA plot to assess batch effects and outliers.\n\n`;

  // Sample Name Dictionary (if available)
  if (parsedOutput.sample_dictionary) {
    formatted += `## SAMPLE NAME DICTIONARY\n\n`;
    formatted += `**CRITICAL:** The PCA plot uses SHORT LABELS (A1, A2, B1, B2, etc.) for clarity.\n`;
    formatted += `Use this dictionary to identify the FULL sample names:\n\n`;
    formatted += `\`\`\`\n`;
    formatted += parsedOutput.sample_dictionary;
    formatted += `\`\`\`\n\n`;
    formatted += `When reporting outliers, you MUST specify BOTH the plot label (e.g., "B1") AND the full sample name from the dictionary.\n\n`;
  }

  // Dataset overview
  formatted += `## Dataset Overview\n`;
  formatted += `- Organism: ${dataInfo.organism || 'unknown'}\n`;
  formatted += `- Samples Analyzed: ${parsedOutput.samples_analyzed || samples.length}\n`;
  formatted += `- Genes Detected: ${parsedOutput.genes_detected?.toLocaleString() || 'unknown'}\n`;
  formatted += `- Sequencing Type: ${dataInfo.pairedEnd ? 'Paired-end' : 'Single-end'}\n`;
  formatted += `\n`;

  // Sample Groups (keep for reference, but note they're in the dictionary now)
  formatted += `## Sample Groups\n`;
  for (const [group, sampleNames] of Object.entries(groups)) {
    formatted += `- ${group}: ${sampleNames.length} samples\n`;
    if (!parsedOutput.sample_dictionary) {
      // Only show full names if no dictionary (fallback)
      formatted += `  (${sampleNames.join(', ')})\n`;
    }
  }
  formatted += `\n`;

  // PCA Variance
  formatted += `## PCA Variance Explained\n`;
  if (parsedOutput.pc1_variance !== null) {
    formatted += `- PC1: ${(parsedOutput.pc1_variance * 100).toFixed(1)}%\n`;
  }
  if (parsedOutput.pc2_variance !== null) {
    formatted += `- PC2: ${(parsedOutput.pc2_variance * 100).toFixed(1)}%\n`;
  }
  formatted += `\n`;

  // PCA Plot Path (for reference)
  if (parsedOutput.pca_plot_path) {
    formatted += `## PCA Plot Location\n`;
    formatted += `\`${parsedOutput.pca_plot_path}\`\n`;
    formatted += `\n`;
  }

  // Warnings and Errors (if any)
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
  formatted += `Based on your VISUAL INSPECTION of the PCA plot, determine:\n`;
  formatted += `1. Whether batch effects are present (look for sub-clusters within groups)\n`;
  formatted += `2. Whether any samples are outliers (samples far from their group cluster)\n`;
  formatted += `3. Which DE analysis method to use (simpleEdger or batch_effect_edger)\n`;
  formatted += `4. Whether to keep all samples or remove outliers\n\n`;

  if (parsedOutput.sample_dictionary) {
    formatted += `**REMINDER:** When identifying outliers, use the Sample Name Dictionary above to provide:\n`;
    formatted += `- The plot label (e.g., "A1", "B2")\n`;
    formatted += `- The full sample name (e.g., "GSM1275870_N052611_untreated")\n\n`;
  }

  return formatted;
}

export default {
  generateStage3Script,
  parseStage3Output,
  formatStage3ForAgents
};
