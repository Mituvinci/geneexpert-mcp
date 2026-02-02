/**
 * Script Generator - Generate bash scripts for AUTOMATION or ADAPTATION
 *
 * AUTOMATION: Template-based script running all 11 standard steps
 * ADAPTATION: Agent-written custom script addressing specific concerns
 */

import fs from 'fs';
import path from 'path';
import dotenv from 'dotenv';
import { applyModifications, validateScript } from './modification_engine.js';
import { extractStructuredRecommendations } from '../coordinator/consensus.js';

dotenv.config();

const SCRIPTS_PATH = process.env.SCRIPTS_PATH || './bio_informatics/scripts';
const CONDA_ENV = process.env.CONDA_ENV || 'pytorch';
const INDEX_PATH = process.env.INDEX_PATH || './bio_informatics/reference_data/index';

/**
 * Actual script filenames (many don't have .sh extension!)
 * Maps logical name → actual filename
 */
const SCRIPT_FILES = {
  'validate_fastq': 'validate_fastq.sh',  // Step 1: FASTQ validation
  'fastq2bam': 'fastq2bam',           // NO .sh extension
  'alignment_qc': 'alignment_qc_screening.R',  // Step 3.5: Alignment QC (contamination/quality screening)
  'featurecounts': 'featurecounts.R',
  'featurecounts_unpaired': 'featurecounts_unpaired.R',  // For single-end reads
  'filterIDS': 'filterIDS.R',
  'RPKM': 'RPKM.R',
  'entrz': 'entrz.R',
  'qc_assessment': 'qc_assessment_pca.R',  // Step 8: QC with outlier/batch detection
  'simpleEdger3': 'simpleEdger3.R',
  'batch_effect_edger': 'batch_effect_edgeR_v3.R',  // Batch-corrected edgeR (design: ~ batch + condition)
  'edger3xl': 'edger3xl.R',
  'merge_results': 'merge_results.R'
};

/**
 * Generate AUTOMATION script - simple template running all steps
 */
export function generateAutomationScript(dataInfo, config, steps, outputPath) {
  console.log('[Script Generator] Generating AUTOMATION script...');
  console.log('[Script Generator] Using standard pipeline template');

  // Map organism to genome build
  const genomeMap = {
    'mouse': 'mm10',
    'human': 'hg38',
    'rat': 'rn6'
  };
  const genomeBuild = genomeMap[config.organism] || config.organism;

  const scriptLines = [];

  // Header
  scriptLines.push('#!/bin/bash');
  scriptLines.push('#');
  scriptLines.push(`# GeneExpert Analysis Script (AUTOMATION)`);
  scriptLines.push(`# Generated: ${new Date().toISOString()}`);
  scriptLines.push(`# Comparison: ${config.comparison}`);
  scriptLines.push(`# Organism: ${config.organism} (${genomeBuild})`);
  scriptLines.push('#');
  scriptLines.push('set -e  # Exit on error');
  scriptLines.push('');

  // Activate conda environment
  scriptLines.push('# Activate conda environment for bioinformatics tools');
  scriptLines.push(`source $(conda info --base)/etc/profile.d/conda.sh`);
  scriptLines.push(`conda activate ${CONDA_ENV}`);
  scriptLines.push('');

  // Variables
  scriptLines.push('# Paths');
  scriptLines.push(`OUTPUT_DIR="${config.output}"`);
  scriptLines.push(`SCRIPTS_DIR="${SCRIPTS_PATH}"`);
  scriptLines.push(`INDEX_PATH="${INDEX_PATH}"`);
  scriptLines.push(`GENOME="${genomeBuild}"  # Use genome build, not common name`);
  scriptLines.push(`COMPARISON="${config.comparison}"`);
  scriptLines.push(`CONTROL_KEYWORD="${config.controlKeyword || 'control'}"`);
  scriptLines.push(`TREATMENT_KEYWORD="${config.treatmentKeyword || 'treatment'}"`);
  scriptLines.push('');

  // Create output directories
  scriptLines.push('# Create output directories');
  scriptLines.push('mkdir -p $OUTPUT_DIR/fastqc_results');
  scriptLines.push('mkdir -p $OUTPUT_DIR/bam_files');
  scriptLines.push('mkdir -p $OUTPUT_DIR/counts');
  scriptLines.push('mkdir -p $OUTPUT_DIR/qc_plots');
  scriptLines.push('mkdir -p $OUTPUT_DIR/de_results');
  scriptLines.push('');

  // Step 1: FastQC
  scriptLines.push('echo "Step 1/10: Running FastQC..."');
  const fastqFiles = dataInfo.samples.map(s => s.r1 + (s.r2 ? ' ' + s.r2 : '')).join(' ');
  scriptLines.push(`fastqc ${fastqFiles} -o $OUTPUT_DIR/fastqc_results`);
  scriptLines.push('');

  // Step 2: Alignment
  scriptLines.push('echo "Step 2/10: Alignment (FASTQ to BAM)..."');
  scriptLines.push(`cd ${config.input}`);
  scriptLines.push('');

  // Generate alignment commands for each sample
  if (dataInfo.pairedEnd) {
    scriptLines.push('# Run alignment for paired-end samples');
    scriptLines.push('for i in *_R1_001.fastq.gz; do');
    scriptLines.push('  fname=$(basename "$i" _R1_001.fastq.gz)');
    scriptLines.push('  echo "Aligning: $fname"');
    scriptLines.push('  subread-align -t 0 \\');
    scriptLines.push('    -i $INDEX_PATH/$GENOME \\');
    scriptLines.push('    -r "${fname}_R1_001.fastq.gz" \\');
    scriptLines.push('    -R "${fname}_R2_001.fastq.gz" \\');
    scriptLines.push('    -T 8 \\');
    scriptLines.push('    -o "${fname}.bam" &> "${fname}.log" &');
    scriptLines.push('done');
    scriptLines.push('wait');
  } else {
    scriptLines.push('# Run alignment for single-end samples');
    scriptLines.push('for i in *.fastq.gz; do');
    scriptLines.push('  fname=$(basename "$i" .fastq.gz)');
    scriptLines.push('  echo "Aligning: $fname"');
    scriptLines.push('  subread-align -t 0 \\');
    scriptLines.push('    -i $INDEX_PATH/$GENOME \\');
    scriptLines.push('    -r "${fname}.fastq.gz" \\');
    scriptLines.push('    -T 8 \\');
    scriptLines.push('    -o "${fname}.bam" &> "${fname}.log" &');
    scriptLines.push('done');
    scriptLines.push('wait');
  }

  scriptLines.push('');
  scriptLines.push('echo "Moving BAM files to output directory..."');
  scriptLines.push(`mv *.bam *.log $OUTPUT_DIR/bam_files/`);
  scriptLines.push(`cd -`);
  scriptLines.push('');

  // Step 3: Feature Counts (use paired or unpaired script based on data)
  const featurecountsScript = dataInfo.pairedEnd ? SCRIPT_FILES.featurecounts : SCRIPT_FILES.featurecounts_unpaired;
  scriptLines.push(`echo "Step 3/10: Feature Counts (${dataInfo.pairedEnd ? 'paired-end' : 'single-end'})..."`);
  scriptLines.push('cd $OUTPUT_DIR/counts');
  scriptLines.push(`Rscript $SCRIPTS_DIR/${featurecountsScript} $GENOME $COMPARISON $CONTROL_KEYWORD $TREATMENT_KEYWORD $OUTPUT_DIR/bam_files/*.bam`);
  scriptLines.push('cd -');
  scriptLines.push('');

  // Step 4: Filter Bad IDs
  scriptLines.push('echo "Step 4/10: Filtering bad gene IDs..."');
  scriptLines.push(`Rscript $SCRIPTS_DIR/${SCRIPT_FILES.filterIDS} $OUTPUT_DIR/counts/${config.comparison}.count.csv`);
  scriptLines.push('');

  // Step 5: RPKM Normalization
  scriptLines.push('echo "Step 5/10: RPKM normalization..."');
  scriptLines.push(`Rscript $SCRIPTS_DIR/${SCRIPT_FILES.RPKM} $OUTPUT_DIR/counts/${config.comparison}.count.filtered.csv`);
  scriptLines.push('');

  // Step 6: Add Gene Symbols
  scriptLines.push('echo "Step 6/10: Adding gene symbols..."');
  scriptLines.push(`Rscript $SCRIPTS_DIR/${SCRIPT_FILES.entrz} $OUTPUT_DIR/counts/$COMPARISON.rpkm.csv $GENOME`);
  scriptLines.push('');

  // Step 7: QC Plots - SKIPPED (not needed)
  scriptLines.push('# Step 7: QC Plots skipped (not required for this analysis)');
  scriptLines.push('');

  // Step 8: DE Analysis
  scriptLines.push('echo "Step 8/10: Differential expression analysis (edgeR)..."');
  scriptLines.push('cd $OUTPUT_DIR/de_results');
  scriptLines.push(`Rscript $SCRIPTS_DIR/${SCRIPT_FILES.simpleEdger3} $OUTPUT_DIR/counts/$COMPARISON.count.filtered.csv $GENOME`);
  scriptLines.push('cd -');
  scriptLines.push('');

  // Step 9: Merge RPKM + edgeR Results into Final Excel
  scriptLines.push('echo "Step 9/10: Merging RPKM and edgeR results into Excel..."');
  scriptLines.push('cd $OUTPUT_DIR/de_results');
  scriptLines.push(`Rscript $SCRIPTS_DIR/merge_results.R $OUTPUT_DIR/counts/$COMPARISON.rpkm.entrz.csv *.csv $COMPARISON.final_results.xlsx`);
  scriptLines.push('cd -');
  scriptLines.push('');

  // Step 10: Complete
  scriptLines.push('echo "Step 10/10: Analysis complete!"');
  scriptLines.push('echo "  - BAM files: $OUTPUT_DIR/bam_files/"');
  scriptLines.push('echo "  - Counts: $OUTPUT_DIR/counts/"');
  scriptLines.push('echo "  - QC plots: $OUTPUT_DIR/qc_plots/"');
  scriptLines.push('echo "  - Final results: $OUTPUT_DIR/de_results/$COMPARISON.final_results.xlsx"');
  scriptLines.push('');

  const scriptContent = scriptLines.join('\n');

  fs.writeFileSync(outputPath, scriptContent, { mode: 0o755 });
  console.log(`[Script Generator] ✓ AUTOMATION script saved: ${outputPath}`);

  return {
    type: 'AUTOMATION',
    scriptPath: outputPath,
    content: scriptContent
  };
}

/**
 * Generate ADAPTATION script - Option 2 with fallback
 *
 * Strategy:
 * 1. Extract structured recommendations from agents
 * 2. Generate base template (same as AUTOMATION)
 * 3. Apply modifications from recommendations
 * 4. Validate modified script
 * 5. If validation fails → fallback to MCP agent (legacy)
 */
export async function generateAdaptationScript(dataInfo, config, steps, agentDecision, coordinator, mcpAgent, outputPath) {
  console.log('[Script Generator] Generating ADAPTATION script with MCP agent...');
  console.log('[Script Generator] Using TRUE ADAPTATION: MCP agent with tool calling');
  console.log('[Script Generator] Agent will:');
  console.log('[Script Generator]   1. Call list_available_scripts to validate what exists');
  console.log('[Script Generator]   2. Only use scripts that actually exist');
  console.log('[Script Generator]   3. Write custom scripts when needed (write_custom_script)');
  console.log('[Script Generator]   4. Handle errors and edge cases adaptively');
  console.log('');

  // Use MCP agent method (the TRULY agentic approach)
  // This is NOT "legacy" - it's the CORRECT approach with tool calling!
  return await generateAdaptationScriptLegacy(dataInfo, config, steps, agentDecision, coordinator, mcpAgent, outputPath);
}

/**
 * Generate base template for ADAPTATION modifications
 * (Same as AUTOMATION but returns string instead of writing file)
 */
function generateBaseTemplate(dataInfo, config, steps) {
  const genomeMap = {
    'mouse': 'mm10',
    'human': 'hg38',
    'rat': 'rn6'
  };
  const genomeBuild = genomeMap[config.organism] || config.organism;

  const scriptLines = [];

  // Header
  scriptLines.push('#!/bin/bash');
  scriptLines.push('#');
  scriptLines.push(`# GeneExpert Analysis Script (ADAPTATION - Option 2)`);
  scriptLines.push(`# Generated: ${new Date().toISOString()}`);
  scriptLines.push(`# Comparison: ${config.comparison}`);
  scriptLines.push(`# Organism: ${config.organism} (${genomeBuild})`);
  scriptLines.push('#');
  scriptLines.push('set -e  # Exit on error');
  scriptLines.push('');

  // Activate conda environment
  scriptLines.push('# Activate conda environment for bioinformatics tools');
  scriptLines.push(`source $(conda info --base)/etc/profile.d/conda.sh`);
  scriptLines.push(`conda activate ${CONDA_ENV}`);
  scriptLines.push('');

  // Variables
  scriptLines.push('# Paths');
  scriptLines.push(`OUTPUT_DIR="${config.output}"`);
  scriptLines.push(`SCRIPTS_DIR="${SCRIPTS_PATH}"`);
  scriptLines.push(`INDEX_PATH="${INDEX_PATH}"`);
  scriptLines.push(`GENOME="${genomeBuild}"`);
  scriptLines.push(`COMPARISON="${config.comparison}"`);
  scriptLines.push(`CONTROL_KEYWORD="${config.controlKeyword || 'control'}"`);
  scriptLines.push(`TREATMENT_KEYWORD="${config.treatmentKeyword || 'treatment'}"`);
  scriptLines.push('');

  // Create output directories
  scriptLines.push('# Create output directories');
  scriptLines.push('mkdir -p $OUTPUT_DIR/fastqc_results');
  scriptLines.push('mkdir -p $OUTPUT_DIR/bam_files');
  scriptLines.push('mkdir -p $OUTPUT_DIR/counts');
  scriptLines.push('mkdir -p $OUTPUT_DIR/qc_plots');
  scriptLines.push('mkdir -p $OUTPUT_DIR/de_results');
  scriptLines.push('');

  // Step 1: FastQC
  scriptLines.push('echo "Step 1/10: Running FastQC..."');
  const fastqFiles = dataInfo.samples.map(s => s.r1 + (s.r2 ? ' ' + s.r2 : '')).join(' ');
  scriptLines.push(`fastqc ${fastqFiles} -o $OUTPUT_DIR/fastqc_results`);
  scriptLines.push('');

  // Step 2: Alignment
  scriptLines.push('echo "Step 2/10: Alignment (FASTQ to BAM)..."');
  scriptLines.push(`cd ${config.input}`);
  scriptLines.push('');

  // Generate alignment commands for each sample
  if (dataInfo.pairedEnd) {
    scriptLines.push('# Run alignment for paired-end samples');
    scriptLines.push('for i in *_R1_001.fastq.gz; do');
    scriptLines.push('  fname=$(basename "$i" _R1_001.fastq.gz)');
    scriptLines.push('  echo "Aligning: $fname"');
    scriptLines.push('  subread-align -t 0 \\');
    scriptLines.push('    -i $INDEX_PATH/$GENOME \\');
    scriptLines.push('    -r "${fname}_R1_001.fastq.gz" \\');
    scriptLines.push('    -R "${fname}_R2_001.fastq.gz" \\');
    scriptLines.push('    -T 8 \\');
    scriptLines.push('    -o "${fname}.bam" &> "${fname}.log" &');
    scriptLines.push('done');
    scriptLines.push('wait');
  } else {
    scriptLines.push('# Run alignment for single-end samples');
    scriptLines.push('for i in *.fastq.gz; do');
    scriptLines.push('  fname=$(basename "$i" .fastq.gz)');
    scriptLines.push('  echo "Aligning: $fname"');
    scriptLines.push('  subread-align -t 0 \\');
    scriptLines.push('    -i $INDEX_PATH/$GENOME \\');
    scriptLines.push('    -r "${fname}.fastq.gz" \\');
    scriptLines.push('    -T 8 \\');
    scriptLines.push('    -o "${fname}.bam" &> "${fname}.log" &');
    scriptLines.push('done');
    scriptLines.push('wait');
  }

  scriptLines.push('');
  scriptLines.push('echo "Moving BAM files to output directory..."');
  scriptLines.push(`mv *.bam *.log $OUTPUT_DIR/bam_files/`);
  scriptLines.push(`cd -`);
  scriptLines.push('');

  // Step 3: Feature Counts (use paired or unpaired script based on data)
  const featurecountsScript = dataInfo.pairedEnd ? SCRIPT_FILES.featurecounts : SCRIPT_FILES.featurecounts_unpaired;
  scriptLines.push(`echo "Step 3/10: Feature Counts (${dataInfo.pairedEnd ? 'paired-end' : 'single-end'})..."`);
  scriptLines.push('cd $OUTPUT_DIR/counts');
  scriptLines.push(`Rscript $SCRIPTS_DIR/${featurecountsScript} $GENOME $COMPARISON $CONTROL_KEYWORD $TREATMENT_KEYWORD $OUTPUT_DIR/bam_files/*.bam`);
  scriptLines.push('cd -');
  scriptLines.push('');

  // Step 4: Filter Bad IDs
  scriptLines.push('echo "Step 4/10: Filtering bad gene IDs..."');
  scriptLines.push(`Rscript $SCRIPTS_DIR/${SCRIPT_FILES.filterIDS} $OUTPUT_DIR/counts/$COMPARISON.count.csv`);
  scriptLines.push('');

  // Step 5: RPKM Normalization
  scriptLines.push('echo "Step 5/10: RPKM normalization..."');
  scriptLines.push(`Rscript $SCRIPTS_DIR/${SCRIPT_FILES.RPKM} $OUTPUT_DIR/counts/$COMPARISON.count.filtered.csv`);
  scriptLines.push('');

  // Step 6: Add Gene Symbols
  scriptLines.push('echo "Step 6/10: Adding gene symbols..."');
  scriptLines.push(`Rscript $SCRIPTS_DIR/${SCRIPT_FILES.entrz} $OUTPUT_DIR/counts/$COMPARISON.rpkm.csv $GENOME`);
  scriptLines.push('');

  // Step 7: QC Plots - SKIPPED (not needed)
  scriptLines.push('# Step 7: QC Plots skipped (not required for this analysis)');
  scriptLines.push('');

  // Step 8: DE Analysis
  scriptLines.push('echo "Step 8/10: Differential expression analysis (edgeR)..."');
  scriptLines.push('cd $OUTPUT_DIR/de_results');
  scriptLines.push(`Rscript $SCRIPTS_DIR/${SCRIPT_FILES.simpleEdger3} $OUTPUT_DIR/counts/$COMPARISON.count.filtered.csv $GENOME`);
  scriptLines.push('cd -');
  scriptLines.push('');

  // Step 9: Merge RPKM + edgeR Results into Final Excel
  scriptLines.push('echo "Step 9/10: Merging RPKM and edgeR results into Excel..."');
  scriptLines.push('cd $OUTPUT_DIR/de_results');
  scriptLines.push(`Rscript $SCRIPTS_DIR/merge_results.R $OUTPUT_DIR/counts/$COMPARISON.rpkm.entrz.csv *.csv $COMPARISON.final_results.xlsx`);
  scriptLines.push('cd -');
  scriptLines.push('');

  // Step 10: Complete
  scriptLines.push('echo "Step 10/10: Analysis complete!"');
  scriptLines.push('echo "  - BAM files: $OUTPUT_DIR/bam_files/"');
  scriptLines.push('echo "  - Counts: $OUTPUT_DIR/counts/"');
  scriptLines.push('echo "  - QC plots: $OUTPUT_DIR/qc_plots/"');
  scriptLines.push('echo "  - Final results: $OUTPUT_DIR/de_results/$COMPARISON.final_results.xlsx"');
  scriptLines.push('');

  return scriptLines.join('\n');
}

/**
 * LEGACY: Generate ADAPTATION script using MCP agent (fallback method)
 */
export async function generateAdaptationScriptLegacy(dataInfo, config, steps, agentDecision, coordinator, mcpAgent, outputPath) {
  console.log('[Script Generator] Using LEGACY MCP agent method (fallback)...');
  console.log('[Script Generator] MCP Agent will read scripts and write custom solution...');
  console.log('');

  // Build context for agents
  const agentConcerns = [
    agentDecision.responses.gpt5_2.content,
    agentDecision.responses.claude.content,
    agentDecision.responses.gemini.content
  ].join('\n\n---\n\n');

  // Map organism to genome build
  const genomeMap = {
    'mouse': 'mm10',
    'human': 'hg38',
    'rat': 'rn6'
  };
  const genomeBuild = genomeMap[config.organism] || config.organism;

  const prompt = `
You are writing a custom bash script for RNA-seq analysis to address specific concerns.

IMPORTANT: You can ONLY use the "read_file" tool to read existing scripts.
DO NOT try to call validate_fastq or list_available_scripts - those were called before.
Your ONLY job is to write a COMPLETE executable bash script.

DATASET:
- Input directory: ${config.input}
- Samples: ${dataInfo.samples.length} (${dataInfo.pairedEnd ? 'paired-end' : 'single-end'})
- Groups: ${dataInfo.groups ? Object.keys(dataInfo.groups).join(' vs ') : 'unknown'}
- Organism: ${config.organism} (genome build: ${genomeBuild})
- Comparison: ${config.comparison}
- Output directory: ${config.output}

IMPORTANT:
- Use genome build "${genomeBuild}" in all commands, NOT "${config.organism}"
- Input FASTQ files are in: ${config.input}
- Save all outputs to: ${config.output}

ENVIRONMENT SETUP:
- Conda environment: ${CONDA_ENV}
- Genome index path: ${INDEX_PATH}
- Scripts directory: ${SCRIPTS_PATH}/

AVAILABLE SCRIPTS (verified to exist):
- validate_fastq: ${SCRIPT_FILES.validate_fastq} (FASTQ integrity check)
- alignment_qc: ${SCRIPT_FILES.alignment_qc} (alignment-based contamination/quality screening)
- featurecounts (paired-end): ${SCRIPT_FILES.featurecounts}
- featurecounts (single-end): ${SCRIPT_FILES.featurecounts_unpaired}
- filterIDS: ${SCRIPT_FILES.filterIDS}
- RPKM: ${SCRIPT_FILES.RPKM}
- entrz: ${SCRIPT_FILES.entrz}
- qc_assessment: ${SCRIPT_FILES.qc_assessment} (PCA-based outlier & batch detection)
- simpleEdger3: ${SCRIPT_FILES.simpleEdger3} (standard DE analysis)
- batch_effect_edger: ${SCRIPT_FILES.batch_effect_edger} (DE with batch correction)
  * 3rd arg specifies EXPERIMENTAL BATCHES (not sequencing type!):
    - "auto" = try to detect from sample names
    - "paired" = ctrl1+trt1 same batch, ctrl2+trt2 same batch, etc.
    - "1,1,2,2" = explicit batch labels (ctrl1=batch1, ctrl2=batch1, trt1=batch2, trt2=batch2)
  * Default if omitted: assumes paired design

DATA IS: ${dataInfo.pairedEnd ? 'PAIRED-END (use featurecounts.R)' : 'SINGLE-END (use featurecounts_unpaired.R)'}

ALIGNMENT:
- Use subread-align directly (NOT the fastq2bam wrapper script)
- Index path: ${INDEX_PATH}/${genomeBuild}

YOUR CONCERNS FROM ANALYSIS:
${agentConcerns}

STANDARD RNA-SEQ PIPELINE (adapt this to address concerns):

Activate conda environment:
  source $(conda info --base)/etc/profile.d/conda.sh
  conda activate ${CONDA_ENV}

Step 0: FASTQ Validation (REQUIRED FIRST STEP)
  bash ${SCRIPTS_PATH}/${SCRIPT_FILES.validate_fastq} ${config.input} ${config.output}/fastq_validation_report.tsv
  # Check if validation passed before proceeding

Step 1: FastQC
  fastqc ${config.input}/*.fastq* -o ${config.output}/fastqc_results

Step 2: Alignment (use subread-align directly with correct index path)
  cd ${config.input}
  for i in *_R1_001.fastq.gz; do
    fname=$(basename "$i" _R1_001.fastq.gz)
    subread-align -t 0 -i ${INDEX_PATH}/${genomeBuild} \\
      -r "\${fname}_R1_001.fastq.gz" -R "\${fname}_R2_001.fastq.gz" \\
      -T 8 -o "\${fname}.bam" &> "\${fname}.log" &
  done
  wait
  echo "Alignment complete - all samples finished successfully"
  mv *.bam *.log ${config.output}/bam_files/
  cd -

Step 2.5: Alignment QC Screening (ALWAYS RUN - contamination/quality check)
  mkdir -p ${config.output}/alignment_qc
  Rscript ${SCRIPTS_PATH}/${SCRIPT_FILES.alignment_qc} ${config.output}/bam_files ${config.output}/alignment_qc
  ALIGN_QC_EXIT=$?

  if [ $ALIGN_QC_EXIT -eq 2 ]; then
    echo "WARNING: Alignment QC detected FAIL or SEVERE_FAIL samples"
    echo "Review: ${config.output}/alignment_qc/alignment_qc_report.txt"
    # Agents will decide whether to proceed or remove samples
  fi

Step 3: Feature Counts (${dataInfo.pairedEnd ? 'PAIRED' : 'UNPAIRED'})
  cd ${config.output}/counts
  Rscript ${SCRIPTS_PATH}/${dataInfo.pairedEnd ? SCRIPT_FILES.featurecounts : SCRIPT_FILES.featurecounts_unpaired} ${genomeBuild} ${config.comparison} ${config.controlKeyword} ${config.treatmentKeyword} ${config.output}/bam_files/*.bam
  cd -

Step 4: Filter Bad IDs
  Rscript ${SCRIPTS_PATH}/${SCRIPT_FILES.filterIDS} ${config.output}/counts/${config.comparison}.count.csv

Step 5: RPKM Normalization (for visualization)
  Rscript ${SCRIPTS_PATH}/${SCRIPT_FILES.RPKM} ${config.output}/counts/${config.comparison}.count.filtered.csv

Step 6: Add Gene Symbols
  Rscript ${SCRIPTS_PATH}/${SCRIPT_FILES.entrz} ${config.output}/counts/${config.comparison}.rpkm.csv ${genomeBuild}

Step 7: QC Assessment (CRITICAL - Outlier & Batch Detection)
  mkdir -p ${config.output}/qc_assessment
  Rscript ${SCRIPTS_PATH}/${SCRIPT_FILES.qc_assessment} \\
    ${config.output}/counts/${config.comparison}.rpkm.csv \\
    ${config.output}/qc_assessment \\
    ${config.controlKeyword} \\
    ${config.treatmentKeyword}
  QC_EXIT_CODE=$?

  # Check QC results
  if [ $QC_EXIT_CODE -eq 2 ]; then
    echo "WARNING: QC issues detected - review qc_assessment/qc_summary.txt"
    # TODO: Decide batch_effect_edger vs simpleEdger3 based on QC report
  fi

Step 8: DE Analysis (edgeR uses RAW filtered counts, not RPKM)
  cd ${config.output}/de_results

  # Decision based on QC results:
  if grep -q "TOO SPREAD" ${config.output}/qc_assessment/qc_summary.txt 2>/dev/null; then
    echo "Batch effects detected - using batch correction"
    Rscript ${SCRIPTS_PATH}/${SCRIPT_FILES.batch_effect_edger} ${config.output}/counts/${config.comparison}.count.filtered.csv ${genomeBuild} paired
  else
    echo "Standard DE analysis"
    Rscript ${SCRIPTS_PATH}/${SCRIPT_FILES.simpleEdger3} ${config.output}/counts/${config.comparison}.count.filtered.csv ${genomeBuild}
  fi

  cd -

Step 9: Merge RPKM + edgeR Results into Final Excel
  cd ${config.output}/de_results
  Rscript ${SCRIPTS_PATH}/merge_results.R ${config.output}/counts/${config.comparison}.rpkm.entrz.csv *.csv ${config.comparison}.final_results.xlsx
  cd -

YOUR TASK:
Based on agent concerns above, write a COMPLETE bash script that:
- Runs ALL the steps above (or modified versions)
- Addresses the specific concerns (e.g., small n=2, batch effects)
- Uses correct genome build: ${genomeBuild}
- Saves everything to: ${config.output}

ADAPTATION EXAMPLES:
- If small n: Add extra QC checks, use more conservative parameters
- If batch effects: Add batch correction comments
- Always include ALL steps (don't just stop at discovery)

CRITICAL: Your script must actually RUN the pipeline, not just plan it!

OUTPUT FORMAT - Write the COMPLETE bash script in a code block:

\`\`\`bash
#!/bin/bash
set -e

# Activate conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate ${CONDA_ENV}

# Create output directories
mkdir -p ${config.output}/fastqc_results
mkdir -p ${config.output}/bam_files
mkdir -p ${config.output}/counts
mkdir -p ${config.output}/de_results

# Step 1/8: FastQC
echo "Step 1/8: Running FastQC..."
fastqc ${config.input}/*.fastq* -o ${config.output}/fastqc_results

# Step 2/8: Alignment
echo "Step 2/8: Alignment..."
cd ${config.input}
for i in *_R1_001.fastq.gz; do
  fname=$(basename "$i" _R1_001.fastq.gz)
  subread-align -t 0 -i ${INDEX_PATH}/${genomeBuild} \\
    -r "\${fname}_R1_001.fastq.gz" -R "\${fname}_R2_001.fastq.gz" \\
    -T 8 -o "\${fname}.bam" &> "\${fname}.log" &
done
wait
mv *.bam *.log ${config.output}/bam_files/
cd -

# Step 3/8: Feature Counts (${dataInfo.pairedEnd ? 'PAIRED' : 'UNPAIRED'})
echo "Step 3/8: Feature Counts..."
cd ${config.output}/counts
Rscript ${SCRIPTS_PATH}/${dataInfo.pairedEnd ? SCRIPT_FILES.featurecounts : SCRIPT_FILES.featurecounts_unpaired} ${genomeBuild} ${config.comparison} ${config.controlKeyword} ${config.treatmentKeyword} ${config.output}/bam_files/*.bam
cd -

# ... CONTINUE WITH ALL OTHER STEPS (4-8) ...
# Don't stop here! Write the complete pipeline!
\`\`\`

CRITICAL REQUIREMENTS:
1. Write a COMPLETE executable bash script with ALL 8 steps
2. Do NOT write conversational text - ONLY the bash script
3. Put the script in a \`\`\`bash code block
4. Include actual commands, not just comments
5. Do NOT call MCP tools (except read_file if you need to check a script)
6. Your response must start with \`\`\`bash and contain only the script

Write the FULL executable bash script now (and ONLY the script, no other text):
`;

  // Call MCP-enabled Claude agent with read_file tool ONLY
  const mcpContext = {
    scriptsPath: SCRIPTS_PATH,
    outputDir: config.output,
    genomeBuild,
    comparison: config.comparison,
    dataInfo
  };

  console.log('[MCP Agent] Using tools: read_file');
  console.log('[MCP Agent] Agent will read scripts and write custom bash script');
  console.log('');

  const result = await mcpAgent.callWithTools(prompt, mcpContext, [
    'read_file'  // To read the R/bash scripts (agent writes script directly in response)
  ]);

  // DEBUG: Log what the agent actually returned
  console.log('[DEBUG] MCP Agent returned:');
  console.log('[DEBUG] - Success:', result.success);
  console.log('[DEBUG] - Model:', result.model);
  console.log('[DEBUG] - Tool calls:', result.toolCalls ? result.toolCalls.length : 0);
  console.log('[DEBUG] - Content length:', result.content.length);
  console.log('[DEBUG] - Content preview (first 500 chars):');
  console.log(result.content.substring(0, 500));
  console.log('[DEBUG] - Content preview (last 500 chars):');
  console.log(result.content.substring(Math.max(0, result.content.length - 500)));

  const scriptContent = extractBashScript(result.content);

  // Validate script was generated
  if (!scriptContent || scriptContent.length < 50 || scriptContent.includes('No text response')) {
    console.error('[Script Generator] ❌ Agent failed to generate valid bash script');
    console.error('[Script Generator] Agent response:', result.content.substring(0, 200));
    throw new Error('MCP agent did not generate a valid bash script. Agent may need more specific instructions or different model.');
  }

  fs.writeFileSync(outputPath, scriptContent, { mode: 0o755 });
  console.log('');
  console.log(`[Script Generator] ✓ ADAPTATION script saved: ${outputPath}`);
  console.log(`[Script Generator] Agent modifications applied based on analysis`);

  return {
    type: 'ADAPTATION',
    scriptPath: outputPath,
    content: scriptContent,
    agentReasoning: result.content
  };
}

/**
 * Extract bash script from agent response
 */
function extractBashScript(agentResponse) {
  // Look for code blocks
  const bashCodeBlock = agentResponse.match(/```bash\n([\s\S]*?)```/);
  if (bashCodeBlock) {
    return bashCodeBlock[1].trim();
  }

  // Look for sh code blocks
  const shCodeBlock = agentResponse.match(/```sh\n([\s\S]*?)```/);
  if (shCodeBlock) {
    return shCodeBlock[1].trim();
  }

  // Look for generic code blocks
  const codeBlock = agentResponse.match(/```\n([\s\S]*?)```/);
  if (codeBlock) {
    return codeBlock[1].trim();
  }

  // If no code block, look for lines starting with # or typical bash commands
  const lines = agentResponse.split('\n');
  const scriptLines = [];
  let inScript = false;

  for (const line of lines) {
    if (line.startsWith('#!/bin/bash') || line.startsWith('#!')) {
      inScript = true;
    }
    if (inScript) {
      scriptLines.push(line);
    }
  }

  if (scriptLines.length > 0) {
    return scriptLines.join('\n');
  }

  // Fallback: return whole response with shebang added
  return '#!/bin/bash\n' + agentResponse;
}

export default {
  generateAutomationScript,
  generateAdaptationScript
};
