#!/usr/bin/env node
/**
 * Preprocess Dataset - Run Stage 1 (FastQC) and Stage 2 (Alignment) ONCE
 *
 * This script runs the deterministic preprocessing steps (FastQC validation and alignment)
 * once per dataset. The outputs can then be reused across all 5 ICML experimental systems
 * (no-agent, single-GPT, single-Claude, parallel, sequential) to save time.
 *
 * Usage:
 *   node bin/preprocess_dataset.js <input_dir> <output_dir> --organism <organism> [options]
 *
 * Example:
 *   node bin/preprocess_dataset.js \
 *     data/1_GSE52778_pe_clean \
 *     experiments/preprocessing/1_GSE52778_pe_clean \
 *     --organism human
 *
 * Outputs:
 *   - output_dir/stage1_fastq_validation/   (FastQC reports, validation CSV)
 *   - output_dir/stage2_alignment/bam_files/  (BAM files)
 *   - output_dir/stage2_alignment/alignment_qc/  (Alignment QC summary)
 */

import path from 'path';
import fs from 'fs';
import { execSync } from 'child_process';
import dotenv from 'dotenv';

dotenv.config();

// Import stage modules
import { generateStage1Script, parseStage1Output } from '../src/stages/stage1_fastq_validation.js';
import { generateStage2Script, parseStage2Output } from '../src/stages/stage2_alignment.js';

// Parse CLI arguments
const args = process.argv.slice(2);

if (args.length < 2) {
  console.log('Usage: node bin/preprocess_dataset.js <input_dir> <output_dir> --organism <organism> [options]');
  console.log('');
  console.log('Required:');
  console.log('  input_dir       Path to FASTQ files');
  console.log('  output_dir      Path to save preprocessing outputs');
  console.log('  --organism      Organism (human, mouse, rat)');
  console.log('');
  console.log('Example:');
  console.log('  node bin/preprocess_dataset.js \\');
  console.log('    data/1_GSE52778_pe_clean \\');
  console.log('    experiments/preprocessing/1_GSE52778_pe_clean \\');
  console.log('    --organism human');
  console.log('');
  process.exit(1);
}

const inputDir = path.resolve(args[0]);
const outputDir = path.resolve(args[1]);
const organism = args.find(a => a.startsWith('--organism='))?.split('=')[1] ||
                 (args.indexOf('--organism') !== -1 ? args[args.indexOf('--organism') + 1] : null);

if (!organism) {
  console.error('ERROR: --organism is required');
  process.exit(1);
}

if (!fs.existsSync(inputDir)) {
  console.error(`ERROR: Input directory not found: ${inputDir}`);
  process.exit(1);
}

// Detect dataset info
function detectDataInfo(inputDir) {
  const files = fs.readdirSync(inputDir).filter(f => f.endsWith('.fastq.gz'));
  const samples = [];
  const pairedEnd = files.some(f => f.includes('_R2_'));

  const sampleMap = {};
  for (const file of files) {
    const sampleName = file.replace(/_R[12]_001\.fastq\.gz$/, '').replace(/\.fastq\.gz$/, '');
    if (!sampleMap[sampleName]) {
      sampleMap[sampleName] = { name: sampleName };
    }
  }

  for (const [name, sample] of Object.entries(sampleMap)) {
    samples.push(sample);
  }

  return {
    samples,
    pairedEnd,
    organism,
    comparison: path.basename(outputDir)
  };
}

// Execute script
function executeScript(scriptPath, stageName) {
  console.log(`[Preprocessing] Executing ${stageName}...`);
  console.log(`[Preprocessing] This may take 20-30 minutes...`);
  console.log('');

  try {
    execSync(`bash ${scriptPath}`, {
      stdio: 'inherit',
      cwd: path.dirname(scriptPath),
      timeout: 3600000  // 1 hour
    });
    console.log('');
    console.log(`[Preprocessing] ${stageName} completed successfully`);
    return { success: true };
  } catch (error) {
    console.error(`[Preprocessing] ${stageName} failed: ${error.message}`);
    return { success: false, error: error.message };
  }
}

async function main() {
  console.log('');
  console.log('='.repeat(70));
  console.log('  PREPROCESSING DATASET (Stage 1 + Stage 2)');
  console.log('='.repeat(70));
  console.log('');
  console.log(`Input:    ${inputDir}`);
  console.log(`Output:   ${outputDir}`);
  console.log(`Organism: ${organism}`);
  console.log('');

  // Detect dataset
  const dataInfo = detectDataInfo(inputDir);
  console.log(`Detected: ${dataInfo.samples.length} samples (${dataInfo.pairedEnd ? 'paired-end' : 'single-end'})`);
  console.log('');

  // Create output directory
  if (!fs.existsSync(outputDir)) {
    fs.mkdirSync(outputDir, { recursive: true });
  }

  const config = {
    input: inputDir,
    output: outputDir,
    organism,
    comparison: dataInfo.comparison,
    controlKeyword: 'control',  // dummy for preprocessing
    treatmentKeyword: 'treatment'  // dummy for preprocessing
  };

  // ===================================================================
  // STAGE 1: FastQC Validation
  // ===================================================================

  console.log('='.repeat(70));
  console.log('  STAGE 1: FastQC Validation');
  console.log('='.repeat(70));
  console.log('');

  const stage1ScriptPath = path.join(outputDir, 'stage_1_preprocessing.sh');
  console.log('[Stage 1] Generating script...');
  generateStage1Script(dataInfo, config, stage1ScriptPath);
  console.log(`[Stage 1] Script: ${stage1ScriptPath}`);
  console.log('');

  const stage1Result = executeScript(stage1ScriptPath, 'Stage 1 (FastQC)');
  if (!stage1Result.success) {
    console.error('ERROR: Stage 1 failed. Cannot proceed to Stage 2.');
    process.exit(1);
  }

  console.log('[Stage 1] Parsing output...');
  const stage1Output = parseStage1Output(outputDir);
  console.log(`[Stage 1] Status: ${stage1Output.overall_status}`);
  console.log(`[Stage 1] Samples: ${stage1Output.samples_validated}`);
  console.log('');

  // ===================================================================
  // STAGE 2: Alignment
  // ===================================================================

  console.log('='.repeat(70));
  console.log('  STAGE 2: Alignment');
  console.log('='.repeat(70));
  console.log('');

  const stage2ScriptPath = path.join(outputDir, 'stage_2_preprocessing.sh');
  console.log('[Stage 2] Generating script...');
  generateStage2Script(dataInfo, config, stage2ScriptPath);
  console.log(`[Stage 2] Script: ${stage2ScriptPath}`);
  console.log('');

  const stage2Result = executeScript(stage2ScriptPath, 'Stage 2 (Alignment)');
  if (!stage2Result.success) {
    console.error('ERROR: Stage 2 failed.');
    process.exit(1);
  }

  console.log('[Stage 2] Parsing output...');
  const stage2Output = parseStage2Output(outputDir);
  console.log(`[Stage 2] Status: ${stage2Output.overall_status}`);
  const bamCount = stage2Output.bam_files ? stage2Output.bam_files.length :
                   (stage2Output.per_sample_alignment ? Object.keys(stage2Output.per_sample_alignment).length : 0);
  console.log(`[Stage 2] BAM files: ${bamCount}`);
  console.log('');

  // ===================================================================
  // SUMMARY
  // ===================================================================

  console.log('='.repeat(70));
  console.log('  PREPROCESSING COMPLETE!');
  console.log('='.repeat(70));
  console.log('');
  console.log('Output directories:');
  console.log(`  Stage 1 (FastQC): ${outputDir}/stage1_fastq_validation/`);
  console.log(`  Stage 2 (BAM):    ${outputDir}/stage2_alignment/bam_files/`);
  console.log('');
  console.log('To use these in ICML experiments, run:');
  console.log('');
  console.log(`  node bin/geneexpert.js analyze ${inputDir} \\`);
  console.log(`    --staged \\`);
  console.log(`    --use-existing-fastqc ${outputDir}/stage1_fastq_validation \\`);
  console.log(`    --use-existing-bam ${outputDir}/stage2_alignment/bam_files \\`);
  console.log(`    --organism ${organism} \\`);
  console.log(`    --output results/your_experiment`);
  console.log('');
}

main().catch(err => {
  console.error('FATAL ERROR:', err);
  process.exit(1);
});
