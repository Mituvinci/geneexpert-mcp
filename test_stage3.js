#!/usr/bin/env node
/**
 * Test Script for Stage 3: Quantification + QC Assessment
 *
 * This script tests the staged architecture by running Stage 3 in isolation:
 * 1. Generate the Stage 3 script (quantification + PCA/QC)
 * 2. Execute it
 * 3. Parse the output
 * 4. Send to agents for review
 * 5. Log the decision
 *
 * Usage:
 *   node test_stage3.js <input_dir> <output_dir> [--dry-run] [--skip-execution]
 *
 * Example:
 *   node test_stage3.js data/DA0036 results/stage3_test
 *   node test_stage3.js data/DA0036 results/stage3_test --dry-run
 *
 * Note: Stage 3 requires Stage 2 to have run (BAM files needed)
 */

import path from 'path';
import fs from 'fs';
import { execSync } from 'child_process';
import dotenv from 'dotenv';

// Load environment variables
dotenv.config();

// Import stage 3 module
import {
  generateStage3Script,
  parseStage3Output,
  formatStage3ForAgents
} from './src/stages/stage3_quantification_qc.js';

// Import coordinator and logger
import { Coordinator } from './src/coordinator/orchestrator.js';
import { Logger } from './src/utils/logger.js';
import { handleStage3UserDecision } from './src/utils/user_input.js';

// Parse command line arguments
const args = process.argv.slice(2);
const inputDir = args[0];
const outputDir = args[1];
const dryRun = args.includes('--dry-run');
const skipExecution = args.includes('--skip-execution');
const organism = args.find(a => a.startsWith('--organism='))?.split('=')[1] || 'mouse';

if (!inputDir || !outputDir) {
  console.log('Usage: node test_stage3.js <input_dir> <output_dir> [options]');
  console.log('');
  console.log('Options:');
  console.log('  --dry-run         Generate script but don\'t execute or call agents');
  console.log('  --skip-execution  Skip script execution, only call agents (if output exists)');
  console.log('  --organism=mouse  Specify organism (mouse, human, rat)');
  console.log('');
  console.log('Examples:');
  console.log('  node test_stage3.js data/DA0036 results/stage3_test');
  console.log('  node test_stage3.js data/DA0036 results/stage3_test --dry-run');
  console.log('  node test_stage3.js data/DA0036 results/stage3_test --organism=human');
  process.exit(1);
}

// Detect data info from input directory
function detectDataInfo(inputDir) {
  console.log(`[Test] Detecting data from: ${inputDir}`);

  const files = fs.readdirSync(inputDir).filter(f => f.endsWith('.fastq.gz'));
  const samples = [];
  const pairedEnd = files.some(f => f.includes('_R2_'));

  // Group files by sample
  const sampleMap = {};
  for (const file of files) {
    const sampleName = file.replace(/_R[12]_001\.fastq\.gz$/, '').replace(/\.fastq\.gz$/, '');

    if (!sampleMap[sampleName]) {
      sampleMap[sampleName] = { name: sampleName };
    }

    if (file.includes('_R1_')) {
      sampleMap[sampleName].r1 = path.join(inputDir, file);
    } else if (file.includes('_R2_')) {
      sampleMap[sampleName].r2 = path.join(inputDir, file);
    } else {
      sampleMap[sampleName].r1 = path.join(inputDir, file);
    }
  }

  for (const [name, sample] of Object.entries(sampleMap)) {
    samples.push(sample);
  }

  // Detect groups
  const groups = {};
  for (const sample of samples) {
    const name = sample.name.toLowerCase();
    let group = 'unknown';

    if (name.includes('cont') || name.includes('ctrl') || name.includes('wt') || name.includes('control')) {
      group = 'control';
    } else if (name.includes('treat') || name.includes('ips') || name.includes('ko') || name.includes('mut')) {
      group = 'treatment';
    }

    if (!groups[group]) {
      groups[group] = [];
    }
    groups[group].push(sample.name);
  }

  const dataInfo = {
    type: 'fastq',
    samples,
    pairedEnd,
    groups,
    organism,
    comparison: path.basename(outputDir)
  };

  console.log(`[Test] Detected ${samples.length} samples (${pairedEnd ? 'paired-end' : 'single-end'})`);
  console.log(`[Test] Groups: ${Object.keys(groups).map(g => `${g}(n=${groups[g].length})`).join(', ')}`);

  return dataInfo;
}

async function main() {
  console.log('');
  console.log('='.repeat(60));
  console.log('  Stage 3 Test: Quantification + QC Assessment');
  console.log('='.repeat(60));
  console.log('');

  // Create output directory
  if (!fs.existsSync(outputDir)) {
    fs.mkdirSync(outputDir, { recursive: true });
  }

  // Check if Stage 2 BAM files exist
  const stage2BamDir = path.join(outputDir, 'stage2_alignment', 'bam_files');
  if (!fs.existsSync(stage2BamDir)) {
    console.error('ERROR: Stage 2 BAM files not found!');
    console.error(`Expected: ${stage2BamDir}`);
    console.error('');
    console.error('You must run Stage 2 (alignment) before Stage 3.');
    console.error('Run: node test_stage2.js data/DA0036 results/stage3_test');
    process.exit(1);
  }

  // Detect data info
  const dataInfo = detectDataInfo(inputDir);

  // Config
  const config = {
    input: path.resolve(inputDir),
    output: path.resolve(outputDir),
    organism,
    comparison: path.basename(outputDir),
    controlKeyword: 'cont',
    treatmentKeyword: 'ips'
  };

  // Initialize logger (use absolute path)
  const absOutputDir = path.resolve(outputDir);
  const logger = new Logger(absOutputDir, 'stage3_test', config);

  // Initialize coordinator
  const coordinator = new Coordinator({ verbose: true });
  coordinator.setDataset(config.comparison);

  // Step 1: Generate Stage 3 script
  console.log('[Test] Step 1: Generating Stage 3 script...');
  const scriptPath = path.resolve(outputDir, `stage_3_${config.comparison}.sh`);
  const scriptResult = generateStage3Script(dataInfo, config, scriptPath);
  console.log(`[Test] Script generated: ${scriptPath}`);

  if (dryRun) {
    console.log('');
    console.log('[Test] DRY RUN - Script generated but not executed');
    console.log('[Test] To run the script manually:');
    console.log(`       bash ${scriptPath}`);
    console.log('');
    console.log('Script content (first 70 lines):');
    console.log('-'.repeat(60));
    const scriptContent = fs.readFileSync(scriptPath, 'utf-8');
    console.log(scriptContent.split('\n').slice(0, 70).join('\n'));
    console.log('...');
    console.log('-'.repeat(60));
    process.exit(0);
  }

  // Step 2: Execute Stage 3 script
  if (!skipExecution) {
    console.log('');
    console.log('[Test] Step 2: Executing Stage 3 script...');
    console.log('[Test] WARNING: This will run quantification which can take 10-20+ minutes!');
    console.log('[Test] Press Ctrl+C to cancel...');
    console.log('');

    try {
      // Run with output streaming
      const result = execSync(`bash ${scriptPath}`, {
        stdio: 'inherit',
        cwd: outputDir,
        timeout: 3600000  // 1 hour timeout
      });
      console.log('');
      console.log('[Test] Stage 3 script completed');
    } catch (error) {
      // Check if it's just a non-zero exit code (issues detected in QC)
      if (error.status === 2) {
        console.log('');
        console.log(`[Test] Stage 3 script completed with QC issues (exit code: 2)`);
      } else {
        console.error('[Test] Error executing Stage 3 script:', error.message);
        process.exit(1);
      }
    }
  } else {
    console.log('[Test] Step 2: Skipping execution (--skip-execution flag)');
  }

  // Step 3: Parse Stage 3 output
  console.log('');
  console.log('[Test] Step 3: Parsing Stage 3 output...');
  const stage3Output = parseStage3Output(absOutputDir);
  console.log(`[Test] Parsed output: ${stage3Output.overall_status}`);

  // Show parsed output summary (BEFORE agent review - only basic metadata)
  console.log('');
  console.log('--- Stage 3 Output Summary ---');
  console.log(`Samples Analyzed: ${stage3Output.samples_analyzed}`);
  console.log(`Genes Detected: ${stage3Output.genes_detected?.toLocaleString()}`);
  console.log(`PC1 Variance: ${stage3Output.pc1_variance ? (stage3Output.pc1_variance * 100).toFixed(1) + '%' : 'N/A'}`);
  console.log(`PC2 Variance: ${stage3Output.pc2_variance ? (stage3Output.pc2_variance * 100).toFixed(1) + '%' : 'N/A'}`);
  console.log(`PCA Plot: ${stage3Output.pca_plot_path || 'N/A'}`);
  console.log(`QC Exit Code: ${stage3Output.exit_code} (0=PASS, 2=Issues detected by R script)`);
  if (stage3Output.warnings.length > 0) {
    console.log(`Warnings: ${stage3Output.warnings.join('; ')}`);
  }
  if (stage3Output.errors.length > 0) {
    console.log(`Errors: ${stage3Output.errors.join('; ')}`);
  }
  console.log('');

  // Step 4: Format output for agents
  console.log('[Test] Step 4: Formatting output for agents...');
  const formattedOutput = formatStage3ForAgents(stage3Output, dataInfo);

  // Note: formattedOutput is sent to agents, not shown to user (contains internal prompts)

  // Step 5: Call agents for review
  console.log('[Test] Step 5: Calling agents for Stage 3 review...');
  console.log('[Test] This will call GPT-5.2, Claude, and Gemini...');
  console.log('');

  try {
    const reviewResult = await coordinator.reviewStage3Output(stage3Output, dataInfo);

    console.log('');
    console.log('='.repeat(60));
    console.log('  Stage 3 Review Results');
    console.log('='.repeat(60));
    console.log('');
    console.log(`Decision: ${reviewResult.consensus.decision.toUpperCase()}`);
    console.log(`Confidence: ${(reviewResult.consensus.confidence * 100).toFixed(0)}%`);
    console.log(`Proceed to Stage 4: ${reviewResult.proceed ? 'YES' : 'NO'}`);
    console.log(`DE Method: ${reviewResult.deMethod}`);
    if (reviewResult.batchSpecification) {
      console.log(`Batch Specification: ${reviewResult.batchSpecification}`);
    }
    console.log(`Outlier Action: ${reviewResult.outlierAction}`);
    console.log(`Outliers to Remove: ${reviewResult.outliersToRemove.length > 0 ? reviewResult.outliersToRemove.join(', ') : 'None'}`);
    console.log('');
    console.log('Agent Votes:');
    console.log(`  GPT-5.2 (Stats): ${reviewResult.responses.gpt5_2.success ? 'Responded' : 'Failed'}`);
    console.log(`  Claude (Pipeline): ${reviewResult.responses.claude.success ? 'Responded' : 'Failed'}`);
    console.log(`  Gemini (Biology): ${reviewResult.responses.gemini.success ? 'Responded' : 'Failed'}`);
    console.log('');

    // Show agent decisions about batch effects and outliers
    console.log('Agent Decisions (from visual PCA inspection):');
    console.log(`  Batch Effects: ${reviewResult.deMethod === 'batch_effect_edger' ? 'YES (agents detected batch effects)' : 'NO (agents saw no batch effects)'}`);
    console.log(`  Outliers: ${reviewResult.outlierAction === 'REMOVE_OUTLIERS' ? `YES (agents recommend removing: ${reviewResult.outliersToRemove.join(', ')})` : 'NO (agents recommend keeping all samples)'}`);
    console.log('');

    // Check if user decision is required
    let finalProceed = reviewResult.proceed;
    let finalDeMethod = reviewResult.deMethod;
    let finalBatchSpecification = reviewResult.batchSpecification;
    let finalOutlierAction = reviewResult.outlierAction;
    let finalOutliersToRemove = reviewResult.outliersToRemove || [];
    let userDecision = null;

    if (reviewResult.consensus.decision.toLowerCase().includes('user_decision')) {
      // Agents couldn't decide - ask user
      userDecision = await handleStage3UserDecision(reviewResult, stage3Output);
      finalProceed = userDecision.proceed;
      finalDeMethod = userDecision.deMethod || finalDeMethod;
      finalBatchSpecification = userDecision.batchSpecification || finalBatchSpecification;
      finalOutlierAction = userDecision.outlierAction || finalOutlierAction;
      finalOutliersToRemove = userDecision.outliersToRemove || [];

      console.log('');
      console.log(`[Test] User decided: ${finalProceed ? 'PROCEED' : 'ABORT'}`);
      console.log(`[Test] DE Method: ${finalDeMethod}`);
      if (finalOutliersToRemove.length > 0) {
        console.log(`[Test] Outliers to remove: ${finalOutliersToRemove.join(', ')}`);
      }
      console.log(`[Test] Reason: ${userDecision.reason}`);
      console.log('');
    }

    // Step 6: Log the decision
    console.log('[Test] Step 6: Logging decision...');
    logger.logStageDecision(
      3,
      'QC Assessment',
      dataInfo,
      stage3Output,
      reviewResult.responses,
      reviewResult.consensus,
      {
        proceed: finalProceed,
        deMethod: finalDeMethod,
        batchSpecification: finalBatchSpecification,
        outlierAction: finalOutlierAction,
        outliersToRemove: finalOutliersToRemove,
        userDecision: userDecision ? {
          madeBy: 'user',
          proceed: userDecision.proceed,
          deMethod: userDecision.deMethod,
          outliersToRemove: userDecision.outliersToRemove,
          reason: userDecision.reason
        } : null
      }
    );

    // Finalize logger
    logger.finalize({
      dataInfo,
      config,
      success: finalProceed
    });

    console.log('');
    console.log('='.repeat(60));
    console.log('  Stage 3 Test Complete!');
    console.log('='.repeat(60));
    console.log('');
    console.log('Output files:');
    console.log(`  Script: ${scriptPath}`);
    console.log(`  Quantification: ${outputDir}/stage3_quantification/`);
    console.log(`  Logs: ${outputDir}/stage3_test_*.json`);
    console.log('');

    if (finalProceed) {
      console.log(`RESULT: Stage 3 PASSED - Ready for Stage 4 (DE Analysis)`);
      console.log(`        DE Method: ${finalDeMethod}`);
      if (finalOutliersToRemove.length > 0) {
        console.log(`        Remove outliers: ${finalOutliersToRemove.join(', ')}`);
      }
    } else {
      console.log('RESULT: Stage 3 FAILED - Analysis aborted');
    }

  } catch (error) {
    console.error('[Test] Error during agent review:', error);
    process.exit(1);
  }
}

main().catch(console.error);
