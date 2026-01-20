#!/usr/bin/env node
/**
 * Test Script for Stage 2: Alignment + Alignment QC
 *
 * This script tests the staged architecture by running Stage 2 in isolation:
 * 1. Generate the Stage 2 script (alignment + QC)
 * 2. Execute it (if not --dry-run)
 * 3. Parse the output
 * 4. Send to agents for review
 * 5. Log the decision
 *
 * Usage:
 *   node test_stage2.js <input_dir> <output_dir> [--dry-run] [--skip-execution]
 *
 * Example:
 *   node test_stage2.js data/DA0036 results/stage2_test
 *   node test_stage2.js data/DA0036 results/stage2_test --dry-run
 *
 * Note: Stage 2 requires Stage 1 to have run (FASTQ files validated)
 */

import path from 'path';
import fs from 'fs';
import { execSync, spawn } from 'child_process';
import dotenv from 'dotenv';

// Load environment variables
dotenv.config();

// Import stage 2 module
import {
  generateStage2Script,
  parseStage2Output,
  formatStage2ForAgents
} from './src/stages/stage2_alignment.js';

// Import coordinator and logger
import { Coordinator } from './src/coordinator/orchestrator.js';
import { Logger } from './src/utils/logger.js';
import { handleStage2UserDecision } from './src/utils/user_input.js';

// Parse command line arguments
const args = process.argv.slice(2);
const inputDir = args[0];
const outputDir = args[1];
const dryRun = args.includes('--dry-run');
const skipExecution = args.includes('--skip-execution');
const organism = args.find(a => a.startsWith('--organism='))?.split('=')[1] || 'mouse';

if (!inputDir || !outputDir) {
  console.log('Usage: node test_stage2.js <input_dir> <output_dir> [options]');
  console.log('');
  console.log('Options:');
  console.log('  --dry-run         Generate script but don\'t execute or call agents');
  console.log('  --skip-execution  Skip script execution, only call agents (if output exists)');
  console.log('  --organism=mouse  Specify organism (mouse, human, rat)');
  console.log('');
  console.log('Examples:');
  console.log('  node test_stage2.js data/DA0036 results/stage2_test');
  console.log('  node test_stage2.js data/DA0036 results/stage2_test --dry-run');
  console.log('  node test_stage2.js data/DA0036 results/stage2_test --organism=human');
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
  console.log('  Stage 2 Test: Alignment + Alignment QC');
  console.log('='.repeat(60));
  console.log('');

  // Create output directory
  if (!fs.existsSync(outputDir)) {
    fs.mkdirSync(outputDir, { recursive: true });
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
  const logger = new Logger(absOutputDir, 'stage2_test', config);

  // Initialize coordinator
  const coordinator = new Coordinator({ verbose: true });
  coordinator.setDataset(config.comparison);

  // Step 1: Generate Stage 2 script
  console.log('[Test] Step 1: Generating Stage 2 script...');
  const scriptPath = path.resolve(outputDir, `stage_2_${config.comparison}.sh`);
  const scriptResult = generateStage2Script(dataInfo, config, scriptPath);
  console.log(`[Test] Script generated: ${scriptPath}`);

  if (dryRun) {
    console.log('');
    console.log('[Test] DRY RUN - Script generated but not executed');
    console.log('[Test] To run the script manually:');
    console.log(`       bash ${scriptPath}`);
    console.log('');
    console.log('Script content (first 60 lines):');
    console.log('-'.repeat(60));
    const scriptContent = fs.readFileSync(scriptPath, 'utf-8');
    console.log(scriptContent.split('\n').slice(0, 60).join('\n'));
    console.log('...');
    console.log('-'.repeat(60));
    process.exit(0);
  }

  // Step 2: Execute Stage 2 script
  if (!skipExecution) {
    console.log('');
    console.log('[Test] Step 2: Executing Stage 2 script...');
    console.log('[Test] WARNING: This will run alignment which can take 10-30+ minutes!');
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
      console.log('[Test] Stage 2 script completed');
    } catch (error) {
      // Check if it's just a non-zero exit code (WARN or FAIL from QC)
      if (error.status === 1 || error.status === 2) {
        console.log('');
        console.log(`[Test] Stage 2 script completed with QC status code: ${error.status}`);
      } else {
        console.error('[Test] Error executing Stage 2 script:', error.message);
        process.exit(1);
      }
    }
  } else {
    console.log('[Test] Step 2: Skipping execution (--skip-execution flag)');
  }

  // Step 3: Parse Stage 2 output
  console.log('');
  console.log('[Test] Step 3: Parsing Stage 2 output...');
  const stage2Output = parseStage2Output(absOutputDir);
  console.log(`[Test] Parsed output: ${stage2Output.overall_status}`);

  // Show parsed output summary
  console.log('');
  console.log('--- Stage 2 Output Summary ---');
  console.log(`Samples Aligned: ${stage2Output.samples_aligned}`);
  console.log(`Mean Mapping Rate: ${stage2Output.mean_mapping_rate ? (stage2Output.mean_mapping_rate * 100).toFixed(1) + '%' : 'N/A'}`);
  console.log(`Median Mapping Rate: ${stage2Output.median_mapping_rate ? (stage2Output.median_mapping_rate * 100).toFixed(1) + '%' : 'N/A'}`);
  console.log(`Range: ${stage2Output.min_mapping_rate ? (stage2Output.min_mapping_rate * 100).toFixed(1) : '?'}% - ${stage2Output.max_mapping_rate ? (stage2Output.max_mapping_rate * 100).toFixed(1) : '?'}%`);
  console.log(`QC Exit Code: ${stage2Output.exit_code}`);
  console.log(`Overall Status: ${stage2Output.overall_status}`);
  console.log(`Passed: ${stage2Output.passed_samples.length}, Warned: ${stage2Output.warned_samples.length}, Failed: ${stage2Output.failed_samples.length}`);
  console.log('');

  // Step 4: Format output for agents
  console.log('[Test] Step 4: Formatting output for agents...');
  const formattedOutput = formatStage2ForAgents(stage2Output, dataInfo);

  // Note: formattedOutput is sent to agents, not shown to user (contains internal prompts)

  // Step 5: Call agents for review
  console.log('[Test] Step 5: Calling agents for Stage 2 review...');
  console.log('[Test] This will call GPT-5.2, Claude, and Gemini...');
  console.log('');

  try {
    const reviewResult = await coordinator.reviewStage2Output(stage2Output, dataInfo);

    console.log('');
    console.log('='.repeat(60));
    console.log('  Stage 2 Review Results');
    console.log('='.repeat(60));
    console.log('');
    console.log(`Decision: ${reviewResult.consensus.decision.toUpperCase()}`);
    console.log(`Confidence: ${(reviewResult.consensus.confidence * 100).toFixed(0)}%`);
    console.log(`Proceed to Stage 3: ${reviewResult.proceed ? 'YES' : 'NO'}`);
    console.log(`Samples to Remove: ${reviewResult.samplesToRemove.length > 0 ? reviewResult.samplesToRemove.join(', ') : 'None'}`);
    console.log('');
    console.log('Agent Votes:');
    console.log(`  GPT-5.2 (Stats): ${reviewResult.responses.gpt5_2.success ? 'Responded' : 'Failed'}`);
    console.log(`  Claude (Pipeline): ${reviewResult.responses.claude.success ? 'Responded' : 'Failed'}`);
    console.log(`  Gemini (Biology): ${reviewResult.responses.gemini.success ? 'Responded' : 'Failed'}`);
    console.log('');

    // Check if user decision is required
    let finalProceed = reviewResult.proceed;
    let finalSamplesToRemove = reviewResult.samplesToRemove || [];
    let userDecision = null;

    if (reviewResult.consensus.decision.toLowerCase().includes('user_decision')) {
      // Agents couldn't decide - ask user
      userDecision = await handleStage2UserDecision(reviewResult, stage2Output);
      finalProceed = userDecision.proceed;
      finalSamplesToRemove = userDecision.samplesToRemove || [];

      console.log('');
      console.log(`[Test] User decided: ${finalProceed ? 'PROCEED' : 'ABORT'}`);
      if (finalSamplesToRemove.length > 0) {
        console.log(`[Test] Samples to remove: ${finalSamplesToRemove.join(', ')}`);
      }
      console.log(`[Test] Reason: ${userDecision.reason}`);
      console.log('');
    }

    // Step 6: Log the decision
    console.log('[Test] Step 6: Logging decision...');
    logger.logStageDecision(
      2,
      'Alignment QC',
      dataInfo,
      stage2Output,
      reviewResult.responses,
      reviewResult.consensus,
      {
        proceed: finalProceed,
        samplesToRemove: finalSamplesToRemove,
        userDecision: userDecision ? {
          madeBy: 'user',
          proceed: userDecision.proceed,
          samplesToRemove: userDecision.samplesToRemove,
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
    console.log('  Stage 2 Test Complete!');
    console.log('='.repeat(60));
    console.log('');
    console.log('Output files:');
    console.log(`  Script: ${scriptPath}`);
    console.log(`  BAM files: ${outputDir}/stage2_alignment/bam_files/`);
    console.log(`  QC: ${outputDir}/stage2_alignment/alignment_qc/`);
    console.log(`  Logs: ${outputDir}/stage2_test_*.json`);
    console.log('');

    if (finalProceed) {
      if (finalSamplesToRemove.length > 0) {
        console.log(`RESULT: Stage 2 PASSED with sample removal - Ready for Stage 3`);
        console.log(`        Remove: ${finalSamplesToRemove.join(', ')}`);
      } else {
        console.log('RESULT: Stage 2 PASSED - Ready for Stage 3 (Quantification + QC)');
      }
    } else {
      console.log('RESULT: Stage 2 FAILED - Analysis aborted');
    }

  } catch (error) {
    console.error('[Test] Error during agent review:', error);
    process.exit(1);
  }
}

main().catch(console.error);
