#!/usr/bin/env node
/**
 * Test Script for Stage 1: FASTQ Validation
 *
 * This script tests the staged architecture by running Stage 1 in isolation:
 * 1. Generate the Stage 1 script (FASTQ validation + FastQC)
 * 2. Execute it
 * 3. Parse the output
 * 4. Send to agents for review
 * 5. Log the decision
 *
 * Usage:
 *   node test_stage1.js <input_dir> <output_dir> [--dry-run] [--skip-execution]
 *
 * Example:
 *   node test_stage1.js data/DA0036 results/stage1_test
 *   node test_stage1.js data/DA0036 results/stage1_test --dry-run
 */

import path from 'path';
import fs from 'fs';
import { execSync } from 'child_process';
import dotenv from 'dotenv';

// Load environment variables
dotenv.config();

// Import stage 1 module
import {
  generateStage1Script,
  parseStage1Output,
  formatStage1ForAgents
} from './src/stages/stage1_fastq_validation.js';

// Import coordinator and logger
import { Coordinator } from './src/coordinator/orchestrator.js';
import { Logger } from './src/utils/logger.js';
import { handleStage1UserDecision } from './src/utils/user_input.js';

// Parse command line arguments
const args = process.argv.slice(2);
const inputDir = args[0];
const outputDir = args[1];
const dryRun = args.includes('--dry-run');
const skipExecution = args.includes('--skip-execution');

if (!inputDir || !outputDir) {
  console.log('Usage: node test_stage1.js <input_dir> <output_dir> [--dry-run] [--skip-execution]');
  console.log('');
  console.log('Options:');
  console.log('  --dry-run         Generate script but don\'t execute or call agents');
  console.log('  --skip-execution  Skip script execution, only call agents (if output exists)');
  console.log('');
  console.log('Examples:');
  console.log('  node test_stage1.js data/DA0036 results/stage1_test');
  console.log('  node test_stage1.js data/DA0036 results/stage1_test --dry-run');
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
    // Extract sample name (remove _R1_001.fastq.gz or _R2_001.fastq.gz)
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

  // Convert to array
  for (const [name, sample] of Object.entries(sampleMap)) {
    samples.push(sample);
  }

  // Detect groups (simple heuristic based on common keywords)
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
    organism: 'mouse',  // Default, can be overridden
    comparison: path.basename(outputDir)
  };

  console.log(`[Test] Detected ${samples.length} samples (${pairedEnd ? 'paired-end' : 'single-end'})`);
  console.log(`[Test] Groups: ${Object.keys(groups).map(g => `${g}(n=${groups[g].length})`).join(', ')}`);

  return dataInfo;
}

async function main() {
  console.log('');
  console.log('='.repeat(60));
  console.log('  Stage 1 Test: FASTQ Validation');
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
    organism: 'mouse',
    comparison: path.basename(outputDir),
    controlKeyword: 'cont',
    treatmentKeyword: 'ips'
  };

  // Initialize logger (use absolute path)
  const absOutputDir = path.resolve(outputDir);
  const logger = new Logger(absOutputDir, 'stage1_test', config);

  // Initialize coordinator
  const coordinator = new Coordinator({ verbose: true });
  coordinator.setDataset(config.comparison);

  // Step 1: Generate Stage 1 script
  console.log('[Test] Step 1: Generating Stage 1 script...');
  const scriptPath = path.resolve(outputDir, `stage_1_${config.comparison}.sh`);
  const scriptResult = generateStage1Script(dataInfo, config, scriptPath);
  console.log(`[Test] Script generated: ${scriptPath}`);

  if (dryRun) {
    console.log('');
    console.log('[Test] DRY RUN - Script generated but not executed');
    console.log('[Test] To run the script manually:');
    console.log(`       bash ${scriptPath}`);
    console.log('');
    console.log('Script content (first 50 lines):');
    console.log('-'.repeat(60));
    const scriptContent = fs.readFileSync(scriptPath, 'utf-8');
    console.log(scriptContent.split('\n').slice(0, 50).join('\n'));
    console.log('...');
    console.log('-'.repeat(60));
    process.exit(0);
  }

  // Step 2: Execute Stage 1 script
  if (!skipExecution) {
    console.log('');
    console.log('[Test] Step 2: Executing Stage 1 script...');
    console.log('[Test] This may take a few minutes for FastQC...');
    console.log('');

    try {
      execSync(`bash ${scriptPath}`, {
        stdio: 'inherit',
        cwd: outputDir
      });
      console.log('');
      console.log('[Test] Stage 1 script executed successfully');
    } catch (error) {
      console.error('[Test] Error executing Stage 1 script:', error.message);
      process.exit(1);
    }
  } else {
    console.log('[Test] Step 2: Skipping execution (--skip-execution flag)');
  }

  // Step 3: Parse Stage 1 output
  console.log('');
  console.log('[Test] Step 3: Parsing Stage 1 output...');
  const stage1Output = parseStage1Output(absOutputDir);
  console.log(`[Test] Parsed output: ${stage1Output.overall_status}`);

  // Show parsed output summary
  console.log('');
  console.log('--- Stage 1 Output Summary ---');
  console.log(`Validation Status: ${stage1Output.validation_status}`);
  console.log(`FastQC Status: ${stage1Output.fastqc_status}`);
  console.log(`Overall Status: ${stage1Output.overall_status}`);
  if (stage1Output.read_statistics) {
    console.log(`Total Reads: ${stage1Output.read_statistics.total_reads?.toLocaleString()}`);
    console.log(`Avg Reads/Sample: ${stage1Output.read_statistics.avg_reads?.toLocaleString()}`);
  }
  console.log(`Warnings: ${stage1Output.warnings.length > 0 ? stage1Output.warnings.join('; ') : 'None'}`);
  console.log(`Errors: ${stage1Output.errors.length > 0 ? stage1Output.errors.join('; ') : 'None'}`);
  console.log('');

  // Step 4: Format output for agents
  console.log('[Test] Step 4: Formatting output for agents...');
  const formattedOutput = formatStage1ForAgents(stage1Output, dataInfo);

  // Note: formattedOutput is sent to agents, not shown to user (contains internal prompts)

  // Step 5: Call agents for review
  console.log('[Test] Step 5: Calling agents for Stage 1 review...');
  console.log('[Test] This will call GPT-5.2, Claude, and Gemini...');
  console.log('');

  try {
    const reviewResult = await coordinator.reviewStage1Output(stage1Output, dataInfo);

    console.log('');
    console.log('='.repeat(60));
    console.log('  Stage 1 Review Results');
    console.log('='.repeat(60));
    console.log('');
    console.log(`Decision: ${reviewResult.consensus.decision.toUpperCase()}`);
    console.log(`Confidence: ${(reviewResult.consensus.confidence * 100).toFixed(0)}%`);
    console.log(`Proceed to Stage 2: ${reviewResult.proceed ? 'YES' : 'NO'}`);
    console.log('');
    console.log('Agent Votes:');
    console.log(`  GPT-5.2 (Stats): ${reviewResult.responses.gpt5_2.success ? 'Responded' : 'Failed'}`);
    console.log(`  Claude (Pipeline): ${reviewResult.responses.claude.success ? 'Responded' : 'Failed'}`);
    console.log(`  Gemini (Biology): ${reviewResult.responses.gemini.success ? 'Responded' : 'Failed'}`);
    console.log('');

    // Check if user decision is required
    let finalProceed = reviewResult.proceed;
    let userDecision = null;

    if (reviewResult.consensus.decision.toLowerCase().includes('user_decision')) {
      // Agents couldn't decide - ask user
      userDecision = await handleStage1UserDecision(reviewResult, stage1Output);
      finalProceed = userDecision.proceed;

      console.log('');
      console.log(`[Test] User decided: ${finalProceed ? 'PROCEED' : 'ABORT'}`);
      console.log(`[Test] Reason: ${userDecision.reason}`);
      console.log('');
    }

    // Step 6: Log the decision
    console.log('[Test] Step 6: Logging decision...');
    logger.logStageDecision(
      1,
      'FASTQ Validation',
      dataInfo,
      stage1Output,
      reviewResult.responses,
      reviewResult.consensus,
      {
        proceed: finalProceed,
        userDecision: userDecision ? {
          madeBy: 'user',
          proceed: userDecision.proceed,
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
    console.log('  Stage 1 Test Complete!');
    console.log('='.repeat(60));
    console.log('');
    console.log('Output files:');
    console.log(`  Script: ${scriptPath}`);
    console.log(`  Validation: ${outputDir}/stage1_validation/`);
    console.log(`  Logs: ${outputDir}/stage1_test_*.json`);
    console.log('');

    if (finalProceed) {
      console.log('RESULT: Stage 1 PASSED - Ready for Stage 2 (Alignment)');
    } else {
      console.log('RESULT: Stage 1 FAILED - Check errors before proceeding');
    }

  } catch (error) {
    console.error('[Test] Error during agent review:', error);
    process.exit(1);
  }
}

main().catch(console.error);
