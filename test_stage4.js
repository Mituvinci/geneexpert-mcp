#!/usr/bin/env node
/**
 * Test Script for Stage 4: Differential Expression Analysis
 *
 * This script tests the staged architecture by running Stage 4 in isolation:
 * 1. Generate the Stage 4 script (DE analysis)
 * 2. Execute it
 * 3. Parse the output
 * 4. Send to agents for review
 * 5. Log the decision
 *
 * Usage:
 *   node test_stage4.js <input_dir> <output_dir> [--dry-run] [--skip-execution] [--de-method=simpleEdger]
 *
 * Example:
 *   node test_stage4.js data/DA0036 results/stage4_test
 *   node test_stage4.js data/DA0036 results/stage4_test --de-method=batch_effect_edger
 *
 * Note: Stage 4 requires Stage 3 to have run (count files needed)
 */

import path from 'path';
import fs from 'fs';
import { execSync } from 'child_process';
import dotenv from 'dotenv';

// Load environment variables
dotenv.config();

// Import stage 4 module
import {
  generateStage4Script,
  parseStage4Output,
  formatStage4ForAgents
} from './src/stages/stage4_de_analysis.js';

// Import coordinator and logger
import { Coordinator } from './src/coordinator/orchestrator.js';
import { Logger } from './src/utils/logger.js';
import { askYesNo } from './src/utils/user_input.js';

// Parse command line arguments
const args = process.argv.slice(2);
const inputDir = args[0];
const outputDir = args[1];
const dryRun = args.includes('--dry-run');
const skipExecution = args.includes('--skip-execution');
const organism = args.find(a => a.startsWith('--organism='))?.split('=')[1] || 'mouse';
const deMethod = args.find(a => a.startsWith('--de-method='))?.split('=')[1] || 'simpleEdger';
const batchSpec = args.find(a => a.startsWith('--batch-spec='))?.split('=')[1] || 'auto';

if (!inputDir || !outputDir) {
  console.log('Usage: node test_stage4.js <input_dir> <output_dir> [options]');
  console.log('');
  console.log('Options:');
  console.log('  --dry-run                      Generate script but don\'t execute or call agents');
  console.log('  --skip-execution               Skip script execution, only call agents (if output exists)');
  console.log('  --organism=mouse               Specify organism (mouse, human, rat)');
  console.log('  --de-method=simpleEdger        DE method (simpleEdger or batch_effect_edger)');
  console.log('  --batch-spec=auto              Batch specification (auto, paired, or explicit)');
  console.log('');
  console.log('Examples:');
  console.log('  node test_stage4.js data/DA0036 results/stage4_test');
  console.log('  node test_stage4.js data/DA0036 results/stage4_test --de-method=batch_effect_edger');
  console.log('  node test_stage4.js data/DA0036 results/stage4_test --dry-run');
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
  console.log('  Stage 4 Test: Differential Expression Analysis');
  console.log('='.repeat(60));
  console.log('');

  // Create output directory
  if (!fs.existsSync(outputDir)) {
    fs.mkdirSync(outputDir, { recursive: true });
  }

  // Check if Stage 3 count files exist
  const stage3Dir = path.join(outputDir, 'stage3_quantification');
  if (!fs.existsSync(stage3Dir)) {
    console.error('ERROR: Stage 3 quantification files not found!');
    console.error(`Expected: ${stage3Dir}`);
    console.error('');
    console.error('You must run Stage 3 (quantification) before Stage 4.');
    console.error('Run: node test_stage3.js data/DA0036 results/stage4_test');
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

  // Mock Stage 3 result (for script generation)
  const stage3Result = {
    deMethod,
    batchSpecification: batchSpec
  };

  // Initialize logger (use absolute path)
  const absOutputDir = path.resolve(outputDir);
  const logger = new Logger(absOutputDir, 'stage4_test', config);

  // Initialize coordinator
  const coordinator = new Coordinator({ verbose: true });
  coordinator.setDataset(config.comparison);

  // Step 1: Generate Stage 4 script
  console.log('[Test] Step 1: Generating Stage 4 script...');
  console.log(`[Test] DE Method: ${deMethod}`);
  if (deMethod === 'batch_effect_edger') {
    console.log(`[Test] Batch Specification: ${batchSpec}`);
  }
  const scriptPath = path.resolve(outputDir, `stage_4_${config.comparison}.sh`);
  const scriptResult = generateStage4Script(dataInfo, config, scriptPath, stage3Result);
  console.log(`[Test] Script generated: ${scriptPath}`);

  if (dryRun) {
    console.log('');
    console.log('[Test] DRY RUN - Script generated but not executed');
    console.log('[Test] To run the script manually:');
    console.log(`       bash ${scriptPath}`);
    console.log('');
    console.log('Script content:');
    console.log('-'.repeat(60));
    const scriptContent = fs.readFileSync(scriptPath, 'utf-8');
    console.log(scriptContent);
    console.log('-'.repeat(60));
    process.exit(0);
  }

  // Step 2: Execute Stage 4 script
  if (!skipExecution) {
    console.log('');
    console.log('[Test] Step 2: Executing Stage 4 script...');
    console.log('[Test] WARNING: This will run DE analysis which can take 5-15+ minutes!');
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
      console.log('[Test] Stage 4 script completed successfully');
    } catch (error) {
      console.error('[Test] Error executing Stage 4 script:', error.message);
      process.exit(1);
    }
  } else {
    console.log('[Test] Step 2: Skipping execution (--skip-execution flag)');
  }

  // Step 3: Parse Stage 4 output
  console.log('');
  console.log('[Test] Step 3: Parsing Stage 4 output...');
  const stage4Output = parseStage4Output(absOutputDir);
  console.log(`[Test] Parsed output: ${stage4Output.overall_status}`);

  // Show parsed output summary
  console.log('');
  console.log('--- Stage 4 Output Summary ---');
  console.log(`Total Genes Tested: ${stage4Output.total_genes_tested?.toLocaleString()}`);
  console.log(`Total DEGs (FDR < ${stage4Output.fdr_threshold}): ${stage4Output.total_degs}`);
  console.log(`Up-regulated: ${stage4Output.num_degs_up}`);
  console.log(`Down-regulated: ${stage4Output.num_degs_down}`);
  if (stage4Output.fc_range[0] !== null) {
    console.log(`LogFC Range: ${stage4Output.fc_range[0].toFixed(2)} to ${stage4Output.fc_range[1].toFixed(2)}`);
  }
  console.log(`Overall Status: ${stage4Output.overall_status}`);
  console.log(`Warnings: ${stage4Output.warnings.length > 0 ? stage4Output.warnings.join('; ') : 'None'}`);
  console.log(`Errors: ${stage4Output.errors.length > 0 ? stage4Output.errors.join('; ') : 'None'}`);

  // Show top genes
  if (stage4Output.top_up_genes.length > 0) {
    console.log('');
    console.log('Top 5 Up-regulated Genes:');
    stage4Output.top_up_genes.slice(0, 5).forEach(g => {
      console.log(`  ${g.name}: LogFC=${g.logFC}, FDR=${g.fdr}`);
    });
  }

  if (stage4Output.top_down_genes.length > 0) {
    console.log('');
    console.log('Top 5 Down-regulated Genes:');
    stage4Output.top_down_genes.slice(0, 5).forEach(g => {
      console.log(`  ${g.name}: LogFC=${g.logFC}, FDR=${g.fdr}`);
    });
  }
  console.log('');

  // Step 4: Format output for agents
  console.log('[Test] Step 4: Formatting output for agents...');
  const formattedOutput = formatStage4ForAgents(stage4Output, dataInfo);

  // Note: formattedOutput is sent to agents, not shown to user (contains internal prompts)

  // Step 5: Call agents for review
  console.log('[Test] Step 5: Calling agents for Stage 4 review...');
  console.log('[Test] This will call GPT-5.2, Claude, and Gemini...');
  console.log('');

  try {
    const reviewResult = await coordinator.reviewStage4Output(stage4Output, dataInfo);

    console.log('');
    console.log('='.repeat(60));
    console.log('  Stage 4 Review Results');
    console.log('='.repeat(60));
    console.log('');
    console.log(`Decision: ${reviewResult.consensus.decision.toUpperCase()}`);
    console.log(`Confidence: ${(reviewResult.consensus.confidence * 100).toFixed(0)}%`);
    console.log(`Final Approval: ${reviewResult.approve ? 'APPROVED' : 'REQUEST REANALYSIS'}`);
    console.log('');
    console.log('Agent Votes:');
    console.log(`  GPT-5.2 (Stats): ${reviewResult.responses.gpt5_2.success ? 'Responded' : 'Failed'}`);
    console.log(`  Claude (Pipeline): ${reviewResult.responses.claude.success ? 'Responded' : 'Failed'}`);
    console.log(`  Gemini (Biology): ${reviewResult.responses.gemini.success ? 'Responded' : 'Failed'}`);
    console.log('');

    // Check if user decision is required
    let finalApprove = reviewResult.approve;
    let userDecision = null;

    if (reviewResult.consensus.decision.toLowerCase().includes('user_decision')) {
      // Agents couldn't decide - ask user
      const userChoice = await askYesNo(
        'Agents could not reach consensus on DE results. Do you approve the results?',
        true
      );
      finalApprove = userChoice;
      userDecision = {
        madeBy: 'user',
        approve: userChoice,
        reason: userChoice ? 'User approved DE results' : 'User rejected DE results'
      };

      console.log('');
      console.log(`[Test] User decided: ${finalApprove ? 'APPROVE' : 'REJECT'}`);
      console.log(`[Test] Reason: ${userDecision.reason}`);
      console.log('');
    }

    // Step 6: Log the decision
    console.log('[Test] Step 6: Logging decision...');
    logger.logStageDecision(
      4,
      'DE Analysis',
      dataInfo,
      stage4Output,
      reviewResult.responses,
      reviewResult.consensus,
      {
        proceed: finalApprove,
        approve: finalApprove,
        userDecision: userDecision
      }
    );

    // Finalize logger
    logger.finalize({
      dataInfo,
      config,
      success: finalApprove
    });

    console.log('');
    console.log('='.repeat(60));
    console.log('  Stage 4 Test Complete!');
    console.log('='.repeat(60));
    console.log('');
    console.log('Output files:');
    console.log(`  Script: ${scriptPath}`);
    console.log(`  DE Results: ${outputDir}/stage4_de_analysis/`);
    console.log(`  Logs: ${outputDir}/stage4_test_*.json`);
    console.log('');

    if (finalApprove) {
      console.log('RESULT: Stage 4 APPROVED - Analysis Complete!');
      console.log(`        Total DEGs: ${stage4Output.total_degs} (FDR < ${stage4Output.fdr_threshold})`);
    } else {
      console.log('RESULT: Stage 4 REJECTED - Reanalysis needed');
    }

  } catch (error) {
    console.error('[Test] Error during agent review:', error);
    process.exit(1);
  }
}

main().catch(console.error);
