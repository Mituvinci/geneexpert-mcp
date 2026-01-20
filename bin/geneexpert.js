#!/usr/bin/env node

/**
 * GeneExpert CLI - Multi-Agent RNA-seq Analysis
 *
 * Usage:
 *   geneexpert analyze <input_dir> --output <output_dir> --organism <species> [options]
 *
 * Example:
 *   geneexpert analyze /destiny/fastqs/DA0036 \
 *     --output /data/halimaakhter/results/DA0036 \
 *     --organism mouse \
 *     --comparison "stroke_vs_control"
 */

import { Command } from 'commander';
import path from 'path';
import fs from 'fs';
import { executeAnalysis } from '../src/pipeline/executor.js';
import { StagedExecutor } from '../src/executor/staged_executor.js';

const program = new Command();

program
  .name('geneexpert')
  .description('Multi-agent RNA-seq analysis system')
  .version('0.1.0');

program
  .command('analyze')
  .description(`Run RNA-seq analysis with multi-agent orchestration

ICML 2026 EXPERIMENTAL MODES:
  1. No-Agent (Baseline):       --staged --force-automation
  2. Single-Agent GPT-5.2:      --staged --single-agent gpt5.2
  3. Single-Agent Claude:       --staged --single-agent claude
  4. Multi-Agent Parallel:      --staged (default, all agents vote independently)
  5. Multi-Agent Sequential:    --staged --sequential-chain (GPT→Gemini→Claude)

EXAMPLES:
  # Multi-agent parallel (default - recommended)
  geneexpert analyze data/DA0036 -o results/test --organism mouse --staged

  # Sequential chain mode (experimental)
  geneexpert analyze data/DA0036 -o results/test --organism mouse --staged --sequential-chain

  # Single-agent baseline (for comparison)
  geneexpert analyze data/DA0036 -o results/test --organism mouse --staged --single-agent gpt5.2

  # No-agent baseline (template only)
  geneexpert analyze data/DA0036 -o results/test --organism mouse --staged --force-automation
`)
  .argument('<input>', 'Input directory (FASTQ files, BAM files, or count matrix)')
  .requiredOption('-o, --output <dir>', 'Output directory for results')
  .option('--organism <species>', 'Organism (mouse, human, rat)', 'mouse')
  .option('--comparison <name>', 'Comparison name (e.g., "treatment_vs_control")', 'comparison')
  .option('--control-keyword <keyword>', 'Keyword in sample names for control group (e.g., "ctrl", "wt", "cont")')
  .option('--treatment-keyword <keyword>', 'Keyword in sample names for treatment group (e.g., "ko", "treat", "ips")')
  .option('--genome <path>', 'Path to genome reference (for alignment)')
  .option('--gtf <path>', 'Path to GTF annotation file (for quantification)')
  .option('--aligner <tool>', 'Alignment tool (subread only)', 'subread')
  .option('--de-tool <tool>', 'DE analysis tool: edger or deseq2', 'edger')
  .option('--threads <n>', 'Number of threads', '4')
  .option('--staged', 'Use staged architecture (4-stage pipeline with agent checkpoints) - REQUIRED for ICML experiments')
  .option('--single-agent <name>', 'Use only one agent: gpt5.2, claude, or gemini (ICML baseline #2-3)')
  .option('--force-automation', 'Skip all agents, use template-based decisions only (ICML baseline #1)')
  .option('--sequential-chain', 'Use sequential chain mode: GPT-5.2 → Gemini → Claude (ICML experiment #5)')
  .option('--verbose', 'Verbose output', false)
  .action(async (input, options) => {
    console.log('🧬 GeneExpert Multi-Agent RNA-seq Analysis');
    console.log('='.repeat(60));
    console.log('');

    // Validate input directory
    if (!fs.existsSync(input)) {
      console.error(`❌ Error: Input directory not found: ${input}`);
      process.exit(1);
    }

    // Validate required flags for staged mode
    if (options.staged && (!options.controlKeyword || !options.treatmentKeyword)) {
      console.error('');
      console.error('❌ Error: --staged mode requires both --control-keyword and --treatment-keyword');
      console.error('');
      console.error('Example:');
      console.error('  geneexpert analyze data/DA0036 -o results/test \\');
      console.error('    --organism mouse --staged \\');
      console.error('    --control-keyword "cont" \\');
      console.error('    --treatment-keyword "ips"');
      console.error('');
      console.error('The keywords are used to identify control vs treatment samples from filenames.');
      console.error('');
      process.exit(1);
    }

    // Create output directory
    const outputDir = path.resolve(options.output);
    if (!fs.existsSync(outputDir)) {
      fs.mkdirSync(outputDir, { recursive: true });
      console.log(`✓ Created output directory: ${outputDir}`);
    }

    // Build analysis configuration
    const config = {
      input: path.resolve(input),
      output: outputDir,
      organism: options.organism,
      comparison: options.comparison,
      controlKeyword: options.controlKeyword,
      treatmentKeyword: options.treatmentKeyword,
      genome: options.genome,
      gtf: options.gtf,
      aligner: options.aligner,
      deTool: options.deTool,
      threads: parseInt(options.threads),
      singleAgent: options.singleAgent, // For experimental comparisons
      forceAutomation: options.forceAutomation, // Skip agents, use template only
      sequentialChain: options.sequentialChain, // NEW: Sequential chain mode
      verbose: options.verbose
    };

    // Detect experimental mode
    let experimentalMode = 'Unknown';
    let modeNumber = 0;

    if (!options.staged) {
      experimentalMode = 'Legacy Monolithic (not for ICML experiments)';
      console.log('⚠️  WARNING: You are using the old monolithic architecture.');
      console.log('   For ICML experiments, use --staged flag!');
      console.log('');
    } else {
      // Detect which ICML experimental mode
      if (options.forceAutomation) {
        modeNumber = 1;
        experimentalMode = 'No-Agent Baseline (template-based decisions)';
        if (options.singleAgent || options.sequentialChain) {
          console.log('⚠️  WARNING: --force-automation conflicts with --single-agent or --sequential-chain');
          console.log('   Using: No-Agent mode (agents will be skipped)');
          console.log('');
        }
      } else if (options.singleAgent) {
        if (options.singleAgent === 'gpt5.2' || options.singleAgent === 'gpt4') {
          modeNumber = 2;
          experimentalMode = `Single-Agent: GPT-5.2 (all roles)`;
        } else if (options.singleAgent === 'claude') {
          modeNumber = 3;
          experimentalMode = `Single-Agent: Claude (all roles)`;
        } else if (options.singleAgent === 'gemini') {
          modeNumber = 3;
          experimentalMode = `Single-Agent: Gemini (all roles)`;
        } else {
          console.error(`❌ Error: Invalid --single-agent value: ${options.singleAgent}`);
          console.error('   Valid options: gpt5.2, claude, gemini');
          process.exit(1);
        }
        if (options.sequentialChain) {
          console.log('⚠️  WARNING: --single-agent conflicts with --sequential-chain');
          console.log('   Using: Single-agent mode (sequential chain will be ignored)');
          console.log('');
        }
      } else if (options.sequentialChain) {
        modeNumber = 5;
        experimentalMode = 'Multi-Agent Sequential Chain (GPT-5.2 → Gemini → Claude)';
      } else {
        modeNumber = 4;
        experimentalMode = 'Multi-Agent Parallel (default - independent voting)';
      }
    }

    console.log('📋 Configuration:');
    console.log(`   Input:      ${config.input}`);
    console.log(`   Output:     ${config.output}`);
    console.log(`   Organism:   ${config.organism}`);
    console.log(`   Comparison: ${config.comparison}`);
    if (config.controlKeyword) {
      console.log(`   Control:    ${config.controlKeyword}`);
    }
    if (config.treatmentKeyword) {
      console.log(`   Treatment:  ${config.treatmentKeyword}`);
    }
    console.log(`   Aligner:    ${config.aligner}`);
    console.log(`   DE Tool:    ${config.deTool}`);
    console.log(`   Architecture: ${options.staged ? 'Staged (4-stage checkpoints)' : 'Monolithic'}`);
    console.log('');
    console.log(`🔬 ICML Experiment Mode #${modeNumber}: ${experimentalMode}`);
    console.log('');

    try {
      // Execute analysis - use staged or monolithic architecture
      if (options.staged) {
        // Use new staged architecture (4-stage pipeline with agent checkpoints)
        const executor = new StagedExecutor({
          inputDir: config.input,
          outputDir: config.output,
          organism: config.organism,
          comparison: config.comparison,
          controlKeyword: config.controlKeyword,
          treatmentKeyword: config.treatmentKeyword,
          verbose: config.verbose,
          singleAgent: config.singleAgent,
          forceAutomation: config.forceAutomation,
          sequentialChain: config.sequentialChain  // NEW: Pass sequential chain flag
        });
        await executor.run();
      } else {
        // Use old monolithic architecture
        await executeAnalysis(config);
      }

      console.log('');
      console.log('='.repeat(60));
      console.log('✅ Analysis complete!');
      console.log(`📊 Results: ${config.output}`);

    } catch (error) {
      console.error('');
      console.error('='.repeat(60));
      console.error('❌ Analysis failed:', error.message);
      if (config.verbose) {
        console.error(error.stack);
      }
      process.exit(1);
    }
  });

program
  .command('config')
  .description('Configure GeneExpert (API keys, default settings)')
  .option('--setup', 'Interactive setup wizard')
  .action((options) => {
    if (options.setup) {
      console.log('🔧 GeneExpert Configuration Setup');
      console.log('='.repeat(60));
      console.log('');
      console.log('Please configure your API keys in .env file:');
      console.log('  OPENAI_API_KEY=sk-...');
      console.log('  ANTHROPIC_API_KEY=sk-ant-...');
      console.log('  GOOGLE_API_KEY=AIza...');
      console.log('');
      console.log('See .env.example for reference');
    }
  });

program.parse();
