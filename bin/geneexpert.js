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

const program = new Command();

program
  .name('geneexpert')
  .description('Multi-agent RNA-seq analysis system')
  .version('0.1.0');

program
  .command('analyze')
  .description('Run RNA-seq analysis with multi-agent orchestration')
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
  .option('--single-agent <name>', 'Use only one agent: gpt5.2, claude, or gemini (for experiments)')
  .option('--force-automation', 'Skip all agents, use template-based AUTOMATION only (no-agent baseline)')
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
      verbose: options.verbose
    };

    console.log('📋 Configuration:');
    console.log(`   Input:      ${config.input}`);
    console.log(`   Output:     ${config.output}`);
    console.log(`   Organism:   ${config.organism}`);
    console.log(`   Comparison: ${config.comparison}`);
    console.log(`   Aligner:    ${config.aligner}`);
    console.log(`   DE Tool:    ${config.deTool}`);
    console.log('');

    try {
      // Execute multi-agent orchestrated analysis
      await executeAnalysis(config);

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
