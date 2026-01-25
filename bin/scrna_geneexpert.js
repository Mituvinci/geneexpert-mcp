#!/usr/bin/env node

/**
 * scRNA-seq GeneExpert CLI - Multi-Agent scRNA-seq Analysis
 *
 * Usage:
 *   scrna_geneexpert analyze <input_dir> --output <output_dir> --organism <species> [options]
 *
 * Example:
 *   node bin/scrna_geneexpert.js analyze data/scRNA_data/10-k-brain-cells_healthy_mouse \
 *     --output results/scrna_test \
 *     --organism mouse
 */

import { Command } from 'commander';
import path from 'path';
import fs from 'fs';
import { ScRNAExecutor } from '../src/scrna_executor/scrna_executor.js';

const program = new Command();

program
  .name('scrna_geneexpert')
  .description('Multi-agent scRNA-seq analysis system (Seurat)')
  .version('0.1.0');

program
  .command('analyze')
  .description(`Run scRNA-seq analysis with multi-agent orchestration

ICML 2026 EXPERIMENTAL MODES:
  1. No-Agent (Baseline):       --force-automation
  2. Single-Agent GPT-5.2:      --single-agent gpt5.2
  3. Single-Agent Claude:       --single-agent claude
  4. Multi-Agent Parallel:      (default, all agents vote independently)
  5. Multi-Agent Sequential:    --sequential-chain (GPT→Gemini→Claude)

EXAMPLES:
  # Multi-agent parallel (default - recommended)
  node bin/scrna_geneexpert.js analyze data/scRNA_data/10-k-brain-cells_healthy_mouse \\
    --output results/scrna_test --organism mouse

  # Sequential chain mode (experimental)
  node bin/scrna_geneexpert.js analyze data/scRNA_data/10-k-brain-cells_healthy_mouse \\
    --output results/scrna_test --organism mouse --sequential-chain

  # Single-agent baseline (for comparison)
  node bin/scrna_geneexpert.js analyze data/scRNA_data/10-k-brain-cells_healthy_mouse \\
    --output results/scrna_test --organism mouse --single-agent gpt5.2

  # No-agent baseline (template only)
  node bin/scrna_geneexpert.js analyze data/scRNA_data/10-k-brain-cells_healthy_mouse \\
    --output results/scrna_test --organism mouse --force-automation
`)
  .argument('<input>', 'Input directory containing filtered_feature_bc_matrix.h5 or matrix files')
  .requiredOption('-o, --output <dir>', 'Output directory for results')
  .option('--organism <species>', 'Organism (mouse, human)', 'mouse')
  .option('--dataset-name <name>', 'Dataset name (default: input directory basename)')
  .option('--single-agent <name>', 'Use only one agent: gpt5.2, claude, or gemini (ICML baseline #2-3)')
  .option('--force-automation', 'Skip all agents, use template-based decisions only (ICML baseline #1)')
  .option('--sequential-chain', 'Use sequential chain mode: GPT-5.2 → Gemini → Claude (ICML experiment #5)')
  .option('--gpt-role <role>', 'Role for GPT-5.2: stats, pipeline, biology (default: stats)', 'stats')
  .option('--claude-role <role>', 'Role for Claude: stats, pipeline, biology (default: pipeline)', 'pipeline')
  .option('--gemini-role <role>', 'Role for Gemini: stats, pipeline, biology (default: biology)', 'biology')
  .option('--auto-resolve <mode>', 'Auto-resolution mode: auto, median, confidence, user (default: auto)', 'auto')
  .option('--escalation-threshold <value>', 'Disagreement score threshold for escalation 0-1 (default: 0.6)', '0.6')
  .option('--verbose', 'Verbose output', false)
  .action(async (input, options) => {
    console.log('🧬 scRNA-seq GeneExpert Multi-Agent Analysis');
    console.log('='.repeat(60));
    console.log('');

    // Validate input directory
    if (!fs.existsSync(input)) {
      console.error(`❌ Error: Input directory not found: ${input}`);
      process.exit(1);
    }

    // Check for data files
    const files = fs.readdirSync(input);
    const hasH5 = files.some(f => f.endsWith('.h5'));
    const hasMatrixDir = fs.existsSync(path.join(input, 'filtered_feature_bc_matrix'));
    const hasCSV = files.some(f => /\.(csv|txt|tsv)$/i.test(f));

    if (!hasH5 && !hasMatrixDir && !hasCSV) {
      console.error('');
      console.error('❌ Error: No scRNA-seq data found in input directory');
      console.error('');
      console.error('Expected one of:');
      console.error('  - filtered_feature_bc_matrix.h5 (10x HDF5 format)');
      console.error('  - filtered_feature_bc_matrix/ directory (10x format) with:');
      console.error('      - barcodes.tsv.gz');
      console.error('      - features.tsv.gz');
      console.error('      - matrix.mtx.gz');
      console.error('  - .csv/.txt/.tsv file (genes as rows, cells as columns)');
      console.error('');
      process.exit(1);
    }

    // Validate role assignments
    const validRoles = ['stats', 'pipeline', 'biology'];
    if (!validRoles.includes(options.gptRole)) {
      console.error(`❌ Error: Invalid --gpt-role: ${options.gptRole}`);
      console.error(`   Valid options: ${validRoles.join(', ')}`);
      process.exit(1);
    }
    if (!validRoles.includes(options.claudeRole)) {
      console.error(`❌ Error: Invalid --claude-role: ${options.claudeRole}`);
      console.error(`   Valid options: ${validRoles.join(', ')}`);
      process.exit(1);
    }
    if (!validRoles.includes(options.geminiRole)) {
      console.error(`❌ Error: Invalid --gemini-role: ${options.geminiRole}`);
      console.error(`   Valid options: ${validRoles.join(', ')}`);
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
      inputDir: path.resolve(input),
      outputDir: outputDir,
      organism: options.organism,
      datasetName: options.datasetName || path.basename(input),
      singleAgent: options.singleAgent,
      forceAutomation: options.forceAutomation,
      sequentialChain: options.sequentialChain,
      roleAssignments: {  // NEW: Role swapping for ablation study
        gptRole: options.gptRole,
        claudeRole: options.claudeRole,
        geminiRole: options.geminiRole
      },
      autoResolveMode: options.autoResolve,
      escalationThreshold: parseFloat(options.escalationThreshold),
      verbose: options.verbose
    };

    // Detect experimental mode
    let experimentalMode = 'Unknown';
    let modeNumber = 0;

    if (options.forceAutomation) {
      modeNumber = 1;
      experimentalMode = 'No-Agent Baseline (template-based decisions)';
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
    } else if (options.sequentialChain) {
      modeNumber = 5;
      experimentalMode = 'Multi-Agent Sequential Chain (GPT-5.2 → Gemini → Claude)';
    } else {
      modeNumber = 4;
      experimentalMode = 'Multi-Agent Parallel (default - independent voting)';
    }

    console.log('📋 Configuration:');
    console.log(`   Input:      ${config.inputDir}`);
    console.log(`   Output:     ${config.outputDir}`);
    console.log(`   Organism:   ${config.organism}`);
    console.log(`   Dataset:    ${config.datasetName}`);
    console.log(`   Data Type:  scRNA-seq (10x Genomics)`);
    console.log('');
    console.log(`🔬 ICML Experiment Mode #${modeNumber}: ${experimentalMode}`);

    // Display role assignments if in multi-agent mode
    if (!options.forceAutomation && !options.singleAgent) {
      const isDefaultRoles = options.gptRole === 'stats' &&
                              options.claudeRole === 'pipeline' &&
                              options.geminiRole === 'biology';
      if (isDefaultRoles) {
        console.log('   Agent Roles: Default (GPT=stats, Claude=pipeline, Gemini=biology)');
      } else {
        console.log('   Agent Roles: Custom Configuration');
        console.log(`   - GPT-5.2:  ${options.gptRole.toUpperCase()} agent`);
        console.log(`   - Claude:   ${options.claudeRole.toUpperCase()} agent`);
        console.log(`   - Gemini:   ${options.geminiRole.toUpperCase()} agent`);
      }
    }
    console.log('');

    try {
      // Execute scRNA-seq analysis
      const executor = new ScRNAExecutor(config);
      await executor.run();

      console.log('');
      console.log('='.repeat(60));
      console.log('✅ Analysis complete!');
      console.log(`📊 Results: ${config.outputDir}`);

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
  .description('Configure scRNA GeneExpert (API keys, default settings)')
  .option('--setup', 'Interactive setup wizard')
  .action((options) => {
    if (options.setup) {
      console.log('🔧 scRNA GeneExpert Configuration Setup');
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
