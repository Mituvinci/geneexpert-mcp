#!/usr/bin/env node

/**
 * Convert scRNA-seq Agent Decision JSONL to CSV
 *
 * Extracts scRNA-specific fields:
 * - Stage 2: nFeature_min, nFeature_max, percent_mt_max
 * - Stage 3A: Cell cycle decision
 * - Stage 4: min_pc, max_pc
 * - Stage 5: Clustering decision
 *
 * Usage:
 *   node bin/json_to_csv_scrna.js convert --dir experiments/scrna_results/
 */

import { Command } from 'commander';
import fs from 'fs';
import path from 'path';
import { globSync } from 'glob';

const program = new Command();

program
  .name('json-to-csv-scrna')
  .description('Convert scRNA-seq agent decision JSONL to CSV')
  .version('1.0.0');

program
  .command('convert')
  .description('Convert JSONL file(s) to CSV')
  .option('--input <file>', 'Input JSONL file')
  .option('--dir <directory>', 'Convert all JSONL files in directory')
  .option('--output <file>', 'Output CSV file (optional, auto-generated if not provided)')
  .action((options) => {
    if (!options.input && !options.dir) {
      console.error('❌ Error: Must provide either --input or --dir');
      process.exit(1);
    }

    let files = [];

    if (options.input) {
      files = [options.input];
    } else if (options.dir) {
      // Find all scRNA JSONL files
      const jsonlPattern = path.join(options.dir, '**/*_analysis_agent_decisions.jsonl');
      files = globSync(jsonlPattern);
    }

    if (files.length === 0) {
      console.log('⚠️  No scRNA JSONL files found.');
      console.log('   Looking for: *_analysis_agent_decisions.jsonl');
      process.exit(0);
    }

    console.log('');
    console.log('='.repeat(70));
    console.log('scRNA-seq JSON → CSV CONVERTER');
    console.log('='.repeat(70));
    console.log('');
    console.log(`Found ${files.length} scRNA JSONL file(s)`);
    console.log('');

    files.forEach(file => {
      try {
        const fileContent = fs.readFileSync(file, 'utf-8');
        const lines = fileContent.trim().split('\n').filter(line => line.trim());
        const decisions = lines.map(line => JSON.parse(line));

        const outputFile = options.output || file.replace('.jsonl', '_metrics.csv');

        const csv = convertToCSV(decisions, file);
        fs.writeFileSync(outputFile, csv);

        console.log(`✓ Converted: ${path.basename(file)}`);
        console.log(`  Output: ${outputFile}`);
        console.log(`  Stages: ${decisions.length}`);
        console.log('');
      } catch (error) {
        console.error(`❌ Error processing ${file}: ${error.message}`);
        console.error(error.stack);
      }
    });

    console.log('='.repeat(70));
    console.log('');
  });

program.parse();

/**
 * Extract QC thresholds from agent content (Stage 2 specific)
 */
function extractThresholds(content) {
  if (!content) return { nFeature_min: 'N/A', nFeature_max: 'N/A', percent_mt_max: 'N/A' };

  const nFeatureMinMatch = content.match(/nFeature_min:\s*(\d+)/i);
  const nFeatureMaxMatch = content.match(/nFeature_max:\s*(\d+)/i);
  const percentMtMaxMatch = content.match(/percent_mt_max:\s*(\d+)/i);

  return {
    nFeature_min: nFeatureMinMatch ? nFeatureMinMatch[1] : 'N/A',
    nFeature_max: nFeatureMaxMatch ? nFeatureMaxMatch[1] : 'N/A',
    percent_mt_max: percentMtMaxMatch ? percentMtMaxMatch[1] : 'N/A'
  };
}

/**
 * Extract PC range from agent content (Stage 4 specific)
 */
function extractPCRange(content) {
  if (!content) return { min_pc: 'N/A', max_pc: 'N/A' };

  const minPcMatch = content.match(/min_pc:\s*(\d+)/i);
  const maxPcMatch = content.match(/max_pc:\s*(\d+)/i);

  // Also try patterns like "PCs 1-30" or "use 1-25"
  if (!maxPcMatch) {
    const rangeMatch = content.match(/(?:PCs?|use|dims?)\s*(?:1-)?(\d+)/i);
    if (rangeMatch) {
      return {
        min_pc: '1',
        max_pc: rangeMatch[1]
      };
    }
  }

  return {
    min_pc: minPcMatch ? minPcMatch[1] : '1',
    max_pc: maxPcMatch ? maxPcMatch[1] : 'N/A'
  };
}

/**
 * Convert JSONL decisions to flat CSV format
 */
function convertToCSV(decisions, file) {
  const rows = [];

  // Extract session-level info
  const firstDecision = decisions[0] || {};
  const dataDir = path.dirname(file);
  const metadataFile = path.join(dataDir, 'scrna_analysis_session_metadata.json');

  let metadataConfig = {};
  if (fs.existsSync(metadataFile)) {
    try {
      const metadata = JSON.parse(fs.readFileSync(metadataFile, 'utf8'));
      metadataConfig = metadata.config || {};
    } catch (e) {
      // Fallback
    }
  }

  const datasetName = extractDatasetName(firstDecision, metadataConfig, file);
  const organism = metadataConfig?.organism || 'N/A';
  const system = extractSystemConfig(file, metadataConfig);

  // CSV Headers
  const headers = [
    // Session info
    'dataset_name',
    'full_dataset',
    'system_config',
    'organism',

    // Stage info
    'stage',
    'stage_name',
    'decision_id',
    'timestamp',

    // Stage output
    'stage_status',

    // Agent responses - GPT-5.2
    'gpt5_2_decision',
    'gpt5_2_confidence',
    'gpt5_2_nFeature_min',  // Stage 2
    'gpt5_2_nFeature_max',  // Stage 2
    'gpt5_2_percent_mt_max', // Stage 2
    'gpt5_2_min_pc',  // Stage 4
    'gpt5_2_max_pc',  // Stage 4
    'gpt5_2_cost_usd',

    // Agent responses - Claude
    'claude_decision',
    'claude_confidence',
    'claude_nFeature_min',  // Stage 2
    'claude_nFeature_max',  // Stage 2
    'claude_percent_mt_max', // Stage 2
    'claude_min_pc',  // Stage 4
    'claude_max_pc',  // Stage 4
    'claude_cost_usd',

    // Agent responses - Gemini
    'gemini_decision',
    'gemini_confidence',
    'gemini_nFeature_min',  // Stage 2
    'gemini_nFeature_max',  // Stage 2
    'gemini_percent_mt_max', // Stage 2
    'gemini_min_pc',  // Stage 4
    'gemini_max_pc',  // Stage 4
    'gemini_cost_usd',

    // Consensus
    'consensus_decision',
    'consensus_confidence_score',
    'votes_approve',
    'votes_reject',
    'votes_uncertain',
    'proceed_to_next_stage',

    // Stage decision
    'final_proceed',
    'user_input_required',

    // Costs
    'stage_total_cost_usd'
  ];

  rows.push(headers.join(','));

  // Process each stage decision
  decisions.forEach(decision => {
    const row = [];

    // Session info
    row.push(quote(datasetName));
    row.push(quote(datasetName));
    row.push(quote(system));
    row.push(organism);

    // Stage info
    row.push(decision.stage || 'N/A');
    row.push(quote(decision.stage_name || 'N/A'));
    row.push(decision.decision_id || 'N/A');
    row.push(decision.timestamp || 'N/A');

    // Stage output
    row.push(decision.stage_output?.overall_status || 'N/A');

    // Determine if this is Stage 2 (QC Filtering) or Stage 4 (PCA)
    const isStage2 = decision.stage === 2 || decision.stage_name?.includes('QC Filtering');
    const isStage4 = decision.stage === 4 || decision.stage_name?.includes('PCA');

    // Agent: GPT-5.2
    const gpt = decision.agent_responses?.gpt5_2 || {};
    row.push(gpt.extracted_decision || 'N/A');
    row.push(gpt.confidence_label || 'N/A');
    if (isStage2) {
      const thresholds = extractThresholds(gpt.content);
      row.push(thresholds.nFeature_min);
      row.push(thresholds.nFeature_max);
      row.push(thresholds.percent_mt_max);
    } else {
      row.push('N/A');
      row.push('N/A');
      row.push('N/A');
    }
    if (isStage4) {
      const pcRange = extractPCRange(gpt.content);
      row.push(pcRange.min_pc);
      row.push(pcRange.max_pc);
    } else {
      row.push('N/A');
      row.push('N/A');
    }
    row.push(gpt.cost_usd || 0);

    // Agent: Claude
    const claude = decision.agent_responses?.claude || {};
    row.push(claude.extracted_decision || 'N/A');
    row.push(claude.confidence_label || 'N/A');
    if (isStage2) {
      const thresholds = extractThresholds(claude.content);
      row.push(thresholds.nFeature_min);
      row.push(thresholds.nFeature_max);
      row.push(thresholds.percent_mt_max);
    } else {
      row.push('N/A');
      row.push('N/A');
      row.push('N/A');
    }
    if (isStage4) {
      const pcRange = extractPCRange(claude.content);
      row.push(pcRange.min_pc);
      row.push(pcRange.max_pc);
    } else {
      row.push('N/A');
      row.push('N/A');
    }
    row.push(claude.cost_usd || 0);

    // Agent: Gemini
    const gemini = decision.agent_responses?.gemini || {};
    row.push(gemini.extracted_decision || 'N/A');
    row.push(gemini.confidence_label || 'N/A');
    if (isStage2) {
      const thresholds = extractThresholds(gemini.content);
      row.push(thresholds.nFeature_min);
      row.push(thresholds.nFeature_max);
      row.push(thresholds.percent_mt_max);
    } else {
      row.push('N/A');
      row.push('N/A');
      row.push('N/A');
    }
    if (isStage4) {
      const pcRange = extractPCRange(gemini.content);
      row.push(pcRange.min_pc);
      row.push(pcRange.max_pc);
    } else {
      row.push('N/A');
      row.push('N/A');
    }
    row.push(gemini.cost_usd || 0);

    // Consensus
    const consensus = decision.consensus || {};
    row.push(consensus.decision || 'N/A');
    row.push(consensus.confidence_score || 0);
    row.push(consensus.votes?.approve || 0);
    row.push(consensus.votes?.reject || 0);
    row.push(consensus.votes?.uncertain || 0);
    row.push(consensus.proceed_to_next_stage || false);

    // Stage decision
    const stageDecision = decision.stage_decision || {};
    row.push(stageDecision.proceed || false);
    row.push(stageDecision.user_input_required || false);

    // Costs
    row.push(decision.costs?.total_usd || 0);

    rows.push(row.join(','));
  });

  return rows.join('\n');
}

/**
 * Extract dataset name
 */
function extractDatasetName(decision, metadataConfig, file) {
  // Try decision_id first
  if (decision.decision_id) {
    const parts = decision.decision_id.split('_stage');
    if (parts[0]) return parts[0];
  }

  // Try metadata
  if (metadataConfig?.datasetName) return metadataConfig.datasetName;

  // Try folder name
  const folderName = path.basename(path.dirname(file));
  if (folderName && folderName !== 'scrna_results') {
    // Remove system suffix like _parallel_gp_bl_cl_pl_gm_st
    return folderName.replace(/_(parallel|sequential|single_claude|single_gpt|single_gemini|no_agent).*$/, '');
  }

  return 'unknown';
}

/**
 * Extract system config from folder name or metadata
 */
function extractSystemConfig(file, metadataConfig) {
  const folderName = path.basename(path.dirname(file));

  // Detect system from folder name
  if (folderName.includes('_parallel_')) {
    // Role swapping system: extract roles from folder name
    // Format: {dataset}_parallel_{gpt_role}_{claude_role}_{gemini_role}
    // Example: 10-k-brain-cells_healthy_mouse_parallel_gp_bl_cl_pl_gm_st
    const match = folderName.match(/_parallel_([a-z]{2})_([a-z]{2})_([a-z]{2})_([a-z]{2})_([a-z]{2})_([a-z]{2})$/);
    if (match) {
      const roleMap = { gp: 'gpt', cl: 'claude', gm: 'gemini', st: 'stats', pl: 'pipeline', bl: 'biology' };
      const gptRole = roleMap[match[1]] + '_' + roleMap[match[2]];
      const claudeRole = roleMap[match[3]] + '_' + roleMap[match[4]];
      const geminiRole = roleMap[match[5]] + '_' + roleMap[match[6]];
      return `parallel_${gptRole}_${claudeRole}_${geminiRole}`;
    }
    return 'parallel_default';
  } else if (folderName.includes('_sequential_')) {
    return 'sequential_chain';
  } else if (folderName.includes('_single_claude')) {
    return 'single_agent_claude';
  } else if (folderName.includes('_single_gpt')) {
    return 'single_agent_gpt';
  } else if (folderName.includes('_single_gemini')) {
    return 'single_agent_gemini';
  } else if (folderName.includes('_no_agent')) {
    return 'no_agent';
  }

  // Fallback to metadata
  if (metadataConfig?.singleAgent) {
    return `single_agent_${metadataConfig.singleAgent}`;
  }
  if (metadataConfig?.sequentialChain) {
    return 'sequential_chain';
  }
  if (metadataConfig?.forceAutomation) {
    return 'no_agent';
  }

  return 'parallel_default';
}

/**
 * Quote a value for CSV (handles commas, quotes, newlines)
 */
function quote(value) {
  if (value === null || value === undefined) return 'N/A';
  const str = String(value);
  // If contains comma, quote, or newline, wrap in quotes and escape quotes
  if (str.includes(',') || str.includes('"') || str.includes('\n')) {
    return '"' + str.replace(/"/g, '""') + '"';
  }
  return str;
}
