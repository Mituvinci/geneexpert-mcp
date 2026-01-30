#!/usr/bin/env node

/**
 * Convert Agent Decision JSON to Human-Readable CSV
 *
 * Takes *_agent_decisions.json and creates a flat CSV with all metrics
 * One row per stage, easy to view in Excel/Google Sheets
 *
 * Usage:
 *   node bin/json_to_csv.js --input results/GSE52778_agent_decisions.json --output results/GSE52778_metrics.csv
 *   node bin/json_to_csv.js --dir results/ (converts all JSON files in directory)
 */

import { Command } from 'commander';
import fs from 'fs';
import path from 'path';
import { globSync } from 'glob';

const program = new Command();

program
  .name('json-to-csv')
  .description('Convert agent decision JSON to CSV for easy viewing')
  .version('1.0.0');

program
  .command('convert')
  .description('Convert JSON file(s) to CSV')
  .option('--input <file>', 'Input JSON file')
  .option('--dir <directory>', 'Convert all JSON files in directory')
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
      // Find both .json (bulk RNA-seq) and .jsonl (scRNA-seq) files
      const jsonPattern = path.join(options.dir, '**/*_agent_decisions.json');
      const jsonlPattern = path.join(options.dir, '**/*_agent_decisions.jsonl');
      files = [...globSync(jsonPattern), ...globSync(jsonlPattern)];
    }

    if (files.length === 0) {
      console.log('⚠️  No JSON/JSONL files found.');
      process.exit(0);
    }

    console.log('');
    console.log('='.repeat(70));
    console.log('JSON → CSV CONVERTER');
    console.log('='.repeat(70));
    console.log('');
    console.log(`Found ${files.length} JSON file(s)`);
    console.log('');

    files.forEach(file => {
      try {
        let data;
        const fileContent = fs.readFileSync(file, 'utf-8');

        // Handle both .json (bulk RNA-seq) and .jsonl (scRNA-seq) formats
        if (file.endsWith('.jsonl')) {
          // JSONL format: one JSON object per line
          const lines = fileContent.trim().split('\n').filter(line => line.trim());
          const decisions = lines.map(line => JSON.parse(line));

          // Convert to same structure as .json format
          data = {
            session_id: decisions[0]?.decision_id?.split('_stage')[0] || 'unknown',
            decisions: decisions,
            config: decisions[0]?.stage_input || {}
          };
        } else {
          // JSON format: single object with decisions array
          data = JSON.parse(fileContent);
        }

        const outputFile = options.output || file.replace(/\.jsonl?$/, '_metrics.csv');

        const csv = convertToCSV(data, file);
        fs.writeFileSync(outputFile, csv);

        console.log(`✓ Converted: ${path.basename(file)}`);
        console.log(`  Output: ${outputFile}`);
        console.log('');
      } catch (error) {
        console.error(`❌ Error processing ${file}: ${error.message}`);
      }
    });

    console.log('='.repeat(70));
    console.log('');
  });

program.parse();

/**
 * Convert JSON to flat CSV format
 */
function convertToCSV(data, file) {
  const rows = [];

  // CSV Headers
  const headers = [
    // Session info
    'session_id',
    'dataset',
    'system',
    'organism',
    'comparison',
    'control_keyword',
    'treatment_keyword',
    'duration_seconds',

    // Stage info
    'stage',
    'stage_name',
    'decision_id',
    'timestamp',

    // Stage output
    'stage_status',
    'warnings',
    'errors',

    // Agent responses
    'gpt5_2_decision',
    'gpt5_2_confidence',
    'gpt5_2_reasoning',
    'gpt5_2_cost_usd',
    'gpt5_2_tokens',

    'claude_decision',
    'claude_confidence',
    'claude_reasoning',
    'claude_cost_usd',
    'claude_tokens',

    'gemini_decision',
    'gemini_confidence',
    'gemini_reasoning',
    'gemini_cost_usd',
    'gemini_tokens',

    // Consensus
    'consensus_decision',
    'consensus_confidence_score',
    'consensus_confidence_label',
    'votes_approve',
    'votes_reject',
    'votes_uncertain',
    'voting_method',
    'proceed_to_next_stage',

    // Stage decision
    'final_proceed',
    'user_input_required',
    'user_made_decision',
    'user_decision_reason',
    'samples_removed',
    'de_method',
    'outlier_action',

    // Costs
    'stage_total_cost_usd'
  ];

  rows.push(headers.join(','));

  // Extract session-level info from session_metadata.json if available
  let metadataConfig = data.config || {};

  // Try to load session metadata from same folder
  const dataDir = path.dirname(file);
  const metadataFile = path.join(dataDir, 'staged_analysis_session_metadata.json');

  if (fs.existsSync(metadataFile)) {
    try {
      const metadata = JSON.parse(fs.readFileSync(metadataFile, 'utf8'));
      metadataConfig = metadata.config || metadataConfig;
    } catch (e) {
      // Fallback to data.config
    }
  }

  // Format comparison as "Treatment vs Control"
  const treatmentKeyword = metadataConfig?.treatmentKeyword || 'N/A';
  const controlKeyword = metadataConfig?.controlKeyword || 'N/A';
  const comparisonFormatted = (treatmentKeyword !== 'N/A' && controlKeyword !== 'N/A')
    ? `${treatmentKeyword} vs ${controlKeyword}`
    : (metadataConfig?.comparison || 'N/A');

  const sessionInfo = {
    session_id: data.session_id || 'N/A',
    dataset: extractDatasetName(data),
    system: data.system || 'unknown',
    organism: metadataConfig?.organism || 'N/A',
    comparison: comparisonFormatted,
    control_keyword: controlKeyword,
    treatment_keyword: treatmentKeyword,
    duration_seconds: data.duration_seconds || 0
  };

  // Process each stage decision
  if (data.decisions && Array.isArray(data.decisions)) {
    data.decisions.forEach(decision => {
      const row = [];

      // Session info
      row.push(sessionInfo.session_id);
      row.push(sessionInfo.dataset);
      row.push(sessionInfo.system);
      row.push(sessionInfo.organism);
      row.push(sessionInfo.comparison);
      row.push(sessionInfo.control_keyword);
      row.push(sessionInfo.treatment_keyword);
      row.push(sessionInfo.duration_seconds);

      // Stage info
      row.push(decision.stage || 'N/A');
      row.push(quote(decision.stage_name || 'N/A'));
      row.push(decision.decision_id || 'N/A');
      row.push(decision.timestamp || 'N/A');

      // Stage output
      row.push(decision.stage_output?.overall_status || 'N/A');
      row.push(quote(JSON.stringify(decision.stage_output?.warnings || [])));
      row.push(quote(JSON.stringify(decision.stage_output?.errors || [])));

      // Agent: GPT-5.2
      const gpt = decision.agent_responses?.gpt5_2 || {};
      row.push(gpt.extracted_decision || 'N/A');
      row.push(gpt.confidence_label || 'N/A');
      row.push(quote(truncate(gpt.content, 200)));
      row.push(gpt.cost_usd || 0);
      row.push(gpt.tokens?.total_tokens || 0);

      // Agent: Claude
      const claude = decision.agent_responses?.claude || {};
      row.push(claude.extracted_decision || 'N/A');
      row.push(claude.confidence_label || 'N/A');
      row.push(quote(truncate(claude.content, 200)));
      row.push(claude.cost_usd || 0);
      row.push(claude.tokens?.output_tokens || 0);

      // Agent: Gemini
      const gemini = decision.agent_responses?.gemini || {};
      row.push(gemini.extracted_decision || 'N/A');
      row.push(gemini.confidence_label || 'N/A');
      row.push(quote(truncate(gemini.content, 200)));
      row.push(gemini.cost_usd || 0);
      row.push(gemini.tokens?.totalTokenCount || 0);

      // Consensus
      const consensus = decision.consensus || {};
      row.push(consensus.decision || 'N/A');
      row.push(consensus.confidence_score || 0);
      row.push(consensus.confidence_label || 'N/A');
      row.push(consensus.votes?.approve || 0);
      row.push(consensus.votes?.reject || 0);
      row.push(consensus.votes?.uncertain || 0);
      row.push(consensus.voting_method || 'N/A');
      row.push(consensus.proceed_to_next_stage || false);

      // Stage decision
      const stageDecision = decision.stage_decision || {};
      row.push(stageDecision.proceed || false);
      row.push(stageDecision.user_input_required || false);
      row.push(stageDecision.user_decision?.madeBy || 'N/A');
      row.push(quote(stageDecision.user_decision?.reason || 'N/A'));
      row.push(quote(JSON.stringify(stageDecision.samples_to_remove || [])));
      row.push(stageDecision.de_method || 'N/A');
      row.push(stageDecision.outlier_action || 'N/A');

      // Costs
      row.push(decision.costs?.total_usd || 0);

      rows.push(row.join(','));
    });
  }

  return rows.join('\n');
}

/**
 * Extract dataset name from JSON data
 */
function extractDatasetName(data) {
  if (data.config?.dataset) return data.config.dataset;
  if (data.config?.comparison) return data.config.comparison;
  if (data.config?.input) {
    const parts = data.config.input.split('/');
    return parts[parts.length - 1];
  }
  return 'unknown';
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

/**
 * Truncate long text for CSV readability
 */
function truncate(text, maxLength) {
  if (!text) return 'N/A';
  const str = String(text).replace(/\n/g, ' ');
  if (str.length <= maxLength) return str;
  return str.substring(0, maxLength) + '...';
}
