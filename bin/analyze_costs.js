#!/usr/bin/env node

/**
 * Cost Analysis Script for ICML Experiments
 *
 * Extracts cost data from agent decision logs and generates CSV report
 *
 * Usage:
 *   node bin/analyze_costs.js --results experiments/results --output costs.csv
 *   node bin/analyze_costs.js --results experiments/results --system multi-agent-parallel
 */

import { Command } from 'commander';
import fs from 'fs';
import path from 'path';
import { globSync } from 'glob';

const program = new Command();

program
  .name('analyze-costs')
  .description('Generate cost analysis CSV from experiment logs')
  .version('1.0.0');

program
  .command('generate')
  .description('Generate cost CSV report')
  .option('--results <dir>', 'Results directory', 'experiments/results')
  .option('--output <file>', 'Output CSV file', 'cost_analysis.csv')
  .option('--system <name>', 'Filter by system (optional)')
  .action((options) => {
    console.log('');
    console.log('='.repeat(70));
    console.log('COST ANALYSIS - ICML 2026 EXPERIMENTS');
    console.log('='.repeat(70));
    console.log('');

    // Find all agent decision JSON files
    const pattern = path.join(options.results, '**/*_agent_decisions.json');
    const files = globSync(pattern);

    if (files.length === 0) {
      console.log('⚠️  No experiment files found.');
      console.log(`   Searched: ${pattern}`);
      console.log('');
      process.exit(0);
    }

    console.log(`Found ${files.length} experiment log(s)`);
    console.log('');

    // Extract cost data from all files
    const allCosts = [];

    files.forEach(file => {
      try {
        const data = JSON.parse(fs.readFileSync(file, 'utf-8'));
        const costs = extractCostData(data, file, options.system);
        allCosts.push(...costs);
      } catch (error) {
        console.warn(`⚠️  Error reading ${file}: ${error.message}`);
      }
    });

    if (allCosts.length === 0) {
      console.log('⚠️  No cost data extracted from experiments.');
      console.log('');
      process.exit(0);
    }

    console.log(`Extracted cost data: ${allCosts.length} stage(s)`);
    console.log('');

    // Generate CSV
    const csv = generateCostCSV(allCosts);
    fs.writeFileSync(options.output, csv);

    console.log(`✓ Cost analysis saved: ${options.output}`);
    console.log('');

    // Display summary
    displayCostSummary(allCosts);

    console.log('='.repeat(70));
    console.log('');
  });

program.parse();

/**
 * Extract cost data from an experiment log
 */
function extractCostData(data, filepath, systemFilter) {
  const costs = [];

  // Detect system type
  let system = data.system || 'unknown';

  // More specific system detection
  if (system === 'multi-agent') {
    // Check if sequential
    if (data.config?.sequentialChain || data.sequential) {
      system = 'multi-agent-sequential';
    } else {
      system = 'multi-agent-parallel';
    }
  } else if (system === 'single-agent') {
    // Detect which single agent
    if (data.config?.singleAgent) {
      system = `single-agent-${data.config.singleAgent}`;
    }
  } else if (data.config?.forceAutomation) {
    system = 'no-agent';
  }

  // Apply filter
  if (systemFilter && system !== systemFilter) {
    return costs;
  }

  // Extract dataset name from config or filepath
  const dataset = extractDatasetName(data, filepath);

  // Process each decision (stage)
  if (data.decisions && Array.isArray(data.decisions)) {
    data.decisions.forEach(decision => {
      const stageCost = {
        dataset,
        system,
        stage: decision.stage,
        stage_name: decision.stage_name || `Stage ${decision.stage}`,

        // Individual agent costs
        gpt5_2_cost: 0,
        claude_cost: 0,
        gemini_cost: 0,

        // Stage totals
        stage_total_cost: decision.costs?.total_usd || 0,

        // Consensus metrics
        decision: decision.consensus?.decision || decision.stage_decision?.proceed ? 'proceed' : 'abort',
        confidence: decision.consensus?.confidence_score || 0,
        user_input_required: decision.stage_decision?.user_input_required || false,

        // Vote breakdown (multi-agent only)
        votes_approve: decision.consensus?.votes?.approve || 0,
        votes_reject: decision.consensus?.votes?.reject || 0,
        votes_uncertain: decision.consensus?.votes?.uncertain || 0,

        // Timestamp
        timestamp: decision.timestamp || data.timestamp_start
      };

      // Extract individual agent costs
      if (decision.costs?.breakdown) {
        stageCost.gpt5_2_cost = decision.costs.breakdown.gpt5_2 || 0;
        stageCost.claude_cost = decision.costs.breakdown.claude || 0;
        stageCost.gemini_cost = decision.costs.breakdown.gemini || 0;
      }

      costs.push(stageCost);
    });
  }

  // Add session-level summary
  if (data.costs?.total_usd) {
    costs.push({
      dataset,
      system,
      stage: 'TOTAL',
      stage_name: 'Session Total',
      gpt5_2_cost: data.costs.breakdown?.gpt5_2 || 0,
      claude_cost: data.costs.breakdown?.claude || 0,
      gemini_cost: data.costs.breakdown?.gemini || 0,
      stage_total_cost: data.costs.total_usd,
      decision: 'N/A',
      confidence: 0,
      user_input_required: false,
      votes_approve: 0,
      votes_reject: 0,
      votes_uncertain: 0,
      timestamp: data.timestamp_end || data.timestamp_start
    });
  }

  return costs;
}

/**
 * Extract dataset name from log data or filepath
 */
function extractDatasetName(data, filepath) {
  // Try config first
  if (data.config?.dataset) {
    return data.config.dataset;
  }

  // Try input directory
  if (data.config?.input) {
    const parts = data.config.input.split('/');
    return parts[parts.length - 1];
  }

  // Try comparison name
  if (data.config?.comparison) {
    return data.config.comparison;
  }

  // Fallback to filename parsing
  const filename = path.basename(filepath, '.json');
  const match = filename.match(/^(.+?)_/);
  return match ? match[1] : 'unknown';
}

/**
 * Generate CSV from cost data
 */
function generateCostCSV(allCosts) {
  const headers = [
    'dataset',
    'system',
    'stage',
    'stage_name',
    'gpt5_2_cost_usd',
    'claude_cost_usd',
    'gemini_cost_usd',
    'stage_total_cost_usd',
    'decision',
    'confidence',
    'user_input_required',
    'votes_approve',
    'votes_reject',
    'votes_uncertain',
    'timestamp'
  ];

  const rows = [headers.join(',')];

  allCosts.forEach(cost => {
    rows.push([
      cost.dataset,
      cost.system,
      cost.stage,
      `"${cost.stage_name}"`, // Quote in case of commas
      cost.gpt5_2_cost.toFixed(6),
      cost.claude_cost.toFixed(6),
      cost.gemini_cost.toFixed(6),
      cost.stage_total_cost.toFixed(6),
      cost.decision,
      cost.confidence.toFixed(2),
      cost.user_input_required,
      cost.votes_approve,
      cost.votes_reject,
      cost.votes_uncertain,
      cost.timestamp
    ].join(','));
  });

  return rows.join('\n');
}

/**
 * Display cost summary to console
 */
function displayCostSummary(allCosts) {
  console.log('COST SUMMARY');
  console.log('-'.repeat(70));

  // Group by system
  const bySystem = {};
  allCosts.forEach(cost => {
    if (!bySystem[cost.system]) {
      bySystem[cost.system] = {
        total: 0,
        gpt5_2: 0,
        claude: 0,
        gemini: 0,
        stages: 0
      };
    }

    // Only count non-TOTAL rows for stage count
    if (cost.stage !== 'TOTAL') {
      bySystem[cost.system].stages++;
    }

    bySystem[cost.system].total += cost.stage_total_cost;
    bySystem[cost.system].gpt5_2 += cost.gpt5_2_cost;
    bySystem[cost.system].claude += cost.claude_cost;
    bySystem[cost.system].gemini += cost.gemini_cost;
  });

  Object.entries(bySystem).forEach(([system, costs]) => {
    console.log(`\n${system.toUpperCase()}`);
    console.log(`  Stages:      ${costs.stages}`);
    console.log(`  Total Cost:  $${costs.total.toFixed(6)}`);
    console.log(`  GPT-5.2:     $${costs.gpt5_2.toFixed(6)}`);
    console.log(`  Claude:      $${costs.claude.toFixed(6)}`);
    console.log(`  Gemini:      $${costs.gemini.toFixed(6)}`);

    if (costs.stages > 0) {
      console.log(`  Avg/Stage:   $${(costs.total / costs.stages).toFixed(6)}`);
    }
  });

  console.log('');

  // Overall total
  const grandTotal = allCosts.reduce((sum, cost) => sum + cost.stage_total_cost, 0);
  console.log(`GRAND TOTAL: $${grandTotal.toFixed(6)}`);
  console.log('');
}
