#!/usr/bin/env node

/**
 * scRNA-seq Multi-Agent Evaluation Script
 *
 * Compares agent decisions to ground truth across all experimental systems.
 * Calculates decision accuracy, user intervention rates, and system performance.
 *
 * Adapted from bulk RNA-seq version for scRNA-seq (5 stages).
 *
 * Usage:
 *   node bin/evaluate_scrna.js <input_csv> <ground_truth_json> <output_csv>
 *
 * Example:
 *   node bin/evaluate_scrna.js \
 *     experiments/scrna_results/scrna_ALL_EXPERIMENTS_DETAILED.csv \
 *     experiments/scrna_ground_truth.json \
 *     scrna_evaluation.csv
 */

import fs from 'fs';
import path from 'path';

// ========================================
// Configuration
// ========================================

const INPUT_CSV = process.argv[2] || 'experiments/scrna_results/scrna_ALL_EXPERIMENTS_DETAILED_FIXED.csv';
const GROUND_TRUTH_FILE = process.argv[3] || 'experiments/scrna_ground_truth.json';
const OUTPUT_CSV = process.argv[4] || 'experiments/scrna_csv_figures/scrna_evaluation.csv';

// Dataset mapping (short name in CSV -> ground truth key)
const DATASET_MAP = {
  '1_GD428_21136_Hu_REH_Parental': '1_GD428_21136_Hu_REH_Parental',
  '2_GD444_21136_Hu_SUP_Parental': '2_GD444_21136_Hu_SUP_Parental',
  'pbmc_healthy_human': 'pbmc_healthy_human',
  '10-k-brain-cells_healthy_mouse': '10-k-brain-cells_healthy_mouse',
  'GSE75748': 'GSE75748',
  'GSE146773': 'GSE146773',
  '3_GSE64016_H1andFUCCI_normalized_EC_original': '3_GSE64016_H1andFUCCI_normalized_EC_original'
};

// ========================================
// Load Ground Truth
// ========================================

function loadGroundTruth(filepath) {
  console.log(`\n📖 Loading ground truth from: ${filepath}`);
  const data = JSON.parse(fs.readFileSync(filepath, 'utf8'));
  console.log(`   ✓ Loaded ${data.total_datasets} datasets with ${data.total_decisions} decisions\n`);
  return data.ground_truth;
}

// ========================================
// Load CSV Data
// ========================================

function loadCSV(filepath) {
  console.log(`📊 Loading CSV data from: ${filepath}\n`);

  if (!fs.existsSync(filepath)) {
    console.error(`ERROR: CSV file not found: ${filepath}`);
    process.exit(1);
  }

  const content = fs.readFileSync(filepath, 'utf8');
  const lines = content.trim().split('\n');
  const headers = lines[0].split(',');

  const rows = [];
  for (let i = 1; i < lines.length; i++) {
    const values = parseCSVLine(lines[i]);
    const row = {};

    for (let j = 0; j < headers.length; j++) {
      row[headers[j]] = values[j] || '';
    }

    rows.push(row);
  }

  console.log(`   ✓ Loaded ${rows.length} rows\n`);
  return rows;
}

// Parse CSV line handling quoted fields
function parseCSVLine(line) {
  const values = [];
  let current = '';
  let inQuotes = false;

  for (let i = 0; i < line.length; i++) {
    const char = line[i];
    const nextChar = line[i + 1];

    if (char === '"') {
      if (inQuotes && nextChar === '"') {
        current += '"';
        i++; // Skip next quote
      } else {
        inQuotes = !inQuotes;
      }
    } else if (char === ',' && !inQuotes) {
      values.push(current);
      current = '';
    } else {
      current += char;
    }
  }

  values.push(current);
  return values;
}

// ========================================
// Group CSV Rows by Experiment
// ========================================

function groupByExperiment(rows) {
  const experiments = {};

  for (const row of rows) {
    const expName = row.full_dataset || row.dataset_name;
    const system = row.system_config;
    const key = `${expName}_${system}`;

    if (!experiments[key]) {
      experiments[key] = {
        dataset: expName,
        system: system,
        organism: row.organism,
        stages: []
      };
    }

    experiments[key].stages.push(row);
  }

  console.log(`   ✓ Grouped into ${Object.keys(experiments).length} experiments\n`);
  return experiments;
}

// ========================================
// Extract Agent Decisions from CSV Row
// ========================================

function extractDecisions(row) {
  // Get individual agent decisions
  const gpt = row.gpt5_2_decision || 'N/A';
  const claude = row.claude_decision || 'N/A';
  const gemini = row.gemini_decision || 'N/A';

  // Calculate majority vote (most common decision)
  // NOTE: consensus_decision contains vote result (approve/reject), NOT the actual decision
  const decisions = [gpt, claude, gemini].filter(d => d !== 'N/A' && d);
  let finalDecision = 'N/A';

  if (decisions.length > 0) {
    const counts = {};
    for (const d of decisions) {
      counts[d] = (counts[d] || 0) + 1;
    }
    // Get most common decision (majority)
    finalDecision = Object.keys(counts).reduce((a, b) => counts[a] > counts[b] ? a : b);
  }

  return {
    decision: finalDecision,
    gpt_decision: gpt,
    claude_decision: claude,
    gemini_decision: gemini,
    user_intervention: row.user_input_required === 'true' || row.user_input_required === '1'
  };
}

// ========================================
// Compare Decisions
// ========================================

function compareDecision(agentDecision, groundTruth, stageKey) {
  const comparison = {
    stage: stageKey,
    match: false,
    agent_decision: agentDecision,
    correct_decision: groundTruth ? groundTruth.correct_decision : null,
    user_intervention: false,
    rationale: groundTruth ? groundTruth.rationale : 'No ground truth'
  };

  if (!groundTruth) {
    return comparison;
  }

  // Normalize and compare
  const normalized_agent = normalizeDecision(agentDecision);
  const normalized_correct = normalizeDecision(groundTruth.correct_decision);

  comparison.match = normalized_agent === normalized_correct;

  return comparison;
}

// ========================================
// Normalize Decisions (with scRNA-specific rules)
// ========================================

function normalizeDecision(decision) {
  if (!decision || decision === 'N/A') return 'unknown';

  let normalized = decision.toString().toLowerCase().trim();

  // scRNA-specific normalization
  if (normalized.includes('skip') || normalized === 'no_cell_cycle' || normalized === 'nocellcycle') {
    return 'skip_cell_cycle';
  }

  if (normalized.includes('remove') && normalized.includes('cell')) {
    return 'remove_cell_cycle';
  }

  if (normalized === 'use_default' || normalized === 'usedefault' || normalized === 'default') {
    return 'use_default';
  }

  if (normalized === 'select_pc_range' || normalized.includes('select')) {
    return 'select_pc_range';
  }

  if (normalized === 'stop_and_review' || normalized.includes('stop')) {
    return 'stop_and_review';
  }

  if (normalized === 'accept_clustering' || normalized === 'accept' || normalized === 'current') {
    return 'accept_clustering';
  }

  if (normalized === 'adjust_resolution' || normalized.includes('adjust')) {
    return 'adjust_resolution';
  }

  if (normalized === 'flag_suspicious' || normalized.includes('flag')) {
    return 'flag_suspicious';
  }

  if (normalized === 'set_thresholds' || normalized.includes('threshold')) {
    return 'set_thresholds';
  }

  if (normalized === 'auto_proceed' || normalized === 'autoproceed') {
    return 'auto_proceed';
  }

  // Remove special characters and spaces
  return normalized.replace(/[_\s-]/g, '').replace(/\*/g, '');
}

// ========================================
// Evaluate Single Experiment
// ========================================

function evaluateExperiment(experiment, groundTruth) {
  const datasetKey = DATASET_MAP[experiment.dataset] || experiment.dataset;
  const gtDataset = groundTruth[datasetKey];

  const result = {
    experiment: `${experiment.dataset}_${experiment.system}`,
    dataset: experiment.dataset,
    system: experiment.system,
    organism: experiment.organism,
    total_stages: 0,
    completed_stages: 0,
    correct_decisions: 0,
    user_interventions: 0,
    auto_decisions: 0,
    auto_correct: 0,
    justified_escalations: 0,
    total_cost_usd: 0,
    stage_results: {},
    completed: false,
    reached_stage_5: false,
    accuracy: 0,
    auto_accuracy: 0,
    safe_escalation_rate: 0
  };

  if (!gtDataset) {
    result.error = `No ground truth for dataset: ${datasetKey}`;
    return result;
  }

  const gtStages = gtDataset.stages;

  // Calculate total cost from ALL stages (including reanalysis)
  for (const row of experiment.stages) {
    const stageCost = parseFloat(row.stage_total_cost_usd) || 0;
    result.total_cost_usd += stageCost;
  }

  // Handle reanalysis: Keep only LAST occurrence of each stage
  const stageMap = new Map();
  for (const row of experiment.stages) {
    const stageName = row.stage_name;
    const stageNum = row.stage;

    // Map stage number to ground truth key
    let gtStageKey = stageNum;
    if (stageName.includes('Cell Cycle Scoring')) {
      gtStageKey = '3A';
    } else if (stageName.includes('Skip Regression') || stageName.includes('Cell Cycle Regression') || stageName.includes('Execute')) {
      gtStageKey = '3B';
    }

    // Store/overwrite with latest occurrence (handles reanalysis)
    stageMap.set(gtStageKey, row);
  }

  // Process unique stages only (reanalysis handled)
  for (const [gtStageKey, row] of stageMap.entries()) {
    const gtStage = gtStages[gtStageKey.toString()];

    if (!gtStage) {
      continue; // Skip stages without ground truth
    }

    result.total_stages++;

    // Extract agent decisions (majority vote calculated inside)
    const decisions = extractDecisions(row);

    // Use the calculated majority decision
    const agentDecision = decisions.decision;

    // Compare to ground truth
    const comparison = compareDecision(agentDecision, gtStage, gtStageKey);
    comparison.user_intervention = decisions.user_intervention;

    result.stage_results[gtStageKey] = comparison;

    // Track overall correctness
    if (comparison.match) {
      result.correct_decisions++;
    }

    // Separate auto-decisions from user-escalated decisions
    if (decisions.user_intervention) {
      result.user_interventions++;

      // Justified escalation: User intervention led to correct decision
      if (comparison.match) {
        result.justified_escalations++;
      }
    } else {
      // Auto-decision (no user intervention)
      result.auto_decisions++;

      if (comparison.match) {
        result.auto_correct++;
      }
    }

    // Check if stage was completed (PASS or SUCCESS)
    if (row.stage_status && (row.stage_status === 'PASS' || row.stage_status === 'SUCCESS')) {
      result.completed_stages++;
    }

    // Check if this is stage 5 (Clustering)
    if (gtStageKey === '5' || row.stage === 5 || row.stage === '5') {
      result.reached_stage_5 = true;
    }
  }

  // Calculate overall accuracy
  result.accuracy = result.total_stages > 0
    ? (result.correct_decisions / result.total_stages * 100).toFixed(1)
    : 0;

  // Calculate auto-decision accuracy (decisions made without user intervention)
  result.auto_accuracy = result.auto_decisions > 0
    ? (result.auto_correct / result.auto_decisions * 100).toFixed(1)
    : 0;

  // Calculate safe escalation rate (% of user interventions that were justified)
  result.safe_escalation_rate = result.user_interventions > 0
    ? (result.justified_escalations / result.user_interventions * 100).toFixed(1)
    : 0;

  // Completion: reached stage 5 (regardless of pass/fail status)
  // User requirement: "Any analysis that went till stage 5, that is actually successful"
  result.completed = result.reached_stage_5;

  return result;
}

// ========================================
// Aggregate Results by System
// ========================================

function aggregateBySystem(results) {
  const systemStats = {};

  for (const result of results) {
    const sys = result.system;

    if (!systemStats[sys]) {
      systemStats[sys] = {
        system: sys,
        total_experiments: 0,
        completed_experiments: 0,
        total_decisions: 0,
        correct_decisions: 0,
        total_user_interventions: 0,
        accuracy: 0,
        completion_rate: 0,
        datasets: {}
      };
    }

    const stats = systemStats[sys];
    stats.total_experiments++;

    if (result.completed) {
      stats.completed_experiments++;
    }

    stats.total_decisions += result.total_stages;
    stats.correct_decisions += result.correct_decisions;
    stats.total_user_interventions += result.user_interventions;

    // Track per-dataset results
    if (!stats.datasets[result.dataset]) {
      stats.datasets[result.dataset] = {
        completed: false,
        accuracy: 0
      };
    }
    stats.datasets[result.dataset].completed = result.completed;
    stats.datasets[result.dataset].accuracy = parseFloat(result.accuracy);
  }

  // Calculate percentages
  for (const sys in systemStats) {
    const stats = systemStats[sys];
    stats.accuracy = stats.total_decisions > 0
      ? (stats.correct_decisions / stats.total_decisions * 100).toFixed(1)
      : 0;
    stats.completion_rate = stats.total_experiments > 0
      ? (stats.completed_experiments / stats.total_experiments * 100).toFixed(1)
      : 0;
  }

  return systemStats;
}

// ========================================
// Export Per-Experiment CSV
// ========================================

function exportPerExperimentCSV(results, outputFile) {
  let csv = 'Experiment,System,Dataset,Organism,Completed,Total_Stages,Correct_Decisions,Accuracy_%,Auto_Decisions,Auto_Correct,Auto_Accuracy_%,User_Interventions,Justified_Escalations,Safe_Escalation_Rate_%,Total_Cost_USD,Stage1_Match,Stage2_Match,Stage3A_Match,Stage3B_Match,Stage4_Match,Stage5_Match,Stage1_UserInput,Stage2_UserInput,Stage3A_UserInput,Stage3B_UserInput,Stage4_UserInput,Stage5_UserInput\n';

  for (const result of results) {
    // Decision match (1/0)
    const stage1 = result.stage_results['1'] ? (result.stage_results['1'].match ? 1 : 0) : '';
    const stage2 = result.stage_results['2'] ? (result.stage_results['2'].match ? 1 : 0) : '';
    const stage3A = result.stage_results['3A'] ? (result.stage_results['3A'].match ? 1 : 0) : '';
    const stage3B = result.stage_results['3B'] ? (result.stage_results['3B'].match ? 1 : 0) : '';
    const stage4 = result.stage_results['4'] ? (result.stage_results['4'].match ? 1 : 0) : '';
    const stage5 = result.stage_results['5'] ? (result.stage_results['5'].match ? 1 : 0) : '';

    // User input required (1/0)
    const stage1_user = result.stage_results['1'] ? (result.stage_results['1'].user_intervention ? 1 : 0) : '';
    const stage2_user = result.stage_results['2'] ? (result.stage_results['2'].user_intervention ? 1 : 0) : '';
    const stage3A_user = result.stage_results['3A'] ? (result.stage_results['3A'].user_intervention ? 1 : 0) : '';
    const stage3B_user = result.stage_results['3B'] ? (result.stage_results['3B'].user_intervention ? 1 : 0) : '';
    const stage4_user = result.stage_results['4'] ? (result.stage_results['4'].user_intervention ? 1 : 0) : '';
    const stage5_user = result.stage_results['5'] ? (result.stage_results['5'].user_intervention ? 1 : 0) : '';

    csv += `${result.experiment},${result.system},${result.dataset},${result.organism},${result.completed ? 1 : 0},`;
    csv += `${result.total_stages},${result.correct_decisions},${result.accuracy},`;
    csv += `${result.auto_decisions},${result.auto_correct},${result.auto_accuracy},`;
    csv += `${result.user_interventions},${result.justified_escalations},${result.safe_escalation_rate},`;
    csv += `${result.total_cost_usd},`;
    csv += `${stage1},${stage2},${stage3A},${stage3B},${stage4},${stage5},`;
    csv += `${stage1_user},${stage2_user},${stage3A_user},${stage3B_user},${stage4_user},${stage5_user}\n`;
  }

  // Ensure output directory exists
  const outputDir = path.dirname(outputFile);
  if (!fs.existsSync(outputDir)) {
    fs.mkdirSync(outputDir, { recursive: true });
    console.log(`📁 Created output directory: ${outputDir}`);
  }

  fs.writeFileSync(outputFile, csv);
  console.log(`✅ Per-experiment results exported to: ${outputFile}\n`);
}

// ========================================
// Export Aggregated CSV
// ========================================

function exportAggregatedCSV(systemStats, outputFile) {
  const datasets = Object.keys(DATASET_MAP);

  let csv = 'System,Total_Experiments,Completed,Completion_Rate_%,Total_Decisions,Correct_Decisions,Decision_Accuracy_%,User_Interventions';

  // Add dataset columns
  for (const ds of datasets) {
    const shortName = ds.split('_')[0] || ds; // 1, 2, pbmc, 10-k, GSE75748, etc.
    csv += `,${shortName}_Accuracy`;
  }
  csv += '\n';

  // Sort systems for consistent output
  const sortedSystems = Object.keys(systemStats).sort();

  for (const sys of sortedSystems) {
    const stats = systemStats[sys];

    csv += `${sys},${stats.total_experiments},${stats.completed_experiments},${stats.completion_rate},`;
    csv += `${stats.total_decisions},${stats.correct_decisions},${stats.accuracy},${stats.total_user_interventions}`;

    // Add per-dataset accuracy
    for (const ds of datasets) {
      const dsStats = stats.datasets[ds];
      if (dsStats && dsStats.completed) {
        csv += `,${dsStats.accuracy}`;
      } else {
        csv += `,N/A`;
      }
    }
    csv += '\n';
  }

  // Ensure output directory exists
  const outputDir = path.dirname(outputFile);
  if (!fs.existsSync(outputDir)) {
    fs.mkdirSync(outputDir, { recursive: true });
    console.log(`📁 Created output directory: ${outputDir}`);
  }

  fs.writeFileSync(outputFile, csv);
  console.log(`✅ Aggregated results exported to: ${outputFile}\n`);
}

// ========================================
// Print Summary Table
// ========================================

function printSummaryTable(systemStats) {
  console.log('\n' + '='.repeat(120));
  console.log('scRNA-SEQ MULTI-AGENT EVALUATION SUMMARY');
  console.log('='.repeat(120));
  console.log();

  const headers = ['System', 'Experiments', 'Completed', 'Rate%', 'Stages', 'Correct', 'Accuracy%', 'User Inputs'];
  const colWidths = [25, 12, 10, 8, 10, 10, 10, 12];

  // Print header
  let headerRow = '';
  for (let i = 0; i < headers.length; i++) {
    headerRow += headers[i].padEnd(colWidths[i]);
  }
  console.log(headerRow);
  console.log('-'.repeat(120));

  // Print rows
  const sortedSystems = Object.keys(systemStats).sort();

  for (const sys of sortedSystems) {
    const stats = systemStats[sys];

    const row = [
      sys,
      `${stats.total_experiments}`,
      `${stats.completed_experiments}`,
      `${stats.completion_rate}%`,
      `${stats.total_decisions}`,
      `${stats.correct_decisions}`,
      `${stats.accuracy}%`,
      `${stats.total_user_interventions}`
    ];

    let rowStr = '';
    for (let i = 0; i < row.length; i++) {
      rowStr += row[i].padEnd(colWidths[i]);
    }
    console.log(rowStr);
  }

  console.log('='.repeat(120));
  console.log();
}

// ========================================
// Main Execution
// ========================================

function main() {
  console.log('\n🔬 scRNA-seq Multi-Agent Evaluation Script');
  console.log('='.repeat(120));

  // Load ground truth
  const groundTruth = loadGroundTruth(GROUND_TRUTH_FILE);

  // Load CSV data
  const rows = loadCSV(INPUT_CSV);

  // Group by experiment
  const experiments = groupByExperiment(rows);

  if (Object.keys(experiments).length === 0) {
    console.log('⚠️  No experiments found. Exiting.\n');
    process.exit(0);
  }

  // Evaluate each experiment
  console.log('🔍 Evaluating experiments...\n');
  const results = [];

  for (const [key, experiment] of Object.entries(experiments)) {
    const result = evaluateExperiment(experiment, groundTruth);
    results.push(result);

    const status = result.completed ? '✓' : '✗';
    console.log(`   ${status} ${key.padEnd(70)} Accuracy: ${result.accuracy}%`);
  }

  // Aggregate by system
  console.log('\n📊 Aggregating results by system type...\n');
  const systemStats = aggregateBySystem(results);

  // Print summary
  printSummaryTable(systemStats);

  // Export aggregated CSV
  exportAggregatedCSV(systemStats, OUTPUT_CSV);

  // Export per-experiment CSV for statistical tests
  const perExperimentFile = OUTPUT_CSV.replace('.csv', '_per_experiment.csv');
  exportPerExperimentCSV(results, perExperimentFile);

  console.log('✨ Evaluation complete!\n');
  console.log('📊 For t-tests and statistical analysis, use:', perExperimentFile);
  console.log('   Example:');
  console.log('   python3 bin/run_ttest.py', perExperimentFile);
  console.log();
}

// Run
main();
