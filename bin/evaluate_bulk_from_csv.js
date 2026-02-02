#!/usr/bin/env node

/**
 * Bulk RNA-seq Multi-Agent Evaluation Script (CSV-Based)
 *
 * Uses pre-generated CSV files instead of reading individual JSONL files.
 * This approach is faster and more convenient.
 *
 * Usage:
 *   node bin/evaluate_bulk_from_csv.js experiments/results/bulk_rna_ALL_EXPERIMENTS_DETAILED.csv experiments/bulk_rna_ground_truth.json bulk_evaluation.csv
 */

import fs from 'fs';
import path from 'path';

// ========================================
// Configuration
// ========================================

const DETAILED_CSV = process.argv[2] || 'experiments/results/bulk_rna_ALL_EXPERIMENTS_DETAILED.csv';
const GROUND_TRUTH_FILE = process.argv[3] || 'experiments/bulk_rna_ground_truth.json';
const OUTPUT_CSV = process.argv[4] || 'bulk_evaluation.csv';

// Dataset mapping (folder prefix -> ground truth key)
const DATASET_MAP = {
  '1_GSE52778': '1_GSE52778_pe_clean',
  '2_GSE114845': '2_GSE114845_se_clean',
  '3_GSE113754': '3_GSE113754_pe_clean',
  '4_GSE141496': '4_GSE141496_batch_effect',
  '5_GSE47774': '5_GSE47774_batch_effect',
  '6_GSE193658': '6_GSE193658_Lab_data',
  '7_GSE114845_se_clean_CONTAM70': '7_GSE114845_se_clean_CONTAM70_with_E.coli_GSE48151'
};

// System classification patterns (order matters - more specific first)
const SYSTEM_PATTERNS = {
  'no_agent': /no_agent/,
  'single_agent_claude': /single_claude/,
  'single_agent_gpt': /single_gpt/,
  'single_agent_gemini': /single_gemini/,
  'default_parallel': /^parallel_default$/,
  'default_sequential': /^sequential_default$/,
  // Parallel role permutations
  'parallel_perm_1': /parallel.*gp_st_cl_bl_gm_pl/,
  'parallel_perm_2': /parallel.*gp_st_cl_pl_gm_bl/,
  'parallel_perm_3': /parallel.*gp_pl_cl_st_gm_bl/,
  'parallel_perm_4': /parallel.*gp_pl_cl_bl_gm_st/,
  'parallel_perm_5': /parallel.*gp_bl_cl_st_gm_pl/,
  'parallel_perm_6': /parallel.*gp_bl_cl_pl_gm_st/,
  // Sequential role permutations
  'sequential_perm_1': /sequential.*gp_st_cl_bl_gm_pl/,
  'sequential_perm_2': /sequential.*gp_st_cl_pl_gm_bl/,
  'sequential_perm_3': /sequential.*gp_pl_cl_st_gm_bl/,
  'sequential_perm_4': /sequential.*gp_pl_cl_bl_gm_st/,
  'sequential_perm_5': /sequential.*gp_bl_cl_st_gm_pl/,
  'sequential_perm_6': /sequential.*gp_bl_cl_pl_gm_st/
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
// Parse CSV
// ========================================

function parseCSV(filepath) {
  console.log(`📁 Loading detailed CSV: ${filepath}\n`);
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

  console.log(`   ✓ Loaded ${rows.length} stage decisions\n`);
  return rows;
}

function parseCSVLine(line) {
  const values = [];
  let current = '';
  let inQuotes = false;

  for (let i = 0; i < line.length; i++) {
    const char = line[i];

    if (char === '"') {
      inQuotes = !inQuotes;
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
// Extract Dataset Name
// ========================================

function extractDatasetName(datasetString) {
  for (const [prefix, gtKey] of Object.entries(DATASET_MAP)) {
    if (datasetString.startsWith(prefix)) {
      return gtKey;
    }
  }
  return null;
}

// ========================================
// Classify System Type
// ========================================

function classifySystem(systemConfig) {
  for (const [systemName, pattern] of Object.entries(SYSTEM_PATTERNS)) {
    if (pattern.test(systemConfig)) {
      return systemName;
    }
  }
  return 'unknown';
}

// ========================================
// Compare Decisions
// ========================================

function compareDecision(agentDecision, groundTruth, stage) {
  const comparison = {
    stage,
    match: false,
    agent_decision: null,
    correct_decision: null,
    details: {}
  };

  if (!agentDecision) {
    comparison.details.error = 'No agent decision found';
    return comparison;
  }

  if (!groundTruth) {
    comparison.details.error = 'No ground truth available';
    return comparison;
  }

  comparison.agent_decision = agentDecision;
  comparison.correct_decision = groundTruth.correct_decision;

  // Normalize and compare
  comparison.match = normalizeDecision(agentDecision) === normalizeDecision(groundTruth.correct_decision);

  return comparison;
}

function normalizeDecision(decision) {
  if (!decision) return 'unknown';
  return decision.toString().toLowerCase().trim().replace(/[_\s-]/g, '');
}

// ========================================
// Process Experiments
// ========================================

function processExperiments(csvRows, groundTruth) {
  const experiments = {};

  for (const row of csvRows) {
    const expKey = `${row.full_dataset}_${row.system_config}`;
    const datasetKey = extractDatasetName(row.full_dataset);
    // Keep original system_config name instead of classifying
    const systemName = row.system_config;

    if (!experiments[expKey]) {
      experiments[expKey] = {
        folder: expKey,
        dataset: datasetKey,
        system: systemName,  // Use full original name
        completed: row.Final_Status === 'SUCCESS' || false,
        total_decisions: 0,
        correct_decisions: 0,
        user_interventions: 0,
        total_cost: 0,
        stage_results: {}
      };
    }

    const exp = experiments[expKey];
    const stage = parseInt(row.stage);
    const stageKey = `stage${stage}`;

    // Get consensus decision
    const consensusDecision = row.consensus_decision;
    const userInputRequired = row.user_input_required === 'True' || row.user_input_required === 'true';

    // Extract individual agent decisions
    let finalDecision = null;
    let agentDecisions = {};

    if (stage === 3) {
      // Stage 3: Use individual agent DE method and outlier action
      const gpt_de = row.gpt5_2_de_method || 'N/A';
      const claude_de = row.claude_de_method || 'N/A';
      const gemini_de = row.gemini_de_method || 'N/A';

      const gpt_outlier = row.gpt5_2_outlier_action || 'N/A';
      const claude_outlier = row.claude_outlier_action || 'N/A';
      const gemini_outlier = row.gemini_outlier_action || 'N/A';

      // Calculate majority for DE method
      const de_decisions = [gpt_de, claude_de, gemini_de].filter(d => d && d !== 'N/A');
      if (de_decisions.length > 0) {
        const counts = {};
        for (const d of de_decisions) {
          counts[d] = (counts[d] || 0) + 1;
        }
        finalDecision = Object.keys(counts).reduce((a, b) => counts[a] > counts[b] ? a : b);
      } else {
        finalDecision = row.de_method || consensusDecision;
      }

      // Store individual agent decisions for Stage 3
      agentDecisions = {
        gpt: { de_method: gpt_de, outlier_action: gpt_outlier },
        claude: { de_method: claude_de, outlier_action: claude_outlier },
        gemini: { de_method: gemini_de, outlier_action: gemini_outlier }
      };
    } else {
      // Use majority decision from agents
      const gptDec = row.gpt5_2_decision;
      const claudeDec = row.claude_decision;
      const geminiDec = row.gemini_decision;

      // Get most common
      const decisions = [gptDec, claudeDec, geminiDec].filter(d => d);
      const counts = {};
      for (const d of decisions) {
        counts[d] = (counts[d] || 0) + 1;
      }
      finalDecision = Object.keys(counts).reduce((a, b) => counts[a] > counts[b] ? a : b, '');

      agentDecisions = {
        gpt: gptDec,
        claude: claudeDec,
        gemini: geminiDec
      };
    }

    // Get ground truth for this stage
    const gtDecisions = groundTruth[datasetKey]?.decisions || {};
    const gtStageKey = Object.keys(gtDecisions).find(k => k.includes(stageKey));
    const gtDecision = gtDecisions[gtStageKey];

    // Compare
    const comparison = compareDecision(finalDecision, gtDecision, stageKey);

    // Add user intervention flag and individual agent decisions to stage result
    comparison.user_intervention = userInputRequired;
    comparison.agent_decisions = agentDecisions;

    exp.stage_results[stageKey] = comparison;
    exp.total_decisions++;

    if (comparison.match) {
      exp.correct_decisions++;
    }

    if (userInputRequired) {
      exp.user_interventions++;
    }

    // Add cost
    const stageCost = parseFloat(row.stage_total_cost_usd) || 0;
    exp.total_cost += stageCost;
  }

  // Calculate accuracy
  for (const exp of Object.values(experiments)) {
    exp.accuracy = exp.total_decisions > 0
      ? (exp.correct_decisions / exp.total_decisions * 100).toFixed(1)
      : 0;
  }

  return Object.values(experiments);
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
        total_cost: 0,
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

    stats.total_decisions += result.total_decisions;
    stats.correct_decisions += result.correct_decisions;
    stats.total_user_interventions += result.user_interventions;
    stats.total_cost += result.total_cost;

    // Track per-dataset
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
    stats.avg_cost = stats.total_experiments > 0
      ? (stats.total_cost / stats.total_experiments).toFixed(4)
      : 0;
  }

  return systemStats;
}

// ========================================
// Export to CSV
// ========================================

function exportToCSV(systemStats, outputFile) {
  const datasets = ['1_GSE52778_pe_clean', '2_GSE114845_se_clean', '3_GSE113754_pe_clean',
                    '4_GSE141496_batch_effect', '5_GSE47774_batch_effect',
                    '6_GSE193658_Lab_data', '7_GSE114845_se_clean_CONTAM70_with_E.coli_GSE48151'];

  let csv = 'System,Total_Experiments,Completed,Completion_Rate_%,Total_Decisions,Correct_Decisions,Decision_Accuracy_%,User_Interventions,Avg_Cost_USD';

  for (const ds of datasets) {
    const shortName = ds.split('_')[1] || ds;
    csv += `,${shortName}_Accuracy`;
  }
  csv += '\n';

  const sortedSystems = Object.keys(systemStats).sort();

  for (const sys of sortedSystems) {
    const stats = systemStats[sys];

    csv += `${sys},${stats.total_experiments},${stats.completed_experiments},${stats.completion_rate},`;
    csv += `${stats.total_decisions},${stats.correct_decisions},${stats.accuracy},${stats.total_user_interventions},${stats.avg_cost}`;

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

  fs.writeFileSync(outputFile, csv);
  console.log(`\n✅ Aggregated results exported to: ${outputFile}\n`);
}

function exportPerExperimentCSV(results, outputFile) {
  let csv = 'Experiment,System,Dataset,Completed,Total_Decisions,Correct_Decisions,Accuracy_%,User_Interventions,Total_Cost_USD,Stage1_Match,Stage2_Match,Stage3_Match,Stage4_Match,Stage1_UserInput,Stage2_UserInput,Stage3_UserInput,Stage4_UserInput\n';

  for (const result of results) {
    // Decision match (1/0)
    const stage1 = result.stage_results.stage1 ? (result.stage_results.stage1.match ? 1 : 0) : '';
    const stage2 = result.stage_results.stage2 ? (result.stage_results.stage2.match ? 1 : 0) : '';
    const stage3 = result.stage_results.stage3 ? (result.stage_results.stage3.match ? 1 : 0) : '';
    const stage4 = result.stage_results.stage4 ? (result.stage_results.stage4.match ? 1 : 0) : '';

    // User input required (1/0)
    const stage1_user = result.stage_results.stage1 ? (result.stage_results.stage1.user_intervention ? 1 : 0) : '';
    const stage2_user = result.stage_results.stage2 ? (result.stage_results.stage2.user_intervention ? 1 : 0) : '';
    const stage3_user = result.stage_results.stage3 ? (result.stage_results.stage3.user_intervention ? 1 : 0) : '';
    const stage4_user = result.stage_results.stage4 ? (result.stage_results.stage4.user_intervention ? 1 : 0) : '';

    csv += `${result.folder},${result.system},${result.dataset},${result.completed ? 1 : 0},`;
    csv += `${result.total_decisions},${result.correct_decisions},${result.accuracy},${result.user_interventions},${result.total_cost.toFixed(4)},`;
    csv += `${stage1},${stage2},${stage3},${stage4},`;
    csv += `${stage1_user},${stage2_user},${stage3_user},${stage4_user}\n`;
  }

  fs.writeFileSync(outputFile, csv);
  console.log(`✅ Per-experiment results exported to: ${outputFile}\n`);
}

// ========================================
// Print Summary Table
// ========================================

function printSummaryTable(systemStats) {
  console.log('\n' + '='.repeat(130));
  console.log('BULK RNA-SEQ MULTI-AGENT EVALUATION SUMMARY');
  console.log('='.repeat(130));
  console.log();

  const headers = ['System', 'Experiments', 'Completed', 'Rate%', 'Decisions', 'Correct', 'Accuracy%', 'User Inputs', 'Avg Cost'];
  const colWidths = [25, 12, 10, 8, 10, 10, 10, 12, 10];

  let headerRow = '';
  for (let i = 0; i < headers.length; i++) {
    headerRow += headers[i].padEnd(colWidths[i]);
  }
  console.log(headerRow);
  console.log('-'.repeat(130));

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
      `${stats.total_user_interventions}`,
      `$${stats.avg_cost}`
    ];

    let rowStr = '';
    for (let i = 0; i < row.length; i++) {
      rowStr += row[i].padEnd(colWidths[i]);
    }
    console.log(rowStr);
  }

  console.log('='.repeat(130));
  console.log();
}

// ========================================
// Main Execution
// ========================================

function main() {
  console.log('\n🔬 Bulk RNA-seq Multi-Agent Evaluation Script (CSV-Based)');
  console.log('=' .repeat(130));

  // Load ground truth
  const groundTruth = loadGroundTruth(GROUND_TRUTH_FILE);

  // Load detailed CSV
  const csvRows = parseCSV(DETAILED_CSV);

  // Process experiments
  console.log('🔍 Evaluating experiments...\n');
  const results = processExperiments(csvRows, groundTruth);

  for (const result of results) {
    const status = result.completed ? '✓' : '✗';
    console.log(`   ${status} ${result.folder.padEnd(70)} Accuracy: ${result.accuracy}%`);
  }

  // Aggregate by system
  console.log('\n📊 Aggregating results by system type...\n');
  const systemStats = aggregateBySystem(results);

  // Print summary
  printSummaryTable(systemStats);

  // Export to CSV
  exportToCSV(systemStats, OUTPUT_CSV);

  const perExperimentFile = OUTPUT_CSV.replace('.csv', '_per_experiment.csv');
  exportPerExperimentCSV(results, perExperimentFile);

  console.log('✨ Evaluation complete!\n');
}

main();
