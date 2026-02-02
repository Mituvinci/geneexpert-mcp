#!/usr/bin/env node

/**
 * Bulk RNA-seq Multi-Agent Evaluation Script
 *
 * Compares agent decisions to ground truth across all experimental systems.
 * Reads from aggregated CSV file (NOT individual JSONL files).
 *
 * Usage:
 *   node bin/evaluate_bulk.js <input_csv> <ground_truth_json> <output_csv>
 *
 * Example:
 *   node bin/evaluate_bulk.js \
 *     experiments/results/bulk_rna_ALL_EXPERIMENTS_DETAILED.csv \
 *     experiments/bulk_rna_ground_truth.json \
 *     bulk_evaluation.csv
 */

import fs from 'fs';
import path from 'path';

// ========================================
// Configuration
// ========================================

const INPUT_CSV = process.argv[2] || 'experiments/results/bulk_rna_ALL_EXPERIMENTS_DETAILED.csv';
const GROUND_TRUTH_FILE = process.argv[3] || 'experiments/bulk_rna_ground_truth.json';
const OUTPUT_CSV = process.argv[4] || 'experiments/bulk_rna_csv_figures/bulk_evaluation.csv';

// Dataset mapping (short name in CSV -> ground truth key)
const DATASET_MAP = {
  '1_GSE52778_pe_clean': '1_GSE52778_pe_clean',
  '2_GSE114845_se_clean': '2_GSE114845_se_clean',
  '3_GSE113754_pe_clean': '3_GSE113754_pe_clean',
  '4_GSE141496_batch_effect': '4_GSE141496_batch_effect',
  '5_GSE47774_batch_effect': '5_GSE47774_batch_effect',
  '6_GSE193658_Lab_data': '6_GSE193658_Lab_data',
  '7_GSE114845_se_clean_CONTAM70_with_E.coli_GSE48151': '7_GSE114845_se_clean_CONTAM70_with_E.coli_GSE48151',
  // Handle typo in CSV: "wiht" instead of "with"
  '7_GSE114845_se_clean_CONTAM70_wiht_E.coli_GSE48151': '7_GSE114845_se_clean_CONTAM70_with_E.coli_GSE48151'
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
  // For bulk RNA: 4 stages
  // Stage 1 & 2: Use individual agent decisions (gpt5_2_decision, claude_decision, gemini_decision)
  // Stage 3: Use individual agent DE method and outlier decisions
  // Stage 4: Use individual agent decisions

  const stage = parseInt(row.stage);

  if (stage === 3) {
    // Stage 3: DE method selection and outlier action
    // Read individual agent decisions for both DE_Method and Outlier_Action
    const gpt_de = row.gpt5_2_de_method || 'N/A';
    const claude_de = row.claude_de_method || 'N/A';
    const gemini_de = row.gemini_de_method || 'N/A';

    const gpt_outlier = row.gpt5_2_outlier_action || 'N/A';
    const claude_outlier = row.claude_outlier_action || 'N/A';
    const gemini_outlier = row.gemini_outlier_action || 'N/A';

    // Calculate majority for DE method
    const de_decisions = [gpt_de, claude_de, gemini_de].filter(d => d !== 'N/A' && d);
    let finalDEMethod = row.de_method || 'N/A';
    if (de_decisions.length > 0) {
      const counts = {};
      for (const d of de_decisions) {
        counts[d] = (counts[d] || 0) + 1;
      }
      finalDEMethod = Object.keys(counts).reduce((a, b) => counts[a] > counts[b] ? a : b);
    }

    // Calculate majority for outlier action
    const outlier_decisions = [gpt_outlier, claude_outlier, gemini_outlier].filter(d => d !== 'N/A' && d);
    let finalOutlierAction = row.outlier_action || 'N/A';
    if (outlier_decisions.length > 0) {
      const counts = {};
      for (const d of outlier_decisions) {
        counts[d] = (counts[d] || 0) + 1;
      }
      finalOutlierAction = Object.keys(counts).reduce((a, b) => counts[a] > counts[b] ? a : b);
    }

    return {
      decision: finalDEMethod,
      de_method: finalDEMethod,
      outlier_action: finalOutlierAction,
      batch_specification: row.batch_specification || 'N/A',
      // Individual agent decisions for evaluation
      gpt_de_method: gpt_de,
      claude_de_method: claude_de,
      gemini_de_method: gemini_de,
      gpt_outlier_action: gpt_outlier,
      claude_outlier_action: claude_outlier,
      gemini_outlier_action: gemini_outlier,
      user_intervention: row.user_input_required === 'true' || row.user_input_required === '1'
    };
  } else {
    // Stage 1, 2, 4: Calculate majority from individual agent decisions
    // NOTE: consensus_decision contains vote result (approve/reject), NOT the actual QC decision
    // We need to use gpt5_2_decision, claude_decision, gemini_decision instead
    const gpt = row.gpt5_2_decision || 'N/A';
    const claude = row.claude_decision || 'N/A';
    const gemini = row.gemini_decision || 'N/A';

    // Calculate majority vote from agent decisions
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
}

// ========================================
// Compare Decisions
// ========================================

function compareDecision(agentData, groundTruth, stageKey) {
  const comparison = {
    stage: stageKey,
    match: false,
    agent_decision: agentData.decision,
    correct_decision: groundTruth ? groundTruth.correct_decision : null,
    user_intervention: agentData.user_intervention || false,
    rationale: groundTruth ? groundTruth.rationale : 'No ground truth'
  };

  if (!groundTruth) {
    return comparison;
  }

  // Normalize and compare
  const normalized_agent = normalizeDecision(agentData.decision);
  const normalized_correct = normalizeDecision(groundTruth.correct_decision);

  comparison.match = normalized_agent === normalized_correct;

  // Stage 3 specific: Also check batch and outlier actions if available
  if (stageKey === 'stage3') {
    if (agentData.batch_specification && groundTruth.correct_batch_specification) {
      comparison.batch_match = normalizeBatchSpec(agentData.batch_specification) ===
                                normalizeBatchSpec(groundTruth.correct_batch_specification);
    }
    if (agentData.outlier_action && groundTruth.correct_outlier_action) {
      comparison.outlier_match = normalizeDecision(agentData.outlier_action) ===
                                 normalizeDecision(groundTruth.correct_outlier_action);
    }

    // Individual agent comparisons for Stage 3
    comparison.individual_agents = {
      gpt: {
        de_method: agentData.gpt_de_method || 'N/A',
        de_method_correct: normalizeDecision(agentData.gpt_de_method) ===
                           normalizeDecision(groundTruth.correct_de_method),
        outlier_action: agentData.gpt_outlier_action || 'N/A',
        outlier_correct: normalizeDecision(agentData.gpt_outlier_action) ===
                         normalizeDecision(groundTruth.correct_outlier_action)
      },
      claude: {
        de_method: agentData.claude_de_method || 'N/A',
        de_method_correct: normalizeDecision(agentData.claude_de_method) ===
                           normalizeDecision(groundTruth.correct_de_method),
        outlier_action: agentData.claude_outlier_action || 'N/A',
        outlier_correct: normalizeDecision(agentData.claude_outlier_action) ===
                         normalizeDecision(groundTruth.correct_outlier_action)
      },
      gemini: {
        de_method: agentData.gemini_de_method || 'N/A',
        de_method_correct: normalizeDecision(agentData.gemini_de_method) ===
                           normalizeDecision(groundTruth.correct_de_method),
        outlier_action: agentData.gemini_outlier_action || 'N/A',
        outlier_correct: normalizeDecision(agentData.gemini_outlier_action) ===
                         normalizeDecision(groundTruth.correct_outlier_action)
      }
    };
  }

  return comparison;
}

// ========================================
// Normalize Decisions
// ========================================

function normalizeDecision(decision) {
  if (!decision || decision === 'N/A') return 'unknown';

  let normalized = decision.toString().toLowerCase().trim();

  // Common synonyms
  if (normalized.includes('pass') && normalized.includes('all')) return 'passall';
  if (normalized.includes('pass') && normalized.includes('warning')) return 'passwithwarning';
  if (normalized === 'pass') return 'pass';
  if (normalized === 'approve' || normalized === 'approved') return 'approve';
  if (normalized === 'reject' || normalized === 'rejected') return 'reject';
  if (normalized === 'simpleedger' || normalized === 'simple_edger') return 'simpleedger';
  if (normalized.includes('batch') && normalized.includes('effect')) return 'batcheffectedger';

  // Remove special characters and spaces
  return normalized.replace(/[_\s-]/g, '').replace(/\*/g, '');
}

function normalizeBatchSpec(spec) {
  if (!spec || spec === 'N/A') return 'auto';
  return spec.toString().toLowerCase().trim();
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
    stage_results: {},
    completed: false,
    accuracy: 0,
    auto_accuracy: 0,
    safe_escalation_rate: 0
  };

  if (!gtDataset) {
    result.error = `No ground truth for dataset: ${datasetKey}`;
    return result;
  }

  const gtDecisions = gtDataset.decisions;

  // Process each stage (1, 2, 3, 4)
  const stagesReached = new Set();

  for (const row of experiment.stages) {
    const stageNum = parseInt(row.stage);
    const stageKey = `stage${stageNum}`;

    // Track which stages were reached
    stagesReached.add(stageNum);

    // Find ground truth for this stage by matching decision_id pattern
    let gtStage = null;
    for (const [gtKey, gtValue] of Object.entries(gtDecisions)) {
      if (gtKey.includes(stageKey)) {
        gtStage = gtValue;
        break;
      }
    }

    if (!gtStage) {
      continue; // Skip stages without ground truth
    }

    result.total_stages++;

    // Extract agent decisions from CSV row
    const agentData = extractDecisions(row);

    // Compare to ground truth
    const comparison = compareDecision(agentData, gtStage, stageKey);

    result.stage_results[stageKey] = comparison;

    // Track overall correctness
    if (comparison.match) {
      result.correct_decisions++;
    }

    // Separate auto-decisions from user-escalated decisions
    if (comparison.user_intervention) {
      result.user_interventions++;

      // Justified escalation: User intervention led to correct decision
      // OR agents disagreed (escalation was appropriate)
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

    // Check if stage was completed (PASS, PASS_WITH_WARNING, or SUCCESS)
    const status = row.stage_status ? row.stage_status.toUpperCase() : '';
    if (status === 'PASS' || status === 'SUCCESS' || status === 'PASS_WITH_WARNING') {
      result.completed_stages++;
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

  // Check completion: Either completed all 4 stages OR reached stage 4
  // (Some experiments may have warnings but still completed all stages)
  result.completed = result.completed_stages >= 4 || stagesReached.has(4);

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
  let csv = 'Experiment,System,Dataset,Organism,Completed,Total_Stages,Correct_Decisions,Accuracy_%,Auto_Decisions,Auto_Correct,Auto_Accuracy_%,User_Interventions,Justified_Escalations,Safe_Escalation_Rate_%,Stage1_Match,Stage2_Match,Stage3_Match,Stage4_Match,Stage1_UserInput,Stage2_UserInput,Stage3_UserInput,Stage4_UserInput\n';

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

    csv += `${result.experiment},${result.system},${result.dataset},${result.organism},${result.completed ? 1 : 0},`;
    csv += `${result.total_stages},${result.correct_decisions},${result.accuracy},`;
    csv += `${result.auto_decisions},${result.auto_correct},${result.auto_accuracy},`;
    csv += `${result.user_interventions},${result.justified_escalations},${result.safe_escalation_rate},`;
    csv += `${stage1},${stage2},${stage3},${stage4},`;
    csv += `${stage1_user},${stage2_user},${stage3_user},${stage4_user}\n`;
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
    const shortName = ds.split('_')[1] || ds; // GSE52778, GSE114845, etc.
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
  console.log('BULK RNA-SEQ MULTI-AGENT EVALUATION SUMMARY');
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
  console.log('\n🔬 Bulk RNA-seq Multi-Agent Evaluation Script (CSV-based)');
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
