#!/usr/bin/env node

/**
 * Evaluation Script for ICML 2026 Experiments
 *
 * Usage:
 *   node bin/evaluate.js --results experiments/results --truth experiments/ground_truth.json
 *   node bin/evaluate.js --system multi-agent --dataset GSE52778
 *
 * Metrics Calculated:
 *   1. Decision accuracy (correct decisions / total)
 *   2. Error reduction (baseline errors - system errors) / baseline errors
 *   3. Success rate (successful analyses / total)
 *   4. Cost efficiency (total cost / successes)
 *   5. Inter-agent agreement (Cohen's kappa)
 *   6. User input frequency (% decisions requiring user input)
 *   7. Error propagation rate (sequential mode only)
 */

import { Command } from 'commander';
import fs from 'fs';
import path from 'path';
import { globSync } from 'glob';
import { normalizeDecision as normalizeDecisionCanonical } from '../src/config/decision_vocabulary.js';

const program = new Command();

program
  .name('evaluate')
  .description('Evaluate ICML experiment results and calculate metrics')
  .version('1.0.0');

program
  .command('metrics')
  .description('Calculate all metrics for experiment results')
  .option('--results <dir>', 'Results directory', 'experiments/results')
  .option('--truth <file>', 'Ground truth JSON file', 'experiments/ground_truth.json')
  .option('--system <name>', 'Filter by system (multi-agent, single-agent-gpt5.2, etc.)')
  .option('--dataset <name>', 'Filter by dataset (GSE52778, GSE114845, etc.)')
  .option('--output <file>', 'Output CSV file (optional)')
  .action((options) => {
    console.log('');
    console.log('='.repeat(70));
    console.log('ICML 2026 EXPERIMENT EVALUATION');
    console.log('='.repeat(70));
    console.log('');

    // Load ground truth
    const groundTruth = loadGroundTruth(options.truth);
    console.log(`✓ Loaded ground truth: ${Object.keys(groundTruth).length} labels`);

    // Load all experiment results
    const experiments = loadExperiments(options.results, options.system, options.dataset);
    console.log(`✓ Loaded experiments: ${experiments.length} total`);
    console.log('');

    if (experiments.length === 0) {
      console.log('⚠️  No experiments found. Run some analyses first!');
      console.log('');
      process.exit(0);
    }

    // Get unique systems
    const systems = [...new Set(experiments.map(e => e.system))];

    console.log('Systems found:', systems.join(', '));
    console.log('');

    // Calculate metrics for each system
    const allMetrics = [];

    systems.forEach(system => {
      const systemExps = experiments.filter(e => e.system === system);

      console.log('='.repeat(70));
      console.log(`SYSTEM: ${system.toUpperCase()}`);
      console.log('='.repeat(70));
      console.log('');

      const metrics = {
        system,
        total_experiments: systemExps.length,
        ...calculateDecisionAccuracy(systemExps, groundTruth),
        ...calculateSuccessRate(systemExps),
        ...calculateCostEfficiency(systemExps),
        ...calculateUserInputFrequency(systemExps)
      };

      // Multi-agent specific metrics
      if (system.includes('multi-agent')) {
        const agentMetrics = calculateInterAgentAgreement(systemExps);
        Object.assign(metrics, agentMetrics);
      }

      // Sequential mode specific metrics
      if (system.includes('sequential')) {
        const errorProp = calculateErrorPropagation(systemExps, groundTruth);
        Object.assign(metrics, errorProp);
      }

      // Display metrics
      displayMetrics(metrics);
      allMetrics.push(metrics);
    });

    // Calculate error reduction (multi-agent vs baseline)
    console.log('');
    console.log('='.repeat(70));
    console.log('ERROR REDUCTION ANALYSIS');
    console.log('='.repeat(70));
    console.log('');

    const baseline = experiments.filter(e => e.system === 'no-agent' || e.system === 'single-agent-gpt5.2');
    const multiAgent = experiments.filter(e => e.system === 'multi-agent-parallel' || e.system === 'multi-agent');

    if (baseline.length > 0 && multiAgent.length > 0) {
      const errorReduction = calculateErrorReduction(baseline, multiAgent, groundTruth);
      console.log(`Baseline Errors:      ${errorReduction.baselineErrors}`);
      console.log(`Multi-Agent Errors:   ${errorReduction.multiAgentErrors}`);
      console.log(`Error Reduction:      ${errorReduction.reductionPercent} (${errorReduction.errorReduction.toFixed(3)})`);
      console.log('');
    } else {
      console.log('⚠️  Need both baseline and multi-agent experiments for error reduction');
      console.log('');
    }

    // Export to CSV if requested
    if (options.output) {
      exportMetricsToCSV(allMetrics, options.output);
      console.log(`✓ Metrics exported to: ${options.output}`);
      console.log('');
    }

    console.log('='.repeat(70));
    console.log('');
  });

program
  .command('compare')
  .description('Compare parallel vs sequential multi-agent modes')
  .option('--results <dir>', 'Results directory', 'experiments/results')
  .option('--truth <file>', 'Ground truth JSON file', 'experiments/ground_truth.json')
  .action((options) => {
    console.log('');
    console.log('='.repeat(70));
    console.log('PARALLEL vs SEQUENTIAL COMPARISON');
    console.log('='.repeat(70));
    console.log('');

    const groundTruth = loadGroundTruth(options.truth);
    const experiments = loadExperiments(options.results);

    const parallel = experiments.filter(e => e.system === 'multi-agent-parallel' || (e.system === 'multi-agent' && !e.sequential));
    const sequential = experiments.filter(e => e.system === 'multi-agent-sequential' || (e.system === 'multi-agent' && e.sequential));

    if (parallel.length === 0 || sequential.length === 0) {
      console.log('⚠️  Need both parallel and sequential experiments for comparison');
      console.log('');
      process.exit(0);
    }

    console.log(`Parallel experiments:   ${parallel.length}`);
    console.log(`Sequential experiments: ${sequential.length}`);
    console.log('');

    // Calculate metrics for both
    const parallelMetrics = {
      ...calculateDecisionAccuracy(parallel, groundTruth),
      ...calculateSuccessRate(parallel),
      ...calculateCostEfficiency(parallel),
      ...calculateInterAgentAgreement(parallel)
    };

    const sequentialMetrics = {
      ...calculateDecisionAccuracy(sequential, groundTruth),
      ...calculateSuccessRate(sequential),
      ...calculateCostEfficiency(sequential),
      ...calculateInterAgentAgreement(sequential),
      ...calculateErrorPropagation(sequential, groundTruth)
    };

    // Display comparison
    console.log('METRIC                          PARALLEL        SEQUENTIAL');
    console.log('-'.repeat(70));
    console.log(`Decision Accuracy               ${(parallelMetrics.accuracy * 100).toFixed(1)}%          ${(sequentialMetrics.accuracy * 100).toFixed(1)}%`);
    console.log(`Success Rate                    ${(parallelMetrics.successRate * 100).toFixed(1)}%          ${(sequentialMetrics.successRate * 100).toFixed(1)}%`);
    console.log(`Cost per Success                $${parallelMetrics.costPerSuccess}         $${sequentialMetrics.costPerSuccess}`);
    console.log(`Cohen's Kappa                   ${parallelMetrics.cohensKappa || 'N/A'}           ${sequentialMetrics.cohensKappa || 'N/A'}`);
    if (sequentialMetrics.errorPropagationRate !== undefined) {
      console.log(`Error Propagation Rate          N/A             ${(sequentialMetrics.errorPropagationRate * 100).toFixed(1)}%`);
    }
    console.log('');
    console.log('='.repeat(70));
    console.log('');
  });

program.parse();

// ============================================================================
// Helper Functions
// ============================================================================

/**
 * Load ground truth labels from JSON file
 */
function loadGroundTruth(filepath) {
  if (!fs.existsSync(filepath)) {
    console.warn(`⚠️  Ground truth file not found: ${filepath}`);
    return {};
  }

  try {
    const data = JSON.parse(fs.readFileSync(filepath, 'utf-8'));
    const groundTruth = {};

    data.forEach(item => {
      if (item.dataset && item.stage && item.expected_decision) {
        const key = `${item.dataset}_stage${item.stage}`;
        groundTruth[key] = {
          decision: item.expected_decision.toLowerCase(),
          reason: item.reason || ''
        };
      }
    });

    return groundTruth;
  } catch (error) {
    console.warn(`⚠️  Error loading ground truth: ${error.message}`);
    return {};
  }
}

/**
 * Load all experiments from results directory
 */
function loadExperiments(resultsDir, systemFilter = null, datasetFilter = null) {
  if (!fs.existsSync(resultsDir)) {
    return [];
  }

  // Find all agent decision files
  const pattern = path.join(resultsDir, '**/*_agent_decisions.json');
  const files = globSync(pattern);

  const experiments = [];

  files.forEach(file => {
    try {
      const data = JSON.parse(fs.readFileSync(file, 'utf-8'));

      // Extract system and dataset from path or metadata
      const pathParts = file.split('/');
      const dataset = data.dataset || pathParts[pathParts.length - 2] || 'unknown';
      const system = data.system || detectSystemFromPath(file);

      // Apply filters
      if (systemFilter && system !== systemFilter) return;
      if (datasetFilter && dataset !== datasetFilter) return;

      // Process each stage decision
      if (data.decisions && Array.isArray(data.decisions)) {
        data.decisions.forEach(decision => {
          experiments.push({
            file,
            dataset,
            system,
            stage: decision.stage,
            decision: decision.consensus?.decision || decision.decision,
            confidence: decision.consensus?.confidence || 0,
            agents: decision.agents || [],
            sequential: data.sequential || false,
            timestamp: data.timestamp || decision.timestamp,
            cost: decision.cost || 0,
            success: decision.proceed !== false,
            userInputRequired: decision.user_input_required || false
          });
        });
      }
    } catch (error) {
      console.warn(`⚠️  Error reading ${file}: ${error.message}`);
    }
  });

  return experiments;
}

/**
 * Detect system type from file path
 */
function detectSystemFromPath(filepath) {
  const lower = filepath.toLowerCase();
  if (lower.includes('no-agent') || lower.includes('no_agent')) return 'no-agent';
  if (lower.includes('single-agent-gpt') || lower.includes('single_gpt')) return 'single-agent-gpt5.2';
  if (lower.includes('single-agent-claude') || lower.includes('single_claude')) return 'single-agent-claude';
  if (lower.includes('single-agent-gemini') || lower.includes('single_gemini')) return 'single-agent-gemini';
  if (lower.includes('sequential')) return 'multi-agent-sequential';
  if (lower.includes('parallel')) return 'multi-agent-parallel';
  if (lower.includes('multi-agent') || lower.includes('multi_agent')) return 'multi-agent';
  return 'unknown';
}

/**
 * Calculate decision accuracy vs ground truth
 */
function calculateDecisionAccuracy(experiments, groundTruth) {
  let correct = 0;
  let total = 0;

  experiments.forEach(exp => {
    const key = `${exp.dataset}_stage${exp.stage}`;
    const truth = groundTruth[key];

    if (truth && truth.decision) {
      total++;

      // Normalize both decisions using canonical vocabulary (pass stage number)
      const normalized = normalizeDecision(exp.decision || '', exp.stage);
      const normalizedTruth = normalizeDecision(truth.decision, exp.stage);

      if (normalized && normalizedTruth && normalized === normalizedTruth) {
        correct++;
      }
    }
  });

  return {
    accuracy: total > 0 ? correct / total : 0,
    correct_decisions: correct,
    total_decisions: total
  };
}

/**
 * Normalize decision names for comparison using canonical vocabulary
 */
function normalizeDecision(decision, stage) {
  if (!decision || typeof decision !== 'string') {
    return null;
  }

  // Map stage numbers to vocabulary identifiers
  const stageMap = {
    1: 'stage1',
    2: 'stage2',
    3: 'stage3_de_method',  // Default for stage 3, may need refinement
    4: 'stage4'
  };

  const stageId = stageMap[stage] || `stage${stage}`;

  // Use canonical vocabulary normalization
  const normalized = normalizeDecisionCanonical(decision, stageId);

  // If canonical normalization succeeds, return it
  if (normalized) {
    return normalized;
  }

  // Fallback: return lowercase trimmed version
  return decision.toLowerCase().trim();
}

/**
 * Calculate success rate
 */
function calculateSuccessRate(experiments) {
  const successes = experiments.filter(e => e.success).length;
  const total = experiments.length;

  return {
    successRate: total > 0 ? successes / total : 0,
    successes,
    total_attempts: total
  };
}

/**
 * Calculate cost efficiency
 */
function calculateCostEfficiency(experiments) {
  const totalCost = experiments.reduce((sum, e) => sum + (e.cost || 0), 0);
  const successes = experiments.filter(e => e.success).length;

  return {
    totalCost: totalCost.toFixed(3),
    costPerSuccess: successes > 0 ? (totalCost / successes).toFixed(3) : 'N/A',
    avgCostPerDecision: experiments.length > 0 ? (totalCost / experiments.length).toFixed(3) : '0.000'
  };
}

/**
 * Calculate user input frequency
 */
function calculateUserInputFrequency(experiments) {
  const userInputCount = experiments.filter(e => e.userInputRequired).length;
  const total = experiments.length;

  return {
    userInputFrequency: total > 0 ? userInputCount / total : 0,
    userInputCount,
    autonomyRate: total > 0 ? (total - userInputCount) / total : 0
  };
}

/**
 * Calculate Cohen's Kappa for inter-agent agreement
 */
function calculateInterAgentAgreement(experiments) {
  const multiAgentExps = experiments.filter(e => e.agents && e.agents.length === 3);

  if (multiAgentExps.length === 0) {
    return { cohensKappa: null };
  }

  // Extract agent decisions
  const agent1Decisions = [];
  const agent2Decisions = [];
  const agent3Decisions = [];

  multiAgentExps.forEach(exp => {
    if (exp.agents.length === 3) {
      // Normalize decisions using stage information
      const d1 = normalizeDecision(exp.agents[0].decision || '', exp.stage);
      const d2 = normalizeDecision(exp.agents[1].decision || '', exp.stage);
      const d3 = normalizeDecision(exp.agents[2].decision || '', exp.stage);

      if (d1 && d2 && d3) {
        agent1Decisions.push(d1);
        agent2Decisions.push(d2);
        agent3Decisions.push(d3);
      }
    }
  });

  if (agent1Decisions.length === 0) {
    return { cohensKappa: null };
  }

  // Calculate pairwise kappas
  const kappa12 = calculatePairwiseKappa(agent1Decisions, agent2Decisions);
  const kappa13 = calculatePairwiseKappa(agent1Decisions, agent3Decisions);
  const kappa23 = calculatePairwiseKappa(agent2Decisions, agent3Decisions);

  const avgKappa = (kappa12 + kappa13 + kappa23) / 3;

  return {
    cohensKappa: isNaN(avgKappa) ? null : avgKappa.toFixed(3),
    kappa_agent1_agent2: kappa12.toFixed(3),
    kappa_agent1_agent3: kappa13.toFixed(3),
    kappa_agent2_agent3: kappa23.toFixed(3),
    agreementDataPoints: agent1Decisions.length
  };
}

/**
 * Calculate pairwise Cohen's Kappa
 */
function calculatePairwiseKappa(rater1, rater2) {
  if (rater1.length !== rater2.length || rater1.length === 0) {
    return 0;
  }

  const n = rater1.length;
  let agreed = 0;

  // Observed agreement
  for (let i = 0; i < n; i++) {
    if (rater1[i] === rater2[i]) {
      agreed++;
    }
  }
  const po = agreed / n;

  // Expected agreement by chance
  const categories = [...new Set([...rater1, ...rater2])];
  let pe = 0;

  categories.forEach(cat => {
    const p1 = rater1.filter(r => r === cat).length / n;
    const p2 = rater2.filter(r => r === cat).length / n;
    pe += p1 * p2;
  });

  // Cohen's Kappa
  const kappa = pe < 1 ? (po - pe) / (1 - pe) : 0;
  return isNaN(kappa) ? 0 : kappa;
}

/**
 * Calculate error propagation rate (sequential mode only)
 *
 * Measures: When first agent is wrong, how often do subsequent agents also get it wrong?
 */
function calculateErrorPropagation(experiments, groundTruth) {
  const sequentialExps = experiments.filter(e => e.sequential && e.agents && e.agents.length === 3);

  if (sequentialExps.length === 0) {
    return { errorPropagationRate: null };
  }

  let firstAgentWrong = 0;
  let subsequentAgentsAlsoWrong = 0;

  sequentialExps.forEach(exp => {
    const key = `${exp.dataset}_stage${exp.stage}`;
    const truth = groundTruth[key];

    if (!truth || !truth.decision) return;

    // Normalize all decisions using stage information
    const truthDecision = normalizeDecision(truth.decision, exp.stage);
    const agent1Decision = normalizeDecision(exp.agents[0].decision || '', exp.stage);
    const agent2Decision = normalizeDecision(exp.agents[1].decision || '', exp.stage);
    const agent3Decision = normalizeDecision(exp.agents[2].decision || '', exp.stage);

    if (!truthDecision || !agent1Decision || !agent2Decision || !agent3Decision) return;

    // Check if first agent (GPT-5.2) was wrong
    if (agent1Decision !== truthDecision) {
      firstAgentWrong++;

      // Check if subsequent agents (Gemini, Claude) also got it wrong
      if (agent2Decision !== truthDecision && agent3Decision !== truthDecision) {
        subsequentAgentsAlsoWrong++;
      }
    }
  });

  return {
    errorPropagationRate: firstAgentWrong > 0 ? subsequentAgentsAlsoWrong / firstAgentWrong : 0,
    firstAgentErrors: firstAgentWrong,
    propagatedErrors: subsequentAgentsAlsoWrong,
    sequentialDecisions: sequentialExps.length
  };
}

/**
 * Calculate error reduction (baseline vs multi-agent)
 */
function calculateErrorReduction(baseline, multiAgent, groundTruth) {
  const baselineAccuracy = calculateDecisionAccuracy(baseline, groundTruth);
  const multiAgentAccuracy = calculateDecisionAccuracy(multiAgent, groundTruth);

  const baselineErrors = baselineAccuracy.total_decisions - baselineAccuracy.correct_decisions;
  const multiAgentErrors = multiAgentAccuracy.total_decisions - multiAgentAccuracy.correct_decisions;

  const errorReduction = baselineErrors > 0 ? (baselineErrors - multiAgentErrors) / baselineErrors : 0;

  return {
    baselineErrors,
    multiAgentErrors,
    errorReduction,
    reductionPercent: (errorReduction * 100).toFixed(1) + '%'
  };
}

/**
 * Display metrics in readable format
 */
function displayMetrics(metrics) {
  console.log(`Experiments:          ${metrics.total_experiments}`);
  console.log(`Decision Accuracy:    ${(metrics.accuracy * 100).toFixed(1)}% (${metrics.correct_decisions}/${metrics.total_decisions})`);
  console.log(`Success Rate:         ${(metrics.successRate * 100).toFixed(1)}% (${metrics.successes}/${metrics.total_attempts})`);
  console.log(`Total Cost:           $${metrics.totalCost}`);
  console.log(`Cost per Success:     $${metrics.costPerSuccess}`);
  console.log(`Autonomy Rate:        ${(metrics.autonomyRate * 100).toFixed(1)}%`);
  console.log(`User Input Required:  ${metrics.userInputCount} decisions (${(metrics.userInputFrequency * 100).toFixed(1)}%)`);

  if (metrics.cohensKappa !== undefined && metrics.cohensKappa !== null) {
    console.log(`Cohen's Kappa:        ${metrics.cohensKappa}`);
    console.log(`  Agent1-Agent2:      ${metrics.kappa_agent1_agent2}`);
    console.log(`  Agent1-Agent3:      ${metrics.kappa_agent1_agent3}`);
    console.log(`  Agent2-Agent3:      ${metrics.kappa_agent2_agent3}`);
  }

  if (metrics.errorPropagationRate !== undefined && metrics.errorPropagationRate !== null) {
    console.log(`Error Propagation:    ${(metrics.errorPropagationRate * 100).toFixed(1)}%`);
    console.log(`  First Agent Errors: ${metrics.firstAgentErrors}`);
    console.log(`  Propagated Errors:  ${metrics.propagatedErrors}`);
  }

  console.log('');
}

/**
 * Export metrics to CSV
 */
function exportMetricsToCSV(allMetrics, outputPath) {
  const headers = [
    'system',
    'total_experiments',
    'accuracy',
    'correct_decisions',
    'total_decisions',
    'success_rate',
    'successes',
    'total_cost',
    'cost_per_success',
    'autonomy_rate',
    'user_input_frequency',
    'cohens_kappa',
    'error_propagation_rate'
  ];

  const rows = [headers.join(',')];

  allMetrics.forEach(m => {
    rows.push([
      m.system,
      m.total_experiments,
      m.accuracy.toFixed(3),
      m.correct_decisions,
      m.total_decisions,
      m.successRate.toFixed(3),
      m.successes,
      m.totalCost,
      m.costPerSuccess,
      m.autonomyRate.toFixed(3),
      m.userInputFrequency.toFixed(3),
      m.cohensKappa || 'N/A',
      m.errorPropagationRate !== undefined && m.errorPropagationRate !== null
        ? m.errorPropagationRate.toFixed(3)
        : 'N/A'
    ].join(','));
  });

  fs.writeFileSync(outputPath, rows.join('\n'));
}
