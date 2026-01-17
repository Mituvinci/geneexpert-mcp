/**
 * Metrics Recording for ICML Experiments
 *
 * Records experiment results for comparing multi-agent vs single-agent systems
 * Tracks: decision accuracy, error reduction, success rate, cost, time, inter-agent agreement
 */

import fs from 'fs';
import path from 'path';
import { globSync } from 'glob';

/**
 * Record a single experiment result
 *
 * @param {Object} result - Experiment result
 * @param {string} result.dataset - Dataset name (e.g., "DA0036")
 * @param {string} result.system - System used: "multi-agent", "gpt5.2", "claude", "gemini", "no-agent"
 * @param {string} result.decision - Agent decision: "automation" or "adaptation"
 * @param {boolean} result.success - Whether analysis succeeded
 * @param {Array} result.errors - List of errors encountered
 * @param {number} result.cost - Cost in USD
 * @param {number} result.time - Time in seconds
 * @param {Object} result.agentResponses - Individual agent responses (for multi-agent)
 * @param {number} result.attempts - Number of retry attempts
 */
export function recordExperiment(result) {
  const resultsDir = path.join(process.cwd(), 'experiments', 'results');

  // Ensure results directory exists
  if (!fs.existsSync(resultsDir)) {
    fs.mkdirSync(resultsDir, { recursive: true });
  }

  const timestamp = new Date().toISOString();
  const resultFile = path.join(resultsDir, `${result.dataset}_${result.system}_${Date.now()}.json`);

  const experimentRecord = {
    timestamp,
    dataset: result.dataset,
    system: result.system,
    decision: result.decision,
    success: result.success,
    errors: result.errors || [],
    cost: result.cost || 0,
    time: result.time || 0,
    attempts: result.attempts || 1,
    agentResponses: result.agentResponses || null,
    outputDir: result.outputDir || null,
    scriptPath: result.scriptPath || null
  };

  // Write individual result
  fs.writeFileSync(resultFile, JSON.stringify(experimentRecord, null, 2));

  // Append to consolidated results file
  const consolidatedFile = path.join(resultsDir, 'all_experiments.jsonl');
  fs.appendFileSync(consolidatedFile, JSON.stringify(experimentRecord) + '\n');

  console.log(`[Metrics] Experiment recorded: ${resultFile}`);

  return experimentRecord;
}

/**
 * Load all experiment decisions from JSON files
 * Reads from: experiments/results/DATASET_NAME/geneexpert_analysis_agent_decisions.json
 *
 * @returns {Array} - Array of all decisions from all experiments
 */
export function loadAllDecisions() {
  const experimentsDir = path.join(process.cwd(), 'experiments', 'results');

  if (!fs.existsSync(experimentsDir)) {
    return [];
  }

  // Find all JSON decision files
  const pattern = path.join(experimentsDir, '**/geneexpert_analysis_agent_decisions.json');
  const files = globSync(pattern);

  const allDecisions = [];

  files.forEach(file => {
    try {
      const data = JSON.parse(fs.readFileSync(file, 'utf-8'));

      // Extract system type from path or data
      const system = data.system || 'unknown';
      const sessionId = data.session_id || path.basename(path.dirname(file));

      // Add each decision with session metadata
      data.decisions.forEach(decision => {
        allDecisions.push({
          ...decision,
          system,
          session_id: sessionId,
          file_path: file
        });
      });
    } catch (error) {
      console.error(`Error reading ${file}:`, error.message);
    }
  });

  return allDecisions;
}

/**
 * Calculate decision accuracy for a system
 *
 * @param {string} system - System name (e.g., 'multi-agent', 'single-agent-claude')
 * @param {Object} groundTruth - Ground truth labels from experiments/ground_truth.json
 * @returns {Object} - Accuracy metrics
 */
export function calculateDecisionAccuracy(system, groundTruth) {
  const allDecisions = loadAllDecisions();

  if (allDecisions.length === 0) {
    return { accuracy: 0, correct: 0, total: 0, message: 'No experiments recorded yet' };
  }

  // Filter by system
  const systemDecisions = allDecisions.filter(d => d.system === system);

  let correct = 0;
  let total = 0;

  systemDecisions.forEach(decision => {
    const truth = groundTruth[decision.decision_id];
    if (truth && truth.expected_decision) {
      total++;
      if (decision.consensus.decision.toLowerCase() === truth.expected_decision.toLowerCase()) {
        correct++;
      }
    }
  });

  return {
    system,
    accuracy: total > 0 ? correct / total : 0,
    correct,
    total,
    total_decisions: systemDecisions.length
  };
}

/**
 * Calculate error reduction: (baseline_errors - system_errors) / baseline_errors
 *
 * @param {string} baselineSystem - Baseline system (e.g., "no-agent")
 * @param {string} comparisonSystem - System to compare (e.g., "multi-agent")
 * @returns {Object} - Error reduction metrics
 */
export function calculateErrorReduction(baselineSystem, comparisonSystem) {
  const resultsFile = path.join(process.cwd(), 'experiments', 'results', 'all_experiments.jsonl');

  if (!fs.existsSync(resultsFile)) {
    return { errorReduction: 0, message: 'No experiments recorded yet' };
  }

  const lines = fs.readFileSync(resultsFile, 'utf-8').split('\n').filter(l => l.trim());
  const allExperiments = lines.map(line => JSON.parse(line));

  const baselineExps = allExperiments.filter(exp => exp.system === baselineSystem);
  const comparisonExps = allExperiments.filter(exp => exp.system === comparisonSystem);

  const baselineErrors = baselineExps.reduce((sum, exp) => sum + exp.errors.length, 0);
  const comparisonErrors = comparisonExps.reduce((sum, exp) => sum + exp.errors.length, 0);

  const errorReduction = baselineErrors > 0 ? (baselineErrors - comparisonErrors) / baselineErrors : 0;

  return {
    baselineSystem,
    comparisonSystem,
    baselineErrors,
    comparisonErrors,
    errorReduction,
    reductionPercent: (errorReduction * 100).toFixed(1) + '%'
  };
}

/**
 * Calculate success rate for a system
 *
 * @param {string} system - System name
 * @returns {Object} - Success rate metrics
 */
export function calculateSuccessRate(system) {
  const resultsFile = path.join(process.cwd(), 'experiments', 'results', 'all_experiments.jsonl');

  if (!fs.existsSync(resultsFile)) {
    return { successRate: 0, successes: 0, total: 0, message: 'No experiments recorded yet' };
  }

  const lines = fs.readFileSync(resultsFile, 'utf-8').split('\n').filter(l => l.trim());
  const experiments = lines.map(line => JSON.parse(line)).filter(exp => exp.system === system);

  const successes = experiments.filter(exp => exp.success).length;
  const total = experiments.length;

  return {
    system,
    successRate: total > 0 ? successes / total : 0,
    successPercent: total > 0 ? ((successes / total) * 100).toFixed(1) + '%' : '0%',
    successes,
    total
  };
}

/**
 * Calculate cost per successful analysis from session metadata
 *
 * @param {string} system - System name
 * @returns {Object} - Cost efficiency metrics
 */
export function calculateCostEfficiency(system) {
  const experimentsDir = path.join(process.cwd(), 'experiments', 'results');

  if (!fs.existsSync(experimentsDir)) {
    return { costPerSuccess: 0, message: 'No experiments recorded yet' };
  }

  // Find all session metadata files
  const pattern = path.join(experimentsDir, '**/geneexpert_analysis_session_metadata.json');
  const files = globSync(pattern);

  const sessions = files
    .map(file => {
      try {
        return JSON.parse(fs.readFileSync(file, 'utf-8'));
      } catch {
        return null;
      }
    })
    .filter(s => s && s.system === system);

  if (sessions.length === 0) {
    return { costPerSuccess: 0, message: 'No experiments found for this system' };
  }

  const successfulSessions = sessions.filter(s => s.execution?.success === true);
  const totalCost = sessions.reduce((sum, s) => sum + (s.costs?.total_usd || 0), 0);
  const successCount = successfulSessions.length;

  return {
    system,
    totalCost: totalCost.toFixed(3),
    totalSessions: sessions.length,
    successCount,
    costPerSuccess: successCount > 0 ? (totalCost / successCount).toFixed(3) : 'N/A',
    averageCost: sessions.length > 0 ? (totalCost / sessions.length).toFixed(3) : '0.000'
  };
}

/**
 * Calculate Cohen's Kappa for inter-agent agreement
 * Measures agreement between agents beyond chance
 *
 * @param {string} dataset - Dataset name to analyze
 * @returns {Object} - Kappa statistics
 */
export function calculateInterAgentAgreement(dataset = null) {
  const resultsFile = path.join(process.cwd(), 'experiments', 'results', 'all_experiments.jsonl');

  if (!fs.existsSync(resultsFile)) {
    return { kappa: 0, message: 'No experiments recorded yet' };
  }

  const lines = fs.readFileSync(resultsFile, 'utf-8').split('\n').filter(l => l.trim());
  let experiments = lines.map(line => JSON.parse(line)).filter(exp => exp.system === 'multi-agent');

  if (dataset) {
    experiments = experiments.filter(exp => exp.dataset === dataset);
  }

  if (experiments.length === 0) {
    return { kappa: 0, message: 'No multi-agent experiments found' };
  }

  // Extract agent decisions from responses
  const agreements = [];
  experiments.forEach(exp => {
    if (exp.agentResponses) {
      const gpt5Decision = extractDecisionFromResponse(exp.agentResponses.gpt5_2?.content);
      const claudeDecision = extractDecisionFromResponse(exp.agentResponses.claude?.content);
      const geminiDecision = extractDecisionFromResponse(exp.agentResponses.gemini?.content);

      if (gpt5Decision && claudeDecision && geminiDecision) {
        agreements.push({
          gpt5: gpt5Decision,
          claude: claudeDecision,
          gemini: geminiDecision,
          consensus: exp.decision
        });
      }
    }
  });

  if (agreements.length === 0) {
    return { kappa: 0, message: 'No agent responses found to calculate agreement' };
  }

  // Calculate pairwise agreement
  const pairwiseKappas = [
    calculatePairwiseKappa(agreements.map(a => a.gpt5), agreements.map(a => a.claude)),
    calculatePairwiseKappa(agreements.map(a => a.gpt5), agreements.map(a => a.gemini)),
    calculatePairwiseKappa(agreements.map(a => a.claude), agreements.map(a => a.gemini))
  ];

  const averageKappa = pairwiseKappas.reduce((sum, k) => sum + k, 0) / pairwiseKappas.length;

  return {
    averageKappa: averageKappa.toFixed(3),
    pairwiseKappas: {
      'gpt5.2-claude': pairwiseKappas[0].toFixed(3),
      'gpt5.2-gemini': pairwiseKappas[1].toFixed(3),
      'claude-gemini': pairwiseKappas[2].toFixed(3)
    },
    dataPoints: agreements.length
  };
}

/**
 * Calculate Cohen's Kappa for two raters
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
  const kappa = (po - pe) / (1 - pe);
  return kappa;
}

/**
 * Extract decision from agent response text
 */
function extractDecisionFromResponse(content) {
  if (!content) return null;

  const lower = content.toLowerCase();

  // Look for explicit decision statements
  if (lower.includes('decision: automation') || lower.includes('recommend automation')) {
    return 'automation';
  }
  if (lower.includes('decision: adaptation') || lower.includes('recommend adaptation')) {
    return 'adaptation';
  }

  // Count mentions
  const automationCount = (lower.match(/automation/g) || []).length;
  const adaptationCount = (lower.match(/adaptation/g) || []).length;

  if (automationCount > adaptationCount) return 'automation';
  if (adaptationCount > automationCount) return 'adaptation';

  return null;
}

/**
 * Generate summary report for all systems
 */
export function generateSummaryReport() {
  const allDecisions = loadAllDecisions();

  if (allDecisions.length === 0) {
    console.log('\n⚠️  No experiment data found. Run some experiments first!');
    return;
  }

  // Get unique systems
  const systems = [...new Set(allDecisions.map(d => d.system))];

  console.log('\n' + '='.repeat(70));
  console.log('ICML EXPERIMENT SUMMARY REPORT');
  console.log('='.repeat(70) + '\n');

  console.log(`Total Experiments: ${allDecisions.length} decisions across ${systems.length} systems\n`);

  systems.forEach(system => {
    console.log(`\n📊 System: ${system.toUpperCase()}`);
    console.log('-'.repeat(70));

    const costEfficiency = calculateCostEfficiency(system);
    const systemDecisions = allDecisions.filter(d => d.system === system);

    console.log(`  Total Decisions:  ${systemDecisions.length}`);
    console.log(`  Total Cost:       $${costEfficiency.totalCost}`);
    console.log(`  Average Cost:     $${costEfficiency.averageCost}`);
    console.log(`  Sessions:         ${costEfficiency.totalSessions || 'N/A'}`);

    // Decision breakdown
    const automationCount = systemDecisions.filter(d => d.consensus.decision === 'automation').length;
    const adaptationCount = systemDecisions.filter(d => d.consensus.decision === 'adaptation').length;
    console.log(`  Decisions:        ${automationCount} automation, ${adaptationCount} adaptation`);

    // Average confidence
    const avgConfidence = systemDecisions.reduce((sum, d) => sum + d.consensus.confidence, 0) / systemDecisions.length;
    console.log(`  Avg Confidence:   ${(avgConfidence * 100).toFixed(1)}%`);
  });

  console.log('\n' + '='.repeat(70));
  console.log('INTER-AGENT AGREEMENT (Multi-Agent Only)');
  console.log('='.repeat(70) + '\n');

  const multiAgentDecisions = allDecisions.filter(d => d.system === 'multi-agent');
  if (multiAgentDecisions.length > 0) {
    const interAgentAgreement = calculateInterAgentAgreement();
    console.log(`Cohen's Kappa: ${interAgentAgreement.averageKappa || 'N/A'}`);
    if (interAgentAgreement.pairwiseKappas) {
      console.log(`  GPT-5.2 vs Claude:  ${interAgentAgreement.pairwiseKappas['gpt5.2-claude']}`);
      console.log(`  GPT-5.2 vs Gemini:  ${interAgentAgreement.pairwiseKappas['gpt5.2-gemini']}`);
      console.log(`  Claude vs Gemini:   ${interAgentAgreement.pairwiseKappas['claude-gemini']}`);
    }
    console.log(`  Data Points:        ${interAgentAgreement.dataPoints || 0}`);
  } else {
    console.log('No multi-agent experiments found.');
  }

  console.log('\n' + '='.repeat(70) + '\n');
}

/**
 * Export results to CSV for statistical analysis
 */
export function exportToCSV(outputPath = null) {
  const resultsFile = path.join(process.cwd(), 'experiments', 'results', 'all_experiments.jsonl');

  if (!fs.existsSync(resultsFile)) {
    console.error('[Metrics] No experiments recorded yet');
    return;
  }

  const lines = fs.readFileSync(resultsFile, 'utf-8').split('\n').filter(l => l.trim());
  const experiments = lines.map(line => JSON.parse(line));

  // CSV header
  const csv = [
    'timestamp,dataset,system,decision,success,num_errors,cost,time,attempts'
  ];

  // CSV rows
  experiments.forEach(exp => {
    csv.push([
      exp.timestamp,
      exp.dataset,
      exp.system,
      exp.decision,
      exp.success ? '1' : '0',
      exp.errors.length,
      exp.cost,
      exp.time,
      exp.attempts
    ].join(','));
  });

  const csvContent = csv.join('\n');
  const csvPath = outputPath || path.join(process.cwd(), 'experiments', 'results', 'experiments.csv');

  fs.writeFileSync(csvPath, csvContent);
  console.log(`[Metrics] CSV exported to: ${csvPath}`);

  return csvPath;
}
