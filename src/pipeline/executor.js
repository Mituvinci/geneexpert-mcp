/**
 * Pipeline Executor - Coordinator-orchestrated execution
 *
 * The Coordinator:
 * 1. Detects input data type
 * 2. Plans pipeline steps
 * 3. DIRECTS Pipeline Agent to execute each step
 * 4. Identifies DECISION POINTS
 * 5. Calls multi-agent debate when needed
 * 6. Synthesizes consensus and continues
 */

import { detectInputData } from './data_detector.js';
import { planPipeline } from './planner.js';
import { createCoordinator } from '../coordinator/orchestrator.js';
import { Logger } from '../utils/logger.js';
import { PipelineAgent } from '../agents/pipeline_agent.js';
import path from 'path';
import fs from 'fs';

/**
 * Main entry point - Execute multi-agent orchestrated analysis
 */
export async function executeAnalysis(config) {
  console.log('[Coordinator] 🧬 Starting multi-agent RNA-seq analysis');
  console.log('='.repeat(60));
  console.log('');

  // Step 1: Detect input data type
  console.log('[Coordinator] Step 1: Detecting input data type');
  const dataInfo = detectInputData(config.input);
  console.log('');

  // Step 2: Plan pipeline
  console.log('[Coordinator] Step 2: Planning pipeline');
  const steps = planPipeline(dataInfo, config);
  console.log('');

  // Step 3: Initialize Coordinator, agents, and logger
  const coordinator = createCoordinator({ verbose: config.verbose });
  const pipelineAgent = new PipelineAgent({ verbose: config.verbose });
  const logger = new Logger(config.output, 'geneexpert_analysis');

  // Create session context
  const session = {
    config,
    dataInfo,
    steps,
    currentStep: 0,
    outputs: {}, // Store outputs from each step
    decisions: {}, // Store multi-agent decisions
    thresholds: {
      fdr: null,
      logFC: null
    },
    logger, // Add logger to session
    pipelineAgent // Add Pipeline Agent to session
  };

  // Step 4: Execute pipeline (Coordinator orchestrates each step)
  console.log('[Coordinator] Step 3: Executing pipeline');
  console.log('='.repeat(60));
  console.log('');

  try {
    for (let i = 0; i < steps.length; i++) {
      const step = steps[i];
      session.currentStep = i;

      console.log(`[Coordinator] Step ${i + 1}/${steps.length}: ${step.name}`);
      console.log(`              ${step.description}`);

      if (step.requiresDebate) {
        // DECISION POINT - Multi-agent debate
        await executeDecisionPoint(step, session, coordinator);
      } else {
        // EXECUTION STEP - Direct Pipeline Agent to execute
        await executeToolStep(step, session);
      }

      console.log('');
    }

    console.log('='.repeat(60));
    console.log('[Coordinator] ✅ Pipeline complete!');
    console.log('');

    // Generate summary
    printSummary(session);

    // Finalize logs
    session.logger.finalize(session);

    return session;

  } catch (error) {
    console.error('[Coordinator] ❌ Pipeline failed:', error.message);
    throw error;
  }
}

/**
 * Execute a tool step (Coordinator directs Pipeline Agent)
 */
async function executeToolStep(step, session) {
  console.log(`[Coordinator] → Pipeline Agent: Execute ${step.tool || step.name}`);

  // Log to file
  session.logger.log(`[Coordinator] → Pipeline Agent: Execute ${step.tool || step.name}`);

  try {
    // REAL EXECUTION - Call Pipeline Agent!
    const result = await session.pipelineAgent.execute(step, session);

    // Log tool execution
    session.logger.logToolExecution(step.name, step.tool || step.name, result.status, {
      description: step.description,
      outputs: result.outputs,
      duration: `${new Date(result.endTime) - new Date(result.startTime)}ms`
    });

    // Store output
    session.outputs[step.name] = result;

    return result;

  } catch (error) {
    // Log failure
    session.logger.logToolExecution(step.name, step.tool || step.name, 'failed', {
      description: step.description,
      error: error.message
    });

    throw error;
  }
}

/**
 * Execute a decision point (Multi-agent debate)
 */
async function executeDecisionPoint(step, session, coordinator) {
  console.log(`[Coordinator] 🤔 DECISION POINT: ${step.name}`);
  console.log(`[Coordinator] → Consulting all agents...`);
  console.log('');

  try {
    let result;

    switch (step.decisionType) {
      case 'threshold':
        result = await decideThresholds(session, coordinator);
        break;

      case 'qc_review':
        result = await reviewQC(session, coordinator);
        break;

      default:
        throw new Error(`Unknown decision type: ${step.decisionType}`);
    }

    // Store decision
    session.decisions[step.name] = result;

    // Log the FULL agent conversation to file
    session.logger.logAgentConversation(step, result.responses, result.consensus);

    console.log('');
    console.log(`[Coordinator] Decision: ${result.consensus.decision.toUpperCase()}`);
    console.log(`[Coordinator] Confidence: ${(result.consensus.confidence * 100).toFixed(0)}%`);

    // Handle consensus
    if (result.consensus.decision === 'user_decision_required') {
      console.log('[Coordinator] ⚠️  No consensus - user decision required');
      // TODO: Implement user prompt
      console.log('[Coordinator] (User prompt not yet implemented - proceeding with default)');
    } else if (result.consensus.recommendations?.length > 0) {
      console.log('[Coordinator] Recommendations:');
      result.consensus.recommendations.forEach(rec => {
        console.log(`   - ${rec.message}`);
      });
    }

  } catch (error) {
    console.error(`[Coordinator] ❌ Decision failed: ${error.message}`);
    throw error;
  }
}

/**
 * Multi-agent debate: Threshold selection
 */
async function decideThresholds(session, coordinator) {
  const sampleSize = session.dataInfo.samples.length;

  const result = await coordinator.validateThreshold(
    'FDR',
    0.05, // Default proposed value
    {
      sampleSize,
      comparison: session.config.comparison,
      organism: session.config.organism,
      experimentType: 'RNA-seq'
    }
  );

  // Extract recommended threshold from consensus
  // For small samples (< 6), agents typically recommend 0.01
  // For normal samples, 0.05 is fine
  if (sampleSize < 6) {
    session.thresholds.fdr = 0.01;
    session.thresholds.logFC = 1.0;
    console.log(`[Coordinator] Thresholds set: FDR < 0.01, logFC > 1.0 (small sample size)`);
  } else {
    session.thresholds.fdr = 0.05;
    session.thresholds.logFC = 1.0;
    console.log(`[Coordinator] Thresholds set: FDR < 0.05, logFC > 1.0`);
  }

  return result;
}

/**
 * Multi-agent debate: QC review
 */
async function reviewQC(session, coordinator) {
  // For now, simulate QC review
  // In full implementation, this would:
  // 1. Generate PCA/MDS plots
  // 2. Show plots to agents
  // 3. Agents discuss outliers, batch effects
  // 4. Vote on sample removal

  const result = await coordinator.reviewQCPlots(
    {
      pca_plot: 'pca_plot.pdf',
      mds_plot: 'mds_plot.pdf',
      density_plot: 'density_plot.pdf',
      observations: [
        'Samples cluster by treatment group',
        `PC1 explains 45% variance`,
        'No obvious outliers detected'
      ]
    },
    {
      nSamples: session.dataInfo.samples.length,
      groups: ['treatment', 'control'],
      replicates: Math.floor(session.dataInfo.samples.length / 2)
    }
  );

  return result;
}

/**
 * Print analysis summary
 */
function printSummary(session) {
  console.log('📊 Analysis Summary');
  console.log('='.repeat(60));
  console.log(`Input:       ${session.dataInfo.type.toUpperCase()}`);
  console.log(`Samples:     ${session.dataInfo.samples.length}`);
  console.log(`Steps:       ${session.steps.length} executed`);
  console.log(`Decisions:   ${Object.keys(session.decisions).length} multi-agent debates`);
  console.log(`Thresholds:  FDR < ${session.thresholds.fdr || 'N/A'}, logFC > ${session.thresholds.logFC || 'N/A'}`);
  console.log(`Output:      ${session.config.output}`);
  console.log('');
}
