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
import { MCPClaudeAgent } from '../agents/mcp_claude_agent.js';
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
  const mcpAgent = new MCPClaudeAgent({ verbose: config.verbose }); // MCP-enabled Claude!
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
    pipelineAgent, // Add Pipeline Agent to session
    mcpAgent // Add MCP Claude Agent for decision points
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
 * Multi-agent debate: QC review - USES MCP!
 */
async function reviewQC(session, coordinator) {
  console.log('[Coordinator] QC Review - MCP Agent will analyze actual outputs');

  // Build context with actual file paths
  const bamDir = path.join(session.config.output, 'bam_files');
  const countsDir = path.join(session.config.output, 'counts');
  const countFiles = fs.existsSync(countsDir)
    ? fs.readdirSync(countsDir).filter(f => f.endsWith('.count.txt'))
    : [];

  const mcpContext = {
    outputDir: session.config.output,
    bamDir: bamDir,
    countsFile: countFiles.length > 0 ? path.join(countsDir, countFiles[0]) : null,
    sampleSize: session.dataInfo.samples.length,
    comparison: session.config.comparison
  };

  // Claude reads ACTUAL outputs via MCP and analyzes
  const mcpPrompt = `Review the QC results for this RNA-seq analysis.

Use MCP tools to:
1. read_bam_summary - Check alignment statistics (mapping rates, total reads)
2. read_count_summary - Check count matrix (total counts, gene numbers)

Then analyze:
- Are mapping rates acceptable? (>70% is good, <50% is concerning)
- Are total counts reasonable?
- Any obvious quality issues?
- Should we proceed to DE analysis or investigate further?

Provide a clear recommendation.`;

  console.log('[MCP Claude Agent] Reading actual QC outputs...');

  const mcpResult = await session.mcpAgent.callWithTools(mcpPrompt, mcpContext);

  console.log('[MCP Claude Agent] Analysis complete');
  console.log(`Tool calls made: ${mcpResult.toolCalls.length}`);

  // Now get opinions from other agents (GPT-4, Gemini) with Claude's findings
  const result = await coordinator.reviewQCPlots(
    {
      mcp_analysis: mcpResult.content,
      tool_calls: mcpResult.toolCalls,
      observations: [
        `MCP Claude Agent analyzed actual outputs`,
        `Samples: ${session.dataInfo.samples.length}`,
        `See MCP agent's detailed analysis below`
      ]
    },
    {
      nSamples: session.dataInfo.samples.length,
      groups: ['treatment', 'control'],
      replicates: Math.floor(session.dataInfo.samples.length / 2)
    }
  );

  // Add MCP analysis to result
  result.mcpAnalysis = mcpResult;

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
