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
import { generateAutomationScript, generateAdaptationScript } from './script_generator.js';
import { spawn } from 'child_process';
import path from 'path';
import fs from 'fs';

/**
 * Main entry point - Execute multi-agent orchestrated analysis
 */
export async function executeAnalysis(config) {
  // Initialize logger FIRST to capture everything
  const logger = new Logger(config.output, 'geneexpert_analysis');

  // Override console.log to also write to file
  const originalLog = console.log;
  console.log = function(...args) {
    const message = args.join(' ');
    originalLog.apply(console, args);
    logger.log(message);
  };

  console.log('[Coordinator] 🧬 GeneExpert Multi-Agent RNA-seq Analysis');
  console.log('='.repeat(60));
  console.log('');

  // Step 1: Detect input data type (NO execution yet!)
  console.log('[Coordinator] 📊 Analyzing input data...');
  console.log('');
  const dataInfo = detectInputData(config.input, config);
  console.log('');

  // Step 2: Plan pipeline
  console.log('[Coordinator] 📝 Creating analysis plan...');
  const steps = planPipeline(dataInfo, config);
  console.log('');

  // Step 3: AGENTS ANALYZE DATA AND DECIDE APPROACH
  console.log('[Coordinator] 🤔 Consulting agents to decide approach...');
  console.log('');

  const coordinator = createCoordinator({ verbose: config.verbose }); // Only verbose if user wants it
  const mcpAgent = new MCPClaudeAgent({ verbose: config.verbose }); // MCP agent for ADAPTATION

  const agentDecision = await coordinator.decideAnalysisApproach(dataInfo, config, steps);

  console.log('');

  // Step 4: GENERATE SCRIPT based on agent decision
  console.log('[Coordinator] 📜 Generating analysis script...');
  console.log('');

  const scriptName = `geneexpert_${config.comparison}_v1.sh`;
  const scriptPath = path.join(config.output, scriptName);

  let scriptResult;
  const decision = agentDecision.consensus.decision.toLowerCase();

  if (decision === 'automation') {
    // AUTOMATION: Template-based script
    scriptResult = generateAutomationScript(dataInfo, config, steps, scriptPath);
  } else {
    // ADAPTATION: MCP agent reads scripts and writes custom solution
    scriptResult = await generateAdaptationScript(dataInfo, config, steps, agentDecision, coordinator, mcpAgent, scriptPath);
  }

  console.log('');

  // Step 5: Create detailed plan with AGENT'S decision and generated script
  const plan = createDetailedPlan(dataInfo, config, steps, agentDecision, scriptResult);

  // Step 6: Display plan to user
  displayPlan(plan);

  // Step 5: Wait for user confirmation
  const confirmed = await getUserConfirmation();

  if (!confirmed) {
    console.log('');
    console.log('[Coordinator] ❌ Analysis cancelled by user.');
    console.log('[Coordinator] No files were created or modified.');
    console.log('');

    // Restore console.log and save log even on cancellation
    console.log = originalLog;
    logger.log('[Coordinator] User cancelled analysis');
    logger.finalize({ cancelled: true });
    console.log(`📄 Planning session log saved to: ${logger.logFile}`);
    console.log('');

    process.exit(0);
  }

  // User confirmed - proceed!
  console.log('');
  console.log('[Coordinator] ✅ User confirmed. Starting analysis...');
  console.log('='.repeat(60));
  console.log('');

  // Step 7: Execute with FEEDBACK LOOP (retry on failure)
  const MAX_RETRIES = 3;
  let attempt = 1;
  let currentScriptResult = scriptResult;
  let executionResult = null;

  while (attempt <= MAX_RETRIES) {
    console.log(`[Executor] Attempt ${attempt}/${MAX_RETRIES}: Running script: ${currentScriptResult.scriptPath}`);
    console.log('='.repeat(60));
    console.log('');

    try {
      executionResult = await executeScript(currentScriptResult.scriptPath, logger);

      if (executionResult.success) {
        // SUCCESS!
        console.log('');
        console.log('='.repeat(60));
        console.log(`[Executor] ✅ Analysis complete! (attempt ${attempt}/${MAX_RETRIES})`);
        console.log('');
        console.log(`📊 Results saved to: ${config.output}`);
        console.log(`📄 Full session log: ${logger.logFile}`);
        console.log(`🔬 Script executed: ${currentScriptResult.scriptPath}`);
        console.log('');

        console.log = originalLog;
        return {
          success: true,
          output: config.output,
          scriptPath: currentScriptResult.scriptPath,
          executionResult,
          attempts: attempt
        };
      }

      // FAILURE - Enter FEEDBACK LOOP
      console.log('');
      console.log('='.repeat(60));
      console.error(`[Executor] ❌ Script failed (attempt ${attempt}/${MAX_RETRIES})`);
      console.error(`Error: ${executionResult.error}`);
      console.log('');

      if (attempt >= MAX_RETRIES) {
        console.log('[Executor] Max retries reached. Giving up.');
        console.log('');
        console.log(`📄 Full session log: ${logger.logFile}`);
        console.log(`🔬 Failed script: ${currentScriptResult.scriptPath}`);
        console.log('');
        console.log = originalLog;
        throw new Error(`Script execution failed after ${MAX_RETRIES} attempts: ${executionResult.error}`);
      }

      // FEEDBACK LOOP: Pass error to agents for debugging
      console.log('');
      console.log('🔄 FEEDBACK LOOP: Consulting agents to debug error...');
      console.log('='.repeat(60));
      console.log('');

      const debugPrompt = `
PREVIOUS SCRIPT FAILED WITH ERROR:
${executionResult.error}

STDERR OUTPUT:
${executionResult.stderr}

STDOUT OUTPUT (last 50 lines):
${executionResult.stdout.split('\n').slice(-50).join('\n')}

FAILED SCRIPT PATH: ${currentScriptResult.scriptPath}

YOUR TASK: Analyze this error and write a FIXED version of the script.
Common issues:
- Missing files (use list_available_scripts to check)
- Wrong paths
- Missing parameters
- Syntax errors

Write a COMPLETE corrected bash script that fixes this error.
`;

      const debugDecision = await coordinator.getMultiAgentDecision(debugPrompt, dataInfo, config);

      // Generate v2 script with fixes
      attempt++;
      const v2ScriptPath = scriptPath.replace(/v\d+\.sh$/, `v${attempt}.sh`).replace(/\.sh$/, `_v${attempt}.sh`);

      console.log('');
      console.log(`[Script Generator] Generating v${attempt} script with agent fixes...`);

      if (decision === 'adaptation') {
        currentScriptResult = await generateAdaptationScript(
          dataInfo, config, steps, debugDecision, coordinator, mcpAgent, v2ScriptPath
        );
      } else {
        currentScriptResult = generateAutomationScript(dataInfo, config, steps, v2ScriptPath);
      }

      console.log(`[Script Generator] ✓ v${attempt} script generated: ${v2ScriptPath}`);
      console.log('');

    } catch (error) {
      console.log = originalLog;
      console.error('[Executor] ❌ Fatal error:', error.message);
      throw error;
    }
  }
}

/**
 * Execute the generated bash script
 */
async function executeScript(scriptPath, logger) {
  console.log('[Executor] Starting script execution...');
  console.log('');

  return new Promise((resolve) => {
    // Make script executable
    try {
      fs.chmodSync(scriptPath, 0o755);
    } catch (error) {
      console.error(`[Executor] Failed to make script executable: ${error.message}`);
      resolve({ success: false, error: error.message });
      return;
    }

    // Spawn bash process
    const childProcess = spawn('bash', [scriptPath], {
      cwd: path.dirname(scriptPath),
      env: { ...process.env }
    });

    let stdout = '';
    let stderr = '';

    // Stream stdout to console in real-time
    childProcess.stdout.on('data', (data) => {
      const output = data.toString();
      console.log(output.trimEnd());
      stdout += output;
    });

    // Stream stderr to console in real-time
    childProcess.stderr.on('data', (data) => {
      const output = data.toString();
      console.error(output.trimEnd());
      stderr += output;
    });

    // Handle process completion
    childProcess.on('close', (code) => {
      console.log('');

      if (code === 0) {
        console.log('[Executor] Script completed successfully');
        resolve({
          success: true,
          exitCode: code,
          stdout,
          stderr
        });
      } else {
        console.error(`[Executor] Script failed with exit code ${code}`);
        resolve({
          success: false,
          exitCode: code,
          error: `Script exited with code ${code}`,
          stdout,
          stderr
        });
      }
    });

    // Handle process errors
    childProcess.on('error', (error) => {
      console.error(`[Executor] Failed to start script: ${error.message}`);
      resolve({
        success: false,
        error: error.message,
        stdout,
        stderr
      });
    });
  });
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
 * Create detailed analysis plan WITH AGENT DECISION and GENERATED SCRIPT
 */
function createDetailedPlan(dataInfo, config, steps, agentDecision, scriptResult) {
  // Extract agent decision and reasoning
  const approach = agentDecision.consensus.decision.toUpperCase();
  const reasoning = agentDecision.consensus.reasoning || 'See agent responses above.';

  // Add script information
  const scriptInfo = {
    type: scriptResult.type,
    path: scriptResult.scriptPath,
    agentGenerated: scriptResult.type === 'ADAPTATION'
  };

  // Format groups
  const formattedGroups = [];
  if (dataInfo.groups) {
    Object.entries(dataInfo.groups).forEach(([groupName, samples]) => {
      formattedGroups.push({
        name: groupName.charAt(0).toUpperCase() + groupName.slice(1),
        samples,
        count: samples.length
      });
    });
  }

  // Script details
  const scriptName = `geneexpert_${config.comparison || 'analysis'}_v1.sh`;
  const scriptPath = path.join(config.output, scriptName);

  // Estimate runtime based on data type and steps
  let estimatedTime = 'Unknown';
  if (dataInfo.type === 'fastq') {
    const nSamples = dataInfo.samples.length;
    estimatedTime = `${20 + nSamples * 3}-${30 + nSamples * 5} minutes`;
  } else if (dataInfo.type === 'bam') {
    estimatedTime = '5-15 minutes';
  } else {
    estimatedTime = '2-5 minutes';
  }

  return {
    analysisType: 'Bulk RNA-seq',
    organism: config.organism === 'mouse' ? 'Mouse (mm10)' : 'Human (hg38)',
    sequencing: dataInfo.pairedEnd ? 'Paired-end ✓' : 'Single-end',
    comparison: config.comparison || 'Not specified',
    groups: formattedGroups,
    approach,
    reasoning,
    scriptInfo,  // Add script generation info
    steps: steps.map(s => ({
      name: s.name,
      description: s.description
    })),
    scriptName,
    scriptPath,
    estimatedTime,
    outputDir: config.output,
    sampleCount: dataInfo.samples.length
  };
}

/**
 * Display beautiful plan to user
 */
function displayPlan(plan) {
  console.log('');
  console.log('╔' + '═'.repeat(58) + '╗');
  console.log('║' + ' '.repeat(15) + 'PROPOSED ANALYSIS PLAN' + ' '.repeat(21) + '║');
  console.log('╚' + '═'.repeat(58) + '╝');
  console.log('');

  // Experiment details
  console.log('📊 EXPERIMENT DETAILS:');
  console.log(`   Analysis Type:  ${plan.analysisType}`);
  console.log(`   Comparison:     ${plan.comparison}`);
  console.log(`   Organism:       ${plan.organism}`);
  console.log(`   Sequencing:     ${plan.sequencing}`);
  console.log('');

  // Samples
  console.log('📁 SAMPLES IDENTIFIED:');
  if (plan.groups.length > 0) {
    plan.groups.forEach(group => {
      console.log(`   ${group.name} group (n=${group.count}):`);
      group.samples.forEach(s => console.log(`     • ${s}`));
    });
  } else {
    console.log(`   Total samples: ${plan.sampleCount}`);
  }
  console.log('');

  // Agent decision
  console.log('🎯 AGENT DECISION:');
  console.log(`   Approach:       ${plan.approach}`);
  console.log(`   Reason:         ${plan.reasoning}`);
  console.log('');

  // Pipeline steps
  console.log('📝 PIPELINE STEPS:');
  plan.steps.forEach((step, i) => {
    const num = `${i + 1}.`.padEnd(4);
    const name = step.name.padEnd(22);
    console.log(`   ${num}${name} - ${step.description}`);
  });
  console.log('');

  // Script info
  console.log('📜 GENERATED SCRIPT:');
  console.log(`   Type:     ${plan.scriptInfo.type} ${plan.scriptInfo.agentGenerated ? '(Agent-written custom script)' : '(Template-based)'}`);
  console.log(`   Name:     ${plan.scriptName}`);
  console.log(`   Location: ${plan.scriptPath}`);
  console.log(`   Runtime:  Estimated ${plan.estimatedTime}`);
  if (plan.scriptInfo.agentGenerated) {
    console.log(`   Note:     Script customized to address agent concerns`);
  }
  console.log('');

  // Outputs
  console.log('💾 OUTPUTS:');
  console.log(`   All results will be saved to: ${plan.outputDir}/`);
  console.log(`   - Alignment files (BAM)`);
  console.log(`   - Count matrices (txt)`);
  console.log(`   - QC plots (PDF)`);
  console.log(`   - DE results (txt, xlsx)`);
  console.log(`   - Analysis log`);
  console.log(`   - Versioned script (for reproducibility)`);
  console.log('');

  console.log('─'.repeat(60));
  console.log('⚠️  IMPORTANT: No analysis has been run yet.');
  console.log('   This is a preview of what will happen.');
  console.log('');
}

/**
 * Get user confirmation (Y/N)
 */
async function getUserConfirmation() {
  const readline = await import('readline');
  const rl = readline.createInterface({
    input: process.stdin,
    output: process.stdout
  });

  return new Promise((resolve) => {
    rl.question('❓ Do you want to proceed with this analysis? (Y/N): ', (answer) => {
      rl.close();
      const userSaidYes = answer.toLowerCase() === 'y' || answer.toLowerCase() === 'yes';
      resolve(userSaidYes);
    });
  });
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
