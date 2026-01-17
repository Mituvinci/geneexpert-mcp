/**
 * Logger - Save all agent interactions and conversations to log file
 * Supports both plain text (human-readable) and JSON (ICML metrics)
 */

import fs from 'fs';
import path from 'path';
import { calculateDecisionCost, calculateAgentCost } from './cost_calculator.js';
import { parseAgentResponse, getConfidenceLabel } from './response_parser.js';

class Logger {
  constructor(outputDir, sessionName = 'analysis', config = {}) {
    this.outputDir = outputDir;
    this.sessionName = sessionName;
    this.config = config;

    // Plain text logs (keep existing for debugging)
    this.logFile = path.join(outputDir, `${sessionName}_log.txt`);
    this.conversationFile = path.join(outputDir, `${sessionName}_agent_conversations.txt`);

    // JSON logs (NEW - for ICML experiments)
    this.sessionMetadataFile = path.join(outputDir, `${sessionName}_session_metadata.json`);
    this.decisionsJsonFile = path.join(outputDir, `${sessionName}_agent_decisions.json`);
    this.decisionsJsonlFile = path.join(outputDir, `${sessionName}_agent_decisions.jsonl`);

    // In-memory storage for JSON export
    this.sessionData = {
      session_id: `${sessionName}_${new Date().toISOString().replace(/[:.]/g, '-')}`,
      timestamp_start: new Date().toISOString(),
      timestamp_end: null,
      duration_seconds: null,
      system: config.singleAgent ? `single-agent-${config.singleAgent}` : 'multi-agent',
      config: config || {},
      decisions: [],
      costs: {
        total_usd: 0,
        breakdown: {
          gpt5_2: 0,
          claude: 0,
          gemini: 0
        }
      }
    };

    // Create output directory if needed
    if (!fs.existsSync(outputDir)) {
      fs.mkdirSync(outputDir, { recursive: true });
    }

    // Initialize log files
    this.initLogFiles();
  }

  initLogFiles() {
    const timestamp = new Date().toISOString();

    // Main log file (explicit UTF-8 encoding)
    fs.writeFileSync(this.logFile, `GeneExpert Multi-Agent Analysis Log\n`, 'utf8');
    fs.appendFileSync(this.logFile, `Session: ${this.sessionName}\n`, 'utf8');
    fs.appendFileSync(this.logFile, `Started: ${timestamp}\n`, 'utf8');
    fs.appendFileSync(this.logFile, `${'='.repeat(80)}\n\n`, 'utf8');

    // Agent conversations file (explicit UTF-8 encoding)
    fs.writeFileSync(this.conversationFile, `GeneExpert Agent Conversations\n`, 'utf8');
    fs.appendFileSync(this.conversationFile, `Session: ${this.sessionName}\n`, 'utf8');
    fs.appendFileSync(this.conversationFile, `Started: ${timestamp}\n`, 'utf8');
    fs.appendFileSync(this.conversationFile, `${'='.repeat(80)}\n\n`, 'utf8');

    console.log(`Logging to:`);
    console.log(`   Main log: ${this.logFile}`);
    console.log(`   Agent conversations: ${this.conversationFile}`);
    console.log('');
  }

  /**
   * Log general pipeline events
   */
  log(message, level = 'INFO') {
    const timestamp = new Date().toISOString();
    const logLine = `[${timestamp}] [${level}] ${message}\n`;

    fs.appendFileSync(this.logFile, logLine, 'utf8');

    // Don't print to console - console.log is already overridden in executor.js
  }

  /**
   * Log multi-agent conversation (the full thinking!)
   * Writes to BOTH plain text (human-readable) and JSON (metrics)
   */
  logAgentConversation(step, agents, consensus, userPrompt = '') {
    const timestamp = new Date().toISOString();

    // 1. WRITE PLAIN TEXT (existing behavior - keep for debugging)
    this.writePlainTextConversation(step, agents, consensus, timestamp);

    // 2. WRITE JSON (NEW - for ICML experiments)
    this.writeJsonConversation(step, agents, consensus, timestamp, userPrompt);
  }

  /**
   * Write plain text conversation (existing format)
   */
  writePlainTextConversation(step, agents, consensus, timestamp) {
    let conversation = `\n${'='.repeat(80)}\n`;
    conversation += `DECISION POINT: ${step.name}\n`;
    conversation += `Timestamp: ${timestamp}\n`;
    conversation += `Description: ${step.description}\n`;
    conversation += `${'='.repeat(80)}\n\n`;

    // Log each agent's full response (using ASCII box characters)
    if (agents.gpt5_2) {
      conversation += `+-- Stats Agent (GPT-5.2) ${'-'.repeat(39)}+\n`;
      conversation += `| Model: ${agents.gpt5_2.model}\n`;
      conversation += `| Status: ${agents.gpt5_2.success ? 'Success' : 'Failed'}\n`;
      conversation += `+${'-'.repeat(64)}+\n\n`;
      conversation += `${agents.gpt5_2.content || 'No response'}\n\n`;
    }

    if (agents.claude) {
      conversation += `+-- Pipeline Agent (Claude) ${'-'.repeat(37)}+\n`;
      conversation += `| Model: ${agents.claude.model}\n`;
      conversation += `| Status: ${agents.claude.success ? 'Success' : 'Failed'}\n`;
      conversation += `+${'-'.repeat(64)}+\n\n`;
      conversation += `${agents.claude.content || 'No response'}\n\n`;
    }

    if (agents.gemini) {
      conversation += `+-- Biology Agent (Gemini) ${'-'.repeat(38)}+\n`;
      conversation += `| Model: ${agents.gemini.model}\n`;
      conversation += `| Status: ${agents.gemini.success ? 'Success' : 'Failed'}\n`;
      conversation += `+${'-'.repeat(64)}+\n\n`;
      conversation += `${agents.gemini.content || 'No response'}\n\n`;
    }

    // Log consensus decision
    conversation += `+-- CONSENSUS DECISION ${'-'.repeat(42)}+\n`;
    conversation += `| Decision: ${consensus.decision.toUpperCase()}\n`;
    conversation += `| Confidence: ${(consensus.confidence * 100).toFixed(0)}%\n`;
    conversation += `| Votes:\n`;
    // Handle both approach decisions (automation/adaptation) and standard decisions (approve/reject)
    if (consensus.votes?.automation !== undefined || consensus.votes?.adaptation !== undefined) {
      conversation += `|   Automation: ${consensus.votes?.automation || 0}/3\n`;
      conversation += `|   Adaptation: ${consensus.votes?.adaptation || 0}/3\n`;
    } else {
      conversation += `|   Approve:   ${consensus.votes?.approve || 0}/3\n`;
      conversation += `|   Reject:    ${consensus.votes?.reject || 0}/3\n`;
      conversation += `|   Uncertain: ${consensus.votes?.uncertain || 0}/3\n`;
    }
    conversation += `+${'-'.repeat(64)}+\n\n`;

    if (consensus.recommendations?.length > 0) {
      conversation += `Recommendations:\n`;
      consensus.recommendations.forEach(rec => {
        conversation += `  - [${rec.priority}] ${rec.action}: ${rec.message}\n`;
      });
      conversation += `\n`;
    }

    // Write to conversation log (explicit UTF-8 encoding)
    fs.appendFileSync(this.conversationFile, conversation, 'utf8');

    // Also log summary to main log
    this.log(`Decision Point: ${step.name} -> ${consensus.decision.toUpperCase()} (${(consensus.confidence * 100).toFixed(0)}% confidence)`);
  }

  /**
   * Write JSON conversation (NEW - for ICML experiments)
   */
  writeJsonConversation(step, agents, consensus, timestamp, userPrompt) {
    // Calculate costs for this decision
    const costs = calculateDecisionCost(agents);

    // Update session totals
    this.sessionData.costs.total_usd += costs.total_usd;
    this.sessionData.costs.breakdown.gpt5_2 += costs.breakdown.gpt5_2;
    this.sessionData.costs.breakdown.claude += costs.breakdown.claude;
    this.sessionData.costs.breakdown.gemini += costs.breakdown.gemini;

    // Build decision object
    const decision = {
      decision_id: step.decision_id || `${this.sessionData.session_id}_step${this.sessionData.decisions.length + 1}`,
      step_name: step.name,
      step_description: step.description || '',
      timestamp: timestamp,
      decision_type: step.decisionType || 'unknown',
      user_prompt: userPrompt || '',

      agent_responses: {
        gpt5_2: agents.gpt5_2 ? (() => {
          const parsed = parseAgentResponse(agents.gpt5_2.content || '');
          return {
            model: agents.gpt5_2.model,
            success: agents.gpt5_2.success,
            content: agents.gpt5_2.content || '',

            // Extracted structured fields
            extracted_decision: parsed.extracted_decision,
            confidence_label: parsed.confidence_label,
            confidence_score: parsed.confidence_score,

            tokens: agents.gpt5_2.usage || null,
            cost_usd: calculateAgentCost(agents.gpt5_2),
            skipped: agents.gpt5_2.skipped || false
          };
        })() : null,

        claude: agents.claude ? (() => {
          const parsed = parseAgentResponse(agents.claude.content || '');
          return {
            model: agents.claude.model,
            success: agents.claude.success,
            content: agents.claude.content || '',

            // Extracted structured fields
            extracted_decision: parsed.extracted_decision,
            confidence_label: parsed.confidence_label,
            confidence_score: parsed.confidence_score,

            tokens: agents.claude.usage || null,
            cost_usd: calculateAgentCost(agents.claude),
            skipped: agents.claude.skipped || false,
            mcp_tools_used: agents.claude.mcp_tools_used || []
          };
        })() : null,

        gemini: agents.gemini ? (() => {
          const parsed = parseAgentResponse(agents.gemini.content || '');
          return {
            model: agents.gemini.model,
            success: agents.gemini.success,
            content: agents.gemini.content || '',

            // Extracted structured fields
            extracted_decision: parsed.extracted_decision,
            confidence_label: parsed.confidence_label,
            confidence_score: parsed.confidence_score,

            tokens: agents.gemini.usage || null,
            cost_usd: calculateAgentCost(agents.gemini),
            skipped: agents.gemini.skipped || false
          };
        })() : null
      },

      consensus: {
        decision: consensus.decision,

        // Consistent confidence format (label + score)
        confidence_label: getConfidenceLabel(consensus.confidence),
        confidence_score: consensus.confidence,

        voting_method: consensus.votingMethod || 'majority',
        votes: consensus.votes || {},
        vote_breakdown: consensus.vote_breakdown || null,
        reasoning: consensus.reasoning || '',
        recommendations: consensus.recommendations || []
      },

      costs: costs
    };

    // Add to in-memory array
    this.sessionData.decisions.push(decision);

    // Append to JSONL file (one decision per line for streaming)
    fs.appendFileSync(
      this.decisionsJsonlFile,
      JSON.stringify(decision) + '\n',
      'utf8'
    );
  }

  /**
   * Log tool execution
   */
  logToolExecution(stepName, tool, status, details = {}) {
    const timestamp = new Date().toISOString();

    let logEntry = `\n--- Tool Execution ---\n`;
    logEntry += `Step: ${stepName}\n`;
    logEntry += `Tool: ${tool}\n`;
    logEntry += `Status: ${status}\n`;
    logEntry += `Timestamp: ${timestamp}\n`;

    if (Object.keys(details).length > 0) {
      logEntry += `Details:\n`;
      Object.entries(details).forEach(([key, value]) => {
        logEntry += `  ${key}: ${value}\n`;
      });
    }

    logEntry += `\n`;

    fs.appendFileSync(this.logFile, logEntry, 'utf8');
  }

  /**
   * Finalize log (write summary + JSON files)
   */
  finalize(session) {
    const timestamp = new Date().toISOString();

    // 1. WRITE PLAIN TEXT SUMMARY (existing)
    this.writePlainTextSummary(session, timestamp);

    // 2. WRITE JSON FILES (NEW)
    this.writeJsonFiles(session, timestamp);

    console.log('');
    console.log(`Logs saved:`);
    console.log(`   Plain text log: ${this.logFile}`);
    console.log(`   Agent conversations: ${this.conversationFile}`);
    console.log(`   JSON decisions: ${this.decisionsJsonFile}`);
    console.log(`   JSONL decisions: ${this.decisionsJsonlFile}`);
    console.log(`   Session metadata: ${this.sessionMetadataFile}`);
  }

  /**
   * Write plain text summary (existing format)
   */
  writePlainTextSummary(session, timestamp) {
    let summary = `\n${'='.repeat(80)}\n`;
    summary += `ANALYSIS COMPLETE\n`;
    summary += `Finished: ${timestamp}\n`;
    summary += `${'='.repeat(80)}\n\n`;
    summary += `Summary:\n`;

    if (session && session.dataInfo) {
      summary += `  Input Type: ${session.dataInfo.type?.toUpperCase() || 'N/A'}\n`;
      summary += `  Samples: ${session.dataInfo.samples?.length || 0}\n`;
    }
    if (session && session.steps) {
      summary += `  Steps Executed: ${session.currentStep + 1}/${session.steps.length}\n`;
    }
    if (session && session.decisions) {
      summary += `  Multi-Agent Decisions: ${Object.keys(session.decisions).length}\n`;
    }
    if (session && session.thresholds) {
      summary += `  Thresholds: FDR < ${session.thresholds.fdr || 'N/A'}, logFC > ${session.thresholds.logFC || 'N/A'}\n`;
    }
    if (session && session.config) {
      summary += `  Output Directory: ${session.config.output}\n`;
    }

    // Add cost summary
    summary += `\n`;
    summary += `  Total API Cost: $${this.sessionData.costs.total_usd.toFixed(3)}\n`;
    summary += `    GPT-5.2: $${this.sessionData.costs.breakdown.gpt5_2.toFixed(3)}\n`;
    summary += `    Claude:  $${this.sessionData.costs.breakdown.claude.toFixed(3)}\n`;
    summary += `    Gemini:  $${this.sessionData.costs.breakdown.gemini.toFixed(3)}\n`;
    summary += `\n`;

    fs.appendFileSync(this.logFile, summary, 'utf8');
    fs.appendFileSync(this.conversationFile, summary, 'utf8');
  }

  /**
   * Write JSON files (NEW - for ICML experiments)
   */
  writeJsonFiles(session, timestamp) {
    // Complete session data
    this.sessionData.timestamp_end = timestamp;
    this.sessionData.duration_seconds =
      (new Date(timestamp) - new Date(this.sessionData.timestamp_start)) / 1000;

    // Add data info if available
    if (session && session.dataInfo) {
      this.sessionData.data_info = {
        type: session.dataInfo.type,
        samples: session.dataInfo.samples,
        sample_count: session.dataInfo.samples?.length || 0,
        paired_end: session.dataInfo.pairedEnd || false,
        groups: session.dataInfo.groups || {}
      };
    }

    // Add execution info
    if (session && session.config) {
      this.sessionData.execution = {
        script_path: session.scriptPath || null,
        script_type: session.scriptType || null,
        success: session.success !== undefined ? session.success : null,
        attempts: session.attempts || 1,
        exit_code: session.exitCode || null
      };
    }

    // Calculate metrics
    this.sessionData.metrics = {
      total_decisions: this.sessionData.decisions.length,
      decisions_automation: this.sessionData.decisions.filter(d => d.consensus.decision === 'automation').length,
      decisions_adaptation: this.sessionData.decisions.filter(d => d.consensus.decision === 'adaptation').length,
      average_confidence: this.sessionData.decisions.length > 0
        ? this.sessionData.decisions.reduce((sum, d) => sum + d.consensus.confidence_score, 0) / this.sessionData.decisions.length
        : 0
    };

    // Write full JSON (all decisions in one file)
    fs.writeFileSync(
      this.decisionsJsonFile,
      JSON.stringify(this.sessionData, null, 2),
      'utf8'
    );

    // Write session metadata (separate file for quick access)
    const metadata = {
      session_id: this.sessionData.session_id,
      timestamp_start: this.sessionData.timestamp_start,
      timestamp_end: this.sessionData.timestamp_end,
      duration_seconds: this.sessionData.duration_seconds,
      system: this.sessionData.system,
      config: this.sessionData.config,
      data_info: this.sessionData.data_info,
      execution: this.sessionData.execution,
      costs: this.sessionData.costs,
      metrics: this.sessionData.metrics
    };

    fs.writeFileSync(
      this.sessionMetadataFile,
      JSON.stringify(metadata, null, 2),
      'utf8'
    );
  }
}

export { Logger };
