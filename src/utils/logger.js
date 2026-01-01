/**
 * Logger - Save all agent interactions and conversations to log file
 */

import fs from 'fs';
import path from 'path';

class Logger {
  constructor(outputDir, sessionName = 'analysis') {
    this.outputDir = outputDir;
    this.sessionName = sessionName;
    this.logFile = path.join(outputDir, `${sessionName}_log.txt`);
    this.conversationFile = path.join(outputDir, `${sessionName}_agent_conversations.txt`);

    // Create output directory if needed
    if (!fs.existsSync(outputDir)) {
      fs.mkdirSync(outputDir, { recursive: true });
    }

    // Initialize log files
    this.initLogFiles();
  }

  initLogFiles() {
    const timestamp = new Date().toISOString();

    // Main log file
    fs.writeFileSync(this.logFile, `GeneExpert Multi-Agent Analysis Log\n`);
    fs.appendFileSync(this.logFile, `Session: ${this.sessionName}\n`);
    fs.appendFileSync(this.logFile, `Started: ${timestamp}\n`);
    fs.appendFileSync(this.logFile, `${'='.repeat(80)}\n\n`);

    // Agent conversations file
    fs.writeFileSync(this.conversationFile, `GeneExpert Agent Conversations\n`);
    fs.appendFileSync(this.conversationFile, `Session: ${this.sessionName}\n`);
    fs.appendFileSync(this.conversationFile, `Started: ${timestamp}\n`);
    fs.appendFileSync(this.conversationFile, `${'='.repeat(80)}\n\n`);

    console.log(`📝 Logging to:`);
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

    fs.appendFileSync(this.logFile, logLine);

    // Also print to console
    console.log(message);
  }

  /**
   * Log multi-agent conversation (the full thinking!)
   */
  logAgentConversation(step, agents, consensus) {
    const timestamp = new Date().toISOString();

    let conversation = `\n${'='.repeat(80)}\n`;
    conversation += `DECISION POINT: ${step.name}\n`;
    conversation += `Timestamp: ${timestamp}\n`;
    conversation += `Description: ${step.description}\n`;
    conversation += `${'='.repeat(80)}\n\n`;

    // Log each agent's full response
    if (agents.gpt4) {
      conversation += `┌─ Stats Agent (GPT-4) ─────────────────────────────────────────┐\n`;
      conversation += `│ Model: ${agents.gpt4.model}\n`;
      conversation += `│ Status: ${agents.gpt4.success ? 'Success' : 'Failed'}\n`;
      conversation += `└────────────────────────────────────────────────────────────────┘\n\n`;
      conversation += `${agents.gpt4.content || 'No response'}\n\n`;
    }

    if (agents.claude) {
      conversation += `┌─ Pipeline Agent (Claude) ──────────────────────────────────────┐\n`;
      conversation += `│ Model: ${agents.claude.model}\n`;
      conversation += `│ Status: ${agents.claude.success ? 'Success' : 'Failed'}\n`;
      conversation += `└────────────────────────────────────────────────────────────────┘\n\n`;
      conversation += `${agents.claude.content || 'No response'}\n\n`;
    }

    if (agents.gemini) {
      conversation += `┌─ Biology Agent (Gemini) ───────────────────────────────────────┐\n`;
      conversation += `│ Model: ${agents.gemini.model}\n`;
      conversation += `│ Status: ${agents.gemini.success ? 'Success' : 'Failed'}\n`;
      conversation += `└────────────────────────────────────────────────────────────────┘\n\n`;
      conversation += `${agents.gemini.content || 'No response'}\n\n`;
    }

    // Log consensus decision
    conversation += `┌─ CONSENSUS DECISION ───────────────────────────────────────────┐\n`;
    conversation += `│ Decision: ${consensus.decision.toUpperCase()}\n`;
    conversation += `│ Confidence: ${(consensus.confidence * 100).toFixed(0)}%\n`;
    conversation += `│ Votes:\n`;
    conversation += `│   Approve:   ${consensus.votes?.approve || 0}/3\n`;
    conversation += `│   Reject:    ${consensus.votes?.reject || 0}/3\n`;
    conversation += `│   Uncertain: ${consensus.votes?.uncertain || 0}/3\n`;
    conversation += `└────────────────────────────────────────────────────────────────┘\n\n`;

    if (consensus.recommendations?.length > 0) {
      conversation += `Recommendations:\n`;
      consensus.recommendations.forEach(rec => {
        conversation += `  - [${rec.priority}] ${rec.action}: ${rec.message}\n`;
      });
      conversation += `\n`;
    }

    // Write to conversation log
    fs.appendFileSync(this.conversationFile, conversation);

    // Also log summary to main log
    this.log(`Decision Point: ${step.name} → ${consensus.decision.toUpperCase()} (${(consensus.confidence * 100).toFixed(0)}% confidence)`);
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

    fs.appendFileSync(this.logFile, logEntry);
  }

  /**
   * Finalize log (write summary)
   */
  finalize(session) {
    const timestamp = new Date().toISOString();

    let summary = `\n${'='.repeat(80)}\n`;
    summary += `ANALYSIS COMPLETE\n`;
    summary += `Finished: ${timestamp}\n`;
    summary += `${'='.repeat(80)}\n\n`;
    summary += `Summary:\n`;
    summary += `  Input Type: ${session.dataInfo.type.toUpperCase()}\n`;
    summary += `  Samples: ${session.dataInfo.samples.length}\n`;
    summary += `  Steps Executed: ${session.currentStep + 1}/${session.steps.length}\n`;
    summary += `  Multi-Agent Decisions: ${Object.keys(session.decisions).length}\n`;
    summary += `  Thresholds: FDR < ${session.thresholds.fdr || 'N/A'}, logFC > ${session.thresholds.logFC || 'N/A'}\n`;
    summary += `  Output Directory: ${session.config.output}\n`;
    summary += `\n`;

    fs.appendFileSync(this.logFile, summary);
    fs.appendFileSync(this.conversationFile, summary);

    console.log('');
    console.log(`✅ Logs saved:`);
    console.log(`   ${this.logFile}`);
    console.log(`   ${this.conversationFile}`);
  }
}

export { Logger };
