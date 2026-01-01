/**
 * Coordinator Agent
 * Orchestrates multi-agent collaboration for bioinformatics analysis
 */

import {
  callGPT4,
  callClaude,
  callGemini,
  callAllAgents,
  askStatsAgent,
  askPipelineAgent,
  askQCAgent,
  askBiologyAgent
} from '../utils/llm_clients.js';

import {
  synthesizeConsensus,
  formatConsensusReport,
  vote,
  categorizeDisagreement
} from './consensus.js';

/**
 * Main Coordinator Class
 */
export class Coordinator {
  constructor(options = {}) {
    this.verbose = options.verbose || false;
    this.sessionHistory = [];
  }

  /**
   * Log message if verbose mode is on
   */
  log(message) {
    if (this.verbose) {
      console.log(`[Coordinator] ${message}`);
    }
  }

  /**
   * Ask all agents a question and synthesize consensus
   */
  async consultAgents(question, context = {}, decisionType = 'threshold') {
    this.log(`Consulting all agents...`);
    this.log(`Question: ${question}`);
    this.log(`Decision type: ${decisionType}`);

    // Call all agents in parallel
    const responses = await callAllAgents(question, {
      gpt4SystemPrompt: 'You are a statistical expert for genomics. Provide statistical validation and recommendations.',
      claudeSystemPrompt: 'You are a bioinformatics pipeline expert. Provide technical and quality control guidance.',
      geminiSystemPrompt: 'You are a molecular biology expert. Provide biological interpretation and validation.'
    });

    // Log responses
    this.log(`GPT-4: ${responses.gpt4.success ? 'responded' : 'failed'}`);
    this.log(`Claude: ${responses.claude.success ? 'responded' : 'failed'}`);
    this.log(`Gemini: ${responses.gemini.success ? 'responded' : 'failed'}`);

    // Synthesize consensus
    const consensus = synthesizeConsensus(responses, decisionType);

    // Store in session history
    this.sessionHistory.push({
      timestamp: new Date().toISOString(),
      question,
      context,
      decisionType,
      responses,
      consensus
    });

    return {
      responses,
      consensus,
      report: formatConsensusReport(consensus, question)
    };
  }

  /**
   * Review QC plots with multi-agent panel
   */
  async reviewQCPlots(plotInfo, sampleInfo = {}) {
    this.log('Reviewing QC plots...');

    const question = `
Please review the following QC analysis:

Sample Information:
${JSON.stringify(sampleInfo, null, 2)}

Plot Information:
${JSON.stringify(plotInfo, null, 2)}

Questions to address:
1. Are there any outlier samples that should be removed?
2. Is there evidence of batch effects?
3. Are the samples clustering as expected by experimental groups?
4. Should we proceed to differential expression analysis?

Provide a clear recommendation: approve to proceed, reject and re-process, or request additional analysis.
`;

    const result = await this.consultAgents(question, { plotInfo, sampleInfo }, 'sample_removal');

    return {
      ...result,
      actionRequired: result.consensus.decision !== 'approve'
    };
  }

  /**
   * Validate differential expression results
   */
  async reviewDEResults(deResults, analysisParams = {}) {
    this.log('Reviewing DE analysis results...');

    const question = `
Please review the following differential expression analysis results:

Analysis Parameters:
${JSON.stringify(analysisParams, null, 2)}

Results Summary:
${JSON.stringify(deResults, null, 2)}

Questions to address:
1. Are the FDR and logFC thresholds appropriate?
2. Is the number of DEGs reasonable given the experimental design?
3. Do the top DEGs make biological sense?
4. Are there any statistical concerns?

Provide a clear assessment and recommendations.
`;

    const result = await this.consultAgents(question, { deResults, analysisParams }, 'threshold');

    return result;
  }

  /**
   * Validate threshold selection
   */
  async validateThreshold(thresholdType, proposedValue, context = {}) {
    this.log(`Validating ${thresholdType} threshold: ${proposedValue}`);

    const question = `
Should we use ${thresholdType} = ${proposedValue} for this analysis?

Context:
${JSON.stringify(context, null, 2)}

Consider:
- Statistical rigor
- Field standards
- Sample size and power
- Biological relevance

Provide a clear yes/no recommendation with reasoning.
`;

    const result = await this.consultAgents(question, context, 'threshold');

    return result;
  }

  /**
   * Decide on sample removal
   */
  async decideSampleRemoval(sampleId, reason, qcMetrics = {}) {
    this.log(`Evaluating sample removal: ${sampleId}`);

    const question = `
Should we remove sample "${sampleId}" from the analysis?

Reason for consideration: ${reason}

QC Metrics:
${JSON.stringify(qcMetrics, null, 2)}

This is a critical decision. Provide your recommendation (yes/no) with detailed reasoning.
`;

    const result = await this.consultAgents(question, { sampleId, reason, qcMetrics }, 'sample_removal');

    return {
      ...result,
      requiresUnanimous: true,
      sampleId
    };
  }

  /**
   * Route task to appropriate agent
   */
  async routeTask(task, taskType = 'general') {
    this.log(`Routing task: ${taskType}`);

    let result;

    switch (taskType) {
      case 'statistical':
        result = await askStatsAgent(task);
        return { agent: 'stats', ...result };

      case 'pipeline':
        result = await askPipelineAgent(task);
        return { agent: 'pipeline', ...result };

      case 'qc':
        result = await askQCAgent(task);
        return { agent: 'qc', ...result };

      case 'biological':
        result = await askBiologyAgent(task);
        return { agent: 'biology', ...result };

      default:
        // Multi-agent consultation for general tasks
        return await this.consultAgents(task);
    }
  }

  /**
   * Interpret biological results
   */
  async interpretBiology(geneList, context = {}) {
    this.log('Requesting biological interpretation...');

    const question = `
Please interpret the biological significance of these genes:

Gene List:
${JSON.stringify(geneList, null, 2)}

Experimental Context:
${JSON.stringify(context, null, 2)}

Provide:
1. Likely biological pathways involved
2. Coherence with experimental hypothesis
3. Literature support (if known)
4. Recommendations for follow-up validation

Focus on biological insight and interpretation.
`;

    // Route to biology agent specifically
    const result = await askBiologyAgent(question, context);

    return {
      agent: 'biology',
      ...result
    };
  }

  /**
   * Get session summary
   */
  getSessionSummary() {
    const totalQueries = this.sessionHistory.length;
    const consensusDecisions = this.sessionHistory.map(h => h.consensus.decision);

    const decisionCounts = {
      approve: consensusDecisions.filter(d => d === 'approve').length,
      reject: consensusDecisions.filter(d => d === 'reject').length,
      user_decision_required: consensusDecisions.filter(d => d === 'user_decision_required').length,
      user_approval_required: consensusDecisions.filter(d => d === 'user_approval_required').length
    };

    const avgConfidence = this.sessionHistory.length > 0
      ? this.sessionHistory.reduce((sum, h) => sum + h.consensus.confidence, 0) / this.sessionHistory.length
      : 0;

    return {
      totalQueries,
      decisionCounts,
      averageConfidence: avgConfidence,
      history: this.sessionHistory
    };
  }

  /**
   * Export session for analysis
   */
  exportSession() {
    return {
      timestamp: new Date().toISOString(),
      summary: this.getSessionSummary(),
      history: this.sessionHistory
    };
  }
}

/**
 * Create a coordinator instance
 */
export function createCoordinator(options = {}) {
  return new Coordinator(options);
}

export default {
  Coordinator,
  createCoordinator
};
