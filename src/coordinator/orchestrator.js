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
    this.currentDataset = null;  // Set when analysis starts, used for decision_id generation
    this.stepCounter = 0;        // Tracks pipeline step for unique decision_ids
  }

  /**
   * Set current dataset for decision_id tracking
   */
  setDataset(datasetName) {
    this.currentDataset = datasetName;
    this.stepCounter = 0;  // Reset step counter for new dataset
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
   *
   * @param {string} question - The question to ask agents
   * @param {Object} context - Context data for the decision
   * @param {string} decisionType - Type of decision (approach_decision, threshold, etc.)
   * @param {string} decision_id - Unique ID for evaluation logging (format: "{dataset}_{step}_{type}")
   */
  async consultAgents(question, context = {}, decisionType = 'threshold', decision_id = null) {
    this.log(`Consulting all agents...`);
    this.log(`Question: ${question}`);
    this.log(`Decision type: ${decisionType}`);
    if (decision_id) this.log(`Decision ID: ${decision_id}`);

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

    // Synthesize consensus (with decision_id for evaluation)
    const consensus = synthesizeConsensus(responses, decisionType, decision_id);

    // Store in session history
    this.sessionHistory.push({
      timestamp: new Date().toISOString(),
      decision_id,
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

    // Generate unique decision_id for QC review
    const decision_id = `${this.currentDataset || 'unknown'}_step7_qc_review`;

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

    const result = await this.consultAgents(question, { plotInfo, sampleInfo }, 'sample_removal', decision_id);

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

    // Generate unique decision_id for DE review
    const decision_id = `${this.currentDataset || 'unknown'}_step8_de_review`;

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

    const result = await this.consultAgents(question, { deResults, analysisParams }, 'threshold', decision_id);

    return result;
  }

  /**
   * AGENTS DECIDE: AUTOMATION or ADAPTATION approach
   */
  async decideAnalysisApproach(dataInfo, config, proposedSteps) {
    console.log('[Coordinator] → Asking all agents to analyze data...');
    console.log('');

    // Set current dataset for all subsequent decision_ids
    const datasetName = config.comparison || 'unknown';
    this.setDataset(datasetName);

    // Generate decision_id for evaluation logging
    const decision_id = `${this.currentDataset}_step0_mode`;

    const question = `
Analyze this RNA-seq dataset and decide whether to use AUTOMATION or ADAPTATION:

DATASET:
- Type: ${dataInfo.type}
- Samples: ${dataInfo.samples.length}
- Sequencing: ${dataInfo.pairedEnd ? 'Paired-end (R1 + R2 files)' : 'Single-end'}
- Organism: ${config.organism}

EXPERIMENTAL DESIGN:
- Comparison: "${config.comparison}"
- Groups: ${dataInfo.groups ? Object.keys(dataInfo.groups).map(g => `${g} (n=${dataInfo.groups[g].length})`).join(' vs ') : 'unknown'}
${dataInfo.groups && Object.keys(dataInfo.groups).length === 2 ? `- Replicates per group: ${Math.min(...Object.values(dataInfo.groups).map(g => g.length))}` : ''}

PROPOSED PIPELINE (Standard Bulk RNA-seq):
${proposedSteps.map((s, i) => `${i + 1}. ${s.name} - ${s.description}`).join('\n')}

PIPELINE CAPABILITIES:
- Alignment (fastq2bam.sh): Automatically handles ${dataInfo.pairedEnd ? 'paired-end (R1 + R2)' : 'single-end'} reads
- Feature counting: Uses Subread featureCounts
- Normalization: RPKM for visualization, edgeR uses raw counts with internal TMM normalization
- Statistics: edgeR for differential expression

DECISION CRITERIA:
- AUTOMATION: Standard pipeline, data looks clean, no edge cases → fast, $0 cost
- ADAPTATION: Issues detected (small n, batch effects, etc.) → agents write custom script

Analyze the data and recommend AUTOMATION or ADAPTATION. Provide clear reasoning (2-3 sentences).
`;

    const result = await this.consultAgents(question, { dataInfo, config }, 'approach_decision', decision_id);

    console.log('[Coordinator] 🤝 Synthesizing consensus from agent responses...');
    console.log(`[Coordinator] Decision ID: ${decision_id}`);
    console.log(`[Coordinator] Decision: ${result.consensus.decision.toUpperCase()}`);
    console.log(`[Coordinator] Confidence: ${(result.consensus.confidence * 100).toFixed(0)}%`);
    console.log(`[Coordinator] Reasoning: ${result.consensus.reasoning}`);

    return result;
  }

  /**
   * Validate threshold selection
   */
  async validateThreshold(thresholdType, proposedValue, context = {}) {
    this.log(`Validating ${thresholdType} threshold: ${proposedValue}`);

    // Generate unique decision_id for threshold validation
    // Use thresholdType to differentiate (e.g., "fdr", "logfc", "cpm")
    this.stepCounter++;
    const decision_id = `${this.currentDataset || 'unknown'}_threshold_${thresholdType}_${this.stepCounter}`;

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

    const result = await this.consultAgents(question, context, 'threshold', decision_id);

    return result;
  }

  /**
   * Decide on sample removal
   */
  async decideSampleRemoval(sampleId, reason, qcMetrics = {}) {
    this.log(`Evaluating sample removal: ${sampleId}`);

    // Generate unique decision_id for sample removal
    // Include sampleId to differentiate multiple removal decisions
    const sanitizedSampleId = sampleId.replace(/[^a-zA-Z0-9]/g, '_');
    const decision_id = `${this.currentDataset || 'unknown'}_outlier_${sanitizedSampleId}`;

    const question = `
Should we remove sample "${sampleId}" from the analysis?

Reason for consideration: ${reason}

QC Metrics:
${JSON.stringify(qcMetrics, null, 2)}

This is a critical decision. Provide your recommendation (yes/no) with detailed reasoning.
`;

    const result = await this.consultAgents(question, { sampleId, reason, qcMetrics }, 'sample_removal', decision_id);

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

  /**
   * Export decisions for evaluation (runtime log)
   * These can be joined with ground_truth offline by decision_id
   */
  exportDecisionsForEvaluation() {
    return this.sessionHistory
      .filter(h => h.decision_id)  // Only include entries with decision_id
      .map(h => ({
        decision_id: h.decision_id,
        timestamp: h.timestamp,
        decision_type: h.decisionType,
        decision: h.consensus.decision,
        confidence: h.consensus.confidence,
        disagreement_score: h.consensus.disagreement_score,
        votes: h.consensus.votes,
        agentDecisions: h.consensus.agentDecisions
      }));
  }
}

/**
 * Generate a decision_id for evaluation logging
 * Format: "{dataset}_{step}_{type}"
 *
 * @param {string} dataset - Dataset name (e.g., "DA0036")
 * @param {string|number} step - Pipeline step (e.g., "step0", "step7")
 * @param {string} type - Decision type short name (e.g., "mode", "qc", "outlier")
 */
export function generateDecisionId(dataset, step, type) {
  return `${dataset}_${step}_${type}`;
}

/**
 * Create a coordinator instance
 */
export function createCoordinator(options = {}) {
  return new Coordinator(options);
}

export default {
  Coordinator,
  createCoordinator,
  generateDecisionId
};
