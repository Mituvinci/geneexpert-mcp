/**
 * Coordinator Agent
 * Orchestrates multi-agent collaboration for bioinformatics analysis
 */

import fs from 'fs';

import {
  callGPT5,
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

import { getMultiAgentPrompts } from '../config/prompts.js';
import { getStagePrompts, formatStageOutputForAgents } from '../config/stage_prompts.js';

/**
 * Main Coordinator Class
 */
export class Coordinator {
  constructor(options = {}) {
    this.verbose = options.verbose || false;
    this.singleAgent = options.singleAgent || null;  // For experimental single-agent mode
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

    // Call all agents in parallel (or single agent if in experimental mode)
    // Note: In single-agent mode, the agent gets a COMBINED prompt (stats + pipeline + bio)
    // In multi-agent mode, each agent gets a domain-separated prompt
    const prompts = getMultiAgentPrompts();
    const responses = await callAllAgents(question, {
      singleAgent: this.singleAgent,
      ...prompts
    });

    // Log responses
    this.log(`GPT-5.2: ${responses.gpt5_2.success ? 'responded' : 'failed'}`);
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

  // ============================================
  // STAGED ARCHITECTURE METHODS
  // ============================================

  /**
   * Consult agents with stage-specific prompts
   * Used in staged architecture where each stage has specialized prompts
   *
   * @param {number} stage - Stage number (1, 2, 3)
   * @param {string} stageOutput - Formatted stage output for agents
   * @param {Object} context - Additional context
   * @param {string} decisionType - Type of decision
   * @param {string} decision_id - Unique ID for evaluation
   */
  async consultAgentsWithStagePrompts(stage, stageOutput, context = {}, decisionType = 'stage_checkpoint', decision_id = null) {
    this.log(`Consulting agents for Stage ${stage} review...`);
    if (decision_id) this.log(`Decision ID: ${decision_id}`);

    // Get stage-specific prompts
    const stagePrompts = getStagePrompts(stage);

    // Build options with images (if provided)
    const agentOptions = {
      singleAgent: this.singleAgent,
      ...stagePrompts
    };

    // Add images to each agent's options if available
    if (context.images && context.images.length > 0) {
      agentOptions.gpt5_2_Options = { ...(agentOptions.gpt5_2_Options || {}), images: context.images };
      agentOptions.claudeOptions = { ...(agentOptions.claudeOptions || {}), images: context.images };
      agentOptions.geminiOptions = { ...(agentOptions.geminiOptions || {}), images: context.images };
    }

    // Call all agents with stage-specific prompts (and images)
    const responses = await callAllAgents(stageOutput, agentOptions);

    // Log responses
    this.log(`GPT-5.2: ${responses.gpt5_2.success ? 'responded' : 'failed'}`);
    this.log(`Claude: ${responses.claude.success ? 'responded' : 'failed'}`);
    this.log(`Gemini: ${responses.gemini.success ? 'responded' : 'failed'}`);

    // Synthesize consensus (with decision_id for evaluation)
    const consensus = synthesizeConsensus(responses, decisionType, decision_id);

    // Store in session history
    this.sessionHistory.push({
      timestamp: new Date().toISOString(),
      decision_id,
      stage,
      stageName: `Stage ${stage}`,
      question: stageOutput,
      context,
      decisionType,
      responses,
      consensus
    });

    return {
      responses,
      consensus,
      report: formatConsensusReport(consensus, `Stage ${stage} Review`)
    };
  }

  /**
   * Review Stage 1 output (FASTQ Validation)
   *
   * @param {Object} stage1Output - Parsed output from stage1_fastq_validation.js
   * @param {Object} dataInfo - Dataset information
   * @returns {Object} - { responses, consensus, proceed }
   */
  async reviewStage1Output(stage1Output, dataInfo) {
    console.log('[Coordinator] Reviewing Stage 1 (FASTQ Validation) output...');
    console.log('');

    // Generate decision_id
    const decision_id = `${this.currentDataset || 'unknown'}_stage1_validation`;

    // Format output for agents
    const formattedOutput = formatStageOutputForAgents(1, stage1Output, dataInfo);

    // Consult agents with stage-specific prompts
    const result = await this.consultAgentsWithStagePrompts(
      1,
      formattedOutput,
      { stage1Output, dataInfo },
      'stage_checkpoint',
      decision_id
    );

    // Determine if we should proceed
    const proceed = this.shouldProceedFromStage1(result.consensus);

    console.log('[Coordinator] Stage 1 Review Complete');
    console.log(`[Coordinator] Decision: ${result.consensus.decision.toUpperCase()}`);
    console.log(`[Coordinator] Confidence: ${(result.consensus.confidence * 100).toFixed(0)}%`);
    console.log(`[Coordinator] Proceed to Stage 2: ${proceed ? 'YES' : 'NO'}`);
    console.log('');

    return {
      ...result,
      proceed,
      stage: 1,
      stageName: 'FASTQ Validation'
    };
  }

  /**
   * Determine if we should proceed from Stage 1 based on consensus
   */
  shouldProceedFromStage1(consensus) {
    const decision = consensus.decision.toLowerCase();

    // PASS or PASS_WITH_WARNING -> proceed
    if (decision.includes('pass')) {
      return true;
    }

    // FAIL -> don't proceed
    if (decision.includes('fail')) {
      return false;
    }

    // Check votes - majority PASS means proceed
    if (consensus.votes) {
      const passVotes = (consensus.votes.pass || 0) + (consensus.votes.pass_with_warning || 0);
      const failVotes = consensus.votes.fail || 0;
      return passVotes > failVotes;
    }

    // Default: proceed if confidence is high enough
    return consensus.confidence >= 0.5;
  }

  /**
   * Review Stage 2 output (Alignment + QC)
   *
   * @param {Object} stage2Output - Parsed output from stage2_alignment.js
   * @param {Object} dataInfo - Dataset information
   * @returns {Object} - { responses, consensus, proceed, samplesToRemove }
   */
  async reviewStage2Output(stage2Output, dataInfo) {
    console.log('[Coordinator] Reviewing Stage 2 (Alignment QC) output...');
    console.log('');

    // Generate decision_id
    const decision_id = `${this.currentDataset || 'unknown'}_stage2_alignment`;

    // Format output for agents
    const formattedOutput = formatStageOutputForAgents(2, stage2Output, dataInfo);

    // Consult agents with stage-specific prompts
    const result = await this.consultAgentsWithStagePrompts(
      2,
      formattedOutput,
      { stage2Output, dataInfo },
      'stage_checkpoint',
      decision_id
    );

    // Determine action based on consensus
    const { proceed, samplesToRemove } = this.parseStage2Decision(result.consensus, result.responses);

    console.log('[Coordinator] Stage 2 Review Complete');
    console.log(`[Coordinator] Decision: ${result.consensus.decision.toUpperCase()}`);
    console.log(`[Coordinator] Samples to Remove: ${samplesToRemove.length > 0 ? samplesToRemove.join(', ') : 'None'}`);
    console.log(`[Coordinator] Proceed to Stage 3: ${proceed ? 'YES' : 'NO'}`);
    console.log('');

    return {
      ...result,
      proceed,
      samplesToRemove,
      stage: 2,
      stageName: 'Alignment QC'
    };
  }

  /**
   * Parse Stage 2 decision to extract samples to remove
   */
  parseStage2Decision(consensus, responses) {
    const decision = consensus.decision.toLowerCase();
    let samplesToRemove = [];
    let proceed = true;

    if (decision.includes('abort')) {
      proceed = false;
    } else if (decision.includes('remove')) {
      // Extract samples to remove from agent responses
      samplesToRemove = this.extractSamplesToRemove(responses);
    }

    return { proceed, samplesToRemove };
  }

  /**
   * Extract samples to remove from agent responses
   */
  extractSamplesToRemove(responses) {
    const samples = new Set();

    for (const [agent, response] of Object.entries(responses)) {
      if (!response.success || !response.content) continue;

      const content = response.content;

      // Look for "Samples_to_Remove:" pattern
      const removeMatch = content.match(/Samples_to_Remove:\s*\[([^\]]+)\]/i);
      if (removeMatch) {
        const sampleList = removeMatch[1].split(',').map(s => s.trim().replace(/["']/g, ''));
        sampleList.forEach(s => {
          if (s && s.toLowerCase() !== 'none') {
            samples.add(s);
          }
        });
      }

      // Also look for sample names mentioned after "remove"
      const removePatterns = content.match(/remove[:\s]+([^\n]+)/gi);
      if (removePatterns) {
        for (const pattern of removePatterns) {
          const sampleMatches = pattern.match(/[A-Za-z0-9_-]+_R[12]/g);
          if (sampleMatches) {
            sampleMatches.forEach(s => samples.add(s.replace(/_R[12]$/, '')));
          }
        }
      }
    }

    return Array.from(samples);
  }

  /**
   * Review Stage 3 output (Quantification + QC Assessment)
   *
   * @param {Object} stage3Output - Parsed output from stage3_quantification_qc.js
   * @param {Object} dataInfo - Dataset information
   * @returns {Object} - { responses, consensus, deMethod, outlierAction, outliersToRemove }
   */
  async reviewStage3Output(stage3Output, dataInfo) {
    console.log('[Coordinator] Reviewing Stage 3 (QC Assessment) output...');
    console.log('');

    // Generate decision_id
    const decision_id = `${this.currentDataset || 'unknown'}_stage3_qc`;

    // Format output for agents
    const formattedOutput = formatStageOutputForAgents(3, stage3Output, dataInfo);

    // Read PCA plot image (if available) to pass to agents
    let pcaImages = [];
    if (stage3Output.pca_plot_path && fs.existsSync(stage3Output.pca_plot_path)) {
      console.log('[Coordinator] Loading PCA plot for agent review...');
      console.log(`[Coordinator] Plot: ${stage3Output.pca_plot_path}`);
      try {
        const imageBuffer = fs.readFileSync(stage3Output.pca_plot_path);
        const base64Image = imageBuffer.toString('base64');
        pcaImages.push({
          data: base64Image,
          mediaType: 'application/pdf'
        });
        console.log('[Coordinator] PCA plot loaded successfully');
      } catch (error) {
        console.log(`[Coordinator] Warning: Could not read PCA plot: ${error.message}`);
      }
      console.log('');
    }

    // Consult agents with stage-specific prompts (and PCA plot image)
    const result = await this.consultAgentsWithStagePrompts(
      3,
      formattedOutput,
      { stage3Output, dataInfo, images: pcaImages },
      'stage_checkpoint',
      decision_id
    );

    // Extract DE method and outlier decision
    const { deMethod, batchSpecification, outlierAction, outliersToRemove } =
      this.parseStage3Decision(result.consensus, result.responses);

    console.log('[Coordinator] Stage 3 Review Complete');
    console.log(`[Coordinator] DE Method: ${deMethod}`);
    console.log(`[Coordinator] Batch Specification: ${batchSpecification || 'N/A'}`);
    console.log(`[Coordinator] Outlier Action: ${outlierAction}`);
    console.log(`[Coordinator] Outliers to Remove: ${outliersToRemove.length > 0 ? outliersToRemove.join(', ') : 'None'}`);
    console.log('');

    return {
      ...result,
      deMethod,
      batchSpecification,
      outlierAction,
      outliersToRemove,
      stage: 3,
      stageName: 'QC Assessment'
    };
  }

  /**
   * Parse Stage 3 decision to extract DE method and outlier handling
   */
  parseStage3Decision(consensus, responses) {
    let deMethod = 'simpleEdger';  // default
    let batchSpecification = null;
    let outlierAction = 'KEEP_ALL';
    let outliersToRemove = [];

    // Count votes for DE method
    let simpleEdgerVotes = 0;
    let batchEffectVotes = 0;
    let keepAllVotes = 0;
    let removeOutlierVotes = 0;

    for (const [agent, response] of Object.entries(responses)) {
      if (!response.success || !response.content) continue;

      const content = response.content;

      // Extract DE_Method
      const deMatch = content.match(/DE_Method:\s*\[?\s*(simpleEdger|batch_effect_edger)\s*\]?/i);
      if (deMatch) {
        if (deMatch[1].toLowerCase().includes('batch')) {
          batchEffectVotes++;
        } else {
          simpleEdgerVotes++;
        }
      }

      // Extract Batch_Specification
      const batchMatch = content.match(/Batch_Specification:\s*\[?\s*([^\]\n]+)\s*\]?/i);
      if (batchMatch && batchMatch[1].toLowerCase() !== 'n/a') {
        batchSpecification = batchMatch[1].trim().replace(/["']/g, '');
      }

      // Extract Outlier_Action
      const outlierMatch = content.match(/Outlier_Action:\s*\[?\s*(KEEP_ALL|REMOVE_OUTLIERS)\s*\]?/i);
      if (outlierMatch) {
        if (outlierMatch[1].toUpperCase().includes('REMOVE')) {
          removeOutlierVotes++;
        } else {
          keepAllVotes++;
        }
      }

      // Extract Outliers_to_Remove
      const outliersMatch = content.match(/Outliers_to_Remove:\s*\[([^\]]+)\]/i);
      if (outliersMatch) {
        const outlierList = outliersMatch[1].split(',').map(s => s.trim().replace(/["']/g, ''));
        outlierList.forEach(s => {
          if (s && s.toLowerCase() !== 'none') {
            outliersToRemove.push(s);
          }
        });
      }
    }

    // Majority vote for DE method
    deMethod = batchEffectVotes > simpleEdgerVotes ? 'batch_effect_edger' : 'simpleEdger';

    // Majority vote for outlier action
    outlierAction = removeOutlierVotes > keepAllVotes ? 'REMOVE_OUTLIERS' : 'KEEP_ALL';

    // Deduplicate outliers
    outliersToRemove = [...new Set(outliersToRemove)];

    return { deMethod, batchSpecification, outlierAction, outliersToRemove };
  }

  /**
   * Review Stage 4 output (DE Analysis)
   */
  async reviewStage4Output(stage4Output, dataInfo) {
    console.log('[Coordinator] Reviewing Stage 4 (DE Analysis) output...');
    console.log('');

    // Generate decision_id
    const decision_id = `${this.currentDataset || 'unknown'}_stage4_de_analysis`;

    // Format output for agents
    const formattedOutput = formatStageOutputForAgents(4, stage4Output, dataInfo);

    // Consult agents with stage-specific prompts
    const result = await this.consultAgentsWithStagePrompts(
      4,
      formattedOutput,
      { stage4Output, dataInfo },
      'stage_checkpoint',
      decision_id
    );

    // Extract approval decision
    const { approve } = this.parseStage4Decision(result.consensus, result.responses);

    console.log('[Coordinator] Stage 4 Review Complete');
    console.log(`[Coordinator] Final Decision: ${approve ? 'APPROVE' : 'REQUEST_REANALYSIS'}`);
    console.log('');

    return {
      ...result,
      approve,
      proceed: approve,  // proceed if approved
      stage: 4,
      stageName: 'DE Analysis'
    };
  }

  /**
   * Parse Stage 4 decision to extract approval
   */
  parseStage4Decision(consensus, responses) {
    let approveVotes = 0;
    let reanalysisVotes = 0;

    for (const [agent, response] of Object.entries(responses)) {
      if (!response.success || !response.content) continue;

      const content = response.content;

      // Extract Final_Decision
      const decisionMatch = content.match(/Final_Decision:\s*\[?\s*(APPROVE|REQUEST_REANALYSIS)\s*\]?/i);
      if (decisionMatch) {
        if (decisionMatch[1].toUpperCase().includes('APPROVE')) {
          approveVotes++;
        } else {
          reanalysisVotes++;
        }
      }
    }

    // Majority vote for approval
    const approve = approveVotes > reanalysisVotes;

    return { approve };
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
