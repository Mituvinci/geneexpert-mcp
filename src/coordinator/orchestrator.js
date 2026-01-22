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
  callAllAgentsSequential,  // NEW: For sequential chain mode
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
import { getStagePrompts, getCombinedStagePrompt, formatStageOutputForAgents } from '../config/stage_prompts.js';
import {
  SCRNA_STAGE_2_THRESHOLD_PROMPTS,
  SCRNA_STAGE_2_PROMPTS,
  SCRNA_STAGE_4_PROMPTS,
  SCRNA_STAGE_5_PROMPTS
} from '../config/scrna_stage_prompts.js';

/**
 * Main Coordinator Class
 */
export class Coordinator {
  constructor(options = {}) {
    this.verbose = options.verbose || false;
    this.singleAgent = options.singleAgent || null;  // For experimental single-agent mode
    this.sequentialChain = options.sequentialChain || false;  // NEW: For sequential chain mode
    this.sessionHistory = [];
    this.currentDataset = null;  // Set when analysis starts, used for decision_id generation
    this.stepCounter = 0;        // Tracks pipeline step for unique decision_ids
  }

  /**
   * Set current dataset for decision_id tracking
   * Strips system suffixes to keep decision IDs consistent with ground_truth.json
   */
  setDataset(datasetName) {
    // Strip system suffixes: _parallel, _single_claude, _single_gpt, _single_gemini, _no_agent, _sequential
    this.currentDataset = datasetName.replace(/_(parallel|single_claude|single_gpt|single_gemini|no_agent|sequential)$/, '');
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

    // Build options
    const agentOptions = {
      singleAgent: this.singleAgent
    };

    // CRITICAL: For single-agent mode, use COMBINED prompt (all 3 perspectives)
    // For multi-agent mode, use separate prompts (one perspective each)
    if (this.singleAgent) {
      // Single-agent: Get combined prompt that merges all 3 perspectives
      const combinedPrompt = getCombinedStagePrompt(stage);

      // Provide combined prompt for all three agent options
      // (only the selected agent will be called, but we structure it this way for consistency)
      agentOptions.gpt5_2_SystemPrompt = combinedPrompt;
      agentOptions.claudeSystemPrompt = combinedPrompt;
      agentOptions.geminiSystemPrompt = combinedPrompt;
    } else {
      // Multi-agent mode: separate prompts (each agent sees only their role)
      agentOptions.gpt5_2_SystemPrompt = stagePrompts.gpt5_2_SystemPrompt;
      agentOptions.claudeSystemPrompt = stagePrompts.claudeSystemPrompt;
      agentOptions.geminiSystemPrompt = stagePrompts.geminiSystemPrompt;
    }

    // Add images to each agent's options if available
    if (context.images && context.images.length > 0) {
      agentOptions.gpt5_2_Options = { ...(agentOptions.gpt5_2_Options || {}), images: context.images };
      agentOptions.claudeOptions = { ...(agentOptions.claudeOptions || {}), images: context.images };
      agentOptions.geminiOptions = { ...(agentOptions.geminiOptions || {}), images: context.images };
    }

    // Call all agents with stage-specific prompts (and images)
    // Choose parallel or sequential mode based on flag
    const responses = this.sequentialChain
      ? await callAllAgentsSequential(stageOutput, agentOptions)  // Sequential: GPT-5.2 → Gemini → Claude
      : await callAllAgents(stageOutput, agentOptions);            // Parallel: All at once (default)

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
    const decision_id = `${this.currentDataset || 'unknown'}_stage1`;

    // Format output for agents
    const formattedOutput = formatStageOutputForAgents(1, stage1Output, dataInfo);

    // Consult agents with stage-specific prompts
    const result = await this.consultAgentsWithStagePrompts(
      1,
      formattedOutput,
      { stage1Output, dataInfo },
      'stage1',
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

    // Approve = proceed (agents agreed on PASS or PASS_WITH_WARNING)
    if (decision === 'approve') {
      return true;
    }

    // Reject = don't proceed (agents agreed on FAIL)
    if (decision === 'reject') {
      return false;
    }

    // User decision required = don't proceed automatically
    if (decision.includes('user_decision')) {
      return false;
    }

    // Check votes - majority approve means proceed
    if (consensus.votes && consensus.votes.approve >= 2) {
      return true;
    }

    // Default: don't proceed without clear consensus
    return false;
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
    const decision_id = `${this.currentDataset || 'unknown'}_stage2`;

    // Format output for agents
    const formattedOutput = formatStageOutputForAgents(2, stage2Output, dataInfo);

    // Consult agents with stage-specific prompts
    const result = await this.consultAgentsWithStagePrompts(
      2,
      formattedOutput,
      { stage2Output, dataInfo },
      'stage2',
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

    // Approve = proceed with all samples (agents agreed on PASS_ALL)
    if (decision === 'approve') {
      proceed = true;
      samplesToRemove = [];
    }
    // Reject = abort pipeline (agents agreed on ABORT)
    else if (decision === 'reject') {
      proceed = false;
      samplesToRemove = [];
    }
    // User decision required = check what agents recommended
    else if (decision.includes('user_decision')) {
      // Extract samples to remove from agent responses
      samplesToRemove = this.extractSamplesToRemove(responses);
      proceed = true;  // Proceed but may need to remove samples
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
    const decision_id = `${this.currentDataset || 'unknown'}_stage3`;

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

        // Detect media type from file extension
        let mediaType = 'application/pdf';  // Default
        if (stage3Output.pca_plot_path.endsWith('.jpg') || stage3Output.pca_plot_path.endsWith('.jpeg')) {
          mediaType = 'image/jpeg';
        } else if (stage3Output.pca_plot_path.endsWith('.png')) {
          mediaType = 'image/png';
        }

        pcaImages.push({
          data: base64Image,
          mediaType: mediaType
        });
        console.log(`[Coordinator] PCA plot loaded successfully (${mediaType})`);
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
      'stage3_de_method',
      decision_id
    );

    // Extract DE method and outlier decision
    const { deMethod, batchSpecification, outlierAction, outliersToRemove } =
      this.parseStage3Decision(result.consensus, result.responses);

    // Determine if we should proceed to Stage 4
    // For Stage 3: proceed if consensus is 'approve' (agents chose a DE method)
    // Both simpleEdger and batch_effect_edger are valid "proceed" decisions
    const proceed = result.consensus.decision.toLowerCase() === 'approve';

    console.log('[Coordinator] Stage 3 Review Complete');
    console.log(`[Coordinator] DE Method: ${deMethod}`);
    console.log(`[Coordinator] Batch Specification: ${batchSpecification || 'N/A'}`);
    console.log(`[Coordinator] Outlier Action: ${outlierAction}`);
    console.log(`[Coordinator] Outliers to Remove: ${outliersToRemove.length > 0 ? outliersToRemove.join(', ') : 'None'}`);
    console.log(`[Coordinator] Proceed to Stage 4: ${proceed ? 'YES' : 'NO'}`);
    console.log('');

    return {
      ...result,
      proceed,
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
    const decision_id = `${this.currentDataset || 'unknown'}_stage4`;

    // Format output for agents
    const formattedOutput = formatStageOutputForAgents(4, stage4Output, dataInfo);

    // Consult agents with stage-specific prompts
    const result = await this.consultAgentsWithStagePrompts(
      4,
      formattedOutput,
      { stage4Output, dataInfo },
      'stage4',
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

  /**
   * ===================================================================
   * scRNA-seq SPECIFIC REVIEW METHODS
   * ===================================================================
   */

  /**
   * Review scRNA Stage 2: QC Filtering
   * Stats Agent decides: PROCEED / PROCEED_WITH_WARNING / STOP_AND_REVIEW
   */
  /**
   * Review scRNA Stage 1 Output: Recommend QC Thresholds (NEW)
   * ALL 3 AGENTS (Stats, Pipeline, Biology) recommend thresholds BEFORE filtering
   */
  async reviewScRNAStage1ForThresholds(stage1Output, dataInfo) {
    this.log('[Coordinator] Reviewing scRNA Stage 1 (QC Metrics) - Recommending Thresholds...');

    // Format data for agents
    const formattedOutput = `# scRNA-seq QC Metrics (Stage 1 - Before Filtering)

## Dataset Overview
- Organism: ${dataInfo.organism || 'unknown'}
- Tissue: ${dataInfo.tissue || 'unknown (please infer from context if possible)'}
- Total Cells: ${stage1Output.cells_total}
- Total Genes: ${stage1Output.genes_total}

## QC Metric Distributions
- **Median Features/Cell**: ${stage1Output.nFeature_median} genes
- **Median UMI Count/Cell**: ${stage1Output.nCount_median} UMIs
- **Median % Mitochondrial**: ${stage1Output.percent_mt_median.toFixed(2)}%

## Full QC Summary (JSON)
\`\`\`json
${JSON.stringify(stage1Output, null, 2)}
\`\`\`

## Your Task
Based on these QC metrics, recommend filtering thresholds that balance:
1. **Quality control**: Remove low-quality cells (empty droplets, dying cells, doublets)
2. **Cell preservation**: Keep enough high-quality cells for downstream analysis
3. **Biological context**: Consider tissue-specific patterns (e.g., brain has higher MT%)

**Important**: This is a ${stage1Output.cells_total < 1000 ? 'SMALL' : stage1Output.cells_total < 5000 ? 'MEDIUM' : 'LARGE'} dataset. ${stage1Output.cells_total < 1000 ? 'Use LENIENT thresholds to preserve cells.' : ''}

## Threshold Parameters to Recommend
- **nFeature_min**: Minimum genes per cell (typical: 200-1000)
- **nFeature_max**: Maximum genes per cell (typical: 5000-15000)
- **percent_mt_max**: Maximum % mitochondrial (typical: 10-30%, tissue-dependent)`;

    // Build agent options
    const agentOptions = {
      singleAgent: this.singleAgent,
      gpt5_2_SystemPrompt: SCRNA_STAGE_2_THRESHOLD_PROMPTS.gpt5_2,
      claudeSystemPrompt: SCRNA_STAGE_2_THRESHOLD_PROMPTS.claude,
      geminiSystemPrompt: SCRNA_STAGE_2_THRESHOLD_PROMPTS.gemini
    };

    // Call agents
    const responses = this.sequentialChain
      ? await callAllAgentsSequential(formattedOutput, agentOptions)
      : await callAllAgents(formattedOutput, agentOptions);

    this.log(`GPT-5.2 (Stats): ${responses.gpt5_2.success ? 'responded' : 'failed'}`);
    this.log(`Claude (Pipeline): ${responses.claude.success ? 'responded' : 'failed'}`);
    this.log(`Gemini (Biology): ${responses.gemini.success ? 'responded' : 'failed'}`);

    // Parse threshold recommendations from each agent
    const thresholds = {
      gpt5_2: this._parseThresholds(responses.gpt5_2),
      claude: this._parseThresholds(responses.claude),
      gemini: this._parseThresholds(responses.gemini)
    };

    // Use AVERAGE of agent recommendations (similar to PC selection)
    const avgThresholds = {
      nFeature_min: Math.round((thresholds.gpt5_2.nFeature_min + thresholds.claude.nFeature_min + thresholds.gemini.nFeature_min) / 3),
      nFeature_max: Math.round((thresholds.gpt5_2.nFeature_max + thresholds.claude.nFeature_max + thresholds.gemini.nFeature_max) / 3),
      percent_mt_max: Math.round((thresholds.gpt5_2.percent_mt_max + thresholds.claude.percent_mt_max + thresholds.gemini.percent_mt_max) / 3)
    };

    // Synthesize consensus
    const consensus = synthesizeConsensus(responses, 'scrna_stage2_thresholds', null);

    // Map vote category back to canonical decision
    let canonicalDecision = 'SET_THRESHOLDS';  // Default
    if (consensus.decision === 'approve') {
      canonicalDecision = 'SET_THRESHOLDS';  // Either custom or default thresholds
    } else if (consensus.decision === 'reject') {
      canonicalDecision = 'INSUFFICIENT_DATA';
    }

    this.log(`[Coordinator] scRNA Stage 2 Threshold Consensus: ${canonicalDecision}`);
    this.log(`[Coordinator] Average Thresholds: nFeature ${avgThresholds.nFeature_min}-${avgThresholds.nFeature_max}, MT% <${avgThresholds.percent_mt_max}%`);

    return {
      decision: canonicalDecision,
      thresholds: avgThresholds,
      reasoning: consensus.summary,
      confidence: consensus.confidence,
      agentResponses: responses,
      consensus,
      individualThresholds: thresholds
    };
  }

  /**
   * Helper: Parse thresholds from agent response
   */
  _parseThresholds(agentResponse) {
    if (!agentResponse.success || !agentResponse.content) {
      // Return conservative defaults if agent failed
      return { nFeature_min: 200, nFeature_max: 6000, percent_mt_max: 20 };
    }

    const content = agentResponse.content;

    // Extract thresholds using regex
    const nFeatureMinMatch = content.match(/nFeature_min:\s*(\d+)/i);
    const nFeatureMaxMatch = content.match(/nFeature_max:\s*(\d+)/i);
    const percentMtMaxMatch = content.match(/percent_mt_max:\s*(\d+)/i);

    return {
      nFeature_min: nFeatureMinMatch ? parseInt(nFeatureMinMatch[1]) : 200,
      nFeature_max: nFeatureMaxMatch ? parseInt(nFeatureMaxMatch[1]) : 6000,
      percent_mt_max: percentMtMaxMatch ? parseInt(percentMtMaxMatch[1]) : 20
    };
  }

  /**
   * Review scRNA Stage 2: QC Filtering (DEPRECATED - now just applies agent thresholds)
   * ALL 3 AGENTS (Stats, Pipeline, Biology) vote on filtering quality
   */
  async reviewScRNAStage2(stage2Output, dataInfo) {
    this.log('[Coordinator] Reviewing scRNA Stage 2 (QC Filtering) with ALL 3 AGENTS...');

    // Format data for agents
    const context = {
      organism: dataInfo.organism || 'unknown',
      tissue: dataInfo.tissue || 'unknown',
      cells_before: stage2Output.cells_before_filtering,
      cells_after: stage2Output.cells_after_filtering,
      percent_removed: stage2Output.percent_removed
    };

    const formattedOutput = `# scRNA-seq QC Filtering Results

## Dataset Context
- Organism: ${context.organism}
- Tissue: ${context.tissue || 'Not specified'}
- Technology: 10x Genomics Chromium

## Filtering Summary
- Cells Before Filtering: ${context.cells_before?.toLocaleString() || 'N/A'}
- Cells After Filtering: ${context.cells_after?.toLocaleString() || 'N/A'}
- Percent Removed: ${context.percent_removed?.toFixed(1) || 'N/A'}%

## QC Metrics (JSON)
\`\`\`json
${JSON.stringify(stage2Output, null, 2)}
\`\`\`

## Your Task
Evaluate whether the QC filtering is appropriate for this scRNA-seq dataset.

**Allowed Decisions:**
- PROCEED
- PROCEED_WITH_WARNING
- STOP_AND_REVIEW`;

    // Build agent options with scRNA Stage 2 prompts
    const agentOptions = {
      singleAgent: this.singleAgent,
      gpt5_2_SystemPrompt: SCRNA_STAGE_2_PROMPTS.gpt5_2,
      claudeSystemPrompt: SCRNA_STAGE_2_PROMPTS.claude,
      geminiSystemPrompt: SCRNA_STAGE_2_PROMPTS.gemini
    };

    // Call agents (parallel or sequential based on flag)
    const responses = this.sequentialChain
      ? await callAllAgentsSequential(formattedOutput, agentOptions)
      : await callAllAgents(formattedOutput, agentOptions);

    this.log(`GPT-5.2 (Stats): ${responses.gpt5_2.success ? 'responded' : 'failed'}`);
    this.log(`Claude (Pipeline): ${responses.claude.success ? 'responded' : 'failed'}`);
    this.log(`Gemini (Biology): ${responses.gemini.success ? 'responded' : 'failed'}`);

    // Synthesize consensus
    const consensus = synthesizeConsensus(responses, 'scrna_stage2', null);

    // Map vote category back to canonical decision for scRNA
    let canonicalDecision = 'STOP_AND_REVIEW';  // Default
    if (consensus.decision === 'approve') {
      // Check if any agent said PROCEED_WITH_WARNING
      const hasWarning = [responses.gpt5_2, responses.claude, responses.gemini].some(r =>
        r.content && r.content.includes('PROCEED_WITH_WARNING')
      );
      canonicalDecision = hasWarning ? 'PROCEED_WITH_WARNING' : 'PROCEED';
    } else if (consensus.decision === 'reject') {
      canonicalDecision = 'STOP_AND_REVIEW';
    }

    this.log(`[Coordinator] scRNA Stage 2 Consensus: ${canonicalDecision} (vote: ${consensus.decision})`);

    return {
      decision: canonicalDecision,
      reasoning: consensus.summary,
      confidence: consensus.confidence,
      agentResponses: responses,
      consensus
    };
  }

  /**
   * Review scRNA Stage 4: PCA + PC Selection
   * Stats Agent decides: USE_DEFAULT / SELECT_PC_RANGE / STOP_AND_REVIEW
   */
  /**
   * Review scRNA Stage 4: PCA + PC Selection
   * ALL 3 AGENTS vote on optimal PC range
   */
  async reviewScRNAStage4(stage4Output, dataInfo) {
    this.log('[Coordinator] Reviewing scRNA Stage 4 (PCA) with ALL 3 AGENTS...');

    // Calculate cumulative variance
    const pcVariance = stage4Output.pc_variance_percent || [];
    const cumVariance = [];
    let sum = 0;
    for (let i = 0; i < Math.min(pcVariance.length, 30); i++) {
      sum += pcVariance[i];
      cumVariance.push(sum.toFixed(2));
    }

    // Format data for agents
    const formattedOutput = `# scRNA-seq PCA Results

## Dataset Context
- Organism: ${dataInfo.organism || 'unknown'}
- Tissue: ${dataInfo.tissue || 'unknown'}
- Total Cells: ${stage4Output.n_cells || 'N/A'}

## PCA Variance Summary
**Individual PC Variance (first 20 PCs):**
${pcVariance.slice(0, 20).map((v, i) => `PC ${i + 1}: ${v.toFixed(2)}%`).join('\n')}

**Cumulative Variance:**
- PC 1-10: ${cumVariance[9] || 'N/A'}%
- PC 1-20: ${cumVariance[19] || 'N/A'}%
- PC 1-30: ${cumVariance[29] || 'N/A'}%

## Full Data (JSON)
\`\`\`json
${JSON.stringify({ pc_variance_percent: pcVariance.slice(0, 50), cumulative_variance: cumVariance }, null, 2)}
\`\`\`

## Your Task
Based on the variance explained, select an appropriate PC range (min_pc to max_pc) for downstream clustering and UMAP.

**Allowed Decisions:**
- USE_DEFAULT (use PCs 1-20)
- SELECT_PC_RANGE (specify custom min_pc and max_pc)
- STOP_AND_REVIEW (if variance pattern is concerning)`;

    // Build agent options with scRNA Stage 4 prompts
    const agentOptions = {
      singleAgent: this.singleAgent,
      gpt5_2_SystemPrompt: SCRNA_STAGE_4_PROMPTS.gpt5_2,
      claudeSystemPrompt: SCRNA_STAGE_4_PROMPTS.claude,
      geminiSystemPrompt: SCRNA_STAGE_4_PROMPTS.gemini
    };

    // Call agents
    const responses = this.sequentialChain
      ? await callAllAgentsSequential(formattedOutput, agentOptions)
      : await callAllAgents(formattedOutput, agentOptions);

    this.log(`GPT-5.2 (Stats): ${responses.gpt5_2.success ? 'responded' : 'failed'}`);
    this.log(`Claude (Pipeline): ${responses.claude.success ? 'responded' : 'failed'}`);
    this.log(`Gemini (Biology): ${responses.gemini.success ? 'responded' : 'failed'}`);

    // Parse PC selections from each agent
    const pcSelections = {
      gpt5_2: this._parsePCSelection(responses.gpt5_2),
      claude: this._parsePCSelection(responses.claude),
      gemini: this._parsePCSelection(responses.gemini)
    };

    // Use AVERAGE (not median) of max_pc from all 3 agents
    const maxPCs = [pcSelections.gpt5_2.max_pc, pcSelections.claude.max_pc, pcSelections.gemini.max_pc];
    const avgMaxPC = Math.round(maxPCs.reduce((a, b) => a + b, 0) / maxPCs.length);

    const consensus = synthesizeConsensus(responses, 'scrna_stage4', null);

    // Map vote category back to canonical decision
    let canonicalDecision = 'SELECT_PC_RANGE';  // Default to success (all agents likely gave valid PC ranges)
    if (consensus.decision === 'reject') {
      canonicalDecision = 'STOP_AND_REVIEW';
    } else if (consensus.decision === 'approve' || consensus.decision === 'user_decision_required') {
      // Either consensus approve, or no clear consensus but all gave valid PC ranges
      canonicalDecision = 'SELECT_PC_RANGE';
    }

    this.log(`[Coordinator] scRNA Stage 4 Consensus: ${canonicalDecision} - PCs 1-${avgMaxPC} (average of ${maxPCs.join(', ')})`);
    this.log(`[Coordinator] Voting result: ${consensus.decision}`);

    return {
      decision: canonicalDecision,
      min_pc: 1,
      max_pc: avgMaxPC,
      reasoning: consensus.summary,
      confidence: consensus.confidence,
      agentResponses: responses,
      consensus,
      pcSelections
    };
  }

  /**
   * Helper: Parse PC selection from agent response
   */
  _parsePCSelection(response) {
    let min_pc = 1;
    let max_pc = 20;  // Default

    if (!response.success) return { min_pc, max_pc };

    try {
      // Look for patterns like "max_pc: 25" or "PCs 1-30" or "use 1-25"
      const content = response.content;

      // Pattern 1: max_pc: 25
      const maxPcMatch = content.match(/max_pc[:\s]+(\d+)/i);
      if (maxPcMatch) {
        max_pc = parseInt(maxPcMatch[1]);
      }

      // Pattern 2: PCs 1-30 or 1-25
      const rangeMatch = content.match(/(?:PCs?|use|dims?)\s*(?:1-)?(\d+)/i);
      if (rangeMatch && !maxPcMatch) {
        max_pc = parseInt(rangeMatch[1]);
      }

      // Sanity check
      if (max_pc < 10) max_pc = 10;
      if (max_pc > 50) max_pc = 50;

    } catch (error) {
      this.log(`Warning: Failed to parse PC selection from agent: ${error.message}`);
    }

    return { min_pc, max_pc };
  }

  /**
   * Review scRNA Stage 5: Clustering + Markers
   * Pipeline Agent decides: ACCEPT_CLUSTERING / ADJUST_RESOLUTION / FLAG_SUSPICIOUS
   */
  /**
   * Review scRNA Stage 5: Clustering + Markers
   * ALL 3 AGENTS vote on clustering quality
   */
  async reviewScRNAStage5(stage5Output, dataInfo) {
    this.log('[Coordinator] Reviewing scRNA Stage 5 (Clustering) with ALL 3 AGENTS...');

    // Format cluster sizes
    const clusterSizes = stage5Output.cluster_sizes || {};
    const clusterSummary = Object.entries(clusterSizes)
      .map(([cluster, size]) => `  Cluster ${cluster}: ${size} cells`)
      .join('\n');

    // Format data for agents
    const formattedOutput = `# scRNA-seq Clustering Results

## Dataset Context
- Organism: ${dataInfo.organism || 'unknown'}
- Tissue: ${dataInfo.tissue || 'unknown'}
- Total Cells: ${stage5Output.n_cells || 'N/A'}
- Number of Clusters: ${stage5Output.n_clusters || 'N/A'}

## Cluster Size Distribution
${clusterSummary}

## Top Markers Per Cluster
${stage5Output.top_markers_summary || '(See full JSON below)'}

## Full Data (JSON)
\`\`\`json
${JSON.stringify(stage5Output, null, 2)}
\`\`\`

## Your Task
Evaluate whether the clustering is biologically and technically valid.

**Allowed Decisions:**
- ACCEPT_CLUSTERING (clustering is good, proceed to analysis)
- ADJUST_RESOLUTION (suggest increasing or decreasing resolution)
- FLAG_SUSPICIOUS (major concerns about clustering quality)`;

    // Build agent options with scRNA Stage 5 prompts
    const agentOptions = {
      singleAgent: this.singleAgent,
      gpt5_2_SystemPrompt: SCRNA_STAGE_5_PROMPTS.gpt5_2,
      claudeSystemPrompt: SCRNA_STAGE_5_PROMPTS.claude,
      geminiSystemPrompt: SCRNA_STAGE_5_PROMPTS.gemini
    };

    // Call agents
    const responses = this.sequentialChain
      ? await callAllAgentsSequential(formattedOutput, agentOptions)
      : await callAllAgents(formattedOutput, agentOptions);

    this.log(`GPT-5.2 (Stats): ${responses.gpt5_2.success ? 'responded' : 'failed'}`);
    this.log(`Claude (Pipeline): ${responses.claude.success ? 'responded' : 'failed'}`);
    this.log(`Gemini (Biology): ${responses.gemini.success ? 'responded' : 'failed'}`);

    // Synthesize consensus
    const consensus = synthesizeConsensus(responses, 'scrna_stage5', null);

    // DEBUG: Log consensus details
    console.log('\n[DEBUG] Stage 5 Consensus Details:');
    console.log('  consensus.decision:', consensus.decision);
    console.log('  consensus.votes:', JSON.stringify(consensus.votes));
    console.log('  consensus.agentDecisions:', JSON.stringify(consensus.agentDecisions));

    // Map vote category back to canonical decision
    let canonicalDecision = 'FLAG_SUSPICIOUS';  // Default
    if (consensus.decision === 'approve') {
      canonicalDecision = 'ACCEPT_CLUSTERING';
    } else if (consensus.decision === 'uncertain') {
      canonicalDecision = 'ADJUST_RESOLUTION';
    } else if (consensus.decision === 'reject') {
      canonicalDecision = 'FLAG_SUSPICIOUS';
    } else if (consensus.decision === 'user_decision_required') {
      // If 2/3 agents approved, accept it anyway (lenient for API failures)
      if (consensus.votes.approve >= 2) {
        console.log('[DEBUG] Overriding user_decision_required: 2/3 agents approved');
        canonicalDecision = 'ACCEPT_CLUSTERING';
      } else {
        canonicalDecision = 'FLAG_SUSPICIOUS';
      }
    }

    this.log(`[Coordinator] scRNA Stage 5 Consensus: ${canonicalDecision} (vote: ${consensus.decision})`);

    return {
      decision: canonicalDecision,
      reasoning: consensus.summary,
      confidence: consensus.confidence,
      agentResponses: responses,
      consensus
    };
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
