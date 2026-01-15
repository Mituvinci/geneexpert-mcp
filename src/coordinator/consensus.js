/**
 * Consensus Mechanism
 * Handles voting, disagreement analysis, and decision synthesis
 */

// Number of agents in the voting system
const TOTAL_AGENTS = 3;

/**
 * Analyze agreement between agents
 */
export function analyzeAgreement(responses) {
  const { gpt4, claude, gemini } = responses;

  // Count successful responses
  const successCount = [gpt4.success, claude.success, gemini.success].filter(Boolean).length;

  if (successCount < 2) {
    return {
      level: 'insufficient',
      message: 'Not enough agents responded successfully',
      canProceed: false
    };
  }

  // Simple agreement analysis (can be enhanced with NLP similarity)
  return {
    level: 'pending_analysis',
    successCount,
    agents: {
      stats: gpt4,
      pipeline: claude,
      biology: gemini
    }
  };
}

/**
 * Extract decision from agent response
 * Returns: 'approve', 'reject', 'uncertain', 'automation', 'adaptation', or null
 */
export function extractDecision(agentResponse, decisionType = 'default') {
  if (!agentResponse.success || !agentResponse.content) {
    return null;
  }

  const content = agentResponse.content.toLowerCase();

  // For APPROACH decisions (automation vs adaptation)
  if (decisionType === 'approach_decision') {
    // Check for ADAPTATION first (more specific, should take priority)
    if (content.includes('adaptation') || content.includes('custom script') || content.includes('edge case')) {
      return 'adaptation';
    }

    // Then check for AUTOMATION
    if (content.includes('automation') || content.includes('standard pipeline') || content.includes('proceed with automation')) {
      return 'automation';
    }

    // If unclear, default to automation as safe choice
    return 'automation';
  }

  // For standard decisions (approve/reject)
  // Look for clear approval signals
  if (
    content.includes('approve') ||
    content.includes('proceed') ||
    content.includes('looks good') ||
    content.includes('acceptable') ||
    content.match(/\b(yes|correct|fine)\b/)
  ) {
    return 'approve';
  }

  // Look for clear rejection signals
  if (
    content.includes('reject') ||
    content.includes('do not proceed') ||
    content.includes('not acceptable') ||
    content.includes('remove') ||
    content.includes('incorrect') ||
    content.match(/\b(no|wrong|bad)\b/)
  ) {
    return 'reject';
  }

  // Look for uncertainty
  if (
    content.includes('uncertain') ||
    content.includes('unclear') ||
    content.includes('need more') ||
    content.includes('investigate') ||
    content.includes('check')
  ) {
    return 'uncertain';
  }

  // Default to uncertain if no clear signal
  return 'uncertain';
}

/**
 * Vote on a decision based on agent responses
 *
 * Voting rules (from AGENTS.md):
 * - Run standard tool: Pipeline Agent only
 * - Select threshold: Majority vote (2/3)
 * - Remove sample: Unanimous (3/3)
 * - Change methods: Majority + user approval required
 */
export function vote(responses, decisionType = 'threshold') {
  const { gpt4, claude, gemini } = responses;

  // Extract decisions from each agent (pass decisionType for context)
  const decisions = {
    stats: extractDecision(gpt4, decisionType),
    pipeline: extractDecision(claude, decisionType),
    biology: extractDecision(gemini, decisionType)
  };

  // Count votes (handle both standard and approach decision types)
  const votes = {
    approve: 0,
    reject: 0,
    uncertain: 0,
    automation: 0,
    adaptation: 0,
    null: 0
  };

  Object.values(decisions).forEach(decision => {
    if (decision && votes.hasOwnProperty(decision)) {
      votes[decision]++;
    } else if (decision) {
      votes[decision] = (votes[decision] || 0) + 1;
    } else {
      votes.null++;
    }
  });

  // Apply voting rules based on decision type
  let result;

  switch (decisionType) {
    case 'approach_decision':
      // Majority vote for automation vs adaptation
      const majorityApproach = votes.automation >= votes.adaptation ? 'automation' : 'adaptation';
      const confidence = Math.max(votes.automation, votes.adaptation);

      result = {
        decision: majorityApproach,
        votes,
        agentDecisions: decisions,
        reasoning: `${confidence}/3 agents recommend ${majorityApproach.toUpperCase()}`
      };
      break;

    case 'sample_removal':
      // Require unanimous approval
      result = {
        decision: votes.approve === TOTAL_AGENTS ? 'approve' : 'user_decision_required',
        votes,
        agentDecisions: decisions,
        reasoning: votes.approve === TOTAL_AGENTS
          ? 'All 3 agents agree on sample removal'
          : 'Sample removal requires unanimous agreement - escalating to user'
      };
      break;

    case 'method_change':
      // Require majority + user approval
      result = {
        decision: 'user_approval_required',
        votes,
        agentDecisions: decisions,
        reasoning: 'Method changes require user approval regardless of agent votes'
      };
      break;

    case 'threshold':
    case 'parameter':
      // Majority vote (2/3)
      if (votes.approve >= 2) {
        result = {
          decision: 'approve',
          votes,
          agentDecisions: decisions,
          reasoning: `Majority approval: ${votes.approve}/3 agents agree`
        };
      } else if (votes.reject >= 2) {
        result = {
          decision: 'reject',
          votes,
          agentDecisions: decisions,
          reasoning: `Majority rejection: ${votes.reject}/3 agents disagree`
        };
      } else {
        result = {
          decision: 'user_decision_required',
          votes,
          agentDecisions: decisions,
          reasoning: 'No clear majority - presenting options to user'
        };
      }
      break;

    case 'standard_tool':
      // Pipeline agent only
      result = {
        decision: decisions.pipeline === 'approve' ? 'approve' : decisions.pipeline,
        votes: { pipeline: decisions.pipeline },
        agentDecisions: { pipeline: decisions.pipeline },
        reasoning: 'Standard tool execution - pipeline agent decision only'
      };
      break;

    default:
      // Default to majority
      result = {
        decision: votes.approve >= 2 ? 'approve' : 'user_decision_required',
        votes,
        agentDecisions: decisions,
        reasoning: 'Default majority vote'
      };
  }

  return result;
}

/**
 * Categorize disagreement type
 */
export function categorizeDisagreement(responses) {
  const { gpt4, claude, gemini } = responses;

  const decisions = {
    stats: extractDecision(gpt4),
    pipeline: extractDecision(claude),
    biology: extractDecision(gemini)
  };

  const uniqueDecisions = new Set(Object.values(decisions).filter(Boolean));

  if (uniqueDecisions.size === 1) {
    return {
      type: 'full_agreement',
      severity: 'none',
      description: 'All agents agree'
    };
  }

  if (uniqueDecisions.size === 2) {
    // Find dissenting agent
    const decisionCounts = {};
    Object.entries(decisions).forEach(([agent, decision]) => {
      decisionCounts[decision] = decisionCounts[decision] || [];
      decisionCounts[decision].push(agent);
    });

    const minority = Object.entries(decisionCounts).find(([_, agents]) => agents.length === 1);

    return {
      type: 'minor_disagreement',
      severity: 'low',
      description: `One agent (${minority?.[1][0]}) disagrees`,
      minorityAgent: minority?.[1][0],
      minorityDecision: minority?.[0]
    };
  }

  if (uniqueDecisions.size === TOTAL_AGENTS) {
    return {
      type: 'major_disagreement',
      severity: 'high',
      description: 'All three agents have different opinions',
      requiresUserInput: true
    };
  }

  return {
    type: 'unclear',
    severity: 'medium',
    description: 'Unable to categorize disagreement pattern'
  };
}

/**
 * Synthesize a consensus summary from all agent responses
 *
 * @param {Object} responses - Agent responses from GPT-4, Claude, Gemini
 * @param {string} decisionType - Type of decision (approach_decision, threshold, etc.)
 * @param {string} decision_id - Unique identifier for this decision (for evaluation logging)
 *                               Format: "{dataset}_{step}_{type}" e.g., "DA0036_step0_mode"
 */
export function synthesizeConsensus(responses, decisionType = 'threshold', decision_id = null) {
  const voteResult = vote(responses, decisionType);
  const disagreement = categorizeDisagreement(responses);

  // Calculate disagreement score (0 = full agreement, 0.67 = maximum disagreement)
  // Decision-type aware: only consider relevant vote categories
  const disagreementScore = calculateDisagreementScore(voteResult.votes, decisionType);

  const summary = {
    // decision_id for evaluation (join with ground_truth offline)
    decision_id: decision_id,
    timestamp: new Date().toISOString(),
    decision_type: decisionType,
    decision: voteResult.decision,
    confidence: calculateConfidence(voteResult.votes),
    disagreement_score: disagreementScore,
    votes: voteResult.votes,
    agentDecisions: voteResult.agentDecisions,
    disagreement,
    reasoning: voteResult.reasoning,
    recommendations: generateRecommendations(responses, voteResult, disagreement)
  };

  return summary;
}

/**
 * Calculate disagreement score (0 = full agreement, ~0.67 = max disagreement)
 * Decision-type aware: only considers relevant vote categories
 */
function calculateDisagreementScore(votes, decisionType) {
  let relevantVotes = [];

  if (decisionType === 'approach_decision') {
    // For approach decisions: only automation vs adaptation
    relevantVotes = [votes.automation || 0, votes.adaptation || 0];
  } else {
    // For threshold/parameter/sample_removal: approve/reject/uncertain
    relevantVotes = [
      votes.approve || 0,
      votes.reject || 0,
      votes.uncertain || 0
    ];
  }

  const maxVotes = Math.max(...relevantVotes);
  return 1 - maxVotes / TOTAL_AGENTS;
}

/**
 * Calculate confidence score (0-1) based on vote distribution
 */
function calculateConfidence(votes) {
  // For approach decisions (automation vs adaptation)
  if (votes.automation !== undefined || votes.adaptation !== undefined) {
    const total = votes.automation + votes.adaptation;

    if (total === 0) return 0;

    // Unanimous = 1.0
    if (votes.automation === TOTAL_AGENTS || votes.adaptation === TOTAL_AGENTS) return 1.0;

    // 2-1 majority = 0.67
    if (votes.automation === 2 || votes.adaptation === 2) return 0.67;

    // All uncertain or split = 0.33
    return 0.33;
  }

  // For standard decisions (approve/reject/uncertain)
  const total = votes.approve + votes.reject + votes.uncertain;

  if (total === 0) return 0;

  // Unanimous = 1.0
  if (votes.approve === TOTAL_AGENTS || votes.reject === TOTAL_AGENTS) return 1.0;

  // 2-1 majority = 0.67
  if (votes.approve === 2 || votes.reject === 2) return 0.67;

  // All uncertain or split = 0.33
  return 0.33;
}

/**
 * Generate actionable recommendations based on consensus
 */
function generateRecommendations(responses, voteResult, disagreement) {
  const recommendations = [];

  if (disagreement.type === 'full_agreement') {
    recommendations.push({
      priority: 'high',
      action: 'proceed',
      reason: 'All agents agree - safe to proceed'
    });
  }

  if (disagreement.type === 'minor_disagreement') {
    recommendations.push({
      priority: 'medium',
      action: 'note_dissent',
      reason: `Proceed with caution - ${disagreement.minorityAgent} has concerns`,
      note: `${disagreement.minorityAgent} suggests: ${disagreement.minorityDecision}`
    });
  }

  if (disagreement.type === 'major_disagreement') {
    recommendations.push({
      priority: 'critical',
      action: 'user_review_required',
      reason: 'Significant disagreement detected - manual review needed',
      options: Object.entries(voteResult.agentDecisions).map(([agent, decision]) => ({
        agent,
        recommendation: decision
      }))
    });
  }

  if (voteResult.votes.uncertain >= 2) {
    recommendations.push({
      priority: 'high',
      action: 'gather_more_info',
      reason: 'Multiple agents uncertain - may need additional analysis'
    });
  }

  return recommendations;
}

/**
 * Format consensus for display to user
 */
export function formatConsensusReport(consensus, question) {
  let report = `\n${'='.repeat(60)}\n`;
  report += `MULTI-AGENT CONSENSUS REPORT\n`;
  report += `${'='.repeat(60)}\n\n`;

  report += `Question: ${question}\n\n`;

  report += `Decision: ${consensus.decision.toUpperCase()}\n`;
  report += `Confidence: ${(consensus.confidence * 100).toFixed(0)}%\n\n`;

  report += `Votes:\n`;
  report += `  Approve:   ${consensus.votes.approve || 0}/3\n`;
  report += `  Reject:    ${consensus.votes.reject || 0}/3\n`;
  report += `  Uncertain: ${consensus.votes.uncertain || 0}/3\n\n`;

  report += `Agent Decisions:\n`;
  Object.entries(consensus.agentDecisions).forEach(([agent, decision]) => {
    const icon = decision === 'approve' ? '✓' : decision === 'reject' ? '✗' : '?';
    report += `  ${icon} ${agent.padEnd(10)} → ${decision || 'no response'}\n`;
  });

  report += `\nDisagreement Analysis:\n`;
  report += `  Type: ${consensus.disagreement.type}\n`;
  report += `  Severity: ${consensus.disagreement.severity}\n`;
  report += `  ${consensus.disagreement.description}\n\n`;

  if (consensus.recommendations.length > 0) {
    report += `Recommendations:\n`;
    consensus.recommendations.forEach((rec, i) => {
      report += `  ${i + 1}. [${rec.priority.toUpperCase()}] ${rec.action}\n`;
      report += `     ${rec.reason}\n`;
      if (rec.note) {
        report += `     Note: ${rec.note}\n`;
      }
    });
  }

  report += `\n${'='.repeat(60)}\n`;

  return report;
}

/**
 * Extract structured modifications from agent responses (Option 2)
 * Returns specific changes agents recommend for ADAPTATION
 */
export function extractStructuredRecommendations(responses, dataInfo, config) {
  const modifications = {
    statistical_adjustments: [],
    pipeline_modifications: [],
    quality_checks: [],
    parameter_changes: {}
  };

  // Analyze all agent responses for specific recommendations
  const allContent = [
    responses.gpt4?.content || '',
    responses.claude?.content || '',
    responses.gemini?.content || ''
  ].join(' ').toLowerCase();

  // 1. STATISTICAL ADJUSTMENTS (from sample size, power analysis)

  // Small sample size (n<3) → use exact test
  if (dataInfo.groups) {
    const minN = Math.min(...Object.values(dataInfo.groups).map(g => g.length));
    if (minN < 3) {
      modifications.statistical_adjustments.push({
        type: 'use_exact_test',
        reason: `Small sample size (n=${minN}) requires exact test instead of quasi-likelihood`,
        parameters: { method: 'exactTest' }
      });

      // Relax FDR threshold for low power
      modifications.parameter_changes.fdr_threshold = 0.1;  // instead of 0.05
    }
  }

  // Batch effects mentioned
  if (allContent.includes('batch') || allContent.includes('confound')) {
    modifications.quality_checks.push({
      type: 'batch_detection',
      step: 'before_de_analysis',
      command: 'check_batch_effects.R'
    });
  }

  // 2. PIPELINE MODIFICATIONS

  // Outlier detection recommended
  if (allContent.includes('outlier') || allContent.includes('pca') || allContent.includes('mds')) {
    modifications.quality_checks.push({
      type: 'outlier_detection',
      step: 'before_de_analysis',
      command: 'qc_plots.R',
      check: 'visual_inspection'
    });
  }

  // Conservative filtering
  if (allContent.includes('conservative') || allContent.includes('stringent') || allContent.includes('careful')) {
    modifications.parameter_changes.cpm_threshold = 0.5;  // instead of 1.0
    modifications.parameter_changes.min_samples = 2;
  }

  // Relaxed filtering (for discovery)
  if (allContent.includes('relaxed') || allContent.includes('exploratory') || allContent.includes('discovery')) {
    modifications.parameter_changes.cpm_threshold = 0.3;
    modifications.parameter_changes.logfc_threshold = 0.5;  // instead of 1.0
  }

  // 3. QUALITY CHECKS

  // Add extra QC for low sample size
  if (dataInfo.samples && dataInfo.samples.length <= 4) {
    modifications.quality_checks.push({
      type: 'enhanced_qc',
      step: 'after_alignment',
      checks: ['alignment_rate', 'duplication_rate', 'library_complexity']
    });
  }

  return modifications;
}

export default {
  analyzeAgreement,
  extractDecision,
  vote,
  categorizeDisagreement,
  synthesizeConsensus,
  formatConsensusReport,
  extractStructuredRecommendations
};
