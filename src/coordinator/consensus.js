/**
 * Consensus Mechanism
 * Handles voting, disagreement analysis, and decision synthesis
 */

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
 * Returns: 'approve', 'reject', 'uncertain', or null
 */
export function extractDecision(agentResponse) {
  if (!agentResponse.success || !agentResponse.content) {
    return null;
  }

  const content = agentResponse.content.toLowerCase();

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

  // Extract decisions from each agent
  const decisions = {
    stats: extractDecision(gpt4),
    pipeline: extractDecision(claude),
    biology: extractDecision(gemini)
  };

  // Count votes
  const votes = {
    approve: 0,
    reject: 0,
    uncertain: 0,
    null: 0
  };

  Object.values(decisions).forEach(decision => {
    if (decision) {
      votes[decision]++;
    } else {
      votes.null++;
    }
  });

  // Apply voting rules based on decision type
  let result;

  switch (decisionType) {
    case 'sample_removal':
      // Require unanimous approval
      result = {
        decision: votes.approve === 3 ? 'approve' : 'user_decision_required',
        votes,
        agentDecisions: decisions,
        reasoning: votes.approve === 3
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

  if (uniqueDecisions.size === 3) {
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
 */
export function synthesizeConsensus(responses, decisionType = 'threshold') {
  const voteResult = vote(responses, decisionType);
  const disagreement = categorizeDisagreement(responses);

  const summary = {
    decision: voteResult.decision,
    confidence: calculateConfidence(voteResult.votes),
    votes: voteResult.votes,
    agentDecisions: voteResult.agentDecisions,
    disagreement,
    reasoning: voteResult.reasoning,
    recommendations: generateRecommendations(responses, voteResult, disagreement)
  };

  return summary;
}

/**
 * Calculate confidence score (0-1) based on vote distribution
 */
function calculateConfidence(votes) {
  const total = votes.approve + votes.reject + votes.uncertain;

  if (total === 0) return 0;

  // Unanimous = 1.0
  if (votes.approve === 3 || votes.reject === 3) return 1.0;

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

export default {
  analyzeAgreement,
  extractDecision,
  vote,
  categorizeDisagreement,
  synthesizeConsensus,
  formatConsensusReport
};
