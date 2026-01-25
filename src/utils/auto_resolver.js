/**
 * Auto-Resolution Logic for Agent Disagreements
 *
 * Implements tiered escalation:
 * - Minor disagreement (score < 0.3): Auto-resolve with median
 * - Moderate disagreement (0.3 <= score < 0.6): Use highest confidence agent
 * - Major disagreement (score >= 0.6): Escalate to user
 */

import {
  calculateDisagreementScore,
  findHighestConfidenceAgent,
  median
} from './consensus_helper.js';

// ============================================
// MAIN AUTO-RESOLUTION FUNCTION
// ============================================

/**
 * Tiered auto-resolution for agent disagreements
 *
 * @param {Object} agentResponses - All 3 agent responses
 * @param {Object} consensus - Consensus result
 * @param {string} decisionType - 'thresholds', 'binary', 'numeric', 'categorical'
 * @param {string} autoResolveMode - 'auto', 'median', 'confidence', 'user'
 * @param {number} escalationThreshold - Disagreement score threshold for escalation (default 0.6)
 * @returns {Object} { autoResolved: boolean, decision: any, method: string, ... }
 */
export async function autoResolveDecision(
  agentResponses,
  consensus,
  decisionType,
  autoResolveMode = 'auto',
  escalationThreshold = 0.6
) {
  // If user mode, never auto-resolve
  if (autoResolveMode === 'user') {
    return { autoResolved: false, reason: 'User mode selected (--auto-resolve user)' };
  }

  // Calculate disagreement severity
  const disagreementScore = calculateDisagreementScore(agentResponses, decisionType);

  // If mode is 'median', always use median (unless major disagreement)
  if (autoResolveMode === 'median' && disagreementScore < escalationThreshold) {
    const decision = calculateMedianDecision(agentResponses, decisionType);
    return {
      autoResolved: true,
      decision,
      method: 'median',
      disagreementScore,
      reasoning: `Median mode selected - using median/average of all agent recommendations (disagreement: ${disagreementScore.toFixed(2)})`
    };
  }

  // If mode is 'confidence', always use highest confidence (unless major disagreement)
  if (autoResolveMode === 'confidence' && disagreementScore < escalationThreshold) {
    const { agent, confidence, allConfidences } = findHighestConfidenceAgent(agentResponses);
    const decision = extractDecisionFromAgent(agentResponses[agent], decisionType);
    return {
      autoResolved: true,
      decision,
      method: 'confidence',
      selectedAgent: agent,
      agentConfidence: confidence,
      allConfidences,
      disagreementScore,
      reasoning: `Confidence mode selected - using ${agent} (confidence: ${confidence}, disagreement: ${disagreementScore.toFixed(2)})`
    };
  }

  // Auto mode: Tiered escalation based on disagreement severity
  if (autoResolveMode === 'auto') {
    // TIER 1: Minor disagreement → Use median
    if (disagreementScore < 0.3) {
      const decision = calculateMedianDecision(agentResponses, decisionType);
      return {
        autoResolved: true,
        decision,
        method: 'median',
        disagreementScore,
        reasoning: `Minor disagreement (${disagreementScore.toFixed(2)}) - using median of all agent recommendations`
      };
    }

    // TIER 2: Moderate disagreement → Use highest confidence agent
    // For numeric/threshold decisions, skip consensus.confidence check (not applicable)
    const skipConsensusCheck = (decisionType === 'thresholds' || decisionType === 'numeric');
    const consensusOK = skipConsensusCheck || (consensus && consensus.confidence >= 0.5);

    if (disagreementScore >= 0.3 && disagreementScore < escalationThreshold && consensusOK) {
      const { agent, confidence, allConfidences } = findHighestConfidenceAgent(agentResponses);
      const decision = extractDecisionFromAgent(agentResponses[agent], decisionType);
      return {
        autoResolved: true,
        decision,
        method: 'confidence',
        selectedAgent: agent,
        agentConfidence: confidence,
        allConfidences,
        disagreementScore,
        reasoning: `Moderate disagreement (${disagreementScore.toFixed(2)}) - selected ${agent} (highest confidence: ${confidence})`
      };
    }

    // TIER 3: Major disagreement → Escalate to user
    return {
      autoResolved: false,
      disagreementScore,
      reasoning: `Major disagreement (${disagreementScore.toFixed(2)}) or low confidence - requires user decision`
    };
  }

  // Unknown mode, escalate
  return {
    autoResolved: false,
    reason: `Unknown auto-resolve mode: ${autoResolveMode}`
  };
}

// ============================================
// MEDIAN DECISION CALCULATION
// ============================================

/**
 * Calculate median decision from all 3 agents
 */
function calculateMedianDecision(agentResponses, decisionType) {
  switch (decisionType) {
    case 'thresholds':
      return calculateMedianThresholds(agentResponses);
    case 'binary':
      return calculateMajorityBinary(agentResponses);
    case 'numeric':
      return calculateMedianNumeric(agentResponses);
    case 'categorical':
      return calculateMajorityCategorical(agentResponses);
    default:
      throw new Error(`Unknown decision type: ${decisionType}`);
  }
}

/**
 * Calculate median thresholds for Stage 2
 */
function calculateMedianThresholds(agentResponses) {
  const extractThreshold = (text, pattern) => {
    const match = text?.match(pattern);
    return match ? parseInt(match[1]) : null;
  };

  const mins = [
    extractThreshold(agentResponses.gpt5_2?.content, /nFeature[_\s]*min[:\s]*(\d+)/i),
    extractThreshold(agentResponses.claude?.content, /nFeature[_\s]*min[:\s]*(\d+)/i),
    extractThreshold(agentResponses.gemini?.content, /nFeature[_\s]*min[:\s]*(\d+)/i)
  ].filter(v => v !== null);

  const maxs = [
    extractThreshold(agentResponses.gpt5_2?.content, /nFeature[_\s]*max[:\s]*(\d+)/i),
    extractThreshold(agentResponses.claude?.content, /nFeature[_\s]*max[:\s]*(\d+)/i),
    extractThreshold(agentResponses.gemini?.content, /nFeature[_\s]*max[:\s]*(\d+)/i)
  ].filter(v => v !== null);

  const mts = [
    extractThreshold(agentResponses.gpt5_2?.content, /percent[_\s]*mt[_\s]*max[:\s]*(\d+)/i),
    extractThreshold(agentResponses.claude?.content, /percent[_\s]*mt[_\s]*max[:\s]*(\d+)/i),
    extractThreshold(agentResponses.gemini?.content, /percent[_\s]*mt[_\s]*max[:\s]*(\d+)/i)
  ].filter(v => v !== null);

  return {
    nFeature_min: mins.length > 0 ? Math.round(median(mins)) : 200,
    nFeature_max: maxs.length > 0 ? Math.round(median(maxs)) : 6000,
    percent_mt_max: mts.length > 0 ? Math.round(median(mts)) : 10
  };
}

/**
 * Calculate majority for binary decisions (Stage 3A)
 */
function calculateMajorityBinary(agentResponses) {
  const extractBinaryDecision = (text) => {
    if (!text) return null;
    if (/REMOVE_CELL_CYCLE/i.test(text)) return 'REMOVE_CELL_CYCLE';
    if (/SKIP_CELL_CYCLE/i.test(text)) return 'SKIP_CELL_CYCLE';
    return null;
  };

  const decisions = [
    extractBinaryDecision(agentResponses.gpt5_2?.content),
    extractBinaryDecision(agentResponses.claude?.content),
    extractBinaryDecision(agentResponses.gemini?.content)
  ].filter(d => d !== null);

  // Count votes
  const removeVotes = decisions.filter(d => d === 'REMOVE_CELL_CYCLE').length;
  const skipVotes = decisions.filter(d => d === 'SKIP_CELL_CYCLE').length;

  // Return majority
  return removeVotes > skipVotes ? 'REMOVE_CELL_CYCLE' : 'SKIP_CELL_CYCLE';
}

/**
 * Calculate median numeric value (Stage 4: PC range)
 */
function calculateMedianNumeric(agentResponses) {
  const extractMaxPC = (text) => {
    const match = text?.match(/max[_\s]*pc[:\s]*(\d+)/i);
    return match ? parseInt(match[1]) : null;
  };

  const maxPCs = [
    extractMaxPC(agentResponses.gpt5_2?.content),
    extractMaxPC(agentResponses.claude?.content),
    extractMaxPC(agentResponses.gemini?.content)
  ].filter(v => v !== null);

  return {
    min_pc: 1,
    max_pc: maxPCs.length > 0 ? Math.round(median(maxPCs)) : 20
  };
}

/**
 * Calculate majority categorical decision (Stage 5: clustering)
 */
function calculateMajorityCategorical(agentResponses) {
  const extractCategorical = (text) => {
    if (!text) return null;
    if (/ACCEPT_CLUSTERING/i.test(text)) return 'ACCEPT_CLUSTERING';
    if (/ADJUST_RESOLUTION/i.test(text)) return 'ADJUST_RESOLUTION';
    if (/FLAG_SUSPICIOUS/i.test(text)) return 'FLAG_SUSPICIOUS';
    return null;
  };

  const decisions = [
    extractCategorical(agentResponses.gpt5_2?.content),
    extractCategorical(agentResponses.claude?.content),
    extractCategorical(agentResponses.gemini?.content)
  ].filter(d => d !== null);

  // Count votes for each option
  const votes = { ACCEPT_CLUSTERING: 0, ADJUST_RESOLUTION: 0, FLAG_SUSPICIOUS: 0 };
  decisions.forEach(d => votes[d]++);

  // Return majority
  const majority = Object.keys(votes).reduce((a, b) => votes[a] > votes[b] ? a : b);
  return majority;
}

// ============================================
// EXTRACT DECISION FROM SPECIFIC AGENT
// ============================================

/**
 * Extract decision from a specific agent's response
 */
function extractDecisionFromAgent(agentResponse, decisionType) {
  switch (decisionType) {
    case 'thresholds':
      return extractThresholdsFromAgent(agentResponse);
    case 'binary':
      return extractBinaryFromAgent(agentResponse);
    case 'numeric':
      return extractNumericFromAgent(agentResponse);
    case 'categorical':
      return extractCategoricalFromAgent(agentResponse);
    default:
      throw new Error(`Unknown decision type: ${decisionType}`);
  }
}

function extractThresholdsFromAgent(agentResponse) {
  const extractThreshold = (text, pattern) => {
    const match = text?.match(pattern);
    return match ? parseInt(match[1]) : null;
  };

  return {
    nFeature_min: extractThreshold(agentResponse.content, /nFeature[_\s]*min[:\s]*(\d+)/i) || 200,
    nFeature_max: extractThreshold(agentResponse.content, /nFeature[_\s]*max[:\s]*(\d+)/i) || 6000,
    percent_mt_max: extractThreshold(agentResponse.content, /percent[_\s]*mt[_\s]*max[:\s]*(\d+)/i) || 10
  };
}

function extractBinaryFromAgent(agentResponse) {
  if (/REMOVE_CELL_CYCLE/i.test(agentResponse.content)) return 'REMOVE_CELL_CYCLE';
  if (/SKIP_CELL_CYCLE/i.test(agentResponse.content)) return 'SKIP_CELL_CYCLE';
  return 'SKIP_CELL_CYCLE'; // Default
}

function extractNumericFromAgent(agentResponse) {
  const extractMaxPC = (text) => {
    const match = text?.match(/max[_\s]*pc[:\s]*(\d+)/i);
    return match ? parseInt(match[1]) : null;
  };

  return {
    min_pc: 1,
    max_pc: extractMaxPC(agentResponse.content) || 20
  };
}

function extractCategoricalFromAgent(agentResponse) {
  if (/ACCEPT_CLUSTERING/i.test(agentResponse.content)) return 'ACCEPT_CLUSTERING';
  if (/ADJUST_RESOLUTION/i.test(agentResponse.content)) return 'ADJUST_RESOLUTION';
  if (/FLAG_SUSPICIOUS/i.test(agentResponse.content)) return 'FLAG_SUSPICIOUS';
  return 'ACCEPT_CLUSTERING'; // Default
}

// ============================================
// EXPORT
// ============================================

export default {
  autoResolveDecision
};
