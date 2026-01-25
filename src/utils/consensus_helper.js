/**
 * Consensus Helper Utilities
 *
 * Functions for:
 * - Calculating disagreement scores
 * - Extracting agent confidence
 * - Finding highest confidence agent
 * - Statistical utilities (median, std, mean)
 */

// ============================================
// STATISTICAL UTILITIES
// ============================================

export function mean(values) {
  return values.reduce((a, b) => a + b, 0) / values.length;
}

export function median(values) {
  const sorted = [...values].sort((a, b) => a - b);
  const mid = Math.floor(sorted.length / 2);
  return sorted.length % 2 === 0
    ? (sorted[mid - 1] + sorted[mid]) / 2
    : sorted[mid];
}

export function std(values) {
  const avg = mean(values);
  const squareDiffs = values.map(v => Math.pow(v - avg, 2));
  return Math.sqrt(mean(squareDiffs));
}

// ============================================
// CONFIDENCE EXTRACTION
// ============================================

/**
 * Parse confidence from agent response
 * Looks for: "Confidence: HIGH/MEDIUM/LOW"
 */
export function extractConfidence(agentResponse) {
  if (!agentResponse || !agentResponse.content) return 0.5;

  const confidenceMatch = agentResponse.content.match(/Confidence:\s*(HIGH|MEDIUM|LOW)/i);
  if (!confidenceMatch) return 0.5; // Default to medium

  const confidence = confidenceMatch[1].toUpperCase();
  return { HIGH: 1.0, MEDIUM: 0.6, LOW: 0.3 }[confidence] || 0.5;
}

/**
 * Find agent with highest confidence
 */
export function findHighestConfidenceAgent(agentResponses) {
  const confidences = {
    gpt5_2: extractConfidence(agentResponses.gpt5_2),
    claude: extractConfidence(agentResponses.claude),
    gemini: extractConfidence(agentResponses.gemini)
  };

  const highest = Object.keys(confidences).reduce((a, b) =>
    confidences[a] > confidences[b] ? a : b
  );

  return { agent: highest, confidence: confidences[highest], allConfidences: confidences };
}

// ============================================
// DISAGREEMENT SCORING
// ============================================

/**
 * Calculate disagreement severity score (0-1)
 * 0 = perfect agreement, 1 = maximum disagreement
 */
export function calculateDisagreementScore(agentResponses, decisionType) {
  switch (decisionType) {
    case 'thresholds':
      return calculateThresholdDisagreement(agentResponses);
    case 'binary':  // REMOVE_CELL_CYCLE vs SKIP_CELL_CYCLE
      return calculateBinaryDisagreement(agentResponses);
    case 'numeric':  // PC range selection
      return calculateNumericDisagreement(agentResponses);
    case 'categorical':  // ACCEPT/ADJUST/FLAG
      return calculateCategoricalDisagreement(agentResponses);
    default:
      return 0.5; // Unknown type, moderate score
  }
}

/**
 * Calculate disagreement for QC thresholds (Stage 2)
 */
function calculateThresholdDisagreement(agentResponses) {
  // Extract thresholds from all 3 agents
  const extractThreshold = (text, pattern) => {
    const match = text?.match(pattern);
    return match ? parseInt(match[1]) : null;
  };

  const gptMin = extractThreshold(agentResponses.gpt5_2?.content, /nFeature[_\s]*min[:\s]*(\d+)/i);
  const gptMax = extractThreshold(agentResponses.gpt5_2?.content, /nFeature[_\s]*max[:\s]*(\d+)/i);
  const gptMT = extractThreshold(agentResponses.gpt5_2?.content, /percent[_\s]*mt[_\s]*max[:\s]*(\d+)/i);

  const claudeMin = extractThreshold(agentResponses.claude?.content, /nFeature[_\s]*min[:\s]*(\d+)/i);
  const claudeMax = extractThreshold(agentResponses.claude?.content, /nFeature[_\s]*max[:\s]*(\d+)/i);
  const claudeMT = extractThreshold(agentResponses.claude?.content, /percent[_\s]*mt[_\s]*max[:\s]*(\d+)/i);

  const geminiMin = extractThreshold(agentResponses.gemini?.content, /nFeature[_\s]*min[:\s]*(\d+)/i);
  const geminiMax = extractThreshold(agentResponses.gemini?.content, /nFeature[_\s]*max[:\s]*(\d+)/i);
  const geminiMT = extractThreshold(agentResponses.gemini?.content, /percent[_\s]*mt[_\s]*max[:\s]*(\d+)/i);

  // Calculate coefficient of variation (CV) for each parameter
  const cvs = [];

  if (gptMin && claudeMin && geminiMin) {
    const mins = [gptMin, claudeMin, geminiMin];
    const minCV = std(mins) / mean(mins);
    cvs.push(minCV);
  }

  if (gptMax && claudeMax && geminiMax) {
    const maxs = [gptMax, claudeMax, geminiMax];
    const maxCV = std(maxs) / mean(maxs);
    cvs.push(maxCV);
  }

  if (gptMT && claudeMT && geminiMT) {
    const mts = [gptMT, claudeMT, geminiMT];
    const mtCV = std(mts) / mean(mts);
    cvs.push(mtCV);
  }

  if (cvs.length === 0) return 0.5; // No data, default to moderate

  // Average CV across all parameters
  const avgCV = mean(cvs);

  // Normalize to 0-1 scale (CV > 0.4 = high disagreement)
  return Math.min(avgCV / 0.4, 1.0);
}

/**
 * Calculate disagreement for binary decisions (Stage 3A: REMOVE vs SKIP)
 */
function calculateBinaryDisagreement(agentResponses) {
  // Count how many agents chose each option
  const extractBinaryDecision = (text) => {
    if (!text) return null;
    if (/REMOVE_CELL_CYCLE/i.test(text)) return 'REMOVE';
    if (/SKIP_CELL_CYCLE/i.test(text)) return 'SKIP';
    return null;
  };

  const decisions = [
    extractBinaryDecision(agentResponses.gpt5_2?.content),
    extractBinaryDecision(agentResponses.claude?.content),
    extractBinaryDecision(agentResponses.gemini?.content)
  ].filter(d => d !== null);

  if (decisions.length === 0) return 0.5;

  // Count votes
  const removeVotes = decisions.filter(d => d === 'REMOVE').length;
  const skipVotes = decisions.filter(d => d === 'SKIP').length;

  // Perfect agreement (3-0 or 0-3) → 0.0
  // Split decision (2-1) → 0.33
  // No majority (1-1-1 if 3rd is different) → 1.0
  if (removeVotes === 3 || skipVotes === 3) return 0.0;
  if (removeVotes === 2 || skipVotes === 2) return 0.33;
  return 1.0; // Tie or no clear pattern
}

/**
 * Calculate disagreement for numeric ranges (Stage 4: PC selection)
 */
function calculateNumericDisagreement(agentResponses) {
  // Extract max_pc from all 3 agents
  const extractMaxPC = (text) => {
    const match = text?.match(/max[_\s]*pc[:\s]*(\d+)/i);
    return match ? parseInt(match[1]) : null;
  };

  const maxPCs = [
    extractMaxPC(agentResponses.gpt5_2?.content),
    extractMaxPC(agentResponses.claude?.content),
    extractMaxPC(agentResponses.gemini?.content)
  ].filter(v => v !== null);

  if (maxPCs.length === 0) return 0.5;

  // Calculate CV
  const cv = std(maxPCs) / mean(maxPCs);

  // Normalize (CV > 0.2 = high disagreement for PC ranges)
  return Math.min(cv / 0.2, 1.0);
}

/**
 * Calculate disagreement for categorical decisions (Stage 5: ACCEPT/ADJUST/FLAG)
 */
function calculateCategoricalDisagreement(agentResponses) {
  // Extract categorical decision
  const extractCategorical = (text) => {
    if (!text) return null;
    if (/ACCEPT_CLUSTERING/i.test(text)) return 'ACCEPT';
    if (/ADJUST_RESOLUTION/i.test(text)) return 'ADJUST';
    if (/FLAG_SUSPICIOUS/i.test(text)) return 'FLAG';
    return null;
  };

  const decisions = [
    extractCategorical(agentResponses.gpt5_2?.content),
    extractCategorical(agentResponses.claude?.content),
    extractCategorical(agentResponses.gemini?.content)
  ].filter(d => d !== null);

  if (decisions.length === 0) return 0.5;

  // Count unique decisions
  const uniqueDecisions = new Set(decisions);

  // All agree → 0.0
  // 2 agree, 1 disagrees → 0.33
  // All disagree → 1.0
  if (uniqueDecisions.size === 1) return 0.0;
  if (uniqueDecisions.size === 2) return 0.33;
  return 1.0;
}

// ============================================
// EXPORT ALL
// ============================================

export default {
  mean,
  median,
  std,
  extractConfidence,
  findHighestConfidenceAgent,
  calculateDisagreementScore
};
