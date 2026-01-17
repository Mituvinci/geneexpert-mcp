/**
 * Response Parser - Extract structured data from agent text responses
 *
 * Purpose: Extract decision, confidence label, and convert to numeric scores
 * from unstructured text content returned by LLM agents.
 */

/**
 * Confidence label to numeric score mapping
 */
const CONFIDENCE_MAP = {
  'HIGH': 0.9,
  'MEDIUM': 0.7,
  'LOW': 0.5,
  'UNCERTAIN': 0.3
};

/**
 * Convert numeric confidence score to label
 * @param {number} score - Confidence score (0-1)
 * @returns {string} - Confidence label (HIGH, MEDIUM, LOW, UNCERTAIN)
 */
export function getConfidenceLabel(score) {
  if (score >= 0.85) return 'HIGH';
  if (score >= 0.65) return 'MEDIUM';
  if (score >= 0.45) return 'LOW';
  return 'UNCERTAIN';
}

/**
 * Convert confidence label to numeric score
 * @param {string} label - Confidence label (HIGH, MEDIUM, LOW, UNCERTAIN)
 * @returns {number} - Confidence score (0-1)
 */
export function convertConfidenceToScore(label) {
  const normalized = label?.toUpperCase().trim();
  return CONFIDENCE_MAP[normalized] || 0.5; // Default to 0.5 (MEDIUM) if not recognized
}

/**
 * Extract decision from agent content
 * Looks for patterns like:
 * - "Recommendation: ADAPTATION"
 * - "Final Recommendation: AUTOMATION"
 * - "Decision: UNCERTAIN"
 *
 * @param {string} content - Full text response from agent
 * @returns {string|null} - Extracted decision or null if not found
 */
export function extractDecision(content) {
  if (!content) return null;

  // Pattern 1: "Recommendation: AUTOMATION"
  const recommendationMatch = content.match(/(?:Final\s+)?Recommendation:\s*([A-Z_]+)/i);
  if (recommendationMatch) {
    return recommendationMatch[1].toUpperCase();
  }

  // Pattern 2: "Decision: ADAPTATION"
  const decisionMatch = content.match(/(?:Final\s+)?Decision:\s*([A-Z_]+)/i);
  if (decisionMatch) {
    return decisionMatch[1].toUpperCase();
  }

  // Pattern 3: Look for keywords in the text
  const upperContent = content.toUpperCase();
  if (upperContent.includes('AUTOMATION')) return 'AUTOMATION';
  if (upperContent.includes('ADAPTATION')) return 'ADAPTATION';
  if (upperContent.includes('UNCERTAIN')) return 'UNCERTAIN';

  return null;
}

/**
 * Extract confidence label from agent content
 * Looks for patterns like:
 * - "Confidence: HIGH"
 * - "Overall Confidence: MEDIUM"
 * - "Confidence Level: LOW"
 *
 * @param {string} content - Full text response from agent
 * @returns {string|null} - Extracted confidence label or null if not found
 */
export function extractConfidenceLabel(content) {
  if (!content) return null;

  // Pattern: "Confidence: HIGH" or "Overall Confidence: MEDIUM"
  const confidenceMatch = content.match(/(?:Overall\s+)?Confidence(?:\s+Level)?:\s*([A-Z]+)/i);
  if (confidenceMatch) {
    return confidenceMatch[1].toUpperCase();
  }

  return null;
}

/**
 * Parse full agent response and extract all structured fields
 * @param {string} content - Full text response from agent
 * @returns {object} - Extracted structured data
 */
export function parseAgentResponse(content) {
  const decision = extractDecision(content);
  const confidenceLabel = extractConfidenceLabel(content);
  const confidenceScore = confidenceLabel ? convertConfidenceToScore(confidenceLabel) : null;

  return {
    extracted_decision: decision,
    confidence_label: confidenceLabel,
    confidence_score: confidenceScore
  };
}
