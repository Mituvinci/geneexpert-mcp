/**
 * Canonical Decision Vocabulary
 *
 * Single source of truth for all valid decisions at each stage.
 * All parsing/matching/logging MUST use these exact strings.
 */

export const DECISION_VOCABULARY = {
  // Stage 1: FASTQ Validation
  stage1: {
    canonical: ['PASS', 'PASS_WITH_WARNING', 'FAIL'],
    synonyms: {
      PASS_WITH_WARNING: ['pass with warning', 'pass_with_warnings', 'pass_with_warn', 'warning', 'minor issues', 'with warning'],
      FAIL: ['fail', 'failed', 'reject', 'abort', 'unacceptable', 'quality failed'],
      PASS: ['pass', 'proceed', 'looks good', 'acceptable', 'quality ok']
    }
  },

  // Stage 2: Alignment + QC
  stage2: {
    canonical: ['PASS_ALL', 'REMOVE_SAMPLES', 'ABORT'],
    synonyms: {
      PASS_ALL: ['pass_all', 'pass all', 'keep all', 'all samples pass', 'proceed with all'],
      REMOVE_SAMPLES: ['remove_samples', 'remove samples', 'filter samples', 'exclude samples', 'drop samples'],
      ABORT: ['abort', 'stop', 'terminate', 'do not proceed', 'fail analysis']
    }
  },

  // Stage 3: Quantification + QC (DE Method)
  stage3_de_method: {
    canonical: ['simpleEdger', 'batch_effect_edger'],
    synonyms: {
      simpleEdger: ['simpleedger', 'simple edger', 'simple_edger', 'standard edger', 'no batch', 'standard analysis'],
      batch_effect_edger: ['batch_effect_edger', 'batch edger', 'batch_edger', 'batchedger', 'batch correction', 'with batch']
    }
  },

  // Stage 3: Quantification + QC (Outlier Action)
  stage3_outlier_action: {
    canonical: ['KEEP_ALL', 'REMOVE_OUTLIERS'],
    synonyms: {
      KEEP_ALL: ['keep_all', 'keep all', 'no removal', 'retain all', 'all samples'],
      REMOVE_OUTLIERS: ['remove_outliers', 'remove outliers', 'filter outliers', 'exclude outliers', 'drop outliers']
    }
  },

  // Stage 4: DE Analysis
  stage4: {
    canonical: ['APPROVE', 'REQUEST_REANALYSIS'],
    synonyms: {
      APPROVE: ['approve', 'approved', 'accept', 'looks good', 'results acceptable', 'pass'],
      REQUEST_REANALYSIS: ['request_reanalysis', 'request reanalysis', 'reanalyze', 're-analyze', 'redo', 'reject']
    }
  }
};

/**
 * Normalize a decision string to canonical form
 *
 * @param {string} decision - Raw decision string from agent
 * @param {string} stage - Stage identifier (stage1, stage2, stage3_de_method, stage3_outlier_action, stage4)
 * @returns {string|null} - Canonical decision or null if unrecognized
 */
export function normalizeDecision(decision, stage) {
  if (!decision || typeof decision !== 'string') {
    return null;
  }

  const lower = decision.toLowerCase().trim();
  const vocab = DECISION_VOCABULARY[stage];

  if (!vocab) {
    console.warn(`Unknown stage: ${stage}`);
    return null;
  }

  // Check if already canonical (case-insensitive)
  for (const canonical of vocab.canonical) {
    if (canonical.toLowerCase() === lower) {
      return canonical;
    }
  }

  // Check synonyms
  for (const [canonical, synonymList] of Object.entries(vocab.synonyms)) {
    for (const synonym of synonymList) {
      if (lower.includes(synonym)) {
        return canonical;
      }
    }
  }

  // Not recognized
  console.warn(`Unrecognized decision for ${stage}: "${decision}"`);
  return null;
}

/**
 * Validate decision against canonical vocabulary
 *
 * @param {string} decision - Decision to validate
 * @param {string} stage - Stage identifier
 * @returns {boolean} - True if valid canonical decision
 */
export function isValidDecision(decision, stage) {
  const vocab = DECISION_VOCABULARY[stage];
  if (!vocab) return false;
  return vocab.canonical.includes(decision);
}

/**
 * Get all canonical decisions for a stage
 *
 * @param {string} stage - Stage identifier
 * @returns {Array<string>} - List of valid canonical decisions
 */
export function getCanonicalDecisions(stage) {
  const vocab = DECISION_VOCABULARY[stage];
  return vocab ? vocab.canonical : [];
}

/**
 * Extract decision from agent response text
 * Uses ONLY structured format extraction (Decision: X, DE_Method: X, etc.)
 * NO fallback text search - if agents don't follow format, we fail loudly
 *
 * @param {string} responseText - Agent response content
 * @param {string} stage - Stage identifier
 * @returns {string|null} - Canonical decision or null (if no structured format found)
 */
export function extractDecisionFromResponse(responseText, stage) {
  if (!responseText || typeof responseText !== 'string') {
    console.warn(`[Decision Extraction] Empty or invalid response for ${stage}`);
    return null;
  }

  const vocab = DECISION_VOCABULARY[stage];
  if (!vocab) {
    console.warn(`[Decision Extraction] Unknown stage: ${stage}`);
    return null;
  }

  // ONLY extract from structured format - no text search fallback!
  const patterns = [
    /decision:\s*([a-z_]+)/i,      // "Decision: PASS"
    /final_decision:\s*([a-z_]+)/i, // "Final_Decision: APPROVE"
    /de_method:\s*([a-z_]+)/i,     // "DE_Method: simpleEdger"
    /outlier_action:\s*([a-z_]+)/i // "Outlier_Action: KEEP_ALL"
  ];

  for (const pattern of patterns) {
    const match = responseText.match(pattern);
    if (match) {
      const extracted = match[1].trim();
      const normalized = normalizeDecision(extracted, stage);

      if (normalized) {
        return normalized;
      } else {
        // Agent used structured format but invalid decision
        console.error(`[Decision Extraction] Agent used structured format but invalid decision for ${stage}: "${extracted}"`);
        console.error(`[Decision Extraction] Valid options: ${vocab.canonical.join(', ')}`);
        return null;
      }
    }
  }

  // NO FALLBACK! If agents don't follow format, we fail loudly
  console.error(`[Decision Extraction] No structured decision found in agent response for ${stage}`);
  console.error(`[Decision Extraction] Expected format: "Decision: [${vocab.canonical.join(' | ')}]"`);
  console.error(`[Decision Extraction] Agent response preview: ${responseText.substring(0, 200)}...`);

  return null;
}
