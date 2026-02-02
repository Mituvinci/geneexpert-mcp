/**
 * Canonical Decision Vocabulary
 *
 * Single source of truth for all valid decisions at each stage.
 * All parsing/matching/logging MUST use these exact strings.
 *
 * ===== scRNA-seq Decision Ranges (aligned with ground truth) =====
 * Stage 1 (Load + QC): auto_proceed (no agent checkpoint)
 * Stage 2 (QC Filtering): SET_THRESHOLDS | USE_DEFAULT_THRESHOLDS | INSUFFICIENT_DATA
 * Stage 3 (Normalization + HVG): auto_proceed (no agent checkpoint)
 * Stage 3A (Cell Cycle Scoring): REMOVE_CELL_CYCLE | SKIP_CELL_CYCLE | UNCERTAIN
 * Stage 3B (Execute Cell Cycle): Skip Regression | REMOVED_CELL_CYCLE (auto-proceed based on 3A)
 * Stage 4 (PCA + PC Selection): USE_DEFAULT | SELECT_PC_RANGE | STOP_AND_REVIEW
 * Stage 5 (Clustering + Markers): ACCEPT_CLUSTERING | ADJUST_RESOLUTION | FLAG_SUSPICIOUS
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
  },

  // ===== scRNA-seq Stages =====

  // scRNA Stage 1: Load + QC (always auto-proceed, no agent checkpoint)
  scrna_stage1: {
    canonical: ['auto_proceed'],
    synonyms: {
      auto_proceed: ['auto_proceed', 'auto proceed', 'automatic', 'no review', 'skip checkpoint']
    }
  },

  // scRNA Stage 2: QC Filtering
  scrna_stage2: {
    canonical: ['PROCEED', 'PROCEED_WITH_WARNING', 'STOP_AND_REVIEW'],
    synonyms: {
      PROCEED: ['proceed', 'pass', 'approve', 'accept', 'looks good', 'continue'],
      PROCEED_WITH_WARNING: ['proceed_with_warning', 'proceed with warning', 'warning', 'borderline', 'caution'],
      STOP_AND_REVIEW: ['stop_and_review', 'stop and review', 'stop', 'reject', 'abort', 'flag', 'concerning']
    }
  },

  // scRNA Stage 2: QC Threshold Recommendation
  scrna_stage2_thresholds: {
    canonical: ['SET_THRESHOLDS', 'USE_DEFAULT_THRESHOLDS', 'INSUFFICIENT_DATA'],
    synonyms: {
      SET_THRESHOLDS: ['set_thresholds', 'set thresholds', 'custom thresholds', 'recommend thresholds', 'adaptive'],
      USE_DEFAULT_THRESHOLDS: ['use_default', 'default thresholds', 'standard', 'typical'],
      INSUFFICIENT_DATA: ['insufficient', 'too few cells', 'abort', 'stop', 'reject']
    }
  },

  // scRNA Stage 3: Normalization + HVG (always auto-proceed, no agent checkpoint - deterministic)
  scrna_stage3: {
    canonical: ['auto_proceed'],
    synonyms: {
      auto_proceed: ['auto_proceed', 'auto proceed', 'automatic', 'no review', 'skip checkpoint', 'deterministic']
    }
  },

  // scRNA Stage 3A: Cell Cycle Scoring (agent checkpoint)
  scrna_stage3A: {
    canonical: ['REMOVE_CELL_CYCLE', 'SKIP_CELL_CYCLE', 'UNCERTAIN'],
    synonyms: {
      REMOVE_CELL_CYCLE: ['remove_cell_cycle', 'remove cell cycle', 'regress', 'regress cell cycle', 'apply regression', 'remove'],
      SKIP_CELL_CYCLE: ['skip_cell_cycle', 'skip cell cycle', 'skip regression', 'no regression', 'skip', 'no removal'],
      UNCERTAIN: ['uncertain', 'unclear', 'borderline', 'unsure', 'escalate', 'user decision']
    }
  },

  // scRNA Stage 3B: Execute Cell Cycle Decision (auto-proceed based on Stage 3A)
  scrna_stage3B: {
    canonical: ['Skip Regression', 'REMOVED_CELL_CYCLE'],
    synonyms: {
      'Skip Regression': ['skip regression', 'skip_regression', 'skipped', 'no regression', 'regression skipped'],
      REMOVED_CELL_CYCLE: ['removed_cell_cycle', 'removed cell cycle', 'regression applied', 'cell cycle removed', 'regressed']
    }
  },

  // scRNA Stage 4: PCA + PC Selection
  scrna_stage4: {
    canonical: ['USE_DEFAULT', 'SELECT_PC_RANGE', 'STOP_AND_REVIEW'],
    synonyms: {
      USE_DEFAULT: ['use_default', 'use default', 'default', 'standard', '1-20', 'pc 1-20'],
      SELECT_PC_RANGE: ['select_pc_range', 'select pc range', 'custom range', 'specify pcs'],
      STOP_AND_REVIEW: ['stop_and_review', 'stop and review', 'stop', 'reject', 'abort', 'concerning']
    }
  },

  // scRNA Stage 5: Clustering + Markers
  scrna_stage5: {
    canonical: ['ACCEPT_CLUSTERING', 'ADJUST_RESOLUTION', 'FLAG_SUSPICIOUS'],
    synonyms: {
      ACCEPT_CLUSTERING: ['accept_clustering', 'accept clustering', 'accept', 'approve', 'looks good', 'proceed'],
      ADJUST_RESOLUTION: ['adjust_resolution', 'adjust resolution', 'change resolution', 'tune resolution', 'modify'],
      FLAG_SUSPICIOUS: ['flag_suspicious', 'flag suspicious', 'suspicious', 'concerning', 'reject', 'problematic']
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
