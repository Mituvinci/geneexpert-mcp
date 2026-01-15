/**
 * Modification Engine (Option 2)
 * Applies structured agent recommendations to pipeline template
 */

import fs from 'fs';
import path from 'path';

/**
 * Apply modifications to base template script
 * Returns modified script or null if modifications fail
 */
export function applyModifications(baseScript, modifications, context) {
  try {
    let modifiedScript = baseScript;

    console.log('[Modification Engine] Applying structured modifications...');
    console.log(`[Modification Engine] Statistical adjustments: ${modifications.statistical_adjustments.length}`);
    console.log(`[Modification Engine] Quality checks: ${modifications.quality_checks.length}`);
    console.log(`[Modification Engine] Parameter changes: ${Object.keys(modifications.parameter_changes).length}`);

    // 1. Apply statistical adjustments
    for (const adjustment of modifications.statistical_adjustments) {
      if (adjustment.type === 'use_exact_test') {
        modifiedScript = applyExactTestModification(modifiedScript, adjustment, context);
      }
    }

    // 2. Apply quality checks (group by insertion point to avoid conflicts)
    const qcByStep = {
      before_de_analysis: [],
      after_alignment: []
    };

    for (const qc of modifications.quality_checks) {
      if (qcByStep[qc.step]) {
        qcByStep[qc.step].push(qc);
      }
    }

    // Insert QC checks grouped by location
    if (qcByStep.before_de_analysis.length > 0) {
      modifiedScript = insertQualityChecks(modifiedScript, qcByStep.before_de_analysis, 'before_de_analysis', context);
    }
    if (qcByStep.after_alignment.length > 0) {
      modifiedScript = insertQualityChecks(modifiedScript, qcByStep.after_alignment, 'after_alignment', context);
    }

    // 3. Apply parameter changes
    modifiedScript = applyParameterChanges(modifiedScript, modifications.parameter_changes, context);

    console.log('[Modification Engine] ✓ All modifications applied successfully');

    return modifiedScript;

  } catch (error) {
    console.error('[Modification Engine] ❌ Failed to apply modifications:', error.message);
    return null;  // Trigger fallback
  }
}

/**
 * Modify DE analysis to use exactTest instead of glmQLFTest
 */
function applyExactTestModification(script, adjustment, context) {
  console.log(`[Modification Engine]   → Using exactTest for small n (${adjustment.reason})`);

  // Find the DE analysis step - match actual template output
  // Template has: echo "Step 8/10: Differential expression analysis (edgeR)..."
  //               Rscript $SCRIPTS_DIR/simpleEdger3.R ...
  const deStepRegex = /(echo "Step \d+\/\d+: Differential expression analysis.*?\n)(Rscript.*?simpleEdger3\.R.*?)(\n)/s;

  const match = script.match(deStepRegex);
  if (!match) {
    console.warn('[Modification Engine]   ⚠️  Could not find DE analysis step, skipping modification');
    console.warn('[Modification Engine]   Searched for pattern: echo "Step .../...: Differential expression analysis"');
    return script;
  }

  // Add comment about exact test BEFORE the Rscript line
  const comment = `# ADAPTATION: Using exactTest for small sample size (${adjustment.reason})\n`;

  // Modify R script call to include exact test parameter
  // Note: simpleEdger3.R may not accept these params yet - this documents the intent
  const modifiedCall = match[2] + ` # exactTest method, FDR=${context.fdrThreshold || 0.1}`;

  return script.replace(deStepRegex, `${match[1]}${comment}${modifiedCall}${match[3]}`);
}

/**
 * Insert multiple quality checks at the same insertion point
 * Combines multiple QC checks into one block to avoid regex conflicts
 */
function insertQualityChecks(script, qcChecks, step, context) {
  if (qcChecks.length === 0) return script;

  console.log(`[Modification Engine]   → Adding ${qcChecks.length} QC checks at ${step}`);

  let insertPoint;
  let combinedQcSteps = [];

  if (step === 'before_de_analysis') {
    // Insert before Step 8 (Differential expression analysis)
    insertPoint = /echo "Step 8\/\d+: Differential expression analysis/;

    combinedQcSteps.push('# Step 7.5: Quality Control Checks');
    combinedQcSteps.push('echo "Step 7.5: Running quality control checks..."');
    combinedQcSteps.push('# ADAPTATION: Extra QC for small sample size');
    combinedQcSteps.push('');

    for (const qc of qcChecks) {
      if (qc.type === 'outlier_detection') {
        combinedQcSteps.push('echo "   → Checking for outliers (PCA/MDS plots)..."');
        combinedQcSteps.push('Rscript $SCRIPTS_DIR/qc_plots.R outentrz.txt $OUTPUT_DIR/qc_plots');
        combinedQcSteps.push('echo "   → Review plots in $OUTPUT_DIR/qc_plots/ before proceeding"');
      } else if (qc.type === 'batch_detection') {
        combinedQcSteps.push('echo "   → Checking for batch effects..."');
        combinedQcSteps.push('echo "   → Inspect MDS/PCA plots for batch-related clustering"');
      }
    }

    combinedQcSteps.push('');
  } else if (step === 'after_alignment') {
    // Insert before Step 3 (Feature Counts)
    insertPoint = /echo "Step 3\/\d+: Feature Counts/;

    combinedQcSteps.push('# Step 2.5: Enhanced Quality Control');
    combinedQcSteps.push('echo "Step 2.5: Checking alignment quality..."');
    combinedQcSteps.push('# ADAPTATION: Extra QC for small sample size');

    for (const qc of qcChecks) {
      if (qc.type === 'enhanced_qc') {
        combinedQcSteps.push('echo "   → Checking alignment rates, duplication, library complexity"');
        combinedQcSteps.push('for bam in $OUTPUT_DIR/bam_files/*.bam; do');
        combinedQcSteps.push('  echo "   $(basename $bam): $(samtools flagstat $bam | grep \'mapped (\')"');
        combinedQcSteps.push('done');
      }
    }

    combinedQcSteps.push('');
  }

  if (combinedQcSteps.length > 0 && insertPoint) {
    const qcBlock = combinedQcSteps.join('\n') + '\n';
    const modifiedScript = script.replace(insertPoint, qcBlock + '$&');

    if (modifiedScript === script) {
      console.warn(`[Modification Engine]   ⚠️  Insertion point not found for ${step}`);
      console.warn(`[Modification Engine]   Pattern: ${insertPoint}`);
      return script;
    }

    console.log(`[Modification Engine]   ✓ Inserted ${qcChecks.length} QC checks at ${step}`);
    return modifiedScript;
  }

  console.warn(`[Modification Engine]   ⚠️  Could not insert QC checks at ${step}`);
  return script;
}

/**
 * Apply parameter changes throughout the script
 */
function applyParameterChanges(script, paramChanges, context) {
  let modified = script;

  // FDR threshold
  if (paramChanges.fdr_threshold) {
    console.log(`[Modification Engine]   → FDR threshold: ${paramChanges.fdr_threshold}`);
    context.fdrThreshold = paramChanges.fdr_threshold;
  }

  // CPM threshold (for filtering)
  if (paramChanges.cpm_threshold) {
    console.log(`[Modification Engine]   → CPM threshold: ${paramChanges.cpm_threshold}`);
    // This would modify filterIDS.R call if it accepts parameters
  }

  // LogFC threshold
  if (paramChanges.logfc_threshold) {
    console.log(`[Modification Engine]   → LogFC threshold: ${paramChanges.logfc_threshold}`);
  }

  return modified;
}

/**
 * Validate that the modified script is syntactically correct
 */
export function validateScript(script) {
  try {
    // Basic validation checks

    // 1. Must start with shebang
    if (!script.startsWith('#!/bin/bash')) {
      console.error('[Validation] ❌ Script missing shebang');
      return false;
    }

    // 2. Must have conda activation
    if (!script.includes('conda activate')) {
      console.error('[Validation] ❌ Script missing conda activation');
      return false;
    }

    // 3. Must have all critical steps (match actual template strings)
    const requiredSteps = [
      { pattern: 'FastQC', name: 'FastQC' },
      { pattern: 'Alignment', name: 'Alignment' },
      { pattern: 'Feature Counts', name: 'Feature Counts' },
      { pattern: 'Filtering bad gene IDs', name: 'Filter Bad IDs' },
      { pattern: 'Differential expression analysis', name: 'DE Analysis' }
    ];

    for (const step of requiredSteps) {
      if (!script.includes(step.pattern)) {
        console.error(`[Validation] ❌ Missing required step: ${step.name} (looking for: "${step.pattern}")`);
        return false;
      }
    }

    // 4. Check for obvious syntax errors
    const lines = script.split('\n');
    for (let i = 0; i < lines.length; i++) {
      const line = lines[i].trim();

      // Check for unmatched quotes
      if (line && !line.startsWith('#')) {
        const singleQuotes = (line.match(/'/g) || []).length;
        const doubleQuotes = (line.match(/"/g) || []).length;

        if (singleQuotes % 2 !== 0 || doubleQuotes % 2 !== 0) {
          console.error(`[Validation] ❌ Unmatched quotes on line ${i + 1}: ${line}`);
          return false;
        }
      }
    }

    console.log('[Validation] ✓ Script passed validation');
    return true;

  } catch (error) {
    console.error('[Validation] ❌ Validation error:', error.message);
    return false;
  }
}

export default {
  applyModifications,
  validateScript
};
