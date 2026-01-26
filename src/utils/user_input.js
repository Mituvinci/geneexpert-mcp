/**
 * User Input Utility
 *
 * Provides interactive prompts for user decisions when agents can't reach consensus
 * or when USER_DECISION_REQUIRED is returned.
 */

import readline from 'readline';

/**
 * Create a readline interface for user input
 */
function createReadlineInterface() {
  return readline.createInterface({
    input: process.stdin,
    output: process.stdout
  });
}

/**
 * Ask user a yes/no question
 * @param {string} question - The question to ask
 * @param {boolean} defaultAnswer - Default answer if user just presses enter
 * @returns {Promise<boolean>} - User's answer
 */
export async function askYesNo(question, defaultAnswer = true) {
  const rl = createReadlineInterface();
  const defaultStr = defaultAnswer ? 'Y/n' : 'y/N';

  return new Promise((resolve) => {
    const askQuestion = () => {
      rl.question(`${question} [${defaultStr}]: `, (answer) => {
        // Empty input → use default
        if (!answer || answer.trim() === '') {
          rl.close();
          resolve(defaultAnswer);
          return;
        }

        const normalized = answer.trim().toLowerCase();

        // Validate: must be y/yes/n/no
        if (normalized === 'y' || normalized === 'yes') {
          rl.close();
          resolve(true);
          return;
        } else if (normalized === 'n' || normalized === 'no') {
          rl.close();
          resolve(false);
          return;
        } else {
          console.log('');
          console.log('❌ Invalid input! Please enter "y" (yes) or "n" (no).');
          console.log('');
          askQuestion();  // Loop: ask again
          return;
        }
      });
    };

    askQuestion();  // Start the loop
  });
}

/**
 * Ask user to select from a list of options
 * @param {string} question - The question to ask
 * @param {string[]} options - Array of options
 * @param {number} defaultIndex - Default option index (0-based)
 * @returns {Promise<{index: number, value: string}>} - Selected option
 */
export async function askChoice(question, options, defaultIndex = 0) {
  const rl = createReadlineInterface();

  console.log('');
  console.log(question);
  console.log('-'.repeat(40));
  options.forEach((opt, i) => {
    const marker = i === defaultIndex ? '*' : ' ';
    console.log(`  ${marker} ${i + 1}. ${opt}`);
  });
  console.log('');

  return new Promise((resolve) => {
    const askQuestion = () => {
      rl.question(`Enter choice [1-${options.length}] (default: ${defaultIndex + 1}): `, (answer) => {
        // Empty input → use default
        if (!answer || answer.trim() === '') {
          rl.close();
          resolve({ index: defaultIndex, value: options[defaultIndex] });
          return;
        }

        const choice = parseInt(answer.trim(), 10);

        // Validate: must be number within range
        if (isNaN(choice) || choice < 1 || choice > options.length) {
          console.log('');
          console.log(`❌ Invalid input! Please enter a number between 1 and ${options.length}.`);
          console.log('');
          askQuestion();  // Loop: ask again
          return;
        }

        // Valid input
        rl.close();
        resolve({ index: choice - 1, value: options[choice - 1] });
      });
    };

    askQuestion();  // Start the loop
  });
}

/**
 * Ask user for free text input
 * @param {string} question - The question to ask
 * @param {string} defaultAnswer - Default answer
 * @returns {Promise<string>} - User's answer
 */
export async function askText(question, defaultAnswer = '') {
  const rl = createReadlineInterface();
  const defaultStr = defaultAnswer ? ` [${defaultAnswer}]` : '';

  return new Promise((resolve) => {
    rl.question(`${question}${defaultStr}: `, (answer) => {
      rl.close();
      resolve(answer.trim() || defaultAnswer);
    });
  });
}

/**
 * Ask user to select multiple items from a list
 * @param {string} question - The question to ask
 * @param {string[]} items - Array of items to choose from
 * @returns {Promise<string[]>} - Selected items
 */
export async function askMultiSelect(question, items) {
  const rl = createReadlineInterface();

  console.log('');
  console.log(question);
  console.log('-'.repeat(40));
  items.forEach((item, i) => {
    console.log(`  ${i + 1}. ${item}`);
  });
  console.log('');

  return new Promise((resolve) => {
    rl.question(`Enter choices (comma-separated, e.g., "1,3,5" or "none"): `, (answer) => {
      rl.close();

      if (!answer || answer.trim().toLowerCase() === 'none') {
        resolve([]);
        return;
      }

      const selections = answer.split(',')
        .map(s => parseInt(s.trim(), 10))
        .filter(n => !isNaN(n) && n >= 1 && n <= items.length)
        .map(n => items[n - 1]);

      resolve(selections);
    });
  });
}

/**
 * Display agent responses and ask for user decision on Stage 1
 * @param {Object} reviewResult - Result from coordinator.reviewStage1Output()
 * @param {Object} stage1Output - Parsed stage 1 output
 * @param {Object} config - Configuration object with singleAgent info
 * @returns {Promise<{proceed: boolean, reason: string}>}
 */
export async function handleStage1UserDecision(reviewResult, stage1Output, config = {}) {
  console.log('');
  console.log('='.repeat(60));
  console.log('  USER DECISION REQUIRED');
  console.log('='.repeat(60));
  console.log('');
  console.log('The agents could not reach a clear consensus on Stage 1 results.');
  console.log('');

  // Helper to get display name for agent
  const getAgentDisplayName = (agentKey, modelStr) => {
    if (config.singleAgent) {
      // Single-agent mode: show model name with role
      const modelName = getModelDisplayName(modelStr);
      const roleMap = { gpt5_2: 'Statistics', claude: 'Pipeline', gemini: 'Biology' };
      return `${modelName} (Role: ${roleMap[agentKey]})`;
    } else {
      // Multi-agent mode: use traditional names
      const defaultNames = { gpt5_2: 'GPT-5.2 (Stats)', claude: 'Claude (Pipeline)', gemini: 'Gemini (Biology)' };
      return defaultNames[agentKey];
    }
  };

  // Show what each agent said
  console.log('--- Agent Recommendations ---');
  if (reviewResult.responses.gpt5_2.success) {
    const gptResponse = reviewResult.responses.gpt5_2.content;
    console.log(`${getAgentDisplayName('gpt5_2', reviewResult.responses.gpt5_2.model)}: ${extractDecisionSummary(gptResponse)}`);
  }
  if (reviewResult.responses.claude.success) {
    const claudeResponse = reviewResult.responses.claude.content;
    console.log(`${getAgentDisplayName('claude', reviewResult.responses.claude.model)}: ${extractDecisionSummary(claudeResponse)}`);
  }
  if (reviewResult.responses.gemini.success) {
    const geminiResponse = reviewResult.responses.gemini.content;
    console.log(`${getAgentDisplayName('gemini', reviewResult.responses.gemini.model)}: ${extractDecisionSummary(geminiResponse)}`);
  }
  console.log('');

  // Show the issues detected
  if (stage1Output.warnings && stage1Output.warnings.length > 0) {
    console.log('--- Warnings Detected ---');
    stage1Output.warnings.forEach(w => console.log(`  - ${w}`));
    console.log('');
  }

  if (stage1Output.errors && stage1Output.errors.length > 0) {
    console.log('--- Errors Detected ---');
    stage1Output.errors.forEach(e => console.log(`  - ${e}`));
    console.log('');
  }

  // Ask user
  const options = [
    'PROCEED - Continue to Stage 2 (Alignment)',
    'ABORT - Stop analysis, investigate issues first',
    'VIEW DETAILS - Show full agent responses then decide'
  ];

  const choice = await askChoice('What would you like to do?', options, 0);

  if (choice.index === 2) {
    // Show full responses
    console.log('');
    console.log('=== Full Agent Responses ===');
    console.log('');
    console.log('--- GPT-5.2 (Statistics Agent) ---');
    console.log(reviewResult.responses.gpt5_2.content || 'No response');
    console.log('');
    console.log('--- Claude (Pipeline Agent) ---');
    console.log(reviewResult.responses.claude.content || 'No response');
    console.log('');
    console.log('--- Gemini (Biology Agent) ---');
    console.log(reviewResult.responses.gemini.content || 'No response');
    console.log('');

    // Ask again
    const finalChoice = await askChoice('After reviewing, what would you like to do?', options.slice(0, 2), 0);
    return {
      proceed: finalChoice.index === 0,
      reason: `User decision after viewing details: ${finalChoice.value}`
    };
  }

  return {
    proceed: choice.index === 0,
    reason: `User decision: ${choice.value}`
  };
}

/**
 * Display agent responses and ask for user decision on Stage 2
 * @param {Object} reviewResult - Result from coordinator.reviewStage2Output()
 * @param {Object} stage2Output - Parsed stage 2 output
 * @param {Object} config - Configuration object with singleAgent info
 * @returns {Promise<{proceed: boolean, samplesToRemove: string[], reason: string}>}
 */
export async function handleStage2UserDecision(reviewResult, stage2Output, config = {}) {
  console.log('');
  console.log('='.repeat(60));
  console.log('  USER DECISION REQUIRED');
  console.log('='.repeat(60));
  console.log('');
  console.log('The agents could not reach a clear consensus on Stage 2 (Alignment QC) results.');
  console.log('');

  // Show alignment QC summary
  console.log('--- Alignment QC Summary ---');
  console.log(`  Mean mapping rate: ${stage2Output.mean_mapping_rate ? (stage2Output.mean_mapping_rate * 100).toFixed(1) + '%' : 'N/A'}`);
  console.log(`  PASS samples: ${stage2Output.passed_samples?.length || 0}`);
  console.log(`  WARN samples: ${stage2Output.warned_samples?.length || 0}`);
  console.log(`  FAIL samples: ${stage2Output.failed_samples?.length || 0}`);
  console.log('');

  // Show problematic samples - FIX: Extract sample names properly
  const problematicSamples = [...(stage2Output.warned_samples || []), ...(stage2Output.failed_samples || [])];
  if (problematicSamples.length > 0) {
    console.log('--- Problematic Samples ---');
    problematicSamples.forEach(sampleItem => {
      // Handle both string and object formats
      const sampleName = typeof sampleItem === 'string' ? sampleItem : (sampleItem.sample || sampleItem.name || String(sampleItem));
      const rate = stage2Output.sample_details?.[sampleName]?.mapping_rate;
      const status = stage2Output.sample_details?.[sampleName]?.status;
      console.log(`  ${sampleName}: ${rate ? (rate * 100).toFixed(1) + '%' : 'N/A'} (${status || 'unknown'})`);
    });
    console.log('');
  }

  // Helper to get display name for agent
  const getAgentDisplayName = (agentKey, modelStr) => {
    if (config.singleAgent) {
      const modelName = getModelDisplayName(modelStr);
      const roleMap = { gpt5_2: 'Statistics', claude: 'Pipeline', gemini: 'Biology' };
      return `${modelName} (Role: ${roleMap[agentKey]})`;
    } else {
      const defaultNames = { gpt5_2: 'GPT-5.2 (Stats)', claude: 'Claude (Pipeline)', gemini: 'Gemini (Biology)' };
      return defaultNames[agentKey];
    }
  };

  // Show agent recommendations
  console.log('--- Agent Recommendations ---');
  if (reviewResult.responses.gpt5_2.success) {
    console.log(`${getAgentDisplayName('gpt5_2', reviewResult.responses.gpt5_2.model)}: ${extractDecisionSummary(reviewResult.responses.gpt5_2.content)}`);
  }
  if (reviewResult.responses.claude.success) {
    console.log(`${getAgentDisplayName('claude', reviewResult.responses.claude.model)}: ${extractDecisionSummary(reviewResult.responses.claude.content)}`);
  }
  if (reviewResult.responses.gemini.success) {
    console.log(`${getAgentDisplayName('gemini', reviewResult.responses.gemini.model)}: ${extractDecisionSummary(reviewResult.responses.gemini.content)}`);
  }
  console.log('');

  // Ask user
  const options = [
    'PROCEED - Continue with all samples',
    'PROCEED_REMOVE - Continue but remove problematic samples',
    'ABORT - Stop analysis, investigate issues first',
    'VIEW DETAILS - Show full agent responses then decide'
  ];

  const choice = await askChoice('What would you like to do?', options, 0);

  if (choice.index === 3) {
    // Show full responses
    console.log('');
    console.log('=== Full Agent Responses ===');
    console.log('');
    console.log('--- GPT-5.2 (Statistics Agent) ---');
    console.log(reviewResult.responses.gpt5_2.content || 'No response');
    console.log('');
    console.log('--- Claude (Pipeline Agent) ---');
    console.log(reviewResult.responses.claude.content || 'No response');
    console.log('');
    console.log('--- Gemini (Biology Agent) ---');
    console.log(reviewResult.responses.gemini.content || 'No response');
    console.log('');

    // Ask again without VIEW DETAILS option
    const finalChoice = await askChoice('After reviewing, what would you like to do?', options.slice(0, 3), 0);

    let samplesToRemove = [];
    if (finalChoice.index === 1 && problematicSamples.length > 0) {
      samplesToRemove = await askMultiSelect(
        'Select samples to remove:',
        problematicSamples
      );
      if (samplesToRemove.length === 0) {
        samplesToRemove = problematicSamples; // Default to removing all problematic
      }
    }

    return {
      proceed: finalChoice.index < 2,
      samplesToRemove,
      reason: `User decision after viewing details: ${finalChoice.value}`
    };
  }

  let samplesToRemove = [];
  if (choice.index === 1 && problematicSamples.length > 0) {
    samplesToRemove = await askMultiSelect(
      'Select samples to remove (or "none" to remove all flagged):',
      problematicSamples
    );
    if (samplesToRemove.length === 0) {
      samplesToRemove = problematicSamples; // Default to removing all problematic
    }
  }

  return {
    proceed: choice.index < 2,
    samplesToRemove,
    reason: `User decision: ${choice.value}`
  };
}

/**
 * Display agent responses and ask for user decision on Stage 3 (PCA/QC)
 * @param {Object} reviewResult - Result from coordinator.reviewStage3Output()
 * @param {Object} stage3Output - Parsed stage 3 output
 * @param {Object} config - Configuration object with singleAgent info
 * @returns {Promise<{proceed: boolean, samplesToRemove: string[], useBatchCorrection: boolean, reason: string}>}
 */
export async function handleStage3UserDecision(reviewResult, stage3Output, config = {}) {
  console.log('');
  console.log('='.repeat(60));
  console.log('  USER DECISION REQUIRED');
  console.log('='.repeat(60));
  console.log('');
  console.log('The agents could not reach a clear consensus on Stage 3 (PCA/QC) results.');
  console.log('');

  // Show QC summary
  console.log('--- PCA/QC Summary ---');
  if (stage3Output.outliers && stage3Output.outliers.length > 0) {
    console.log(`  Outliers detected: ${stage3Output.outliers.join(', ')}`);
  } else {
    console.log('  No outliers detected');
  }
  if (stage3Output.batch_effect_detected) {
    console.log('  Batch effects: DETECTED');
  } else {
    console.log('  Batch effects: Not detected');
  }
  console.log(`  QC Exit Code: ${stage3Output.exit_code || 'N/A'}`);
  console.log('');

  // Helper to get display name for agent
  const getAgentDisplayName = (agentKey, modelStr) => {
    if (config.singleAgent) {
      const modelName = getModelDisplayName(modelStr);
      const roleMap = { gpt5_2: 'Statistics', claude: 'Pipeline', gemini: 'Biology' };
      return `${modelName} (Role: ${roleMap[agentKey]})`;
    } else {
      const defaultNames = { gpt5_2: 'GPT-5.2 (Stats)', claude: 'Claude (Pipeline)', gemini: 'Gemini (Biology)' };
      return defaultNames[agentKey];
    }
  };

  // Show agent recommendations
  console.log('--- Agent Recommendations ---');
  if (reviewResult.responses.gpt5_2.success) {
    console.log(`${getAgentDisplayName('gpt5_2', reviewResult.responses.gpt5_2.model)}: ${extractDecisionSummary(reviewResult.responses.gpt5_2.content)}`);
  }
  if (reviewResult.responses.claude.success) {
    console.log(`${getAgentDisplayName('claude', reviewResult.responses.claude.model)}: ${extractDecisionSummary(reviewResult.responses.claude.content)}`);
  }
  if (reviewResult.responses.gemini.success) {
    console.log(`${getAgentDisplayName('gemini', reviewResult.responses.gemini.model)}: ${extractDecisionSummary(reviewResult.responses.gemini.content)}`);
  }
  console.log('');

  // Ask user
  const options = [
    'PROCEED - Continue with standard DE analysis',
    'PROCEED_BATCH - Continue with batch effect correction',
    'PROCEED_REMOVE - Continue but remove outlier samples',
    'ABORT - Stop analysis, investigate issues first',
    'VIEW DETAILS - Show full agent responses then decide'
  ];

  const choice = await askChoice('What would you like to do?', options, 0);

  if (choice.index === 4) {
    // Show full responses
    console.log('');
    console.log('=== Full Agent Responses ===');
    console.log('');
    console.log('--- GPT-5.2 (Statistics Agent) ---');
    console.log(reviewResult.responses.gpt5_2.content || 'No response');
    console.log('');
    console.log('--- Claude (Pipeline Agent) ---');
    console.log(reviewResult.responses.claude.content || 'No response');
    console.log('');
    console.log('--- Gemini (Biology Agent) ---');
    console.log(reviewResult.responses.gemini.content || 'No response');
    console.log('');

    // Ask again without VIEW DETAILS option
    const finalChoice = await askChoice('After reviewing, what would you like to do?', options.slice(0, 4), 0);

    let samplesToRemove = [];
    if (finalChoice.index === 2 && stage3Output.outliers?.length > 0) {
      samplesToRemove = await askMultiSelect(
        'Select samples to remove:',
        stage3Output.outliers
      );
      if (samplesToRemove.length === 0) {
        samplesToRemove = stage3Output.outliers; // Default to removing all outliers
      }
    }

    return {
      proceed: finalChoice.index < 3,
      samplesToRemove,
      useBatchCorrection: finalChoice.index === 1,
      reason: `User decision after viewing details: ${finalChoice.value}`
    };
  }

  let samplesToRemove = [];
  if (choice.index === 2 && stage3Output.outliers?.length > 0) {
    samplesToRemove = await askMultiSelect(
      'Select samples to remove:',
      stage3Output.outliers
    );
    if (samplesToRemove.length === 0) {
      samplesToRemove = stage3Output.outliers; // Default to removing all outliers
    }
  }

  return {
    proceed: choice.index < 3,
    samplesToRemove,
    useBatchCorrection: choice.index === 1,
    reason: `User decision: ${choice.value}`
  };
}

/**
 * Get display name for model (for single-agent mode)
 * @param {string} modelStr - Model string (e.g., 'gpt-4o', 'claude-opus-4-5')
 * @returns {string} - Display name
 */
function getModelDisplayName(modelStr) {
  if (!modelStr) return 'Unknown';
  if (modelStr.includes('gpt')) return 'GPT-5.2';
  if (modelStr.includes('claude')) return 'Claude';
  if (modelStr.includes('gemini')) return 'Gemini';
  return modelStr;
}

/**
 * Extract a brief decision summary from an agent's response
 * @param {string} response - Full agent response
 * @returns {string} - Brief summary
 */
function extractDecisionSummary(response) {
  if (!response) return 'No response';

  // Look for decision keywords
  const lower = response.toLowerCase();

  // Check for explicit decision patterns
  const decisionPatterns = [
    /decision:\s*(\w+)/i,
    /verdict:\s*(\w+)/i,
    /recommendation:\s*(\w+)/i,
    /status:\s*(\w+)/i
  ];

  for (const pattern of decisionPatterns) {
    const match = response.match(pattern);
    if (match) {
      return match[1].toUpperCase();
    }
  }

  // Check for pass/fail keywords
  if (lower.includes('pass') && !lower.includes('fail')) {
    if (lower.includes('warning') || lower.includes('caution')) {
      return 'PASS_WITH_WARNING';
    }
    return 'PASS';
  }

  if (lower.includes('fail')) {
    return 'FAIL';
  }

  if (lower.includes('proceed') && lower.includes('remove')) {
    return 'PROCEED_REMOVE_SAMPLES';
  }

  if (lower.includes('proceed')) {
    return 'PROCEED';
  }

  if (lower.includes('abort') || lower.includes('stop')) {
    return 'ABORT';
  }

  // Return full response (no truncation) if no pattern found
  return response;
}

/**
 * Ask user to confirm removal of outlier samples
 * @param {Array<string>} outlierSamples - List of outlier sample names
 * @returns {Promise<boolean>} - True if user confirms removal, false otherwise
 */
export async function confirmOutlierRemoval(outlierSamples) {
  if (!outlierSamples || outlierSamples.length === 0) {
    return false;
  }

  const question = `\nRemove ${outlierSamples.length} outlier sample(s) before DE analysis?`;
  return askYesNo(question, true);  // Default: yes (recommend removal)
}

/**
 * ============================================================================
 * scRNA-seq SPECIFIC USER DECISION HANDLERS
 * ============================================================================
 */

/**
 * Handle user decision for scRNA Stage 2 (QC Thresholds)
 * Called when agents disagree on QC filtering thresholds
 */
export async function handleScRNAStage2UserDecision(reviewResult, stage1Output) {
  console.log('');
  console.log('='.repeat(60));
  console.log('  USER DECISION REQUIRED - QC Thresholds');
  console.log('='.repeat(60));
  console.log('');
  console.log('The agents disagreed on QC filtering thresholds.');
  console.log('');

  // Extract threshold recommendations from each agent
  const gptThresholds = extractThresholds(reviewResult.agentResponses?.gpt5_2?.content);
  const claudeThresholds = extractThresholds(reviewResult.agentResponses?.claude?.content);
  const geminiThresholds = extractThresholds(reviewResult.agentResponses?.gemini?.content);

  console.log('--- Agent Recommendations ---');
  console.log('');
  console.log('GPT-5.2 (Stats Agent):');
  console.log(`  nFeature: ${gptThresholds.nFeature_min} - ${gptThresholds.nFeature_max}`);
  console.log(`  %MT max: ${gptThresholds.percent_mt_max}%`);
  console.log('');
  console.log('Claude (Pipeline Agent):');
  console.log(`  nFeature: ${claudeThresholds.nFeature_min} - ${claudeThresholds.nFeature_max}`);
  console.log(`  %MT max: ${claudeThresholds.percent_mt_max}%`);
  console.log('');
  console.log('Gemini (Biology Agent):');
  console.log(`  nFeature: ${geminiThresholds.nFeature_min} - ${geminiThresholds.nFeature_max}`);
  console.log(`  %MT max: ${geminiThresholds.percent_mt_max}%`);
  console.log('');

  console.log('--- Dataset Summary ---');
  console.log(`Cells: ${stage1Output.cells_total || stage1Output.n_cells || 'unknown'}`);
  console.log(`Median nFeature: ${stage1Output.nFeature_median || stage1Output.median_nFeature || 'unknown'}`);
  console.log(`Median %MT: ${(stage1Output.percent_mt_median !== undefined ? stage1Output.percent_mt_median.toFixed(2) : stage1Output.median_percent_mt) || 'unknown'}%`);
  console.log('');

  const options = [
    'Use GPT-5.2 thresholds (Statistics-focused)',
    'Use Claude thresholds (Pipeline-focused)',
    'Use Gemini thresholds (Biology-focused)',
    'Use default conservative thresholds (200-6000, <10% MT)',
    'ABORT - Review data manually'
  ];

  const choice = await askChoice('Which thresholds would you like to use?', options, 0);

  if (choice.index === 4) {
    return { proceed: false, reason: 'User aborted to review data manually' };
  }

  let selectedThresholds;
  if (choice.index === 0) selectedThresholds = gptThresholds;
  else if (choice.index === 1) selectedThresholds = claudeThresholds;
  else if (choice.index === 2) selectedThresholds = geminiThresholds;
  else selectedThresholds = { nFeature_min: 200, nFeature_max: 6000, percent_mt_max: 10 };

  return {
    proceed: true,
    thresholds: selectedThresholds,
    reason: `User selected: ${choice.value}`
  };
}

/**
 * Handle user decision for scRNA Stage 4 (PC Selection)
 * Called when agents disagree on PC range
 */
export async function handleScRNAStage4UserDecision(reviewResult, stage4Output) {
  console.log('');
  console.log('='.repeat(60));
  console.log('  USER DECISION REQUIRED - PC Range Selection');
  console.log('='.repeat(60));
  console.log('');
  console.log('The agents disagreed on the optimal PC range for clustering.');
  console.log('');

  // Extract PC recommendations
  const gptPC = extractPCRange(reviewResult.agentResponses?.gpt5_2?.content);
  const claudePC = extractPCRange(reviewResult.agentResponses?.claude?.content);
  const geminiPC = extractPCRange(reviewResult.agentResponses?.gemini?.content);

  console.log('--- Agent Recommendations ---');
  console.log('');
  console.log(`GPT-5.2 (Stats): PCs 1-${gptPC.max_pc} (${gptPC.reasoning || 'no reasoning'})`);
  console.log(`Claude (Pipeline): PCs 1-${claudePC.max_pc} (${claudePC.reasoning || 'no reasoning'})`);
  console.log(`Gemini (Biology): PCs 1-${geminiPC.max_pc} (${geminiPC.reasoning || 'no reasoning'})`);
  console.log('');

  console.log('--- PCA Summary ---');
  console.log(`PC1 variance: ${stage4Output.pc1_variance || 'unknown'}%`);
  console.log(`Cumulative variance (PC20): ${stage4Output.cumulative_variance_pc20 || 'unknown'}%`);
  console.log('');

  const options = [
    `Use GPT-5.2 recommendation (PCs 1-${gptPC.max_pc})`,
    `Use Claude recommendation (PCs 1-${claudePC.max_pc})`,
    `Use Gemini recommendation (PCs 1-${geminiPC.max_pc})`,
    'Use default (PCs 1-20)',
    'Enter custom PC range'
  ];

  const choice = await askChoice('Which PC range would you like to use?', options, 0);

  let selectedPC;
  if (choice.index === 0) selectedPC = gptPC;
  else if (choice.index === 1) selectedPC = claudePC;
  else if (choice.index === 2) selectedPC = geminiPC;
  else if (choice.index === 3) selectedPC = { min_pc: 1, max_pc: 20 };
  else {
    // Custom range
    const maxPC = await askText('Enter maximum PC (min will be 1):', '20');
    selectedPC = { min_pc: 1, max_pc: parseInt(maxPC) };
  }

  return {
    proceed: true,
    pcSelection: selectedPC,
    reason: `User selected: ${choice.value}`
  };
}

/**
 * Handle user decision for scRNA Stage 5 (Clustering)
 * Called when agents flag suspicious clustering or disagree
 */
export async function handleScRNAStage5UserDecision(reviewResult, stage5Output) {
  console.log('');
  console.log('='.repeat(60));
  console.log('  USER DECISION REQUIRED - Clustering Validation');
  console.log('='.repeat(60));
  console.log('');

  // Check if agents actually flagged concerns or just low confidence
  const hasRealConcerns = reviewResult.decision === 'FLAG_SUSPICIOUS' ||
                          reviewResult.decision === 'ADJUST_RESOLUTION' ||
                          (reviewResult.consensus && reviewResult.consensus.votes.reject >= 2);

  if (hasRealConcerns) {
    console.log('Agents flagged concerns with the clustering results.');
  } else {
    console.log('Low consensus confidence - requesting user confirmation.');
  }
  console.log('');

  // Show agent summaries
  console.log('--- Agent Assessments ---');
  console.log('');
  if (reviewResult.agentResponses?.gpt5_2?.success) {
    console.log('GPT-5.2 (Stats):');
    console.log(`  ${extractDecisionSummary(reviewResult.agentResponses.gpt5_2.content)}`);
    console.log('');
  }
  if (reviewResult.agentResponses?.claude?.success) {
    console.log('Claude (Pipeline):');
    console.log(`  ${extractDecisionSummary(reviewResult.agentResponses.claude.content)}`);
    console.log('');
  }
  if (reviewResult.agentResponses?.gemini?.success) {
    console.log('Gemini (Biology):');
    console.log(`  ${extractDecisionSummary(reviewResult.agentResponses.gemini.content)}`);
    console.log('');
  }

  console.log('--- Clustering Summary ---');
  console.log(`Number of clusters: ${stage5Output.n_clusters || 'unknown'}`);
  console.log(`Cells clustered: ${stage5Output.n_cells || 'unknown'}`);
  console.log(`Resolution: ${stage5Output.resolution || 'unknown'}`);
  console.log('');

  // Check if cell cycle is the issue
  const cellCycleConcern = checkCellCycleConcern(reviewResult);

  const options = cellCycleConcern
    ? [
        'PROCEED - Accept current clustering (cell cycle may be biologically relevant)',
        'RE-CLUSTER - Run cell cycle correction and re-cluster (recommended)',
        'VIEW DETAILS - Show full agent responses',
        'ABORT - Stop and review manually'
      ]
    : [
        'PROCEED - Accept current clustering',
        'ADJUST RESOLUTION - Try different resolution parameter',
        'VIEW DETAILS - Show full agent responses',
        'ABORT - Stop and review manually'
      ];

  const choice = await askChoice('What would you like to do?', options, cellCycleConcern ? 1 : 0);

  if (choice.index === 2) {
    // Show full responses
    console.log('');
    console.log('=== Full Agent Responses ===');
    console.log('');
    console.log('--- GPT-5.2 ---');
    console.log(reviewResult.agentResponses?.gpt5_2?.content || 'No response');
    console.log('');
    console.log('--- Claude ---');
    console.log(reviewResult.agentResponses?.claude?.content || 'No response');
    console.log('');
    console.log('--- Gemini ---');
    console.log(reviewResult.agentResponses?.gemini?.content || 'No response');
    console.log('');

    const finalChoice = await askChoice('After reviewing, what would you like to do?', options.slice(0, 2), 0);
    return {
      proceed: true,  // Always true - either accept or re-cluster (not abort)
      action: finalChoice.index === 0 ? 'accept' : (cellCycleConcern ? 'recluster_cell_cycle' : 'recluster_resolution'),
      reason: `User decision after viewing details: ${finalChoice.value}`
    };
  }

  if (choice.index === 3) {
    return { proceed: false, action: 'abort', reason: 'User aborted to review manually' };
  }

  // Handle different user choices
  if (choice.index === 0) {
    // PROCEED - Accept clustering
    return { proceed: true, action: 'accept', reason: `User decision: ${choice.value}` };
  } else if (choice.index === 1) {
    // RE-CLUSTER or ADJUST RESOLUTION
    const action = cellCycleConcern ? 'recluster_cell_cycle' : 'recluster_resolution';
    return {
      proceed: true,  // Don't abort - we want to re-run
      action,
      reason: `User decision: ${choice.value}`
    };
  }
}

/**
 * Helper: Extract QC thresholds from agent response text
 */
function extractThresholds(responseText) {
  if (!responseText) return { nFeature_min: 200, nFeature_max: 6000, percent_mt_max: 10 };

  const minMatch = responseText.match(/nFeature[_\s]*min[:\s]*(\d+)/i);
  const maxMatch = responseText.match(/nFeature[_\s]*max[:\s]*(\d+)/i);
  const mtMatch = responseText.match(/percent[_\s]*mt[_\s]*max[:\s]*(\d+)/i);

  return {
    nFeature_min: minMatch ? parseInt(minMatch[1]) : 200,
    nFeature_max: maxMatch ? parseInt(maxMatch[1]) : 6000,
    percent_mt_max: mtMatch ? parseInt(mtMatch[1]) : 10
  };
}

/**
 * Helper: Extract PC range from agent response text
 */
function extractPCRange(responseText) {
  if (!responseText) return { min_pc: 1, max_pc: 20, reasoning: 'default' };

  const maxMatch = responseText.match(/max[_\s]*pc[:\s]*(\d+)/i);
  const reasoningMatch = responseText.match(/Reasoning[:\s]*([^\n]{50,200})/i);

  return {
    min_pc: 1,
    max_pc: maxMatch ? parseInt(maxMatch[1]) : 20,
    reasoning: reasoningMatch ? reasoningMatch[1].trim() : 'See full response'
  };
}

/**
 * Helper: Check if agents flagged cell cycle concerns
 */
function checkCellCycleConcern(reviewResult) {
  const allResponses = [
    reviewResult.agentResponses?.gpt5_2?.content || '',
    reviewResult.agentResponses?.claude?.content || '',
    reviewResult.agentResponses?.gemini?.content || ''
  ].join(' ').toLowerCase();

  return allResponses.includes('cell cycle') ||
         allResponses.includes('hist1') ||
         allResponses.includes('proliferat');
}

/**
 * Handle user decision for bulk RNA Stage 4 (DE Analysis)
 */
export async function handleStage4UserDecision(reviewResult, stage4Output, stage3Results, autoResolution, config = {}) {
  console.log('');
  console.log('='.repeat(60));
  console.log('  USER DECISION REQUIRED - DE Analysis Review');
  console.log('='.repeat(60));
  console.log('');

  if (autoResolution) {
    console.log('Auto-resolution attempted but escalated to user:');
    console.log(`  ${autoResolution.reasoning}`);
    console.log('');
  }

  console.log('Agents flagged concerns with the differential expression results.');
  console.log('');

  // Helper to get display name for agent
  const getAgentDisplayName = (agentKey, modelStr) => {
    if (config.singleAgent) {
      const modelName = getModelDisplayName(modelStr);
      const roleMap = { gpt5_2: 'Statistics', claude: 'Pipeline', gemini: 'Biology' };
      return `${modelName} (Role: ${roleMap[agentKey]})`;
    } else {
      const defaultNames = { gpt5_2: 'GPT-5.2 (Stats)', claude: 'Claude (Pipeline)', gemini: 'Gemini (Biology)' };
      return defaultNames[agentKey];
    }
  };

  // Show agent assessments
  console.log('--- Agent Assessments ---');
  console.log('');
  if (reviewResult.responses?.gpt5_2?.success) {
    console.log(`${getAgentDisplayName('gpt5_2', reviewResult.responses.gpt5_2.model)}:`);
    console.log(`  ${extractDecisionSummary(reviewResult.responses.gpt5_2.content)}`);
    console.log('');
  }
  if (reviewResult.responses?.claude?.success) {
    console.log(`${getAgentDisplayName('claude', reviewResult.responses.claude.model)}:`);
    console.log(`  ${extractDecisionSummary(reviewResult.responses.claude.content)}`);
    console.log('');
  }
  if (reviewResult.responses?.gemini?.success) {
    console.log(`${getAgentDisplayName('gemini', reviewResult.responses.gemini.model)}:`);
    console.log(`  ${extractDecisionSummary(reviewResult.responses.gemini.content)}`);
    console.log('');
  }

  console.log('--- DE Analysis Summary ---');
  console.log(`Total DEGs: ${stage4Output.total_degs || 0} (Up: ${stage4Output.num_degs_up || 0}, Down: ${stage4Output.num_degs_down || 0})`);
  console.log(`FDR threshold: ${stage4Output.fdr_threshold || 0.05}`);
  console.log(`LogFC threshold: ${stage4Output.logfc_threshold || 0.585}`);
  console.log(`Method: ${stage4Output.de_method || 'unknown'}`);
  console.log('');

  // Check Stage 3 context for re-analysis strategy
  const batchEffectDetected = stage3Results?.batchEffectDetected || false;
  const outliersDetected = stage3Results?.outliersDetected || false;
  const currentMethod = stage4Output.de_method || 'simpleEdger';

  const options = [
    'APPROVE - Accept current DE results',
    `RE-ANALYSIS - Try different approach (${batchEffectDetected ? 'batch correction' : 'different thresholds'})`,
    'VIEW DETAILS - Show full agent responses',
    'ABORT - Stop and review manually'
  ];

  const choice = await askChoice('What would you like to do?', options, 0);

  if (choice.index === 2) {
    // Show full responses
    console.log('');
    console.log('=== Full Agent Responses ===');
    console.log('');
    console.log('--- GPT-5.2 ---');
    console.log(reviewResult.responses?.gpt5_2?.content || 'No response');
    console.log('');
    console.log('--- Claude ---');
    console.log(reviewResult.responses?.claude?.content || 'No response');
    console.log('');
    console.log('--- Gemini ---');
    console.log(reviewResult.responses?.gemini?.content || 'No response');
    console.log('');

    const finalChoice = await askChoice('After reviewing, what would you like to do?', options.slice(0, 2).concat(['ABORT']), 0);

    if (finalChoice.index === 2) {
      return { proceed: false, action: 'abort', reason: 'User aborted after viewing details' };
    }

    return {
      proceed: finalChoice.index === 0,
      action: finalChoice.index === 0 ? 'approve' : 'reanalysis',
      reason: `User decision after viewing details: ${finalChoice.value}`,
      batchEffectDetected,
      outliersDetected,
      currentMethod
    };
  }

  if (choice.index === 3) {
    return { proceed: false, action: 'abort', reason: 'User aborted to review manually' };
  }

  if (choice.index === 0) {
    // APPROVE
    return { proceed: true, action: 'approve', reason: `User decision: ${choice.value}` };
  } else {
    // RE-ANALYSIS
    return {
      proceed: true,
      action: 'reanalysis',
      reason: `User decision: ${choice.value}`,
      batchEffectDetected,
      outliersDetected,
      currentMethod
    };
  }
}

export default {
  askYesNo,
  askChoice,
  askText,
  askMultiSelect,
  handleStage1UserDecision,
  handleStage2UserDecision,
  handleStage3UserDecision,
  handleStage4UserDecision,
  confirmOutlierRemoval,
  handleScRNAStage2UserDecision,
  handleScRNAStage4UserDecision,
  handleScRNAStage5UserDecision
};
