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
 * @returns {Promise<{proceed: boolean, reason: string}>}
 */
export async function handleStage1UserDecision(reviewResult, stage1Output) {
  console.log('');
  console.log('='.repeat(60));
  console.log('  USER DECISION REQUIRED');
  console.log('='.repeat(60));
  console.log('');
  console.log('The agents could not reach a clear consensus on Stage 1 results.');
  console.log('');

  // Show what each agent said
  console.log('--- Agent Recommendations ---');
  if (reviewResult.responses.gpt5_2.success) {
    const gptResponse = reviewResult.responses.gpt5_2.content;
    console.log(`GPT-5.2 (Stats): ${extractDecisionSummary(gptResponse)}`);
  }
  if (reviewResult.responses.claude.success) {
    const claudeResponse = reviewResult.responses.claude.content;
    console.log(`Claude (Pipeline): ${extractDecisionSummary(claudeResponse)}`);
  }
  if (reviewResult.responses.gemini.success) {
    const geminiResponse = reviewResult.responses.gemini.content;
    console.log(`Gemini (Biology): ${extractDecisionSummary(geminiResponse)}`);
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
 * @returns {Promise<{proceed: boolean, samplesToRemove: string[], reason: string}>}
 */
export async function handleStage2UserDecision(reviewResult, stage2Output) {
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

  // Show problematic samples
  const problematicSamples = [...(stage2Output.warned_samples || []), ...(stage2Output.failed_samples || [])];
  if (problematicSamples.length > 0) {
    console.log('--- Problematic Samples ---');
    problematicSamples.forEach(s => {
      const rate = stage2Output.sample_details?.[s]?.mapping_rate;
      const status = stage2Output.sample_details?.[s]?.status;
      console.log(`  ${s}: ${rate ? (rate * 100).toFixed(1) + '%' : 'N/A'} (${status || 'unknown'})`);
    });
    console.log('');
  }

  // Show agent recommendations
  console.log('--- Agent Recommendations ---');
  if (reviewResult.responses.gpt5_2.success) {
    console.log(`GPT-5.2 (Stats): ${extractDecisionSummary(reviewResult.responses.gpt5_2.content)}`);
  }
  if (reviewResult.responses.claude.success) {
    console.log(`Claude (Pipeline): ${extractDecisionSummary(reviewResult.responses.claude.content)}`);
  }
  if (reviewResult.responses.gemini.success) {
    console.log(`Gemini (Biology): ${extractDecisionSummary(reviewResult.responses.gemini.content)}`);
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
 * @returns {Promise<{proceed: boolean, samplesToRemove: string[], useBatchCorrection: boolean, reason: string}>}
 */
export async function handleStage3UserDecision(reviewResult, stage3Output) {
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

  // Show agent recommendations
  console.log('--- Agent Recommendations ---');
  if (reviewResult.responses.gpt5_2.success) {
    console.log(`GPT-5.2 (Stats): ${extractDecisionSummary(reviewResult.responses.gpt5_2.content)}`);
  }
  if (reviewResult.responses.claude.success) {
    console.log(`Claude (Pipeline): ${extractDecisionSummary(reviewResult.responses.claude.content)}`);
  }
  if (reviewResult.responses.gemini.success) {
    console.log(`Gemini (Biology): ${extractDecisionSummary(reviewResult.responses.gemini.content)}`);
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

export default {
  askYesNo,
  askChoice,
  askText,
  askMultiSelect,
  handleStage1UserDecision,
  handleStage2UserDecision,
  handleStage3UserDecision
};
