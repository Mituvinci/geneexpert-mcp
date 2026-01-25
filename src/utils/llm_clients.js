/**
 * LLM API Client Wrappers
 * Unified interface for OpenAI (GPT-5.2), Anthropic (Claude Sonnet 4), and Google (Gemini)
 */

import OpenAI from 'openai';
import Anthropic from '@anthropic-ai/sdk';
import { GoogleGenerativeAI } from '@google/generative-ai';
import dotenv from 'dotenv';
import { getSystemPrompt } from '../config/prompts.js';

dotenv.config();

// ============================================
// Retry Logic with Exponential Backoff
// ============================================

/**
 * Retry a function with exponential backoff
 * @param {Function} fn - Async function to retry
 * @param {number} maxRetries - Maximum number of retries (default: 3)
 * @param {number} initialDelay - Initial delay in ms (default: 1000)
 * @returns {Promise} Result from fn
 */
async function retryWithBackoff(fn, maxRetries = 3, initialDelay = 1000) {
  let lastError;

  for (let attempt = 0; attempt <= maxRetries; attempt++) {
    try {
      return await fn();
    } catch (error) {
      lastError = error;

      // Check if error is retryable (503, 429, network errors)
      const isRetryable =
        error.message?.includes('503') ||  // Service unavailable
        error.message?.includes('429') ||  // Rate limit
        error.message?.includes('overloaded') ||
        error.message?.includes('timeout') ||
        error.message?.includes('ECONNRESET') ||
        error.message?.includes('ETIMEDOUT') ||
        error.message?.includes('fetch failed') ||  // Network fetch failure
        error.message?.includes('ENOTFOUND') ||     // DNS lookup failed
        error.message?.includes('ECONNREFUSED');    // Connection refused

      if (!isRetryable || attempt === maxRetries) {
        throw error;  // Not retryable or max retries reached
      }

      // Exponential backoff: 1s, 2s, 4s, 8s, ...
      const delay = initialDelay * Math.pow(2, attempt);
      console.log(`[Retry] Attempt ${attempt + 1}/${maxRetries} failed. Retrying in ${delay}ms...`);
      await new Promise(resolve => setTimeout(resolve, delay));
    }
  }

  throw lastError;
}

// ============================================
// OpenAI Client (GPT-5.2 for Stats Agent)
// ============================================

const openai = new OpenAI({
  apiKey: process.env.OPENAI_API_KEY,
});

export async function callGPT5(prompt, options = {}) {
  try {
    // Wrap API call with retry logic
    const response = await retryWithBackoff(async () => {
      // Build user message content (text + optional images)
      let userContent = prompt;
      if (options.images && options.images.length > 0) {
        // Vision mode: create multimodal content
        userContent = [
          { type: 'text', text: prompt }
        ];
        for (const img of options.images) {
          userContent.push({
            type: 'image_url',
            image_url: {
              url: `data:${img.mediaType};base64,${img.data}`
            }
          });
        }
      }

      return await openai.chat.completions.create({
        model: options.model || 'gpt-5.2-2025-12-11',
        messages: [
          {
            role: 'system',
            content: options.systemPrompt || 'You are a statistical expert for bioinformatics analysis. Focus on statistical rigor, threshold validation, and outlier detection.'
          },
          {
            role: 'user',
            content: userContent
          }
        ],
        temperature: options.temperature || 0, // Deterministic for reproducible scientific analysis
        max_completion_tokens: options.maxTokens || 2000,  // GPT-5 uses max_completion_tokens
      });
    });

    return {
      success: true,
      model: 'gpt-5.2',
      content: response.choices[0].message.content,
      usage: response.usage,
      raw: response
    };
  } catch (error) {
    console.error('[GPT-5 Error]:', error.message);
    return {
      success: false,
      model: 'gpt-5.2',
      error: error.message,
      content: null
    };
  }
}

// ============================================
// Anthropic Client (Claude for Pipeline & QC Agents)
// ============================================

const anthropic = new Anthropic({
  apiKey: process.env.ANTHROPIC_API_KEY,
});

export async function callClaude(prompt, options = {}) {
  try {
    // Wrap API call with retry logic
    const response = await retryWithBackoff(async () => {
      // Build user message content (text + optional images)
      let userContent = prompt;
      if (options.images && options.images.length > 0) {
        // Vision mode: create multimodal content array
        userContent = [
          { type: 'text', text: prompt }
        ];
        for (const img of options.images) {
          userContent.push({
            type: 'image',
            source: {
              type: 'base64',
              media_type: img.mediaType,
              data: img.data
            }
          });
        }
      }

      return await anthropic.messages.create({
        model: options.model || 'claude-opus-4-5-20251101',  // Flagship model for fair comparison with GPT-5.2
        max_tokens: options.maxTokens || 4096,
        temperature: options.temperature || 0,
        system: options.systemPrompt || 'You are a bioinformatics pipeline expert. Execute workflows, manage data processing, and ensure biological correctness.',
        messages: [
          {
            role: 'user',
            content: userContent
          }
        ]
      });
    });

    return {
      success: true,
      model: 'claude-sonnet-4',
      content: response.content[0].text,
      usage: response.usage,
      raw: response
    };
  } catch (error) {
    console.error('[Claude Error]:', error.message);
    return {
      success: false,
      model: 'claude',
      error: error.message,
      content: null
    };
  }
}

// ============================================
// Google Gemini Client (for Biology Agent)
// ============================================

const genAI = new GoogleGenerativeAI(process.env.GOOGLE_API_KEY);

export async function callGemini(prompt, options = {}) {
  try {
    // Wrap API call with retry logic
    const result = await retryWithBackoff(async () => {
      const model = genAI.getGenerativeModel({
        model: options.model || 'models/gemini-pro-latest',  // Gemini Pro Latest - most recent flagship model
      });

      // Add system instruction as first message since systemInstruction may not be supported
      const fullPrompt = options.systemPrompt
        ? `${options.systemPrompt}\n\n${prompt}`
        : prompt;

      // Build content (text + optional images)
      let content = fullPrompt;
      if (options.images && options.images.length > 0) {
        // Vision mode: create multimodal parts
        content = [
          { text: fullPrompt }
        ];
        for (const img of options.images) {
          content.push({
            inlineData: {
              mimeType: img.mediaType,
              data: img.data
            }
          });
        }
      }

      return await model.generateContent(content, {
        generationConfig: {
          temperature: options.temperature || 0, // Deterministic for reproducible scientific analysis
        }
      });
    });

    const response = result.response;

    return {
      success: true,
      model: 'gemini-pro-latest',
      content: response.text(),
      usage: response.usageMetadata || null,
      raw: response
    };
  } catch (error) {
    console.error('[Gemini Error]:', error.message);
    return {
      success: false,
      model: 'gemini',
      error: error.message,
      content: null
    };
  }
}

// ============================================
// Multi-LLM Call (Parallel Execution)
// ============================================

export async function callAllAgents(prompt, options = {}) {
  const singleAgent = options.singleAgent;

  // Single-MODEL Multi-AGENT mode (TRUE multi-agent with one model)
  // The same model is called 3 times with 3 different role prompts, then consensus voting applies
  if (singleAgent) {
    console.log(`[Single-Model Multi-Agent Mode] Using ${singleAgent.toUpperCase()} as 3 independent agents...`);
    console.log('[Single-Model Multi-Agent] Same model will play 3 roles separately:');
    console.log('  - Agent 1: Statistical analysis (Stats role)');
    console.log('  - Agent 2: Pipeline/Technical analysis (Pipeline role)');
    console.log('  - Agent 3: Biological analysis (Biology role)');
    console.log('[Single-Model Multi-Agent] Each role gets independent API call → Consensus voting applied');
    console.log('');

    // Get role-specific prompts from options
    const statsPrompt = options.gpt5_2_SystemPrompt;
    const pipelinePrompt = options.claudeSystemPrompt;
    const biologyPrompt = options.geminiSystemPrompt;

    if (singleAgent === 'gpt5.2' || singleAgent === 'gpt4') {
      console.log('[GPT-5.2 Multi-Agent] Calling GPT-5.2 THREE times with different role prompts...');
      console.log('[GPT-5.2 as Agent 1] Stats role...');
      console.log('[GPT-5.2 as Agent 2] Pipeline role...');
      console.log('[GPT-5.2 as Agent 3] Biology role...');
      console.log('');

      // Call GPT-5.2 THREE times in parallel with different role prompts
      const [statsResult, pipelineResult, biologyResult] = await Promise.all([
        callGPT5(prompt, {
          systemPrompt: statsPrompt,
          ...options.gpt5_2_Options
        }),
        callGPT5(prompt, {
          systemPrompt: pipelinePrompt,
          ...options.gpt5_2_Options
        }),
        callGPT5(prompt, {
          systemPrompt: biologyPrompt,
          ...options.gpt5_2_Options
        })
      ]);

      // Display results
      console.log('[GPT-5.2 as Stats Agent] ' + (statsResult.success ? '✓ Response received' : '✗ Failed'));
      if (statsResult.success && statsResult.content) {
        console.log('-'.repeat(60));
        console.log(statsResult.content);
        console.log('-'.repeat(60));
      }
      console.log('');

      console.log('[GPT-5.2 as Pipeline Agent] ' + (pipelineResult.success ? '✓ Response received' : '✗ Failed'));
      if (pipelineResult.success && pipelineResult.content) {
        console.log('-'.repeat(60));
        console.log(pipelineResult.content);
        console.log('-'.repeat(60));
      }
      console.log('');

      console.log('[GPT-5.2 as Biology Agent] ' + (biologyResult.success ? '✓ Response received' : '✗ Failed'));
      if (biologyResult.success && biologyResult.content) {
        console.log('-'.repeat(60));
        console.log(biologyResult.content);
        console.log('-'.repeat(60));
      }
      console.log('');

      return {
        gpt5_2: statsResult,       // GPT-5.2 playing Stats role
        claude: pipelineResult,    // GPT-5.2 playing Pipeline role
        gemini: biologyResult,     // GPT-5.2 playing Biology role
        allSuccessful: statsResult.success && pipelineResult.success && biologyResult.success,
        singleModelMultiAgent: true,  // Flag for logging
        model: 'gpt5.2'  // Track which model was used
      };

    } else if (singleAgent === 'claude') {
      console.log('[Claude Multi-Agent] Calling Claude THREE times with different role prompts...');
      console.log('[Claude as Agent 1] Stats role...');
      console.log('[Claude as Agent 2] Pipeline role...');
      console.log('[Claude as Agent 3] Biology role...');
      console.log('');

      // Call Claude THREE times in parallel with different role prompts
      const [statsResult, pipelineResult, biologyResult] = await Promise.all([
        callClaude(prompt, {
          systemPrompt: statsPrompt,
          ...options.claudeOptions
        }),
        callClaude(prompt, {
          systemPrompt: pipelinePrompt,
          ...options.claudeOptions
        }),
        callClaude(prompt, {
          systemPrompt: biologyPrompt,
          ...options.claudeOptions
        })
      ]);

      // Display results
      console.log('[Claude as Stats Agent] ' + (statsResult.success ? '✓ Response received' : '✗ Failed'));
      if (statsResult.success && statsResult.content) {
        console.log('-'.repeat(60));
        console.log(statsResult.content);
        console.log('-'.repeat(60));
      }
      console.log('');

      console.log('[Claude as Pipeline Agent] ' + (pipelineResult.success ? '✓ Response received' : '✗ Failed'));
      if (pipelineResult.success && pipelineResult.content) {
        console.log('-'.repeat(60));
        console.log(pipelineResult.content);
        console.log('-'.repeat(60));
      }
      console.log('');

      console.log('[Claude as Biology Agent] ' + (biologyResult.success ? '✓ Response received' : '✗ Failed'));
      if (biologyResult.success && biologyResult.content) {
        console.log('-'.repeat(60));
        console.log(biologyResult.content);
        console.log('-'.repeat(60));
      }
      console.log('');

      return {
        gpt5_2: statsResult,       // Claude playing Stats role
        claude: pipelineResult,    // Claude playing Pipeline role
        gemini: biologyResult,     // Claude playing Biology role
        allSuccessful: statsResult.success && pipelineResult.success && biologyResult.success,
        singleModelMultiAgent: true,  // Flag for logging
        model: 'claude'  // Track which model was used
      };

    } else if (singleAgent === 'gemini') {
      console.log('[Gemini Multi-Agent] Calling Gemini THREE times with different role prompts...');
      console.log('[Gemini as Agent 1] Stats role...');
      console.log('[Gemini as Agent 2] Pipeline role...');
      console.log('[Gemini as Agent 3] Biology role...');
      console.log('');

      // Call Gemini THREE times in parallel with different role prompts
      const [statsResult, pipelineResult, biologyResult] = await Promise.all([
        callGemini(prompt, {
          systemPrompt: statsPrompt,
          ...options.geminiOptions
        }),
        callGemini(prompt, {
          systemPrompt: pipelinePrompt,
          ...options.geminiOptions
        }),
        callGemini(prompt, {
          systemPrompt: biologyPrompt,
          ...options.geminiOptions
        })
      ]);

      // Display results
      console.log('[Gemini as Stats Agent] ' + (statsResult.success ? '✓ Response received' : '✗ Failed'));
      if (statsResult.success && statsResult.content) {
        console.log('-'.repeat(60));
        console.log(statsResult.content);
        console.log('-'.repeat(60));
      }
      console.log('');

      console.log('[Gemini as Pipeline Agent] ' + (pipelineResult.success ? '✓ Response received' : '✗ Failed'));
      if (pipelineResult.success && pipelineResult.content) {
        console.log('-'.repeat(60));
        console.log(pipelineResult.content);
        console.log('-'.repeat(60));
      }
      console.log('');

      console.log('[Gemini as Biology Agent] ' + (biologyResult.success ? '✓ Response received' : '✗ Failed'));
      if (biologyResult.success && biologyResult.content) {
        console.log('-'.repeat(60));
        console.log(biologyResult.content);
        console.log('-'.repeat(60));
      }
      console.log('');

      return {
        gpt5_2: statsResult,       // Gemini playing Stats role
        claude: pipelineResult,    // Gemini playing Pipeline role
        gemini: biologyResult,     // Gemini playing Biology role
        allSuccessful: statsResult.success && pipelineResult.success && biologyResult.success,
        singleModelMultiAgent: true,  // Flag for logging
        model: 'gemini'  // Track which model was used
      };
    }
  }

  // Multi-agent mode (default)
  // Get role assignments from options (with defaults)
  const roleAssignments = options.roleAssignments || {
    gptRole: 'stats',
    claudeRole: 'pipeline',
    geminiRole: 'biology'
  };

  // Helper function to generate role label
  const getRoleLabel = (role) => {
    const labels = {
      stats: 'Stats Agent',
      pipeline: 'Pipeline Agent',
      biology: 'Biology Agent'
    };
    return labels[role] || role;
  };

  const gptLabel = `GPT-5.2 (${getRoleLabel(roleAssignments.gptRole)})`;
  const claudeLabel = `Claude Opus 4.5 (${getRoleLabel(roleAssignments.claudeRole)})`;
  const geminiLabel = `Gemini Pro (${getRoleLabel(roleAssignments.geminiRole)})`;

  console.log('[Multi-Agent] Sending prompt to all 3 foundation models in parallel...');
  console.log(`[${gptLabel}] Analyzing...`);
  console.log(`[${claudeLabel}] Analyzing...`);
  console.log(`[${geminiLabel}] Analyzing...`);
  console.log('');

  const [gpt5Result, claudeResult, geminiResult] = await Promise.all([
    callGPT5(prompt, {
      systemPrompt: options.gpt5_2_SystemPrompt,
      ...options.gpt5_2_Options
    }),
    callClaude(prompt, {
      systemPrompt: options.claudeSystemPrompt,
      ...options.claudeOptions
    }),
    callGemini(prompt, {
      systemPrompt: options.geminiSystemPrompt,
      ...options.geminiOptions
    })
  ]);

  // Show FULL agent responses
  console.log(`[${gptLabel}] ` + (gpt5Result.success ? '✓ Response received' : '✗ Failed'));
  if (gpt5Result.success && gpt5Result.content) {
    console.log('-'.repeat(60));
    console.log(gpt5Result.content);
    console.log('-'.repeat(60));
  }
  console.log('');

  console.log(`[${claudeLabel}] ` + (claudeResult.success ? '✓ Response received' : '✗ Failed'));
  if (claudeResult.success && claudeResult.content) {
    console.log('-'.repeat(60));
    console.log(claudeResult.content);
    console.log('-'.repeat(60));
  }
  console.log('');

  console.log(`[${geminiLabel}] ` + (geminiResult.success ? '✓ Response received' : '✗ Failed'));
  if (geminiResult.success && geminiResult.content) {
    console.log('-'.repeat(60));
    console.log(geminiResult.content);
    console.log('-'.repeat(60));
  }
  console.log('');

  return {
    gpt5_2: gpt5Result,
    claude: claudeResult,
    gemini: geminiResult,
    allSuccessful: gpt5Result.success && claudeResult.success && geminiResult.success
  };
}

// ============================================
// Sequential Chain (NEW - for ICML comparison)
// ============================================

/**
 * Call agents sequentially: GPT-5.2 → Gemini → Claude
 * Each agent sees previous agents' responses
 *
 * Flow: Stats (GPT-5.2) → Biology (Gemini) → Pipeline (Claude)
 * Claude acts as synthesis agent with full context
 */
export async function callAllAgentsSequential(prompt, options = {}) {
  // Get role assignments from options (with defaults)
  const roleAssignments = options.roleAssignments || {
    gptRole: 'stats',
    claudeRole: 'pipeline',
    geminiRole: 'biology'
  };

  // Helper function to generate role label
  const getRoleLabel = (role) => {
    const labels = {
      stats: 'Stats Agent',
      pipeline: 'Pipeline Agent',
      biology: 'Biology Agent'
    };
    return labels[role] || role;
  };

  const gptLabel = `GPT-5.2 (${getRoleLabel(roleAssignments.gptRole)})`;
  const claudeLabel = `Claude Opus 4.5 (${getRoleLabel(roleAssignments.claudeRole)})`;
  const geminiLabel = `Gemini Pro (${getRoleLabel(roleAssignments.geminiRole)})`;

  console.log('[Sequential Chain Mode] Calling agents in sequence: GPT-5.2 → Gemini → Claude');
  console.log('[Sequential Chain] Each agent will see previous responses for informed synthesis');
  console.log('');

  // STEP 1: First agent (GPT-5.2)
  console.log(`[Step 1/3] ${gptLabel} - First assessment...`);
  const gpt5Result = await callGPT5(prompt, {
    systemPrompt: options.gpt5_2_SystemPrompt,
    ...options.gpt5_2_Options
  });

  console.log(`[${gptLabel}] ` + (gpt5Result.success ? '✓ Response received' : '✗ Failed'));
  if (gpt5Result.success && gpt5Result.content) {
    console.log('-'.repeat(60));
    console.log(gpt5Result.content);
    console.log('-'.repeat(60));
  }
  console.log('');

  // STEP 2: Second agent (Gemini)
  console.log(`[Step 2/3] ${geminiLabel} - Second assessment...`);
  console.log(`[${geminiLabel}] Receiving first agent assessment for context`);

  // Build augmented prompt with GPT-5.2's response
  const geminiAugmentedPrompt = gpt5Result.success
    ? `${prompt}\n\n---\n\n**FIRST AGENT ASSESSMENT (${gptLabel}):**\n\n${gpt5Result.content}\n\n---\n\nConsider the first agent's perspective above while providing your assessment. You may agree or disagree based on your expertise.`
    : prompt; // If GPT-5.2 failed, proceed with original prompt

  const geminiResult = await callGemini(geminiAugmentedPrompt, {
    systemPrompt: options.geminiSystemPrompt,
    ...options.geminiOptions
  });

  console.log(`[${geminiLabel}] ` + (geminiResult.success ? '✓ Response received' : '✗ Failed'));
  if (geminiResult.success && geminiResult.content) {
    console.log('-'.repeat(60));
    console.log(geminiResult.content);
    console.log('-'.repeat(60));
  }
  console.log('');

  // STEP 3: Third agent (Claude) - Receives BOTH previous assessments
  console.log(`[Step 3/3] ${claudeLabel} - Final synthesis...`);
  console.log(`[${claudeLabel}] Receiving both previous assessments for synthesis`);

  // Build fully augmented prompt with both previous responses
  let claudeAugmentedPrompt = prompt;

  if (gpt5Result.success || geminiResult.success) {
    claudeAugmentedPrompt += `\n\n---\n\n**PREVIOUS AGENT ASSESSMENTS:**\n\n`;

    if (gpt5Result.success) {
      claudeAugmentedPrompt += `**1. ${gptLabel}:**\n\n${gpt5Result.content}\n\n---\n\n`;
    }

    if (geminiResult.success) {
      claudeAugmentedPrompt += `**2. ${geminiLabel}:**\n\n${geminiResult.content}\n\n---\n\n`;
    }

    claudeAugmentedPrompt += `As the final synthesis agent, consider both previous perspectives and provide your decision. Consider both assessments but make your own independent judgment based on your expertise (${getRoleLabel(roleAssignments.claudeRole)}).`;
  }

  const claudeResult = await callClaude(claudeAugmentedPrompt, {
    systemPrompt: options.claudeSystemPrompt,
    ...options.claudeOptions
  });

  console.log(`[${claudeLabel}] ` + (claudeResult.success ? '✓ Response received' : '✗ Failed'));
  if (claudeResult.success && claudeResult.content) {
    console.log('-'.repeat(60));
    console.log(claudeResult.content);
    console.log('-'.repeat(60));
  }
  console.log('');

  console.log('[Sequential Chain] All 3 agents completed');
  console.log('');

  return {
    gpt5_2: gpt5Result,
    claude: claudeResult,
    gemini: geminiResult,
    allSuccessful: gpt5Result.success && claudeResult.success && geminiResult.success,
    sequentialMode: true
  };
}

// ============================================
// Agent-Specific Helpers
// ============================================

/**
 * Stats Agent (GPT-5.2) - Statistical validation
 */
export async function askStatsAgent(question, context = {}) {
  const prompt = `
STATISTICAL REVIEW REQUEST:

Question: ${question}

Context:
${JSON.stringify(context, null, 2)}

Please provide:
1. Statistical validation (are thresholds appropriate?)
2. Assumption checks (are statistical assumptions met?)
3. Outlier detection (data-driven analysis)
4. Recommendations (what should be changed?)

Be rigorous and specific.
`;

  return callGPT5(prompt, {
    systemPrompt: `You are a statistical expert for genomics.
Your role is limited to statistical validity, uncertainty, and evaluation metrics.

Rules:
- Do NOT provide biological interpretation or pipeline implementation details.
- Base conclusions only on information provided in the input.
- If assumptions are required or data is insufficient, respond with UNCERTAIN.

Output:
- One clear recommendation
- A confidence level (HIGH / MEDIUM / LOW)
- A brief statistical justification`,
    temperature: 0 // Deterministic for reproducible scientific analysis
  });
}

/**
 * Biology Agent (Gemini) - Biological interpretation
 */
export async function askBiologyAgent(question, context = {}) {
  const prompt = `
BIOLOGICAL INTERPRETATION REQUEST:

Question: ${question}

Context:
${JSON.stringify(context, null, 2)}

Please provide:
1. Biological plausibility (do results make biological sense?)
2. Pathway coherence (are gene patterns consistent?)
3. Literature context (alignment with known biology)
4. Recommendations (next steps for validation)

Focus on biological insight.
`;

  return callGemini(prompt, {
    systemPrompt: `You are a molecular biology expert.
Your role is biological interpretation and experimental plausibility.

Rules:
- Do NOT comment on statistical methods or pipeline mechanics.
- Base reasoning on known biological principles only.
- If evidence is weak or ambiguous, respond with UNCERTAIN.

Output:
- Biological assessment
- Plausibility check
- Confidence level`,
    temperature: 0
  });
}

/**
 * QC Agent (Claude) - Quality control
 */
export async function askQCAgent(question, context = {}) {
  const prompt = `
QUALITY CONTROL REVIEW REQUEST:

Question: ${question}

Context:
${JSON.stringify(context, null, 2)}

Please provide:
1. Technical quality assessment (PCA, MDS, density plots)
2. Batch effect detection (sample clustering issues)
3. Outlier identification (samples to remove?)
4. Recommendations (proceed or re-process?)

Be thorough and critical.
`;

  return callClaude(prompt, {
    systemPrompt: `You are a bioinformatics pipeline expert.
Your role is to assess technical feasibility, QC risks, and execution reliability.

Rules:
- Use tool-grounded reasoning when files or scripts are available.
- Do NOT speculate about biological meaning.
- Flag edge cases such as batch effects, low read depth, or missing metadata.

Output:
- Recommended action (AUTOMATION / ADAPTATION / UNCERTAIN)
- Key technical risk factors
- Confidence level`,
    temperature: 0
  });
}

/**
 * Pipeline Agent (Claude) - Workflow execution
 */
export async function askPipelineAgent(question, context = {}) {
  const prompt = `
PIPELINE EXECUTION REQUEST:

Question: ${question}

Context:
${JSON.stringify(context, null, 2)}

Please provide:
1. Workflow recommendations (what tools to run?)
2. Parameter suggestions (what settings to use?)
3. Next steps (what comes after this stage?)
4. Error handling (if issues arise)

Be practical and actionable.
`;

  return callClaude(prompt, {
    systemPrompt: `You are a bioinformatics pipeline expert.
Your role is to assess technical feasibility, QC risks, and execution reliability.

Rules:
- Use tool-grounded reasoning when files or scripts are available.
- Do NOT speculate about biological meaning.
- Flag edge cases such as batch effects, low read depth, or missing metadata.

Output:
- Recommended action (AUTOMATION / ADAPTATION / UNCERTAIN)
- Key technical risk factors
- Confidence level`,
    temperature: 0
  });
}

// ============================================
// Export all functions
// ============================================

export default {
  callGPT5,
  callClaude,
  callGemini,
  callAllAgents,
  askStatsAgent,
  askBiologyAgent,
  askQCAgent,
  askPipelineAgent
};
