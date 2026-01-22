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
        model: options.model || 'claude-sonnet-4-20250514',
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
        model: options.model || 'models/gemini-2.5-flash'
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
      model: 'gemini-2.5-flash',
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

  // Single-agent mode for experiments (agent performs all 3 roles with stage-specific format)
  if (singleAgent) {
    console.log(`[Single-Agent Mode] Using only ${singleAgent.toUpperCase()} agent...`);
    console.log('[Single-Agent] This agent will perform ALL roles: Statistical + Pipeline + Biological analysis');
    console.log('');

    let result;
    if (singleAgent === 'gpt5.2' || singleAgent === 'gpt4') {
      console.log('[GPT-5.2 Full-Stack Agent] Analyzing (Stats + Pipeline + Biology)...');

      // CRITICAL FIX: Use stage-specific prompt if provided (for correct decision format)
      // Orchestrator provides stage prompt as options.gpt5_2_SystemPrompt
      const stagePrompt = options.gpt5_2_SystemPrompt;

      if (stagePrompt) {
        console.log('[GPT-5.2] Using stage-specific prompt (correct decision format)');
        result = await callGPT5(prompt, {
          systemPrompt: stagePrompt,
          ...options.gpt5_2_Options
        });
      } else {
        // Fallback to legacy combined prompt (for non-staged usage)
        console.log('[GPT-5.2] Using legacy combined prompt (fallback)');
        const combinedPrompt = getSystemPrompt('gpt5.2', true);
        result = await callGPT5(prompt, {
          systemPrompt: combinedPrompt,
          ...options.gpt5_2_Options
        });
      }

      console.log('[GPT-5.2 Full-Stack Agent] ' + (result.success ? '✓ Response received' : '✗ Failed'));
      if (result.success && result.content) {
        console.log('-'.repeat(60));
        console.log(result.content);
        console.log('-'.repeat(60));
      }
      console.log('');

      return {
        gpt5_2: result,
        claude: { success: false, model: 'claude', content: null, skipped: true },
        gemini: { success: false, model: 'gemini', content: null, skipped: true },
        allSuccessful: result.success,
        singleAgentMode: true
      };
    } else if (singleAgent === 'claude') {
      console.log('[Claude Full-Stack Agent] Analyzing (Stats + Pipeline + Biology + MCP Tools)...');

      // CRITICAL FIX: Use stage-specific prompt if provided (for correct decision format)
      // Orchestrator provides stage prompt as options.claudeSystemPrompt
      const stagePrompt = options.claudeSystemPrompt;

      if (stagePrompt) {
        console.log('[Claude] Using stage-specific prompt (correct decision format)');
        result = await callClaude(prompt, {
          systemPrompt: stagePrompt,
          ...options.claudeOptions
        });
      } else {
        // Fallback to legacy combined prompt (for non-staged usage)
        console.log('[Claude] Using legacy combined prompt (fallback)');
        const combinedPrompt = getSystemPrompt('claude', true);
        result = await callClaude(prompt, {
          systemPrompt: combinedPrompt,
          ...options.claudeOptions
        });
      }

      console.log('[Claude Full-Stack Agent] ' + (result.success ? '✓ Response received' : '✗ Failed'));
      if (result.success && result.content) {
        console.log('-'.repeat(60));
        console.log(result.content);
        console.log('-'.repeat(60));
      }
      console.log('');

      return {
        gpt5_2: { success: false, model: 'gpt-5.2', content: null, skipped: true },
        claude: result,
        gemini: { success: false, model: 'gemini', content: null, skipped: true },
        allSuccessful: result.success,
        singleAgentMode: true
      };
    } else if (singleAgent === 'gemini') {
      console.log('[Gemini Full-Stack Agent] Analyzing (Stats + Pipeline + Biology)...');

      // CRITICAL FIX: Use stage-specific prompt if provided (for correct decision format)
      // Orchestrator provides stage prompt as options.geminiSystemPrompt
      const stagePrompt = options.geminiSystemPrompt;

      if (stagePrompt) {
        console.log('[Gemini] Using stage-specific prompt (correct decision format)');
        result = await callGemini(prompt, {
          systemPrompt: stagePrompt,
          ...options.geminiOptions
        });
      } else {
        // Fallback to legacy combined prompt (for non-staged usage)
        console.log('[Gemini] Using legacy combined prompt (fallback)');
        const combinedPrompt = getSystemPrompt('gemini', true);
        result = await callGemini(prompt, {
          systemPrompt: combinedPrompt,
          ...options.geminiOptions
        });
      }

      console.log('[Gemini Full-Stack Agent] ' + (result.success ? '✓ Response received' : '✗ Failed'));
      if (result.success && result.content) {
        console.log('-'.repeat(60));
        console.log(result.content);
        console.log('-'.repeat(60));
      }
      console.log('');

      return {
        gpt5_2: { success: false, model: 'gpt-5.2', content: null, skipped: true },
        claude: { success: false, model: 'claude', content: null, skipped: true },
        gemini: result,
        allSuccessful: result.success,
        singleAgentMode: true
      };
    }
  }

  // Multi-agent mode (default)
  console.log('[Multi-Agent] Sending prompt to all 3 foundation models in parallel...');
  console.log('[GPT-5.2 Stats Agent] Analyzing from statistical perspective...');
  console.log('[Claude Sonnet 4 Pipeline Agent] Analyzing pipeline requirements...');
  console.log('[Gemini Biology Agent] Analyzing biological context...');
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
  console.log('[GPT-5.2 Stats Agent] ' + (gpt5Result.success ? '✓ Response received' : '✗ Failed'));
  if (gpt5Result.success && gpt5Result.content) {
    console.log('-'.repeat(60));
    console.log(gpt5Result.content);
    console.log('-'.repeat(60));
  }
  console.log('');

  console.log('[Claude Sonnet 4 Pipeline Agent] ' + (claudeResult.success ? '✓ Response received' : '✗ Failed'));
  if (claudeResult.success && claudeResult.content) {
    console.log('-'.repeat(60));
    console.log(claudeResult.content);
    console.log('-'.repeat(60));
  }
  console.log('');

  console.log('[Gemini Biology Agent] ' + (geminiResult.success ? '✓ Response received' : '✗ Failed'));
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
  console.log('[Sequential Chain Mode] Calling agents in sequence: GPT-5.2 → Gemini → Claude');
  console.log('[Sequential Chain] Each agent will see previous responses for informed synthesis');
  console.log('');

  // STEP 1: Stats Agent (GPT-5.2) - First assessment
  console.log('[Step 1/3] GPT-5.2 (Stats Agent) - Statistical assessment...');
  const gpt5Result = await callGPT5(prompt, {
    systemPrompt: options.gpt5_2_SystemPrompt,
    ...options.gpt5_2_Options
  });

  console.log('[GPT-5.2 Stats Agent] ' + (gpt5Result.success ? '✓ Response received' : '✗ Failed'));
  if (gpt5Result.success && gpt5Result.content) {
    console.log('-'.repeat(60));
    console.log(gpt5Result.content);
    console.log('-'.repeat(60));
  }
  console.log('');

  // STEP 2: Biology Agent (Gemini) - Receives GPT-5.2's assessment
  console.log('[Step 2/3] Gemini (Biology Agent) - Biological interpretation...');
  console.log('[Gemini] Receiving Stats Agent assessment for context');

  // Build augmented prompt with GPT-5.2's response
  const geminiAugmentedPrompt = gpt5Result.success
    ? `${prompt}\n\n---\n\n**STATS AGENT ASSESSMENT (GPT-5.2):**\n\n${gpt5Result.content}\n\n---\n\nConsider the Stats Agent's statistical perspective above while providing your biological interpretation. You may agree or disagree with their assessment based on biological reasoning.`
    : prompt; // If GPT-5.2 failed, proceed with original prompt

  const geminiResult = await callGemini(geminiAugmentedPrompt, {
    systemPrompt: options.geminiSystemPrompt,
    ...options.geminiOptions
  });

  console.log('[Gemini Biology Agent] ' + (geminiResult.success ? '✓ Response received' : '✗ Failed'));
  if (geminiResult.success && geminiResult.content) {
    console.log('-'.repeat(60));
    console.log(geminiResult.content);
    console.log('-'.repeat(60));
  }
  console.log('');

  // STEP 3: Pipeline Agent (Claude) - Receives BOTH GPT-5.2 and Gemini assessments
  console.log('[Step 3/3] Claude (Pipeline Agent) - Technical synthesis...');
  console.log('[Claude] Receiving Stats + Biology assessments for final synthesis');

  // Build fully augmented prompt with both previous responses
  let claudeAugmentedPrompt = prompt;

  if (gpt5Result.success || geminiResult.success) {
    claudeAugmentedPrompt += `\n\n---\n\n**PREVIOUS AGENT ASSESSMENTS:**\n\n`;

    if (gpt5Result.success) {
      claudeAugmentedPrompt += `**1. STATS AGENT (GPT-5.2) - Statistical Perspective:**\n\n${gpt5Result.content}\n\n---\n\n`;
    }

    if (geminiResult.success) {
      claudeAugmentedPrompt += `**2. BIOLOGY AGENT (Gemini) - Biological Perspective:**\n\n${geminiResult.content}\n\n---\n\n`;
    }

    claudeAugmentedPrompt += `As the Pipeline Agent, synthesize the statistical and biological perspectives above and provide your final technical decision. Consider both assessments but make your own independent judgment based on pipeline requirements and technical feasibility.`;
  }

  const claudeResult = await callClaude(claudeAugmentedPrompt, {
    systemPrompt: options.claudeSystemPrompt,
    ...options.claudeOptions
  });

  console.log('[Claude Pipeline Agent] ' + (claudeResult.success ? '✓ Response received' : '✗ Failed'));
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
