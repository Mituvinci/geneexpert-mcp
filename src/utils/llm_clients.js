/**
 * LLM API Client Wrappers
 * Unified interface for OpenAI (GPT-4), Anthropic (Claude), and Google (Gemini)
 */

import OpenAI from 'openai';
import Anthropic from '@anthropic-ai/sdk';
import { GoogleGenerativeAI } from '@google/generative-ai';
import dotenv from 'dotenv';

dotenv.config();

// ============================================
// OpenAI Client (GPT-4 for Stats Agent)
// ============================================

const openai = new OpenAI({
  apiKey: process.env.OPENAI_API_KEY,
});

export async function callGPT4(prompt, options = {}) {
  try {
    const response = await openai.chat.completions.create({
      model: options.model || 'gpt-4-turbo-preview',
      messages: [
        {
          role: 'system',
          content: options.systemPrompt || 'You are a statistical expert for bioinformatics analysis. Focus on statistical rigor, threshold validation, and outlier detection.'
        },
        {
          role: 'user',
          content: prompt
        }
      ],
      temperature: options.temperature || 0.3, // Lower for more deterministic statistical analysis
      max_tokens: options.maxTokens || 2000,
    });

    return {
      success: true,
      model: 'gpt-4',
      content: response.choices[0].message.content,
      usage: response.usage,
      raw: response
    };
  } catch (error) {
    console.error('[GPT-4 Error]:', error.message);
    return {
      success: false,
      model: 'gpt-4',
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
    const response = await anthropic.messages.create({
      model: options.model || 'claude-3-haiku-20240307',
      max_tokens: options.maxTokens || 4096,
      temperature: options.temperature || 0.5,
      system: options.systemPrompt || 'You are a bioinformatics pipeline expert. Execute workflows, manage data processing, and ensure biological correctness.',
      messages: [
        {
          role: 'user',
          content: prompt
        }
      ]
    });

    return {
      success: true,
      model: 'claude-3-haiku',
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
    const model = genAI.getGenerativeModel({
      model: options.model || 'models/gemini-2.5-flash'
    });

    // Add system instruction as first message since systemInstruction may not be supported
    const fullPrompt = options.systemPrompt
      ? `${options.systemPrompt}\n\n${prompt}`
      : prompt;

    const result = await model.generateContent(fullPrompt);
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
  console.log('[Multi-Agent] Sending prompt to all 3 foundation models in parallel...');
  console.log('[GPT-4 Stats Agent] Analyzing from statistical perspective...');
  console.log('[Claude Pipeline Agent] Analyzing pipeline requirements...');
  console.log('[Gemini Biology Agent] Analyzing biological context...');
  console.log('');

  const [gpt4Result, claudeResult, geminiResult] = await Promise.all([
    callGPT4(prompt, {
      systemPrompt: options.gpt4SystemPrompt,
      ...options.gpt4Options
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
  console.log('[GPT-4 Stats Agent] ' + (gpt4Result.success ? '✓ Response received' : '✗ Failed'));
  if (gpt4Result.success && gpt4Result.content) {
    console.log('─'.repeat(60));
    console.log(gpt4Result.content);
    console.log('─'.repeat(60));
  }
  console.log('');

  console.log('[Claude Pipeline Agent] ' + (claudeResult.success ? '✓ Response received' : '✗ Failed'));
  if (claudeResult.success && claudeResult.content) {
    console.log('─'.repeat(60));
    console.log(claudeResult.content);
    console.log('─'.repeat(60));
  }
  console.log('');

  console.log('[Gemini Biology Agent] ' + (geminiResult.success ? '✓ Response received' : '✗ Failed'));
  if (geminiResult.success && geminiResult.content) {
    console.log('─'.repeat(60));
    console.log(geminiResult.content);
    console.log('─'.repeat(60));
  }
  console.log('');

  return {
    gpt4: gpt4Result,
    claude: claudeResult,
    gemini: geminiResult,
    allSuccessful: gpt4Result.success && claudeResult.success && geminiResult.success
  };
}

// ============================================
// Agent-Specific Helpers
// ============================================

/**
 * Stats Agent (GPT-4) - Statistical validation
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

  return callGPT4(prompt, {
    systemPrompt: 'You are a statistical expert for genomics. Focus on: threshold validation (FDR, logFC), multiple testing correction, sample size adequacy, outlier detection, and statistical assumptions.',
    temperature: 0.2 // Very deterministic for stats
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
    systemPrompt: 'You are a molecular biology and genomics expert. Interpret results through biological lens: pathway analysis, gene function, regulatory networks, disease mechanisms, and experimental validation.',
    temperature: 0.5
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
    systemPrompt: 'You are a quality control expert for RNA-seq/ATAC-seq. Monitor technical artifacts: batch effects, outliers, low-quality samples, normalization issues, and sequencing depth problems.',
    temperature: 0.3
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
    systemPrompt: 'You are a bioinformatics pipeline expert. Guide RNA-seq/ATAC-seq workflows: tool selection, parameter optimization, quality checks, and troubleshooting. Ensure biological best practices.',
    temperature: 0.4
  });
}

// ============================================
// Export all functions
// ============================================

export default {
  callGPT4,
  callClaude,
  callGemini,
  callAllAgents,
  askStatsAgent,
  askBiologyAgent,
  askQCAgent,
  askPipelineAgent
};
