/**
 * Cost Calculator for LLM API Usage
 *
 * Pricing (as of 2026-01-17):
 * - GPT-5.2: $2.50/1M input tokens, $10.00/1M output tokens
 * - Claude Sonnet 4: $3.00/1M input tokens, $15.00/1M output tokens
 * - Gemini 2.5 Flash: $0.075/1M input tokens, $0.30/1M output tokens
 */

const PRICING = {
  'gpt-5.2': {
    input: 2.50 / 1_000_000,    // $2.50 per 1M tokens
    output: 10.00 / 1_000_000   // $10.00 per 1M tokens
  },
  'claude-sonnet-4': {
    input: 3.00 / 1_000_000,    // $3.00 per 1M tokens
    output: 15.00 / 1_000_000   // $15.00 per 1M tokens
  },
  'gemini-2.5-flash': {
    input: 0.075 / 1_000_000,   // $0.075 per 1M tokens
    output: 0.30 / 1_000_000    // $0.30 per 1M tokens
  }
};

/**
 * Calculate cost for GPT-5.2 API call
 * @param {Object} usage - Token usage from OpenAI API
 * @returns {number} - Cost in USD
 */
export function calculateGPT5Cost(usage) {
  if (!usage) return 0;

  const inputTokens = usage.prompt_tokens || 0;
  const outputTokens = usage.completion_tokens || 0;

  const cost =
    (inputTokens * PRICING['gpt-5.2'].input) +
    (outputTokens * PRICING['gpt-5.2'].output);

  return parseFloat(cost.toFixed(6));
}

/**
 * Calculate cost for Claude Sonnet 4 API call
 * @param {Object} usage - Token usage from Anthropic API
 * @returns {number} - Cost in USD
 */
export function calculateClaudeCost(usage) {
  if (!usage) return 0;

  const inputTokens = usage.input_tokens || 0;
  const outputTokens = usage.output_tokens || 0;

  const cost =
    (inputTokens * PRICING['claude-sonnet-4'].input) +
    (outputTokens * PRICING['claude-sonnet-4'].output);

  return parseFloat(cost.toFixed(6));
}

/**
 * Calculate cost for Gemini 2.5 Flash API call
 * @param {Object} usage - Token usage from Google API
 * @returns {number} - Cost in USD
 */
export function calculateGeminiCost(usage) {
  if (!usage) return 0;

  const inputTokens = usage.promptTokenCount || 0;
  const outputTokens = usage.candidatesTokenCount || 0;

  const cost =
    (inputTokens * PRICING['gemini-2.5-flash'].input) +
    (outputTokens * PRICING['gemini-2.5-flash'].output);

  return parseFloat(cost.toFixed(6));
}

/**
 * Calculate cost for any agent response
 * @param {Object} agentResponse - Agent response with model and usage info
 * @returns {number} - Cost in USD
 */
export function calculateAgentCost(agentResponse) {
  if (!agentResponse || !agentResponse.success) return 0;

  const model = agentResponse.model;

  if (model === 'gpt-5.2' || model.startsWith('gpt-5')) {
    return calculateGPT5Cost(agentResponse.usage);
  } else if (model === 'claude-sonnet-4' || model.startsWith('claude')) {
    return calculateClaudeCost(agentResponse.usage);
  } else if (model === 'gemini-2.5-flash' || model.startsWith('gemini')) {
    return calculateGeminiCost(agentResponse.usage);
  }

  return 0;
}

/**
 * Calculate total cost for all agents in a decision
 * @param {Object} agents - Object with gpt5_2, claude, gemini responses
 * @returns {Object} - Total cost and breakdown
 */
export function calculateDecisionCost(agents) {
  const gpt5Cost = agents.gpt5_2 ? calculateAgentCost(agents.gpt5_2) : 0;
  const claudeCost = agents.claude ? calculateAgentCost(agents.claude) : 0;
  const geminiCost = agents.gemini ? calculateAgentCost(agents.gemini) : 0;

  return {
    total_usd: parseFloat((gpt5Cost + claudeCost + geminiCost).toFixed(6)),
    breakdown: {
      gpt5_2: gpt5Cost,
      claude: claudeCost,
      gemini: geminiCost
    }
  };
}

/**
 * Format cost for display
 * @param {number} cost - Cost in USD
 * @returns {string} - Formatted cost string
 */
export function formatCost(cost) {
  if (cost === 0) return '$0.000';
  if (cost < 0.001) return `$${(cost * 1000).toFixed(3)}µ`; // micro-dollars
  return `$${cost.toFixed(3)}`;
}
