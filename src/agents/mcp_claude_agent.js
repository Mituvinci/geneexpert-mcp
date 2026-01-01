/**
 * MCP-Enabled Claude Agent
 *
 * This agent wraps Claude API with MCP tool calling capabilities.
 * At decision points, Claude can:
 * 1. READ outputs via MCP (BAM stats, count summaries, files)
 * 2. ANALYZE data and detect problems
 * 3. WRITE custom scripts if needed
 * 4. EXECUTE custom solutions
 *
 * This is what makes the system AGENTIC instead of just automation!
 */

import Anthropic from '@anthropic-ai/sdk';
import dotenv from 'dotenv';
import { tools } from '../mcp/tools.js';

dotenv.config();

const anthropic = new Anthropic({
  apiKey: process.env.ANTHROPIC_API_KEY,
});

export class MCPClaudeAgent {
  constructor(config = {}) {
    this.config = config;
    this.model = config.model || 'claude-3-5-sonnet-20241022'; // Sonnet for tool use
    this.verbose = config.verbose || false;
  }

  /**
   * Call Claude with MCP tool capabilities
   * @param {string} prompt - The question/task for Claude
   * @param {Object} context - Additional context (session data, outputs, etc.)
   * @param {Array} availableTools - MCP tools Claude can use (default: all analysis tools)
   * @returns {Object} - Claude's response with any tool calls made
   */
  async callWithTools(prompt, context = {}, availableTools = null) {
    // Default to analysis tools (read-only tools for decision points)
    const mcpTools = availableTools || [
      'read_bam_summary',
      'read_count_summary',
      'read_file',
      'write_custom_script',
      'execute_custom_script'
    ];

    // Convert MCP tools to Anthropic tool format
    const anthropicTools = mcpTools.map(toolName => {
      const mcpTool = tools.find(t => t.name === toolName);
      if (!mcpTool) return null;

      return {
        name: mcpTool.name,
        description: mcpTool.description,
        input_schema: mcpTool.inputSchema
      };
    }).filter(t => t !== null);

    const messages = [
      {
        role: 'user',
        content: this.buildPrompt(prompt, context)
      }
    ];

    const systemPrompt = `You are an expert bioinformatics pipeline agent with access to MCP tools.

Your role at DECISION POINTS:
1. READ the actual outputs using MCP tools (read_bam_summary, read_count_summary, read_file)
2. ANALYZE the data to detect problems (low mapping rates, batch effects, outliers)
3. RECOMMEND solutions based on what you see
4. If needed, WRITE custom R scripts to handle edge cases (batch correction, special filtering)

Available MCP tools:
- read_bam_summary: Read alignment statistics from BAM log files
- read_count_summary: Read count matrix and get summary stats
- read_file: Read any text file (logs, results, etc.)
- write_custom_script: Write custom R/bash scripts for edge cases
- execute_custom_script: Execute custom scripts you write

IMPORTANT:
- Use tools to SEE the actual data before making recommendations
- Don't guess - read the files and check the numbers
- If you detect problems (low QC, batch effects), explain clearly
- Provide specific, actionable recommendations

Current analysis context:
${JSON.stringify(context, null, 2)}`;

    this.log('[MCP Claude Agent] Calling Claude with tool use enabled...');
    this.log(`Available tools: ${anthropicTools.map(t => t.name).join(', ')}`);

    try {
      const response = await anthropic.messages.create({
        model: this.model,
        max_tokens: 4096,
        system: systemPrompt,
        messages,
        tools: anthropicTools
      });

      // Handle tool calls if Claude wants to use MCP tools
      const toolCalls = [];
      let finalText = '';

      for (const block of response.content) {
        if (block.type === 'text') {
          finalText += block.text;
        } else if (block.type === 'tool_use') {
          this.log(`[MCP Claude Agent] Tool call: ${block.name}`);

          // Execute the MCP tool
          const toolResult = await this.executeMCPTool(block.name, block.input);

          toolCalls.push({
            tool: block.name,
            input: block.input,
            output: toolResult
          });

          // Continue conversation with tool result
          messages.push({
            role: 'assistant',
            content: response.content
          });

          messages.push({
            role: 'user',
            content: [{
              type: 'tool_result',
              tool_use_id: block.id,
              content: JSON.stringify(toolResult)
            }]
          });

          // Get Claude's analysis after seeing tool results
          const followUp = await anthropic.messages.create({
            model: this.model,
            max_tokens: 4096,
            system: systemPrompt,
            messages,
            tools: anthropicTools
          });

          for (const followBlock of followUp.content) {
            if (followBlock.type === 'text') {
              finalText += '\n\n' + followBlock.text;
            }
          }
        }
      }

      this.log('[MCP Claude Agent] Response received');

      return {
        success: true,
        model: this.model,
        content: finalText || 'No text response',
        toolCalls,
        usage: response.usage
      };

    } catch (error) {
      console.error('[MCP Claude Agent] Error:', error.message);
      return {
        success: false,
        model: this.model,
        content: `Error: ${error.message}`,
        toolCalls: [],
        error: error.message
      };
    }
  }

  /**
   * Execute an MCP tool
   */
  async executeMCPTool(toolName, toolInput) {
    const mcpTool = tools.find(t => t.name === toolName);

    if (!mcpTool) {
      return {
        error: `Tool ${toolName} not found`
      };
    }

    this.log(`  Executing MCP tool: ${toolName}`);
    this.log(`  Input: ${JSON.stringify(toolInput, null, 2)}`);

    try {
      const result = await mcpTool.handler(toolInput);
      const output = result.content[0].text;

      this.log(`  Output: ${output.substring(0, 200)}...`);

      return {
        success: true,
        output
      };
    } catch (error) {
      this.log(`  Error: ${error.message}`);
      return {
        success: false,
        error: error.message
      };
    }
  }

  /**
   * Build prompt with context
   */
  buildPrompt(prompt, context) {
    let fullPrompt = prompt;

    // Add relevant file paths if available
    if (context.outputDir) {
      fullPrompt += `\n\nOutput directory: ${context.outputDir}`;
    }

    if (context.bamDir) {
      fullPrompt += `\nBAM files directory: ${context.bamDir}`;
    }

    if (context.countsFile) {
      fullPrompt += `\nCount matrix: ${context.countsFile}`;
    }

    fullPrompt += `\n\nPlease use the MCP tools to READ the actual data before making your recommendation.`;

    return fullPrompt;
  }

  /**
   * Logging helper
   */
  log(message) {
    if (this.verbose || true) {
      console.log(message);
    }
  }
}
