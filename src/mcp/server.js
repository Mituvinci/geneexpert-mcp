#!/usr/bin/env node

/**
 * GeneExpert MCP Server
 * Multi-agent bioinformatics system with Claude, GPT-4, and Gemini collaboration
 */

import { Server } from '@modelcontextprotocol/sdk/server/index.js';
import { StdioServerTransport } from '@modelcontextprotocol/sdk/server/stdio.js';
import {
  CallToolRequestSchema,
  ListToolsRequestSchema,
} from '@modelcontextprotocol/sdk/types.js';
import dotenv from 'dotenv';
import { tools } from './tools.js';

// Load environment variables
dotenv.config();

// Create server instance
const server = new Server(
  {
    name: 'geneexpert-mcp',
    version: '0.1.0',
  },
  {
    capabilities: {
      tools: {},
    },
  }
);

/**
 * List all available tools
 */
server.setRequestHandler(ListToolsRequestSchema, async () => {
  return {
    tools: tools.map(tool => ({
      name: tool.name,
      description: tool.description,
      inputSchema: tool.inputSchema,
    })),
  };
});

/**
 * Handle tool execution
 */
server.setRequestHandler(CallToolRequestSchema, async (request) => {
  const { name, arguments: args } = request.params;

  // Find the tool
  const tool = tools.find(t => t.name === name);

  if (!tool) {
    return {
      content: [{
        type: 'text',
        text: `Error: Tool '${name}' not found. Available tools: ${tools.map(t => t.name).join(', ')}`
      }],
      isError: true,
    };
  }

  try {
    // Execute the tool handler
    console.error(`[GeneExpert] Executing tool: ${name}`);
    console.error(`[GeneExpert] Arguments:`, JSON.stringify(args, null, 2));

    const result = await tool.handler(args);

    console.error(`[GeneExpert] Tool ${name} completed successfully`);

    return result;
  } catch (error) {
    console.error(`[GeneExpert] Tool ${name} failed:`, error);

    return {
      content: [{
        type: 'text',
        text: `Error executing ${name}: ${error.message}\n\nStack trace:\n${error.stack}`
      }],
      isError: true,
    };
  }
});

/**
 * Start the server
 */
async function main() {
  console.error('[GeneExpert] Starting MCP server...');
  console.error(`[GeneExpert] Loaded ${tools.length} tools`);
  console.error('[GeneExpert] Available tools:');
  tools.forEach(tool => {
    console.error(`  - ${tool.name}: ${tool.description}`);
  });

  const transport = new StdioServerTransport();
  await server.connect(transport);

  console.error('[GeneExpert] MCP server running on stdio');
  console.error('[GeneExpert] Ready to receive requests');
}

// Handle errors
process.on('uncaughtException', (error) => {
  console.error('[GeneExpert] Uncaught exception:', error);
  process.exit(1);
});

process.on('unhandledRejection', (error) => {
  console.error('[GeneExpert] Unhandled rejection:', error);
  process.exit(1);
});

// Start the server
main().catch((error) => {
  console.error('[GeneExpert] Failed to start server:', error);
  process.exit(1);
});
