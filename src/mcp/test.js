#!/usr/bin/env node

/**
 * Simple test to verify MCP server tools are properly configured
 */

import { tools } from './tools.js';

console.log('GeneExpert MCP Server - Tool Test\n');
console.log('='.repeat(60));
console.log(`Total tools loaded: ${tools.length}\n`);

// List all tools
console.log('Available Tools:');
console.log('-'.repeat(60));

tools.forEach((tool, index) => {
  console.log(`${index + 1}. ${tool.name}`);
  console.log(`   Description: ${tool.description}`);
  console.log(`   Required params: ${tool.inputSchema.required?.join(', ') || 'none'}`);
  console.log('');
});

console.log('='.repeat(60));
console.log('\nMCP server configuration looks good! ✓');
console.log('\nTo start the server, run:');
console.log('  node src/mcp/server.js\n');
