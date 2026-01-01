#!/usr/bin/env node

/**
 * Test LLM API Connections
 * Verifies OpenAI, Anthropic, and Google API keys are configured correctly
 */

import { callGPT4, callClaude, callGemini } from './llm_clients.js';

console.log('🧪 Testing LLM API Connections...\n');
console.log('='.repeat(60));

const testPrompt = 'Respond with exactly: "API connection successful"';

async function testAPIs() {
  let allPassed = true;

  // Test OpenAI (GPT-4)
  console.log('\n1. Testing OpenAI (GPT-4)...');
  try {
    const gpt4Result = await callGPT4(testPrompt, {
      systemPrompt: 'You are a test assistant.',
      maxTokens: 50
    });

    if (gpt4Result.success) {
      console.log('   ✅ GPT-4 Connected!');
      console.log(`   Model: ${gpt4Result.model}`);
      console.log(`   Response: ${gpt4Result.content.substring(0, 100)}...`);
      console.log(`   Tokens used: ${gpt4Result.usage?.total_tokens || 'N/A'}`);
    } else {
      console.log('   ❌ GPT-4 Failed!');
      console.log(`   Error: ${gpt4Result.error}`);
      allPassed = false;
    }
  } catch (error) {
    console.log('   ❌ GPT-4 Error!');
    console.log(`   ${error.message}`);
    allPassed = false;
  }

  // Test Anthropic (Claude)
  console.log('\n2. Testing Anthropic (Claude)...');
  try {
    const claudeResult = await callClaude(testPrompt, {
      systemPrompt: 'You are a test assistant.',
      maxTokens: 50
    });

    if (claudeResult.success) {
      console.log('   ✅ Claude Connected!');
      console.log(`   Model: ${claudeResult.model}`);
      console.log(`   Response: ${claudeResult.content.substring(0, 100)}...`);
      console.log(`   Tokens used: ${claudeResult.usage?.input_tokens + claudeResult.usage?.output_tokens || 'N/A'}`);
    } else {
      console.log('   ❌ Claude Failed!');
      console.log(`   Error: ${claudeResult.error}`);
      allPassed = false;
    }
  } catch (error) {
    console.log('   ❌ Claude Error!');
    console.log(`   ${error.message}`);
    allPassed = false;
  }

  // Test Google (Gemini)
  console.log('\n3. Testing Google (Gemini)...');
  try {
    const geminiResult = await callGemini(testPrompt, {
      systemPrompt: 'You are a test assistant.',
    });

    if (geminiResult.success) {
      console.log('   ✅ Gemini Connected!');
      console.log(`   Model: ${geminiResult.model}`);
      console.log(`   Response: ${geminiResult.content.substring(0, 100)}...`);
      console.log(`   Tokens used: ${geminiResult.usage?.totalTokenCount || 'N/A'}`);
    } else {
      console.log('   ❌ Gemini Failed!');
      console.log(`   Error: ${geminiResult.error}`);
      allPassed = false;
    }
  } catch (error) {
    console.log('   ❌ Gemini Error!');
    console.log(`   ${error.message}`);
    allPassed = false;
  }

  // Summary
  console.log('\n' + '='.repeat(60));
  if (allPassed) {
    console.log('✅ All API connections successful!\n');
    console.log('You are ready to use the multi-agent system.');
    console.log('\nNext steps:');
    console.log('  1. Build Coordinator Agent');
    console.log('  2. Test multi-agent collaboration');
    console.log('  3. Run RNA-seq pipeline with agent reviews\n');
  } else {
    console.log('❌ Some API connections failed.\n');
    console.log('Please check:');
    console.log('  1. API keys in .env file are correct');
    console.log('  2. API keys have not expired');
    console.log('  3. You have API credits/quota available\n');
  }
}

// Run tests
testAPIs().catch(error => {
  console.error('\n💥 Test script crashed:', error);
  process.exit(1);
});
