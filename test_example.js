#!/usr/bin/env node

/**
 * Simple Example: Test the Coordinator Yourself
 *
 * This shows you how to:
 * 1. Import the coordinator
 * 2. Ask a question to all 3 agents
 * 3. See their responses and consensus
 */

import { createCoordinator } from './src/coordinator/orchestrator.js';
import dotenv from 'dotenv';

dotenv.config();

async function myTest() {
  console.log('🧪 My Custom Coordinator Test\n');

  // Step 1: Create a coordinator
  const coordinator = createCoordinator({ verbose: true });

  // Step 2: Ask a question
  const result = await coordinator.consultAgents(
    "Is it safe to use an FDR threshold of 0.01 for a small RNA-seq study with 2 vs 2 samples?",
    {
      sampleSize: 4,
      design: "2 control vs 2 treatment"
    },
    'threshold' // Decision type: 'threshold', 'sample_removal', 'method_change'
  );

  // Step 3: See the results
  console.log(result.report);

  // Step 4: See individual agent responses
  console.log('\n\n📝 Individual Agent Responses:\n');

  console.log('Stats Agent (GPT-4) says:');
  console.log(result.responses.gpt4.content);
  console.log('\n' + '='.repeat(60) + '\n');

  console.log('Pipeline Agent (Claude) says:');
  console.log(result.responses.claude.content);
  console.log('\n' + '='.repeat(60) + '\n');

  console.log('Biology Agent (Gemini) says:');
  console.log(result.responses.gemini.content);
  console.log('\n' + '='.repeat(60) + '\n');

  // Step 5: Check the decision
  if (result.consensus.decision === 'approve') {
    console.log('✅ Agents approved your approach!');
  } else {
    console.log('⚠️  Agents recommend caution or have concerns');
  }
}

// Run the test
myTest().catch(error => {
  console.error('Error:', error.message);
});
