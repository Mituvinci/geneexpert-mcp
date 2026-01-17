/**
 * Test JSON Logging System
 *
 * Creates a mock agent conversation and saves it in both plain text and JSON
 */

import { Logger } from './src/utils/logger.js';
import { generateSummaryReport } from './src/utils/metrics.js';
import fs from 'fs';
import path from 'path';

console.log('🧪 Testing JSON Logging System\n');

// Create test output directory
const testOutputDir = './test_output';
if (!fs.existsSync(testOutputDir)) {
  fs.mkdirSync(testOutputDir, { recursive: true });
}

// Initialize logger with mock config
const mockConfig = {
  input: '/data/test_dataset',
  output: testOutputDir,
  organism: 'mouse',
  comparison: 'test_comparison',
  singleAgent: null  // Multi-agent mode
};

const logger = new Logger(testOutputDir, 'test_analysis', mockConfig);

console.log('✅ Logger initialized');
console.log(`   Output directory: ${testOutputDir}\n`);

// Mock agent responses (simulating GPT-5.2, Claude, Gemini)
const mockAgents = {
  gpt5_2: {
    model: 'gpt-5.2-2025-12-11',
    success: true,
    content: `Statistical Assessment:
- Sample size n=2 is below recommended minimum (n≥3)
- Statistical power will be limited
- Recommendation: ADAPTATION with robust methods

Confidence: HIGH
Reasoning: Small sample size requires specialized statistical approach`,
    usage: {
      prompt_tokens: 450,
      completion_tokens: 230,
      total_tokens: 680
    }
  },
  claude: {
    model: 'claude-sonnet-4-20250514',
    success: true,
    content: `Pipeline/QC Assessment:
- Dataset: 2 samples per group (edge case)
- Risk: Low statistical power, batch effects could dominate
- Recommendation: ADAPTATION (custom QC pipeline)

Confidence: MEDIUM
Technical risks: Limited replication, sensitivity to outliers`,
    usage: {
      input_tokens: 450,
      output_tokens: 310
    },
    mcp_tools_used: [
      {
        tool: 'read_file',
        file: '/data/test_dataset/metadata.txt',
        success: true
      }
    ]
  },
  gemini: {
    model: 'gemini-2.5-flash',
    success: true,
    content: `Biological Assessment:
- Mouse stroke model is well-established
- n=2 per group is acceptable for pilot studies
- Recommendation: AUTOMATION (standard pipeline is fine)

Confidence: MEDIUM
Plausibility: Biological signal likely detectable despite low n`,
    usage: {
      promptTokenCount: 450,
      candidatesTokenCount: 280
    }
  }
};

// Mock consensus
const mockConsensus = {
  decision: 'adaptation',
  confidence: 0.72,
  votingMethod: 'majority',
  votes: {
    automation: 1,
    adaptation: 2,
    uncertain: 0
  },
  reasoning: 'Majority of agents (2/3) recommend ADAPTATION due to small sample size',
  recommendations: [
    {
      priority: 'HIGH',
      action: 'normalization',
      message: 'Use robust normalization (TMM or DESeq2)'
    },
    {
      priority: 'MEDIUM',
      action: 'qc',
      message: 'Manual QC inspection recommended'
    }
  ]
};

// Mock step info
const mockStep = {
  name: 'Automation vs Adaptation Decision',
  description: 'Agents decide analysis approach',
  decision_id: 'test_dataset_step1_approach',
  decisionType: 'approach_decision'
};

const mockUserPrompt = `Dataset: test_dataset
Organism: mouse
Samples: n=2 per group (control: cont1, cont2; treatment: ips1, ips2)
Data type: FASTQ paired-end

Question: Should we use AUTOMATION (template pipeline) or ADAPTATION (custom solution)?`;

// Log the conversation
console.log('📝 Logging mock agent conversation...\n');
logger.logAgentConversation(mockStep, mockAgents, mockConsensus, mockUserPrompt);
console.log('✅ Agent conversation logged\n');

// Simulate another decision
const mockStep2 = {
  name: 'FDR Threshold Selection',
  description: 'Agents decide statistical threshold',
  decision_id: 'test_dataset_step2_threshold',
  decisionType: 'threshold'
};

const mockAgents2 = {
  gpt5_2: {
    model: 'gpt-5.2-2025-12-11',
    success: true,
    content: 'Recommend FDR < 0.05 (standard threshold)\nConfidence: HIGH',
    usage: { prompt_tokens: 200, completion_tokens: 50, total_tokens: 250 }
  },
  claude: {
    model: 'claude-sonnet-4-20250514',
    success: true,
    content: 'FDR < 0.05 is appropriate\nConfidence: HIGH',
    usage: { input_tokens: 200, output_tokens: 50 }
  },
  gemini: {
    model: 'gemini-2.5-flash',
    success: true,
    content: 'Standard FDR < 0.05 acceptable\nConfidence: MEDIUM',
    usage: { promptTokenCount: 200, candidatesTokenCount: 50 }
  }
};

const mockConsensus2 = {
  decision: 'FDR_0.05',
  confidence: 0.90,
  votingMethod: 'unanimous',
  votes: { 'FDR_0.05': 3 },
  reasoning: 'All agents agree on standard FDR threshold',
  recommendations: []
};

logger.logAgentConversation(mockStep2, mockAgents2, mockConsensus2, 'What FDR threshold should we use?');
console.log('✅ Second decision logged\n');

// Finalize (write JSON files)
const mockSession = {
  dataInfo: {
    type: 'fastq',
    samples: ['cont1', 'cont2', 'ips1', 'ips2'],
    pairedEnd: true,
    groups: {
      control: ['cont1', 'cont2'],
      treatment: ['ips1', 'ips2']
    }
  },
  config: mockConfig,
  scriptPath: path.join(testOutputDir, 'test_script_v1.sh'),
  scriptType: 'ADAPTATION',
  success: true,
  attempts: 1,
  exitCode: 0,
  currentStep: 1,
  steps: [mockStep, mockStep2],
  decisions: {},
  thresholds: { fdr: 0.05, logFC: 1.0 }
};

console.log('💾 Finalizing logs (writing JSON files)...\n');
logger.finalize(mockSession);

// Check what files were created
console.log('\n📂 Files created:');
const files = fs.readdirSync(testOutputDir);
files.forEach(file => {
  const filePath = path.join(testOutputDir, file);
  const stats = fs.statSync(filePath);
  console.log(`   ${file} (${(stats.size / 1024).toFixed(2)} KB)`);
});

// Read and display JSON summary
console.log('\n📊 JSON Decision Summary:');
const decisionsJson = JSON.parse(
  fs.readFileSync(path.join(testOutputDir, 'test_analysis_agent_decisions.json'), 'utf-8')
);

console.log(`   Session ID: ${decisionsJson.session_id}`);
console.log(`   System: ${decisionsJson.system}`);
console.log(`   Total Cost: $${decisionsJson.costs.total_usd.toFixed(3)}`);
console.log(`     - GPT-5.2: $${decisionsJson.costs.breakdown.gpt5_2.toFixed(3)}`);
console.log(`     - Claude:  $${decisionsJson.costs.breakdown.claude.toFixed(3)}`);
console.log(`     - Gemini:  $${decisionsJson.costs.breakdown.gemini.toFixed(3)}`);
console.log(`   Decisions Made: ${decisionsJson.decisions.length}`);

decisionsJson.decisions.forEach((decision, i) => {
  console.log(`\n   Decision ${i + 1}: ${decision.step_name}`);
  console.log(`     - Final Decision: ${decision.consensus.decision}`);
  console.log(`     - Confidence: ${decision.consensus.confidence_label} (${(decision.consensus.confidence_score * 100).toFixed(0)}%)`);
  console.log(`     - Cost: $${decision.costs.total_usd.toFixed(4)}`);

  // Show extracted agent decisions
  console.log(`     - Agent Decisions:`);
  if (decision.agent_responses.gpt5_2?.extracted_decision) {
    console.log(`       * GPT-5.2: ${decision.agent_responses.gpt5_2.extracted_decision} (${decision.agent_responses.gpt5_2.confidence_label})`);
  }
  if (decision.agent_responses.claude?.extracted_decision) {
    console.log(`       * Claude:  ${decision.agent_responses.claude.extracted_decision} (${decision.agent_responses.claude.confidence_label})`);
  }
  if (decision.agent_responses.gemini?.extracted_decision) {
    console.log(`       * Gemini:  ${decision.agent_responses.gemini.extracted_decision} (${decision.agent_responses.gemini.confidence_label})`);
  }
});

console.log('\n✅ JSON Logging Test Complete!');
console.log('\nYou can now inspect the files in:', testOutputDir);
console.log('\n📖 Files to check:');
console.log('   - test_analysis_log.txt (plain text log)');
console.log('   - test_analysis_agent_conversations.txt (plain text conversations)');
console.log('   - test_analysis_agent_decisions.json (JSON - all decisions)');
console.log('   - test_analysis_agent_decisions.jsonl (JSONL - streaming format)');
console.log('   - test_analysis_session_metadata.json (JSON - session info)');
