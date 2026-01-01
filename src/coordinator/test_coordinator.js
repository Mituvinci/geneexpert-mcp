#!/usr/bin/env node

/**
 * Test Coordinator Agent
 * Verifies multi-agent collaboration and consensus mechanisms
 */

import { createCoordinator } from './orchestrator.js';
import dotenv from 'dotenv';

dotenv.config();

console.log('🧪 Testing Multi-Agent Coordinator\n');
console.log('='.repeat(60));

async function runTests() {
  const coordinator = createCoordinator({ verbose: true });

  // Test 1: Simple threshold validation
  console.log('\n\n📊 TEST 1: Validate FDR Threshold');
  console.log('-'.repeat(60));

  try {
    const result1 = await coordinator.validateThreshold(
      'FDR',
      0.05,
      {
        sampleSize: 6,
        comparison: '3 vs 3 replicates',
        organism: 'mouse',
        experimentType: 'RNA-seq'
      }
    );

    console.log(result1.report);

    console.log('\n📈 Individual Agent Responses (FULL THINKING):');
    console.log('\n' + '='.repeat(60));
    console.log('Stats Agent (GPT-4):');
    console.log('='.repeat(60));
    console.log(result1.responses.gpt4.content || 'No response');
    console.log('\n' + '='.repeat(60));
    console.log('Pipeline Agent (Claude):');
    console.log('='.repeat(60));
    console.log(result1.responses.claude.content || 'No response');
    console.log('\n' + '='.repeat(60));
    console.log('Biology Agent (Gemini):');
    console.log('='.repeat(60));
    console.log(result1.responses.gemini.content || 'No response');

  } catch (error) {
    console.error('❌ Test 1 failed:', error.message);
  }

  // Test 2: Sample removal decision (should require unanimous)
  console.log('\n\n🔍 TEST 2: Sample Removal Decision (Requires Unanimous Approval)');
  console.log('-'.repeat(60));

  try {
    const result2 = await coordinator.decideSampleRemoval(
      'Sample_WT_3',
      'PCA shows this sample is 3 standard deviations from the group centroid',
      {
        pcDistance: 3.2,
        librarySize: 15000000,
        mappingRate: 0.89,
        rRNA_contamination: 0.02
      }
    );

    console.log(result2.report);

    if (result2.consensus.decision !== 'approve') {
      console.log('\n⚠️  Sample removal NOT unanimous - user decision required');
    } else {
      console.log('\n✅ All agents agree on sample removal');
    }

  } catch (error) {
    console.error('❌ Test 2 failed:', error.message);
  }

  // Test 3: QC plot review
  console.log('\n\n📈 TEST 3: QC Plot Review');
  console.log('-'.repeat(60));

  try {
    const result3 = await coordinator.reviewQCPlots(
      {
        pca_plot: 'pca_plot.pdf',
        mds_plot: 'mds_plot.pdf',
        density_plot: 'density_plot.pdf',
        observations: [
          'Samples cluster by treatment group',
          'PC1 explains 45% variance',
          'One sample (WT_3) appears as outlier'
        ]
      },
      {
        nSamples: 6,
        groups: ['WT_Control', 'WT_Treatment', 'KO_Control', 'KO_Treatment'],
        replicates: 3
      }
    );

    console.log(result3.report);

    if (result3.actionRequired) {
      console.log('\n⚠️  Action required before proceeding');
    } else {
      console.log('\n✅ OK to proceed to DE analysis');
    }

  } catch (error) {
    console.error('❌ Test 3 failed:', error.message);
  }

  // Test 4: Biological interpretation
  console.log('\n\n🧬 TEST 4: Biological Interpretation (Biology Agent Only)');
  console.log('-'.repeat(60));

  try {
    const result4 = await coordinator.interpretBiology(
      ['Tnf', 'Il6', 'Il1b', 'Cxcl10', 'Ccl2'],
      {
        experiment: 'Macrophage response to E. coli',
        timepoint: '6 hours',
        foldChange: 'All upregulated 5-10 fold',
        fdr: '< 0.001'
      }
    );

    console.log('\n🧬 Biology Agent Response:');
    console.log(result4.content);

  } catch (error) {
    console.error('❌ Test 4 failed:', error.message);
  }

  // Session summary
  console.log('\n\n📊 SESSION SUMMARY');
  console.log('='.repeat(60));

  const summary = coordinator.getSessionSummary();
  console.log(`Total queries: ${summary.totalQueries}`);
  console.log(`Average confidence: ${(summary.averageConfidence * 100).toFixed(1)}%`);
  console.log(`\nDecisions:`);
  console.log(`  Approved: ${summary.decisionCounts.approve}`);
  console.log(`  Rejected: ${summary.decisionCounts.reject}`);
  console.log(`  User decision required: ${summary.decisionCounts.user_decision_required}`);
  console.log(`  User approval required: ${summary.decisionCounts.user_approval_required}`);

  // Export session
  const sessionData = coordinator.exportSession();
  console.log(`\n✅ Session exported successfully`);
  console.log(`   Timestamp: ${sessionData.timestamp}`);

  console.log('\n' + '='.repeat(60));
  console.log('✅ All coordinator tests complete!\n');
}

// Run tests
runTests().catch(error => {
  console.error('\n💥 Test suite failed:', error);
  console.error(error.stack);
  process.exit(1);
});
