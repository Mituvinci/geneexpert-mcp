#!/usr/bin/env node

/**
 * Re-extract Decision Fields from Existing JSON Files
 *
 * Fixes missing extracted_decision, confidence_label, confidence_score
 * that were null due to Claude's markdown bold formatting (**Field:** value)
 *
 * Usage:
 *   node bin/reextract_decisions.js experiments/results
 */

import fs from 'fs';
import path from 'path';
import { parseAgentResponse } from '../src/utils/response_parser.js';

const RESULTS_DIR = process.argv[2] || 'experiments/results';

function reextractFile(jsonFile) {
  console.log(`Processing: ${path.basename(path.dirname(jsonFile))}`);

  let data;
  try {
    data = JSON.parse(fs.readFileSync(jsonFile, 'utf8'));
  } catch (e) {
    console.error(`  ✗ Error reading ${jsonFile}: ${e.message}`);
    return;
  }

  // Handle both array and single object formats
  const entries = Array.isArray(data) ? data : [data];
  let fixedCount = 0;

  for (const entry of entries) {
    if (!entry.agent_responses) continue;

    for (const [agentName, response] of Object.entries(entry.agent_responses)) {
      if (!response.content) continue;

      // Re-extract using fixed patterns
      const extracted = parseAgentResponse(response.content);

      // Update if previously null
      if (!response.extracted_decision && extracted.extracted_decision) {
        response.extracted_decision = extracted.extracted_decision;
        fixedCount++;
      }

      if (!response.confidence_label && extracted.confidence_label) {
        response.confidence_label = extracted.confidence_label;
        fixedCount++;
      }

      if (!response.confidence_score && extracted.confidence_score) {
        response.confidence_score = extracted.confidence_score;
        fixedCount++;
      }
    }
  }

  if (fixedCount > 0) {
    // Write back
    fs.writeFileSync(jsonFile, JSON.stringify(Array.isArray(data) ? entries : entries[0], null, 2));
    console.log(`  ✓ Fixed ${fixedCount} missing fields`);
  } else {
    console.log(`  - No missing fields found`);
  }
}

function main() {
  console.log('\n🔧 Re-extracting Decision Fields from JSON Files\n');
  console.log('=' .repeat(80));
  console.log();

  if (!fs.existsSync(RESULTS_DIR)) {
    console.error(`ERROR: Directory not found: ${RESULTS_DIR}`);
    process.exit(1);
  }

  let totalFixed = 0;
  let processedCount = 0;

  const folders = fs.readdirSync(RESULTS_DIR)
    .filter(f => fs.statSync(path.join(RESULTS_DIR, f)).isDirectory());

  for (const folder of folders) {
    const jsonFile = path.join(RESULTS_DIR, folder, 'staged_analysis_agent_decisions.json');

    if (!fs.existsSync(jsonFile)) continue;

    reextractFile(jsonFile);
    processedCount++;
  }

  console.log();
  console.log('=' .repeat(80));
  console.log(`✨ Processed ${processedCount} experiment folders`);
  console.log();
  console.log('📋 Next Steps:');
  console.log('  1. Regenerate CSV files:');
  console.log('     node bin/json_to_csv.js convert --dir experiments/results');
  console.log('     python bin/aggregate_experiments.py');
  console.log();
  console.log('  2. Run evaluation:');
  console.log('     node bin/evaluate_bulk_from_csv.js \\');
  console.log('       experiments/results/bulk_rna_ALL_EXPERIMENTS_DETAILED.csv \\');
  console.log('       experiments/bulk_rna_ground_truth.json \\');
  console.log('       bulk_evaluation_final.csv');
  console.log();
}

main();
