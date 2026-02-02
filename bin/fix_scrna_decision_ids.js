#!/usr/bin/env node

/**
 * Fix "unknown" decision_id in scRNA JSONL files
 *
 * Changes:
 *   FROM: "decision_id": "unknown_stage1_load_+_qc"
 *   TO:   "decision_id": "1_GD428_21136_Hu_REH_Parental_stage1_load_+_qc"
 *
 * Usage:
 *   node bin/fix_scrna_decision_ids.js experiments/scrna_results
 */

import fs from 'fs';
import path from 'path';
import { globSync } from 'glob';

const RESULTS_DIR = process.argv[2] || 'experiments/scrna_results';

console.log('\n' + '='.repeat(80));
console.log('FIX SCRNA DECISION IDs - Replace "unknown" with dataset name');
console.log('='.repeat(80));
console.log('');

if (!fs.existsSync(RESULTS_DIR)) {
  console.error(`ERROR: Directory not found: ${RESULTS_DIR}`);
  process.exit(1);
}

// Find all JSONL files
const jsonlFiles = globSync(path.join(RESULTS_DIR, '**/*_agent_decisions.jsonl'));

console.log(`Found ${jsonlFiles.length} JSONL files to fix\n`);

let totalFixed = 0;
let filesProcessed = 0;

for (const jsonlFile of jsonlFiles) {
  const folderName = path.basename(path.dirname(jsonlFile));

  // Extract dataset name from folder name
  // Example: "1_GD428_21136_Hu_REH_Parental_single_gpt" -> "1_GD428_21136_Hu_REH_Parental"
  // Pattern: remove "_parallel_*", "_sequential_*", "_single_*", "_no_agent"
  let datasetId = folderName
    .replace(/_parallel_.*$/, '')
    .replace(/_sequential_.*$/, '')
    .replace(/_single_.*$/, '')
    .replace(/_no_agent$/, '');

  console.log(`Processing: ${folderName}`);
  console.log(`  Dataset ID: ${datasetId}`);

  try {
    // Read JSONL file
    const content = fs.readFileSync(jsonlFile, 'utf-8');
    const lines = content.trim().split('\n').filter(line => line.trim());

    let fixedCount = 0;
    const fixedLines = [];

    for (const line of lines) {
      const decision = JSON.parse(line);

      // Fix decision_id if it starts with "unknown"
      if (decision.decision_id && decision.decision_id.startsWith('unknown')) {
        const oldId = decision.decision_id;
        decision.decision_id = oldId.replace('unknown', datasetId);
        fixedCount++;
        console.log(`    Fixed: ${oldId} -> ${decision.decision_id}`);
      }

      fixedLines.push(JSON.stringify(decision));
    }

    if (fixedCount > 0) {
      // Write back to file
      fs.writeFileSync(jsonlFile, fixedLines.join('\n') + '\n');
      console.log(`  ✓ Fixed ${fixedCount} decision_ids\n`);
      totalFixed += fixedCount;
      filesProcessed++;
    } else {
      console.log(`  - No fixes needed\n`);
    }

  } catch (error) {
    console.error(`  ✗ Error: ${error.message}\n`);
  }
}

console.log('='.repeat(80));
console.log(`✨ Complete! Fixed ${totalFixed} decision_ids in ${filesProcessed} files`);
console.log('='.repeat(80));
console.log('');
