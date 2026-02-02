#!/usr/bin/env node

/**
 * Add decision_id to scRNA ground truth JSON
 *
 * Adds decision_id field to each stage in scrna_ground_truth.json
 * to match the structure of bulk_rna_ground_truth.json
 *
 * Usage:
 *   node bin/add_decision_ids_to_scrna_ground_truth.js
 */

import fs from 'fs';
import path from 'path';

const GROUND_TRUTH_FILE = 'experiments/scrna_ground_truth.json';

console.log('\n' + '='.repeat(80));
console.log('ADD DECISION IDs TO SCRNA GROUND TRUTH');
console.log('='.repeat(80));
console.log('');

if (!fs.existsSync(GROUND_TRUTH_FILE)) {
  console.error(`ERROR: File not found: ${GROUND_TRUTH_FILE}`);
  process.exit(1);
}

// Read ground truth
const groundTruth = JSON.parse(fs.readFileSync(GROUND_TRUTH_FILE, 'utf-8'));

// Stage name to suffix mapping (for decision_id)
const stageSuffixes = {
  '1': 'load_+_qc',
  '2': 'qc_filtering',
  '3A': 'cell_cycle_scoring',
  '3B': 'execute_cell_cycle_decision',
  '4': 'pca_+_pc_selection',
  '5': 'clustering_+_markers'
};

let totalAdded = 0;

// Process each dataset
for (const [datasetKey, dataset] of Object.entries(groundTruth.ground_truth)) {
  const datasetId = dataset.dataset_id;

  console.log(`Processing: ${datasetId}`);

  if (!dataset.stages) {
    console.log(`  ⚠ No stages found, skipping\n`);
    continue;
  }

  // Add decision_id to each stage
  for (const [stageKey, stage] of Object.entries(dataset.stages)) {
    const suffix = stageSuffixes[stageKey];

    if (!suffix) {
      console.log(`  ⚠ Unknown stage key: ${stageKey}`);
      continue;
    }

    const decisionId = `${datasetId}_stage${stageKey}_${suffix}`;

    // Add decision_id as the first field after stage_name
    const updatedStage = {
      stage_name: stage.stage_name,
      decision_id: decisionId,
      ...stage
    };

    // Remove old stage_name (since we added it at the beginning)
    delete updatedStage.stage_name;
    updatedStage.stage_name = stage.stage_name;

    // Reorder: stage_name, decision_id, then rest
    const orderedStage = {
      stage_name: stage.stage_name,
      decision_id: decisionId
    };

    for (const key of Object.keys(stage)) {
      if (key !== 'stage_name') {
        orderedStage[key] = stage[key];
      }
    }

    dataset.stages[stageKey] = orderedStage;

    console.log(`  ✓ Stage ${stageKey}: ${decisionId}`);
    totalAdded++;
  }

  console.log('');
}

// Update metadata
groundTruth.updated = new Date().toISOString().split('T')[0];
groundTruth.total_decisions = totalAdded;

// Write back to file
fs.writeFileSync(GROUND_TRUTH_FILE, JSON.stringify(groundTruth, null, 2));

console.log('='.repeat(80));
console.log(`✨ Added ${totalAdded} decision_ids to ${GROUND_TRUTH_FILE}`);
console.log('='.repeat(80));
console.log('');
