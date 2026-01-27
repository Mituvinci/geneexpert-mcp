#!/usr/bin/env node
/**
 * Generate Missing dataset_metadata.json Files
 *
 * This script creates dataset_metadata.json files for preprocessing folders
 * that are missing them. It detects paired-end/single-end from FASTQ files.
 *
 * Usage:
 *   # Generate for a single preprocessing folder
 *   node bin/generate_metadata.js \
 *     --preprocessing experiments/preprocessing/6_GSE193658_Lab_data \
 *     --data data/download_new_bulk_RNA_Data/6_GSE193658_Lab_data \
 *     --organism human
 *
 *   # Batch process all preprocessing folders
 *   node bin/generate_metadata.js --batch
 */

import path from 'path';
import fs from 'fs';
import { fileURLToPath } from 'url';

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);
const projectRoot = path.resolve(__dirname, '..');

// Parse CLI arguments
const args = process.argv.slice(2);

function printUsage() {
  console.log('');
  console.log('Usage:');
  console.log('  Single folder:');
  console.log('    node bin/generate_metadata.js \\');
  console.log('      --preprocessing <preprocessing_dir> \\');
  console.log('      --data <original_data_dir> \\');
  console.log('      --organism <organism>');
  console.log('');
  console.log('  Batch process all:');
  console.log('    node bin/generate_metadata.js --batch');
  console.log('');
  console.log('Example:');
  console.log('  node bin/generate_metadata.js \\');
  console.log('    --preprocessing experiments/preprocessing/6_GSE193658_Lab_data \\');
  console.log('    --data data/download_new_bulk_RNA_Data/6_GSE193658_Lab_data \\');
  console.log('    --organism human');
  console.log('');
  process.exit(1);
}

// Detect paired-end/single-end from validation report (fallback)
function detectFromValidationReport(preprocessingDir) {
  const validationPath = path.join(preprocessingDir, 'stage1_validation/validation_report.tsv');

  if (!fs.existsSync(validationPath)) {
    return null;
  }

  try {
    const content = fs.readFileSync(validationPath, 'utf-8');
    const lines = content.trim().split('\n');

    if (lines.length < 2) {
      return null;
    }

    // Skip header, get file names
    const files = lines.slice(1).map(line => line.split('\t')[0]);

    if (files.length === 0) {
      return null;
    }

    const pairedEnd = files.some(f => f.includes('_R2_') || f.includes('_2.fastq.gz'));

    // Extract sample names
    const sampleMap = {};
    for (const file of files) {
      const sampleName = file
        .replace(/_R[12]_001\.fastq\.gz$/, '')
        .replace(/_[12]\.fastq\.gz$/, '')
        .replace(/\.fastq\.gz$/, '');

      if (!sampleMap[sampleName]) {
        sampleMap[sampleName] = true;
      }
    }

    const samples = Object.keys(sampleMap);

    return {
      pairedEnd,
      sequencingType: pairedEnd ? 'paired-end' : 'single-end',
      numSamples: samples.length,
      samples,
      source: 'validation_report'
    };
  } catch (error) {
    console.error(`Warning: Could not parse validation report: ${error.message}`);
    return null;
  }
}

// Detect paired-end/single-end from FASTQ files
function detectSequencingType(dataDir, preprocessingDir) {
  // First try: Look at actual FASTQ files in data directory
  if (fs.existsSync(dataDir)) {
    const files = fs.readdirSync(dataDir).filter(f => f.endsWith('.fastq.gz'));

    if (files.length > 0) {
      const pairedEnd = files.some(f => f.includes('_R2_') || f.includes('_2.fastq.gz'));

      // Extract sample names
      const sampleMap = {};
      for (const file of files) {
        const sampleName = file
          .replace(/_R[12]_001\.fastq\.gz$/, '')
          .replace(/_[12]\.fastq\.gz$/, '')
          .replace(/\.fastq\.gz$/, '');

        if (!sampleMap[sampleName]) {
          sampleMap[sampleName] = true;
        }
      }

      const samples = Object.keys(sampleMap);

      return {
        pairedEnd,
        sequencingType: pairedEnd ? 'paired-end' : 'single-end',
        numSamples: samples.length,
        samples,
        source: 'fastq_files'
      };
    }
  }

  // Fallback: Try to read from validation report
  console.log('  FASTQ files not found, checking validation report...');
  const fromReport = detectFromValidationReport(preprocessingDir);
  if (fromReport) {
    return fromReport;
  }

  console.error(`ERROR: No FASTQ files found in ${dataDir} and no validation report available`);
  return null;
}

// Infer organism from dataset name or ask user
function inferOrganism(datasetName, providedOrganism) {
  if (providedOrganism) {
    return providedOrganism;
  }

  const nameLower = datasetName.toLowerCase();

  if (nameLower.includes('mouse') || nameLower.includes('mm10') || nameLower.includes('mm39')) {
    return 'mouse';
  }
  if (nameLower.includes('human') || nameLower.includes('hg38') || nameLower.includes('hg19')) {
    return 'human';
  }
  if (nameLower.includes('rat') || nameLower.includes('rn6') || nameLower.includes('rn7')) {
    return 'rat';
  }

  // Default to human if can't determine
  console.warn(`Warning: Could not infer organism from dataset name. Defaulting to 'human'.`);
  return 'human';
}

// Get genome build
function getGenomeBuild(organism) {
  const builds = {
    'human': 'hg38',
    'mouse': 'mm10',
    'rat': 'rn6'
  };
  return builds[organism] || 'hg38';
}

// Generate metadata for a single preprocessing folder
function generateMetadata(preprocessingDir, dataDir, organism) {
  console.log('');
  console.log('='.repeat(70));
  console.log('  GENERATING METADATA');
  console.log('='.repeat(70));
  console.log('');
  console.log(`Preprocessing: ${preprocessingDir}`);
  console.log(`Data:          ${dataDir}`);
  console.log(`Organism:      ${organism}`);
  console.log('');

  // Check if metadata already exists
  const metadataPath = path.join(preprocessingDir, 'dataset_metadata.json');
  if (fs.existsSync(metadataPath)) {
    console.log(`⚠️  Metadata already exists: ${metadataPath}`);
    console.log('   Skipping (use --force to overwrite)');
    return false;
  }

  // Detect sequencing type
  const detectedInfo = detectSequencingType(dataDir, preprocessingDir);
  if (!detectedInfo) {
    console.error('ERROR: Could not detect sequencing type');
    return false;
  }

  console.log(`Detected: ${detectedInfo.numSamples} samples (${detectedInfo.sequencingType})`);
  if (detectedInfo.source) {
    console.log(`Source:   ${detectedInfo.source}`);
  }
  console.log('');

  // Create metadata object
  const datasetName = path.basename(preprocessingDir);
  const metadata = {
    dataset_name: datasetName,
    organism: organism,
    genome_build: getGenomeBuild(organism),
    sequencing_type: detectedInfo.sequencingType,
    paired_end: detectedInfo.pairedEnd,
    num_samples: detectedInfo.numSamples,
    preprocessing_date: new Date().toISOString(),
    samples: detectedInfo.samples,
    notes: 'Generated automatically by generate_metadata.js'
  };

  // Save metadata
  fs.writeFileSync(metadataPath, JSON.stringify(metadata, null, 2));
  console.log(`✓ Metadata saved: ${metadataPath}`);
  console.log('');
  console.log('Contents:');
  console.log(JSON.stringify(metadata, null, 2));
  console.log('');

  return true;
}

// Batch process all preprocessing folders
function batchProcess() {
  console.log('');
  console.log('='.repeat(70));
  console.log('  BATCH PROCESSING: Generating Missing Metadata Files');
  console.log('='.repeat(70));
  console.log('');

  const preprocessingRoot = path.join(projectRoot, 'experiments/preprocessing');

  if (!fs.existsSync(preprocessingRoot)) {
    console.error(`ERROR: Preprocessing directory not found: ${preprocessingRoot}`);
    process.exit(1);
  }

  const folders = fs.readdirSync(preprocessingRoot).filter(f => {
    const fullPath = path.join(preprocessingRoot, f);
    return fs.statSync(fullPath).isDirectory();
  });

  console.log(`Found ${folders.length} preprocessing folders`);
  console.log('');

  let processed = 0;
  let skipped = 0;
  let errors = 0;

  for (const folder of folders) {
    const preprocessingDir = path.join(preprocessingRoot, folder);
    const metadataPath = path.join(preprocessingDir, 'dataset_metadata.json');

    // Skip if metadata exists
    if (fs.existsSync(metadataPath)) {
      console.log(`[${folder}] Metadata exists, skipping`);
      skipped++;
      continue;
    }

    // Try to find corresponding data directory
    const possibleDataDirs = [
      path.join(projectRoot, 'data', folder),
      path.join(projectRoot, 'data/download_new_bulk_RNA_Data', folder)
    ];

    let dataDir = null;
    for (const dir of possibleDataDirs) {
      if (fs.existsSync(dir)) {
        dataDir = dir;
        break;
      }
    }

    if (!dataDir) {
      console.log(`[${folder}] ⚠️  Could not find data directory, skipping`);
      errors++;
      continue;
    }

    // Infer organism from folder name
    const organism = inferOrganism(folder, null);

    // Generate metadata
    console.log(`[${folder}] Processing...`);
    const success = generateMetadata(preprocessingDir, dataDir, organism);

    if (success) {
      processed++;
    } else {
      errors++;
    }
  }

  console.log('');
  console.log('='.repeat(70));
  console.log('  BATCH PROCESSING COMPLETE');
  console.log('='.repeat(70));
  console.log('');
  console.log(`Processed: ${processed}`);
  console.log(`Skipped:   ${skipped}`);
  console.log(`Errors:    ${errors}`);
  console.log('');
}

// Main
function main() {
  if (args.length === 0) {
    printUsage();
  }

  // Batch mode
  if (args.includes('--batch')) {
    batchProcess();
    return;
  }

  // Single folder mode
  const preprocessingIndex = args.indexOf('--preprocessing');
  const dataIndex = args.indexOf('--data');
  const organismIndex = args.indexOf('--organism');

  if (preprocessingIndex === -1 || dataIndex === -1) {
    console.error('ERROR: --preprocessing and --data are required');
    printUsage();
  }

  const preprocessingDir = path.resolve(args[preprocessingIndex + 1]);
  const dataDir = path.resolve(args[dataIndex + 1]);
  const organism = organismIndex !== -1 ? args[organismIndex + 1] : null;

  if (!preprocessingDir || !dataDir) {
    console.error('ERROR: Invalid arguments');
    printUsage();
  }

  if (!fs.existsSync(preprocessingDir)) {
    console.error(`ERROR: Preprocessing directory not found: ${preprocessingDir}`);
    process.exit(1);
  }

  const inferredOrganism = inferOrganism(path.basename(preprocessingDir), organism);
  const success = generateMetadata(preprocessingDir, dataDir, inferredOrganism);

  if (!success) {
    process.exit(1);
  }
}

main();
