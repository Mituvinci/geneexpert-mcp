/**
 * Data Detector - Auto-detect input data type
 *
 * Detects:
 * - FASTQ files (paired-end or single-end)
 * - BAM files (aligned reads)
 * - Count matrix (CSV/TSV/TXT)
 */

import fs from 'fs';
import path from 'path';

/**
 * Detect input data type and extract sample information
 * @param {string} inputDir - Input directory path
 * @param {Object} config - Configuration with group keywords (optional)
 * @returns {Object} - { type, samples, pairedEnd, files }
 */
export function detectInputData(inputDir, config = {}) {
  console.log('[Data Detector] Analyzing input directory...');

  if (!fs.existsSync(inputDir)) {
    throw new Error(`Input directory not found: ${inputDir}`);
  }

  const stat = fs.statSync(inputDir);

  // If it's a file (count matrix), not a directory
  if (stat.isFile()) {
    return detectCountMatrix(inputDir);
  }

  // If it's a directory, scan for files
  const files = fs.readdirSync(inputDir);

  // Check for FASTQ files
  const fastqFiles = files.filter(f =>
    f.endsWith('.fastq.gz') ||
    f.endsWith('.fq.gz') ||
    f.endsWith('.fastq') ||
    f.endsWith('.fq')
  );

  if (fastqFiles.length > 0) {
    return detectFASTQ(inputDir, fastqFiles, config);
  }

  // Check for BAM files
  const bamFiles = files.filter(f => f.endsWith('.bam'));

  if (bamFiles.length > 0) {
    return detectBAM(inputDir, bamFiles);
  }

  // Check for count matrix
  const matrixFiles = files.filter(f =>
    f.endsWith('.csv') ||
    f.endsWith('.tsv') ||
    f.endsWith('.txt') ||
    f.endsWith('.counts')
  );

  if (matrixFiles.length > 0) {
    return detectCountMatrix(path.join(inputDir, matrixFiles[0]));
  }

  throw new Error('No recognized data files found. Expected: FASTQ, BAM, or count matrix.');
}

/**
 * Detect FASTQ files and determine if paired-end
 */
function detectFASTQ(inputDir, fastqFiles, config = {}) {
  console.log(`[Data Detector] Found ${fastqFiles.length} FASTQ files`);

  // Check for paired-end pattern (R1/R2 or _1/_2)
  const r1Files = fastqFiles.filter(f => f.includes('_R1_') || f.includes('_1.f'));
  const r2Files = fastqFiles.filter(f => f.includes('_R2_') || f.includes('_2.f'));

  const pairedEnd = r1Files.length > 0 && r2Files.length > 0;

  // Extract sample names
  let samples = [];

  if (pairedEnd) {
    samples = r1Files.map(f => {
      // Extract sample name (remove _R1_ and extensions)
      const baseName = f
        .replace(/_R1_.*/, '')
        .replace(/_1\.f.*/, '')
        .replace(/\.fastq.*/, '');

      const r1Path = path.join(inputDir, f);
      const r2Name = f.replace(/_R1_/, '_R2_').replace(/_1\.f/, '_2.f');
      const r2Path = path.join(inputDir, r2Name);

      return {
        name: baseName,
        r1: r1Path,
        r2: r2Files.includes(r2Name) ? r2Path : null
      };
    });
  } else {
    samples = fastqFiles.map(f => ({
      name: f.replace(/\.fastq.*/, '').replace(/\.fq.*/, ''),
      r1: path.join(inputDir, f),
      r2: null
    }));
  }

  console.log(`[Data Detector] Detected ${samples.length} samples (${pairedEnd ? 'paired-end' : 'single-end'})`);

  samples.forEach((s, i) => {
    console.log(`   ${i + 1}. ${s.name}`);
  });

  // Auto-detect experimental groups from sample names (using user-provided keywords if available)
  const groups = detectGroups(samples, config);

  return {
    type: 'fastq',
    pairedEnd,
    samples,
    files: fastqFiles.map(f => path.join(inputDir, f)),
    groups  // Add detected groups
  };
}

/**
 * Auto-detect experimental groups from sample names
 * Uses user-provided keywords if available, otherwise falls back to common patterns
 */
function detectGroups(samples, config = {}) {
  const groups = {};

  // Get user-provided keywords (convert to lowercase for matching)
  const controlKeyword = config.controlKeyword?.toLowerCase();
  const treatmentKeyword = config.treatmentKeyword?.toLowerCase();

  samples.forEach(sample => {
    const name = sample.name.toLowerCase();

    // If user provided keywords, use those FIRST
    if (controlKeyword && name.includes(controlKeyword)) {
      if (!groups.control) groups.control = [];
      groups.control.push(sample.name);
    } else if (treatmentKeyword && name.includes(treatmentKeyword)) {
      if (!groups.treatment) groups.treatment = [];
      groups.treatment.push(sample.name);
    }
    // Fallback to common patterns if no keywords provided or no match
    else if (name.includes('control') || name.includes('ctrl') || name.includes('wt')) {
      if (!groups.control) groups.control = [];
      groups.control.push(sample.name);
    } else if (name.includes('ko') || name.includes('knockout')) {
      if (!groups.knockout) groups.knockout = [];
      groups.knockout.push(sample.name);
    } else if (name.includes('treat') || name.includes('drug')) {
      if (!groups.treatment) groups.treatment = [];
      groups.treatment.push(sample.name);
    } else {
      // Unknown group
      if (!groups.unknown) groups.unknown = [];
      groups.unknown.push(sample.name);
    }
  });

  // Log what keywords were used
  if (controlKeyword || treatmentKeyword) {
    console.log(`[Data Detector] Using user-provided keywords:`);
    if (controlKeyword) console.log(`   Control: "${controlKeyword}"`);
    if (treatmentKeyword) console.log(`   Treatment: "${treatmentKeyword}"`);
  }

  return groups;
}

/**
 * Detect BAM files
 */
function detectBAM(inputDir, bamFiles) {
  console.log(`[Data Detector] Found ${bamFiles.length} BAM files`);

  const samples = bamFiles.map(f => ({
    name: f.replace('.bam', ''),
    bam: path.join(inputDir, f),
    bai: fs.existsSync(path.join(inputDir, f + '.bai'))
      ? path.join(inputDir, f + '.bai')
      : null
  }));

  console.log(`[Data Detector] Detected ${samples.length} aligned samples`);

  samples.forEach((s, i) => {
    console.log(`   ${i + 1}. ${s.name} ${s.bai ? '(indexed)' : '(no index)'}`);
  });

  return {
    type: 'bam',
    pairedEnd: null,
    samples,
    files: bamFiles.map(f => path.join(inputDir, f))
  };
}

/**
 * Detect count matrix
 */
function detectCountMatrix(filePath) {
  console.log(`[Data Detector] Found count matrix: ${path.basename(filePath)}`);

  // Read first few lines to determine format
  const content = fs.readFileSync(filePath, 'utf-8');
  const lines = content.split('\n').slice(0, 5);

  // Detect delimiter (comma, tab, or space)
  const firstLine = lines[0];
  let delimiter = '\t';
  if (firstLine.includes(',')) delimiter = ',';
  else if (firstLine.includes('\t')) delimiter = '\t';
  else if (firstLine.includes(' ')) delimiter = ' ';

  // Extract sample names from header
  const header = firstLine.split(delimiter);
  const sampleNames = header.slice(1); // First column is gene names

  console.log(`[Data Detector] Detected count matrix with ${sampleNames.length} samples`);
  console.log(`   Samples: ${sampleNames.join(', ')}`);

  return {
    type: 'counts',
    pairedEnd: null,
    samples: sampleNames.map(name => ({ name, counts: filePath })),
    files: [filePath],
    format: {
      delimiter,
      hasHeader: true
    }
  };
}
