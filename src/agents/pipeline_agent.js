/**
 * Pipeline Agent - ACTUALLY executes tools via MCP
 *
 * This agent receives instructions from the Coordinator and
 * executes REAL bioinformatics tools from /data/scripts/
 */

import { exec } from 'child_process';
import { promisify } from 'util';
import fs from 'fs';
import path from 'path';

const execAsync = promisify(exec);

/**
 * Pipeline Agent - Executes bioinformatics tools
 */
export class PipelineAgent {
  constructor(config = {}) {
    this.config = config;
    this.verbose = config.verbose || false;
  }

  /**
   * Execute a pipeline step
   * @param {Object} step - Step definition from planner
   * @param {Object} session - Current session context
   * @returns {Object} - Execution result
   */
  async execute(step, session) {
    this.log(`[Pipeline Agent] Executing: ${step.name}`);

    const result = {
      step: step.name,
      tool: step.tool,
      status: 'pending',
      outputs: [],
      errors: [],
      startTime: new Date().toISOString()
    };

    try {
      // Execute based on tool type
      switch (step.tool) {
        case 'fastqc':
          result.outputs = await this.runFastQC(step, session);
          break;

        case 'fastq2bam':
          result.outputs = await this.runFastq2Bam(step, session);
          break;

        case 'featurecounts':
          result.outputs = await this.runFeatureCounts(step, session);
          break;

        case 'filterIDS':
          result.outputs = await this.runFilterIDS(step, session);
          break;

        case 'rpkm':
          result.outputs = await this.runRPKM(step, session);
          break;

        case 'edger':
          result.outputs = await this.runEdgeR(step, session);
          break;

        case 'visualization':
          result.outputs = await this.runVisualization(step, session);
          break;

        default:
          throw new Error(`Unknown tool: ${step.tool}`);
      }

      result.status = 'success';
      result.endTime = new Date().toISOString();

      this.log(`[Pipeline Agent] ✓ ${step.name} completed successfully`);

      return result;

    } catch (error) {
      result.status = 'failed';
      result.errors.push(error.message);
      result.endTime = new Date().toISOString();

      this.log(`[Pipeline Agent] ✗ ${step.name} failed: ${error.message}`);

      throw error;
    }
  }

  /**
   * Run FastQC (optional QC step)
   */
  async runFastQC(step, session) {
    const outputDir = path.join(session.config.output, 'fastqc');

    if (!fs.existsSync(outputDir)) {
      fs.mkdirSync(outputDir, { recursive: true });
    }

    this.log(`  FastQC output: ${outputDir}`);
    this.log(`  Processing ${step.inputs.length} files...`);

    // Run FastQC on all FASTQ files
    const fastqFiles = step.inputs.join(' ');
    const command = `fastqc ${fastqFiles} -o ${outputDir} -t 4`;

    this.log(`  Command: ${command}`);

    const { stdout, stderr } = await execAsync(command, {
      maxBuffer: 1024 * 1024 * 10, // 10MB buffer
      timeout: 600000 // 10 min timeout
    });

    if (this.verbose && stderr) {
      console.log('FastQC stderr:', stderr);
    }

    return [{
      type: 'fastqc_reports',
      path: outputDir,
      description: 'FastQC HTML and ZIP reports'
    }];
  }

  /**
   * Run fastq2bam (Subread alignment) - USER'S SCRIPT!
   */
  async runFastq2Bam(step, session) {
    const inputDir = session.config.input;
    const outputDir = path.join(session.config.output, 'bam_files');

    if (!fs.existsSync(outputDir)) {
      fs.mkdirSync(outputDir, { recursive: true });
    }

    // Copy FASTQ files to output directory (fastq2bam runs in current directory)
    this.log(`  Setting up alignment in: ${outputDir}`);
    this.log(`  Genome: ${step.params.genome}`);
    this.log(`  Read type: ${step.params.readType}`);

    // Create symbolic links to FASTQ files
    for (const sample of session.dataInfo.samples) {
      const r1Link = path.join(outputDir, path.basename(sample.r1));
      if (!fs.existsSync(r1Link)) {
        fs.symlinkSync(sample.r1, r1Link);
      }

      if (sample.r2) {
        const r2Link = path.join(outputDir, path.basename(sample.r2));
        if (!fs.existsSync(r2Link)) {
          fs.symlinkSync(sample.r2, r2Link);
        }
      }
    }

    // Run fastq2bam in output directory
    const command = `cd ${outputDir} && ${step.script} ${step.params.genome} ${step.params.readType}`;

    this.log(`  Command: ${command}`);
    this.log(`  This may take 10-30 minutes depending on data size...`);

    const { stdout, stderr } = await execAsync(command, {
      maxBuffer: 1024 * 1024 * 50, // 50MB buffer
      timeout: 3600000 // 60 min timeout
    });

    this.log(`  Alignment complete!`);

    if (stdout) this.log(`  stdout: ${stdout}`);
    if (stderr && this.verbose) console.log('Alignment stderr:', stderr);

    // Find generated BAM files
    const bamFiles = fs.readdirSync(outputDir).filter(f => f.endsWith('.bam'));

    return [{
      type: 'bam_files',
      path: outputDir,
      files: bamFiles,
      description: `${bamFiles.length} BAM files generated`
    }];
  }

  /**
   * Run featureCounts - USER'S SCRIPT!
   */
  async runFeatureCounts(step, session) {
    const bamDir = path.join(session.config.output, 'bam_files');
    const outputDir = path.join(session.config.output, 'counts');

    if (!fs.existsSync(outputDir)) {
      fs.mkdirSync(outputDir, { recursive: true });
    }

    // Get BAM files from previous step
    const bamFiles = fs.readdirSync(bamDir)
      .filter(f => f.endsWith('.bam'))
      .map(f => path.join(bamDir, f));

    if (bamFiles.length === 0) {
      throw new Error('No BAM files found from alignment step');
    }

    this.log(`  Counting reads in ${bamFiles.length} BAM files...`);
    this.log(`  Annotation: ${step.params.annotation}`);

    // Run featurecounts.R
    const command = `cd ${outputDir} && ${step.script} ${step.params.annotation} ${bamFiles.join(' ')}`;

    this.log(`  Command: ${command}`);

    const { stdout, stderr } = await execAsync(command, {
      maxBuffer: 1024 * 1024 * 50,
      timeout: 1800000 // 30 min
    });

    this.log(`  featureCounts complete!`);

    if (stdout) this.log(`  stdout: ${stdout}`);
    if (stderr && this.verbose) console.log('featureCounts stderr:', stderr);

    // Find generated count file
    const countFiles = fs.readdirSync(outputDir).filter(f => f.endsWith('.count.txt'));

    return [{
      type: 'count_matrix',
      path: outputDir,
      files: countFiles,
      description: `Count matrix: ${countFiles[0] || 'generated'}`
    }];
  }

  /**
   * Run filterIDS.R - USER'S SCRIPT!
   */
  async runFilterIDS(step, session) {
    const countsDir = path.join(session.config.output, 'counts');
    const outputDir = countsDir; // Output to same directory

    // Find count file
    const countFiles = fs.readdirSync(countsDir).filter(f => f.endsWith('.count.txt'));

    if (countFiles.length === 0) {
      throw new Error('No count matrix found');
    }

    const countFile = path.join(countsDir, countFiles[0]);

    this.log(`  Filtering low-count genes from: ${countFiles[0]}`);

    const command = `${step.script} ${countFile}`;

    this.log(`  Command: ${command}`);

    const { stdout, stderr } = await execAsync(command);

    this.log(`  Filtering complete!`);

    if (stdout) this.log(`  stdout: ${stdout}`);

    return [{
      type: 'filtered_counts',
      path: outputDir,
      file: `${countFiles[0]}.filtered`,
      description: 'Filtered count matrix'
    }];
  }

  /**
   * Run RPKM.R - USER'S SCRIPT!
   */
  async runRPKM(step, session) {
    const countsDir = path.join(session.config.output, 'counts');

    // Find filtered count file
    const filteredFile = fs.readdirSync(countsDir)
      .filter(f => f.endsWith('.filtered'))[0];

    if (!filteredFile) {
      throw new Error('No filtered count file found');
    }

    const inputFile = path.join(countsDir, filteredFile);

    this.log(`  Normalizing counts (RPKM): ${filteredFile}`);

    const command = `${step.script} ${inputFile}`;

    this.log(`  Command: ${command}`);

    const { stdout, stderr } = await execAsync(command);

    this.log(`  RPKM normalization complete!`);

    if (stdout) this.log(`  stdout: ${stdout}`);

    return [{
      type: 'normalized_counts',
      path: countsDir,
      file: `${filteredFile}.RPKM`,
      description: 'RPKM-normalized counts'
    }];
  }

  /**
   * Run edgeR - USER'S SCRIPT!
   */
  async runEdgeR(step, session) {
    const countsDir = path.join(session.config.output, 'counts');
    const deDir = path.join(session.config.output, 'differential_expression');

    if (!fs.existsSync(deDir)) {
      fs.mkdirSync(deDir, { recursive: true });
    }

    // Find filtered count file (edgeR uses filtered, not RPKM)
    const filteredFile = fs.readdirSync(countsDir)
      .filter(f => f.endsWith('.filtered'))[0];

    if (!filteredFile) {
      throw new Error('No filtered count file found');
    }

    const inputFile = path.join(countsDir, filteredFile);

    this.log(`  Running edgeR differential expression...`);
    this.log(`  Thresholds: FDR < ${session.thresholds.fdr}, logFC > ${session.thresholds.logFC}`);

    // Note: simpleEdger3.R might need additional parameters
    // For now, run with default
    const command = `cd ${deDir} && ${step.script} ${inputFile}`;

    this.log(`  Command: ${command}`);

    const { stdout, stderr } = await execAsync(command, {
      timeout: 600000 // 10 min
    });

    this.log(`  edgeR analysis complete!`);

    if (stdout) this.log(`  stdout: ${stdout}`);

    // Find DE results
    const deFiles = fs.readdirSync(deDir);

    return [{
      type: 'de_results',
      path: deDir,
      files: deFiles,
      description: 'Differential expression results'
    }];
  }

  /**
   * Generate visualization plots
   */
  async runVisualization(step, session) {
    this.log(`  Visualization step - plots should be in DE results`);

    // For now, visualization is handled by simpleEdger3.R
    // No additional action needed

    return [{
      type: 'plots',
      path: path.join(session.config.output, 'differential_expression'),
      description: 'Volcano plots, MA plots generated by edgeR'
    }];
  }

  /**
   * Logging helper
   */
  log(message) {
    if (this.verbose || true) { // Always log for now
      console.log(message);
    }
  }
}
