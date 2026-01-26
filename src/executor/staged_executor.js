/**
 * Staged Executor
 *
 * Orchestrates the 4-stage RNA-seq pipeline with agent checkpoints:
 *   Stage 1: FASTQ Validation
 *   Stage 2: Alignment + QC
 *   Stage 3: Quantification + PCA/QC
 *   Stage 4: DE Analysis
 *
 * Each stage: Run script → Parse output → Agents review → Decision → Next stage
 */

import path from 'path';
import fs from 'fs';
import { execSync } from 'child_process';

// Stage modules
import {
  generateStage1Script,
  parseStage1Output,
  formatStage1ForAgents
} from '../stages/stage1_fastq_validation.js';

import {
  generateStage2Script,
  parseStage2Output,
  formatStage2ForAgents
} from '../stages/stage2_alignment.js';

import {
  generateStage3Script,
  parseStage3Output,
  formatStage3ForAgents
} from '../stages/stage3_quantification_qc.js';

import {
  generateStage4Script,
  parseStage4Output,
  formatStage4ForAgents
} from '../stages/stage4_de_analysis.js';

// scRNA-seq stages
import {
  generateStage3aScript,
  parseStage3aOutput,
  formatStage3aForAgents
} from '../stages/stage3a_cell_cycle_scoring.js';

import {
  generateStage3bScript,
  parseStage3bOutput,
  formatStage3bForAgents
} from '../stages/stage3b_cell_cycle_regression.js';

import {
  generateStage3cScript,
  parseStage3cOutput,
  formatStage3cForAgents
} from '../stages/stage3c_clustering_umap.js';

// Coordinator and utilities
import { Coordinator } from '../coordinator/orchestrator.js';
import { Logger } from '../utils/logger.js';
import {
  handleStage1UserDecision,
  handleStage2UserDecision,
  handleStage3UserDecision,
  askYesNo
} from '../utils/user_input.js';

/**
 * Staged Executor Class
 */
export class StagedExecutor {
  constructor(options = {}) {
    this.inputDir = options.inputDir;
    this.outputDir = options.outputDir;
    this.organism = options.organism || 'mouse';
    this.comparison = options.comparison || path.basename(options.outputDir);
    this.controlKeyword = options.controlKeyword || 'cont';
    this.treatmentKeyword = options.treatmentKeyword || 'ips';
    this.verbose = options.verbose !== false;
    this.dryRun = options.dryRun || false;
    this.singleAgent = options.singleAgent || null;
    this.forceAutomation = options.forceAutomation || false;
    this.sequentialChain = options.sequentialChain || false;  // NEW: Sequential chain mode
    this.roleAssignments = options.roleAssignments || {  // NEW: Role assignments for ablation study
      gptRole: 'stats',
      claudeRole: 'pipeline',
      geminiRole: 'biology'
    };
    this.useExistingFastqc = options.useExistingFastqc || null;  // FIX: Use preprocessed FastQC
    this.useExistingBam = options.useExistingBam || null;  // FIX: Use preprocessed BAM
    this.dataType = options.dataType || 'bulk';  // NEW: 'bulk' or 'scrna'

    // State tracking
    this.stageResults = {
      stage1: null,
      stage2: null,
      stage3: null,
      stage3a: null,  // scRNA-seq: Cell cycle scoring
      stage3b: null,  // scRNA-seq: Cell cycle regression
      stage3c: null,  // scRNA-seq: Clustering + UMAP
      stage4: null
    };
    this.dataInfo = null;
    this.config = null;

    // Initialize coordinator
    this.coordinator = new Coordinator({
      verbose: this.verbose,
      singleAgent: this.singleAgent,
      sequentialChain: this.sequentialChain,  // NEW: Pass sequential chain flag
      roleAssignments: this.roleAssignments  // NEW: Pass role assignments for ablation study
    });

    // Logger will be initialized when we have output dir
    this.logger = null;
  }

  /**
   * Log message if verbose
   */
  log(message) {
    if (this.verbose) {
      console.log(`[Executor] ${message}`);
    }
  }

  /**
   * Detect data info from input directory
   */
  detectDataInfo() {
    this.log(`Detecting data from: ${this.inputDir}`);

    // FIRST: Check if we're using preprocessed BAM files (metadata might exist)
    let pairedEndFromMetadata = null;
    let experimentalDesignFromMetadata = null;
    if (this.useExistingBam) {
      // Metadata is in: preprocessing_dir/dataset_metadata.json
      // BAM dir is:     preprocessing_dir/stage2_alignment/bam_files/
      const metadataPath = path.join(path.dirname(path.dirname(this.useExistingBam)), 'dataset_metadata.json');
      if (fs.existsSync(metadataPath)) {
        try {
          const metadata = JSON.parse(fs.readFileSync(metadataPath, 'utf-8'));
          pairedEndFromMetadata = metadata.paired_end || metadata.pairedEnd;
          experimentalDesignFromMetadata = metadata.experimental_design || null;
          this.log(`✓ Loaded metadata from preprocessing: ${metadata.sequencing_type || (pairedEndFromMetadata ? 'paired-end' : 'single-end')}`);
          if (experimentalDesignFromMetadata && experimentalDesignFromMetadata.design) {
            this.log(`✓ Experimental design: ${experimentalDesignFromMetadata.design}`);
          }
        } catch (error) {
          this.log(`Warning: Could not read metadata file: ${error.message}`);
        }
      }
    }

    const files = fs.readdirSync(this.inputDir).filter(f => f.endsWith('.fastq.gz'));
    const samples = [];
    const pairedEnd = pairedEndFromMetadata !== null ? pairedEndFromMetadata : files.some(f => f.includes('_R2_'));

    // Group files by sample
    const sampleMap = {};
    for (const file of files) {
      const sampleName = file.replace(/_R[12]_001\.fastq\.gz$/, '').replace(/\.fastq\.gz$/, '');

      if (!sampleMap[sampleName]) {
        sampleMap[sampleName] = { name: sampleName };
      }

      if (file.includes('_R1_')) {
        sampleMap[sampleName].r1 = path.join(this.inputDir, file);
      } else if (file.includes('_R2_')) {
        sampleMap[sampleName].r2 = path.join(this.inputDir, file);
      } else {
        sampleMap[sampleName].r1 = path.join(this.inputDir, file);
      }
    }

    for (const [name, sample] of Object.entries(sampleMap)) {
      samples.push(sample);
    }

    // Detect groups
    const groups = {};
    for (const sample of samples) {
      const name = sample.name.toLowerCase();
      let group = 'unknown';

      if (name.includes(this.controlKeyword) || name.includes('ctrl') || name.includes('wt') || name.includes('control')) {
        group = 'control';
      } else if (name.includes(this.treatmentKeyword) || name.includes('ko') || name.includes('mut')) {
        group = 'treatment';
      }

      if (!groups[group]) {
        groups[group] = [];
      }
      groups[group].push(sample.name);
    }

    this.dataInfo = {
      type: 'fastq',
      samples,
      pairedEnd,
      groups,
      organism: this.organism,
      comparison: this.comparison,
      experimental_design: experimentalDesignFromMetadata  // NEW: Store experimental design from metadata
    };

    const seqTypeSource = pairedEndFromMetadata !== null ? ' (from metadata)' : ' (from FASTQ detection)';
    this.log(`Detected ${samples.length} samples (${pairedEnd ? 'paired-end' : 'single-end'}${seqTypeSource})`);
    this.log(`Groups: ${Object.keys(groups).map(g => `${g}(n=${groups[g].length})`).join(', ')}`);

    return this.dataInfo;
  }

  /**
   * Initialize the executor
   */
  initialize() {
    // Create output directory
    const absOutputDir = path.resolve(this.outputDir);
    if (!fs.existsSync(absOutputDir)) {
      fs.mkdirSync(absOutputDir, { recursive: true });
    }

    // Build config
    this.config = {
      input: path.resolve(this.inputDir),
      output: absOutputDir,
      organism: this.organism,
      comparison: this.comparison,
      controlKeyword: this.controlKeyword,
      treatmentKeyword: this.treatmentKeyword,
      singleAgent: this.singleAgent,
      forceAutomation: this.forceAutomation,
      sequentialChain: this.sequentialChain,  // CRITICAL: Pass sequential chain flag to logger
      roleAssignments: this.roleAssignments,  // CRITICAL BUG FIX: Pass role assignments to logger!
      useExistingFastqc: this.useExistingFastqc,  // FIX: Use preprocessed FastQC
      useExistingBam: this.useExistingBam  // FIX: Use preprocessed BAM
    };

    // Initialize logger
    this.logger = new Logger(absOutputDir, 'staged_analysis', this.config);

    // Set dataset for coordinator
    this.coordinator.setDataset(this.comparison);

    // Detect data
    this.detectDataInfo();

    return this;
  }

  /**
   * Execute a bash script with output streaming
   */
  executeScript(scriptPath, description, timeout = 3600000) {
    this.log(`Executing: ${description}`);

    try {
      execSync(`bash ${scriptPath}`, {
        stdio: 'inherit',
        cwd: path.dirname(scriptPath),
        timeout
      });
      return { success: true, exitCode: 0 };
    } catch (error) {
      // Check if it's just a non-zero exit code (WARN or FAIL from QC)
      if (error.status === 1 || error.status === 2) {
        this.log(`Script completed with QC status code: ${error.status}`);
        return { success: true, exitCode: error.status };
      }
      return { success: false, exitCode: error.status, error: error.message };
    }
  }

  /**
   * Run Stage 1: FASTQ Validation
   */
  async runStage1() {
    console.log('');
    console.log('='.repeat(60));
    console.log('  STAGE 1: FASTQ Validation');
    console.log('='.repeat(60));
    console.log('');

    const absOutputDir = path.resolve(this.outputDir);
    const scriptPath = path.resolve(absOutputDir, `stage_1_${this.comparison}.sh`);

    // Check if using existing FastQC outputs (preprocessing mode)
    if (this.config.useExistingFastqc) {
      console.log('📋 Using existing FastQC outputs (preprocessing mode)');
      console.log(`   Source: ${this.config.useExistingFastqc}`);
      console.log('');

      const sourceDir = this.config.useExistingFastqc;
      const targetDir = path.join(absOutputDir, 'stage1_fastq_validation');

      // Create target directory
      if (!fs.existsSync(targetDir)) {
        fs.mkdirSync(targetDir, { recursive: true });
      }

      // Copy FastQC outputs
      this.log('Copying FastQC outputs from preprocessing...');
      try {
        const { execSync } = await import('child_process');
        execSync(`cp -r "${sourceDir}"/* "${targetDir}"/`, { stdio: 'inherit' });
        this.log('FastQC outputs copied successfully');
      } catch (error) {
        console.error(`Failed to copy FastQC outputs: ${error.message}`);
        return { proceed: false, error: error.message };
      }
    } else {
      // Normal mode: generate and execute script
      // 1. Generate script
      this.log('Generating Stage 1 script...');
      generateStage1Script(this.dataInfo, this.config, scriptPath);
      this.log(`Script: ${scriptPath}`);

      if (this.dryRun) {
        this.log('DRY RUN - Script generated but not executed');
        return { proceed: true, dryRun: true };
      }

      // 2. Execute script
      this.log('Executing Stage 1 script...');
      const execResult = this.executeScript(scriptPath, 'FASTQ validation + FastQC');

      if (!execResult.success) {
        console.error(`Stage 1 execution failed: ${execResult.error}`);
        return { proceed: false, error: execResult.error };
      }
    }

    // 3. Parse output
    this.log('Parsing Stage 1 output...');
    const stage1Output = parseStage1Output(absOutputDir);
    this.log(`Status: ${stage1Output.overall_status}`);

    // 4. Skip agents if forceAutomation
    if (this.forceAutomation) {
      this.log('FORCE AUTOMATION - Skipping agent review');
      const proceed = stage1Output.overall_status !== 'FAIL';
      this.stageResults.stage1 = { output: stage1Output, proceed, automated: true };
      return { proceed, automated: true };
    }

    // 5. Agent review
    this.log('Calling agents for Stage 1 review...');
    const reviewResult = await this.coordinator.reviewStage1Output(stage1Output, this.dataInfo);

    // 6. Check if user decision needed
    let finalProceed = reviewResult.proceed;
    let userDecision = null;

    if (reviewResult.consensus.decision.toLowerCase().includes('user_decision')) {
      userDecision = await handleStage1UserDecision(reviewResult, stage1Output, this.config);
      finalProceed = userDecision.proceed;
      this.log(`User decided: ${finalProceed ? 'PROCEED' : 'ABORT'}`);
    }

    // 7. Log decision
    this.logger.logStageDecision(
      1, 'FASTQ Validation',
      this.dataInfo, stage1Output,
      reviewResult.responses, reviewResult.consensus,
      { proceed: finalProceed, userDecision }
    );

    // 8. Store result
    this.stageResults.stage1 = {
      output: stage1Output,
      review: reviewResult,
      proceed: finalProceed,
      userDecision
    };

    console.log('');
    console.log(`Stage 1 Result: ${finalProceed ? 'PASS' : 'FAIL'}`);

    return { proceed: finalProceed, output: stage1Output, review: reviewResult };
  }

  /**
   * Run Stage 2: Alignment + QC
   */
  async runStage2() {
    console.log('');
    console.log('='.repeat(60));
    console.log('  STAGE 2: Alignment + Alignment QC');
    console.log('='.repeat(60));
    console.log('');

    const absOutputDir = path.resolve(this.outputDir);
    const scriptPath = path.resolve(absOutputDir, `stage_2_${this.comparison}.sh`);

    // Check if using existing BAM files (preprocessing mode)
    if (this.config.useExistingBam) {
      console.log('📋 Using existing BAM files (preprocessing mode)');
      console.log(`   Source: ${this.config.useExistingBam}`);
      console.log('');

      const sourceDir = this.config.useExistingBam;
      const targetDir = path.join(absOutputDir, 'stage2_alignment', 'bam_files');
      const targetQcDir = path.join(absOutputDir, 'stage2_alignment', 'alignment_qc');

      // Create target directories
      if (!fs.existsSync(targetDir)) {
        fs.mkdirSync(targetDir, { recursive: true });
      }
      if (!fs.existsSync(targetQcDir)) {
        fs.mkdirSync(targetQcDir, { recursive: true });
      }

      // Copy all files from preprocessing (BAM files, count.csv, logs, etc.)
      this.log('Copying Stage 2 outputs from preprocessing...');
      try {
        const { execSync } = await import('child_process');

        // Copy everything from bam_files directory (BAMs, count.csv, logs, etc.)
        execSync(`cp -r "${sourceDir}"/* "${targetDir}"/`, { stdio: 'inherit' });

        // Copy alignment QC
        const sourceQcDir = path.join(path.dirname(sourceDir), 'alignment_qc');
        if (fs.existsSync(sourceQcDir)) {
          execSync(`cp -r "${sourceQcDir}"/* "${targetQcDir}"/`, { stdio: 'inherit' });
        }

        this.log('Stage 2 outputs copied successfully');
      } catch (error) {
        console.error(`Failed to copy Stage 2 outputs: ${error.message}`);
        return { proceed: false, error: error.message };
      }
    } else {
      // Normal mode: generate and execute script
      // 1. Generate script (pass stage1 decision for context)
      this.log('Generating Stage 2 script...');
      generateStage2Script(this.dataInfo, this.config, scriptPath, this.stageResults.stage1);
      this.log(`Script: ${scriptPath}`);

      if (this.dryRun) {
        this.log('DRY RUN - Script generated but not executed');
        return { proceed: true, dryRun: true };
      }

      // 2. Execute script (longer timeout for alignment)
      this.log('Executing Stage 2 script (this may take 10-30+ minutes)...');
      const execResult = this.executeScript(scriptPath, 'Alignment + QC', 7200000); // 2 hour timeout

      if (!execResult.success) {
        console.error(`Stage 2 execution failed: ${execResult.error}`);
        return { proceed: false, error: execResult.error };
      }
    }

    // 3. Parse output
    this.log('Parsing Stage 2 output...');
    const stage2Output = parseStage2Output(absOutputDir);
    this.log(`Status: ${stage2Output.overall_status}`);

    // 4. Skip agents if forceAutomation
    if (this.forceAutomation) {
      this.log('FORCE AUTOMATION - Skipping agent review');
      const proceed = stage2Output.overall_status !== 'FAIL' && stage2Output.overall_status !== 'SEVERE_FAIL';
      this.stageResults.stage2 = { output: stage2Output, proceed, automated: true, samplesToRemove: [] };
      return { proceed, automated: true };
    }

    // 5. Agent review
    this.log('Calling agents for Stage 2 review...');
    const reviewResult = await this.coordinator.reviewStage2Output(stage2Output, this.dataInfo);

    // 6. Check if user decision needed
    let finalProceed = reviewResult.proceed;
    let finalSamplesToRemove = reviewResult.samplesToRemove || [];
    let userDecision = null;

    if (reviewResult.consensus.decision.toLowerCase().includes('user_decision')) {
      userDecision = await handleStage2UserDecision(reviewResult, stage2Output, this.config);
      finalProceed = userDecision.proceed;
      finalSamplesToRemove = userDecision.samplesToRemove || [];
      this.log(`User decided: ${finalProceed ? 'PROCEED' : 'ABORT'}`);
    }

    // 7. Log decision
    this.logger.logStageDecision(
      2, 'Alignment QC',
      this.dataInfo, stage2Output,
      reviewResult.responses, reviewResult.consensus,
      { proceed: finalProceed, samplesToRemove: finalSamplesToRemove, userDecision }
    );

    // 8. Store result
    this.stageResults.stage2 = {
      output: stage2Output,
      review: reviewResult,
      proceed: finalProceed,
      samplesToRemove: finalSamplesToRemove,
      userDecision
    };

    // 9. Update dataInfo if samples removed
    if (finalSamplesToRemove.length > 0) {
      this.dataInfo.samples = this.dataInfo.samples.filter(
        s => !finalSamplesToRemove.includes(s.name)
      );
      this.log(`Removed ${finalSamplesToRemove.length} samples from analysis`);
    }

    console.log('');
    console.log(`Stage 2 Result: ${finalProceed ? 'PASS' : 'FAIL'}`);
    if (finalSamplesToRemove.length > 0) {
      console.log(`Samples removed: ${finalSamplesToRemove.join(', ')}`);
    }

    return { proceed: finalProceed, output: stage2Output, review: reviewResult, samplesToRemove: finalSamplesToRemove };
  }

  /**
   * Run Stage 3: Quantification + PCA/QC
   * (Placeholder - will be implemented with stage3 module)
   */
  async runStage3() {
    console.log('');
    console.log('='.repeat(60));
    console.log('  STAGE 3: Quantification + PCA/QC');
    console.log('='.repeat(60));
    console.log('');

    const absOutputDir = path.resolve(this.outputDir);
    const scriptPath = path.resolve(absOutputDir, `stage_3_${this.comparison}.sh`);

    // 1. Generate script (pass stage2 result for BAM file locations)
    this.log('Generating Stage 3 script...');
    generateStage3Script(this.dataInfo, this.config, scriptPath, this.stageResults.stage2);
    this.log(`Script: ${scriptPath}`);

    if (this.dryRun) {
      this.log('DRY RUN - Script generated but not executed');
      return { proceed: true, dryRun: true };
    }

    // 2. Execute script (longer timeout for quantification)
    this.log('Executing Stage 3 script (this may take 10-20 minutes)...');
    const execResult = this.executeScript(scriptPath, 'Quantification + QC', 3600000); // 1 hour timeout

    if (!execResult.success && execResult.exitCode !== 2) {
      // Exit code 2 is OK (just means issues detected in QC)
      console.error(`Stage 3 execution failed: ${execResult.error}`);
      return { proceed: false, error: execResult.error };
    }

    // 3. Parse output
    this.log('Parsing Stage 3 output...');
    const stage3Output = parseStage3Output(absOutputDir);
    this.log(`Status: ${stage3Output.overall_status}`);

    // Notify user about sample dictionary file (for interpreting PCA plot labels)
    if (stage3Output.sample_dictionary) {
      console.log('');
      console.log('📋 Sample Dictionary Created:');
      const dictPath = path.join(absOutputDir, 'stage3_quantification', `${this.comparison}_sample_dictionary.txt`);
      console.log(`   ${dictPath}`);
      console.log('   (Maps short PCA labels A1, B1... to full sample names)');
      console.log('');
    }

    // 4. Skip agents if forceAutomation
    if (this.forceAutomation) {
      this.log('FORCE AUTOMATION - Skipping agent review');
      const proceed = stage3Output.overall_status !== 'ERROR';
      // Use defaults: simpleEdger if no batch effects, keep all samples
      const deMethod = stage3Output.batch_effect_detected ? 'batch_effect_edger' : 'simpleEdger';

      // FIX: Automatically determine batch specification from metadata
      let batchSpecification = 'auto';
      if (deMethod === 'batch_effect_edger' && this.dataInfo.experimental_design) {
        const design = this.dataInfo.experimental_design.design;
        if (design === 'paired_subjects') {
          batchSpecification = 'paired';
          this.log(`✓ Auto-detected batch specification: paired (from experimental_design.design="${design}")`);
        } else {
          this.log(`⚠️ Unknown experimental design: ${design}, using auto batch detection`);
        }
      }

      this.stageResults.stage3 = {
        output: stage3Output,
        proceed,
        automated: true,
        deMethod,
        batchSpecification,
        outlierAction: 'KEEP_ALL',
        outliersToRemove: []
      };
      return { proceed, automated: true, deMethod };
    }

    // 5. Agent review
    this.log('Calling agents for Stage 3 review...');
    const reviewResult = await this.coordinator.reviewStage3Output(stage3Output, this.dataInfo);

    // 6. Check if user decision needed
    let finalProceed = reviewResult.proceed;
    let finalDeMethod = reviewResult.deMethod || 'simpleEdger';
    let finalBatchSpecification = reviewResult.batchSpecification || null;
    let finalOutlierAction = reviewResult.outlierAction || 'KEEP_ALL';
    let finalOutliersToRemove = reviewResult.outliersToRemove || [];
    let userDecision = null;

    if (reviewResult.consensus.decision.toLowerCase().includes('user_decision')) {
      userDecision = await handleStage3UserDecision(reviewResult, stage3Output, this.config);
      finalProceed = userDecision.proceed;
      finalDeMethod = userDecision.deMethod || finalDeMethod;
      finalBatchSpecification = userDecision.batchSpecification || finalBatchSpecification;
      finalOutlierAction = userDecision.outlierAction || finalOutlierAction;
      finalOutliersToRemove = userDecision.outliersToRemove || [];
      this.log(`User decided: ${finalProceed ? 'PROCEED' : 'ABORT'}`);
    }

    // FIX: Auto-detect batch specification from metadata if agents chose batch_effect_edger
    if (finalDeMethod === 'batch_effect_edger' && (!finalBatchSpecification || finalBatchSpecification === 'auto')) {
      if (this.dataInfo.experimental_design && this.dataInfo.experimental_design.design) {
        const design = this.dataInfo.experimental_design.design;
        if (design === 'paired_subjects') {
          finalBatchSpecification = 'paired';
          this.log(`✓ Auto-detected batch specification: paired (from experimental_design.design="${design}")`);
        } else {
          this.log(`⚠️ Unknown experimental design: ${design}, batch_effect_edgeR will use 'auto' (may fail)`);
          finalBatchSpecification = 'auto';
        }
      } else {
        this.log(`⚠️ No experimental design in metadata, batch_effect_edgeR will use 'auto' (may fail)`);
        finalBatchSpecification = 'auto';
      }
    }

    // 6b. Check if agents recommended removing outliers and ask user for confirmation
    if (finalOutlierAction === 'REMOVE_OUTLIERS' && finalOutliersToRemove.length > 0 && !this.forceAutomation) {
      console.log('');
      console.log('='.repeat(60));
      console.log('  OUTLIER DETECTION');
      console.log('='.repeat(60));
      console.log('');
      console.log(`Agents detected ${finalOutliersToRemove.length} outlier sample(s):`);
      finalOutliersToRemove.forEach((outlier, idx) => {
        console.log(`  ${idx + 1}. ${outlier}`);
      });
      console.log('');
      console.log('Removing outliers will help improve DE analysis quality,');
      console.log('but will reduce sample size and statistical power.');
      console.log('');

      const { confirmOutlierRemoval } = await import('../utils/user_input.js');
      const shouldRemove = await confirmOutlierRemoval(finalOutliersToRemove);

      if (!shouldRemove) {
        this.log('User chose to KEEP all samples (no outlier removal)');
        finalOutlierAction = 'KEEP_ALL';
        finalOutliersToRemove = [];
        // Proceed anyway (keeping all samples is valid)
        finalProceed = true;
      } else {
        this.log(`User confirmed removal of ${finalOutliersToRemove.length} outlier(s)`);
        // CRITICAL FIX: User confirmed outlier removal, so proceed to Stage 4
        finalProceed = true;
      }

      if (!userDecision) {
        userDecision = {
          madeBy: 'user',
          outlierRemovalConfirmed: shouldRemove,
          outliersToRemove: shouldRemove ? finalOutliersToRemove : []
        };
      } else {
        userDecision.outlierRemovalConfirmed = shouldRemove;
        userDecision.outliersToRemove = shouldRemove ? finalOutliersToRemove : [];
      }
    }

    // 7. Log decision
    this.logger.logStageDecision(
      3, 'QC Assessment',
      this.dataInfo, stage3Output,
      reviewResult.responses, reviewResult.consensus,
      {
        proceed: finalProceed,
        deMethod: finalDeMethod,
        batchSpecification: finalBatchSpecification,
        outlierAction: finalOutlierAction,
        outliersToRemove: finalOutliersToRemove,
        userDecision
      }
    );

    // 8. Store result
    this.stageResults.stage3 = {
      output: stage3Output,
      review: reviewResult,
      proceed: finalProceed,
      deMethod: finalDeMethod,
      batchSpecification: finalBatchSpecification,
      outlierAction: finalOutlierAction,
      outliersToRemove: finalOutliersToRemove,
      userDecision
    };

    // 9. Update dataInfo and filter count/RPKM files if outliers removed
    if (finalOutliersToRemove.length > 0) {
      // Parse sample dictionary to convert short labels (A1, B1) to full names
      const stage3Dir = path.join(this.config.output, 'stage3_quantification');
      const dictionaryPath = path.join(stage3Dir, `${this.config.comparison}_sample_dictionary.csv`);

      let shortLabelToFullName = {};
      if (fs.existsSync(dictionaryPath)) {
        const dictData = fs.readFileSync(dictionaryPath, 'utf-8').split('\n');
        // Parse CSV: short_label,full_name,group (skip header row)
        for (let i = 1; i < dictData.length; i++) {
          const line = dictData[i].trim();
          if (!line) continue;
          const [shortLabel, fullName] = line.split(',');
          if (shortLabel && fullName) {
            // Remove quotes and trim
            const cleanShortLabel = shortLabel.trim().replace(/^"|"$/g, '');
            const cleanFullName = fullName.trim().replace(/^"|"$/g, '');
            shortLabelToFullName[cleanShortLabel] = cleanFullName;
            this.log(`  Loaded mapping: ${cleanShortLabel} → ${cleanFullName}`);
          }
        }
      }

      // Convert short labels to full names (if agents used short labels)
      const fullOutlierNames = finalOutliersToRemove.map(outlier => {
        // If outlier looks like a short label (A1, B2, etc.), convert to full name
        if (shortLabelToFullName[outlier]) {
          this.log(`Converting short label "${outlier}" → "${shortLabelToFullName[outlier]}"`);
          return shortLabelToFullName[outlier];
        }
        // Otherwise, assume it's already a full name
        return outlier;
      });

      this.dataInfo.samples = this.dataInfo.samples.filter(
        s => !fullOutlierNames.includes(s.name)
      );
      this.log(`Removed ${fullOutlierNames.length} outlier(s) from analysis`);

      // Helper function to filter CSV files
      const filterCSV = (filePath, fileDesc) => {
        if (fs.existsSync(filePath)) {
          this.log(`Filtering ${fileDesc}...`);
          const data = fs.readFileSync(filePath, 'utf-8').split('\n');

          if (data.length > 0) {
            // Parse header
            const header = data[0].split(',');

            // Find indices of columns to keep
            const keepIndices = header.map((colName, idx) => {
              // Keep first 2 columns (ID, Length or ID, GeneSymbol)
              if (idx < 2) return true;
              // Keep columns that are NOT in outlier list (use full names)
              return !fullOutlierNames.some(outlier => colName.includes(outlier));
            });

            // Filter all rows
            const filteredData = data.map(line => {
              const cols = line.split(',');
              return cols.filter((_, idx) => keepIndices[idx]).join(',');
            });

            // Write filtered file
            fs.writeFileSync(filePath, filteredData.join('\n'));
            this.log(`  ${fileDesc} filtered: ${header.length - 2} → ${this.dataInfo.samples.length} samples`);
          }
        }
      };

      // Filter count and RPKM files (reuse stage3Dir from above)
      filterCSV(path.join(stage3Dir, `${this.config.comparison}.count.filtered.csv`), 'count file');
      filterCSV(path.join(stage3Dir, `${this.config.comparison}.rpkm.csv`), 'RPKM file');
    }

    console.log('');
    console.log('='.repeat(60));
    console.log(`Stage 3 Result: ${finalProceed ? 'PASS' : 'FAIL'}`);
    console.log('='.repeat(60));
    console.log('');
    console.log('Agent Decisions (from visual PCA inspection):');
    console.log(`  Batch Effects: ${finalDeMethod === 'batch_effect_edger' ? 'YES (using batch correction)' : 'NO (standard analysis)'}`);
    console.log(`  DE Method: ${finalDeMethod}`);
    if (finalBatchSpecification) {
      console.log(`  Batch Specification: ${finalBatchSpecification}`);
    }
    console.log(`  Outlier Action: ${finalOutlierAction}`);
    if (finalOutliersToRemove.length > 0) {
      console.log(`  Outliers Removed: ${finalOutliersToRemove.join(', ')}`);
    }

    return {
      proceed: finalProceed,
      output: stage3Output,
      review: reviewResult,
      deMethod: finalDeMethod,
      batchSpecification: finalBatchSpecification,
      outliersToRemove: finalOutliersToRemove
    };
  }

  /**
   * Run Stage 4: DE Analysis
   */
  async runStage4() {
    console.log('');
    console.log('='.repeat(60));
    console.log('  STAGE 4: Differential Expression Analysis');
    console.log('='.repeat(60));
    console.log('');

    const absOutputDir = path.resolve(this.outputDir);
    const scriptPath = path.resolve(absOutputDir, `stage_4_${this.comparison}.sh`);

    // 1. Generate script (pass stage3 result for DE method decision)
    this.log('Generating Stage 4 script...');
    generateStage4Script(this.dataInfo, this.config, scriptPath, this.stageResults.stage3);
    this.log(`Script: ${scriptPath}`);

    if (this.dryRun) {
      this.log('DRY RUN - Script generated but not executed');
      return { proceed: true, dryRun: true };
    }

    // 2. Execute script
    this.log('Executing Stage 4 script (this may take 5-15 minutes)...');
    const execResult = this.executeScript(scriptPath, 'DE Analysis', 3600000); // 1 hour timeout

    if (!execResult.success) {
      console.error(`Stage 4 execution failed: ${execResult.error}`);
      return { proceed: false, error: execResult.error };
    }

    // 3. Parse output
    this.log('Parsing Stage 4 output...');
    const stage4Output = parseStage4Output(absOutputDir);
    this.log(`Status: ${stage4Output.overall_status}`);
    this.log(`Total DEGs: ${stage4Output.total_degs} (Up: ${stage4Output.num_degs_up}, Down: ${stage4Output.num_degs_down})`);

    // 4. Skip agents if forceAutomation
    if (this.forceAutomation) {
      this.log('FORCE AUTOMATION - Skipping agent review');
      const approve = stage4Output.overall_status !== 'ERROR';
      this.stageResults.stage4 = {
        output: stage4Output,
        proceed: approve,
        approve,
        automated: true
      };
      return { proceed: approve, automated: true };
    }

    // 5. Agent review
    this.log('Calling agents for Stage 4 review...');
    const reviewResult = await this.coordinator.reviewStage4Output(stage4Output, this.dataInfo);

    // 6. Check if user decision needed
    let finalApprove = reviewResult.approve;
    let userDecision = null;

    if (reviewResult.consensus.decision.toLowerCase().includes('user_decision')) {
      // For Stage 4, we can create a simple user decision handler
      // For now, just ask yes/no to approve
      const userChoice = await askYesNo(
        'Agents could not reach consensus on DE results. Do you approve the results?',
        true
      );
      finalApprove = userChoice;
      userDecision = {
        madeBy: 'user',
        approve: userChoice,
        reason: userChoice ? 'User approved DE results' : 'User rejected DE results'
      };
      this.log(`User decided: ${finalApprove ? 'APPROVE' : 'REJECT'}`);
    }

    // 7. Log decision
    this.logger.logStageDecision(
      4, 'DE Analysis',
      this.dataInfo, stage4Output,
      reviewResult.responses, reviewResult.consensus,
      {
        proceed: finalApprove,
        approve: finalApprove,
        userDecision
      }
    );

    // 8. Store result
    this.stageResults.stage4 = {
      output: stage4Output,
      review: reviewResult,
      proceed: finalApprove,
      approve: finalApprove,
      userDecision
    };

    console.log('');
    console.log(`Stage 4 Result: ${finalApprove ? 'APPROVED' : 'REJECTED'}`);
    console.log(`Total DEGs: ${stage4Output.total_degs} (FDR < ${stage4Output.fdr_threshold})`);

    return {
      proceed: finalApprove,
      approve: finalApprove,
      output: stage4Output,
      review: reviewResult
    };
  }

  /**
   * Run the full pipeline
   */
  /**
   * Stage 3a: Cell Cycle Scoring (scRNA-seq only)
   */
  async runStage3a() {
    console.log('');
    console.log('='.repeat(60));
    console.log('  STAGE 3a: Cell Cycle Scoring (scRNA-seq)');
    console.log('='.repeat(60));
    console.log('');

    const absOutputDir = path.resolve(this.outputDir);
    const scriptPath = path.resolve(absOutputDir, `stage_3a_${this.comparison}.sh`);

    // 1. Generate script
    this.log('Generating Stage 3a script...');
    generateStage3aScript(this.dataInfo, this.config, scriptPath);
    this.log(`Script: ${scriptPath}`);

    if (this.dryRun) {
      this.log('DRY RUN - Script generated but not executed');
      return { proceed: true, dryRun: true };
    }

    // 2. Execute script
    this.log('Executing Stage 3a script (cell cycle scoring)...');
    const execResult = this.executeScript(scriptPath, 'Cell Cycle Scoring', 900000); // 15 min timeout

    if (!execResult.success) {
      console.error(`Stage 3a execution failed: ${execResult.error}`);
      return { proceed: false, error: execResult.error };
    }

    // 3. Parse output
    this.log('Parsing Stage 3a output...');
    const stage3aOutput = parseStage3aOutput(absOutputDir, this.comparison);
    this.log(`Status: ${stage3aOutput.overall_status}`);

    // 4. Agent review (informational only - always proceeds)
    if (!this.forceAutomation) {
      this.log('Calling agents for Stage 3a review (informational only)...');
      const formattedOutput = formatStage3aForAgents(stage3aOutput, this.dataInfo);

      // Add PCA plot as image
      const images = stage3aOutput.pca_plot_path ? [{ path: stage3aOutput.pca_plot_path }] : [];

      await this.coordinator.consultAgentsWithStagePrompts(
        '3a',
        formattedOutput,
        { stage3aOutput, dataInfo: this.dataInfo, images },
        'stage3a',
        `${this.comparison}_stage3a`
      );
    }

    // 5. Store result (always proceed to Stage 3b)
    this.stageResults.stage3a = {
      output: stage3aOutput,
      proceed: true
    };

    console.log('');
    console.log('✓ Stage 3a complete - Proceeding to Stage 3b (Cell Cycle Regression)');
    return { proceed: true };
  }

  /**
   * Stage 3b: Cell Cycle Regression (scRNA-seq only)
   */
  async runStage3b() {
    console.log('');
    console.log('='.repeat(60));
    console.log('  STAGE 3b: Cell Cycle Regression (scRNA-seq)');
    console.log('='.repeat(60));
    console.log('');

    const absOutputDir = path.resolve(this.outputDir);
    const scriptPath = path.resolve(absOutputDir, `stage_3b_${this.comparison}.sh`);

    // 1. Generate script (pass Stage 3a result for RDS path)
    this.log('Generating Stage 3b script...');
    generateStage3bScript(this.dataInfo, this.config, scriptPath, this.stageResults.stage3a?.output);
    this.log(`Script: ${scriptPath}`);

    if (this.dryRun) {
      this.log('DRY RUN - Script generated but not executed');
      return { proceed: true, dryRun: true };
    }

    // 2. Execute script
    this.log('Executing Stage 3b script (cell cycle regression)...');
    const execResult = this.executeScript(scriptPath, 'Cell Cycle Regression', 1200000); // 20 min timeout

    if (!execResult.success) {
      console.error(`Stage 3b execution failed: ${execResult.error}`);
      return { proceed: false, error: execResult.error };
    }

    // 3. Parse output
    this.log('Parsing Stage 3b output...');
    const stage3bOutput = parseStage3bOutput(absOutputDir, this.comparison);
    this.log(`Status: ${stage3bOutput.overall_status}`);
    this.log(`Regression Status: ${stage3bOutput.regression_status}`);

    // 4. Agent review (decision required)
    let finalProceed = stage3bOutput.regression_status === 'SUCCESS';

    if (!this.forceAutomation) {
      this.log('Calling agents for Stage 3b review...');
      const formattedOutput = formatStage3bForAgents(stage3bOutput, this.dataInfo);

      // Add comparison plot as image
      const images = stage3bOutput.comparison_plot_path ? [{ path: stage3bOutput.comparison_plot_path }] : [];

      const reviewResult = await this.coordinator.consultAgentsWithStagePrompts(
        '3b',
        formattedOutput,
        { stage3bOutput, dataInfo: this.dataInfo, images },
        'stage3b',
        `${this.comparison}_stage3b`
      );

      // Check consensus decision
      const decision = reviewResult.consensus.decision.toUpperCase();
      if (decision.includes('SUCCESS')) {
        finalProceed = true;
      } else if (decision.includes('PARTIAL')) {
        console.log('⚠️  Regression partially successful - proceeding with caution');
        finalProceed = true;
      } else {
        console.log('❌ Regression failed - aborting');
        finalProceed = false;
      }
    }

    // 5. Store result
    this.stageResults.stage3b = {
      output: stage3bOutput,
      proceed: finalProceed
    };

    if (finalProceed) {
      console.log('');
      console.log('✓ Stage 3b complete - Proceeding to Stage 3c (Clustering + UMAP)');
    }

    return { proceed: finalProceed };
  }

  /**
   * Stage 3c: Clustering + UMAP (scRNA-seq only)
   */
  async runStage3c() {
    console.log('');
    console.log('='.repeat(60));
    console.log('  STAGE 3c: Clustering + UMAP (scRNA-seq)');
    console.log('='.repeat(60));
    console.log('');

    const absOutputDir = path.resolve(this.outputDir);
    const scriptPath = path.resolve(absOutputDir, `stage_3c_${this.comparison}.sh`);

    // 1. Generate script (pass Stage 3b result for RDS path)
    this.log('Generating Stage 3c script...');
    generateStage3cScript(this.dataInfo, this.config, scriptPath, this.stageResults.stage3b?.output);
    this.log(`Script: ${scriptPath}`);

    if (this.dryRun) {
      this.log('DRY RUN - Script generated but not executed');
      return { proceed: true, dryRun: true };
    }

    // 2. Execute script
    this.log('Executing Stage 3c script (clustering + UMAP)...');
    const execResult = this.executeScript(scriptPath, 'Clustering + UMAP', 1800000); // 30 min timeout

    if (!execResult.success) {
      console.error(`Stage 3c execution failed: ${execResult.error}`);
      return { proceed: false, error: execResult.error };
    }

    // 3. Parse output
    this.log('Parsing Stage 3c output...');
    const stage3cOutput = parseStage3cOutput(absOutputDir, this.comparison);
    this.log(`Status: ${stage3cOutput.overall_status}`);
    this.log(`Clusters found: ${stage3cOutput.n_clusters}`);

    // 4. Agent review (decision required)
    let finalProceed = stage3cOutput.overall_status === 'SUCCESS';

    if (!this.forceAutomation) {
      this.log('Calling agents for Stage 3c review...');
      const formattedOutput = formatStage3cForAgents(stage3cOutput, this.dataInfo);

      // Add both UMAP plots as images
      const images = [];
      if (stage3cOutput.umap_clusters_plot_path) {
        images.push({ path: stage3cOutput.umap_clusters_plot_path });
      }
      if (stage3cOutput.umap_phase_plot_path) {
        images.push({ path: stage3cOutput.umap_phase_plot_path });
      }

      const reviewResult = await this.coordinator.consultAgentsWithStagePrompts(
        '3c',
        formattedOutput,
        { stage3cOutput, dataInfo: this.dataInfo, images },
        'stage3c',
        `${this.comparison}_stage3c`
      );

      // Check consensus decision
      const decision = reviewResult.consensus.decision.toUpperCase();
      if (decision.includes('APPROVE')) {
        finalProceed = true;
      } else if (decision.includes('ADJUST')) {
        console.log('⚠️  Agents recommend adjusting resolution - proceeding anyway');
        finalProceed = true;
      } else {
        console.log('❌ Clustering failed - aborting');
        finalProceed = false;
      }
    }

    // 5. Store result
    this.stageResults.stage3c = {
      output: stage3cOutput,
      proceed: finalProceed
    };

    if (finalProceed) {
      console.log('');
      console.log('✓ Stage 3c complete - Proceeding to Stage 4 (Differential Expression)');
    }

    return { proceed: finalProceed };
  }

  async run() {
    console.log('');
    console.log('#'.repeat(60));
    console.log('  GeneExpert Staged Pipeline');
    console.log('#'.repeat(60));
    console.log('');
    console.log(`Input:  ${this.inputDir}`);
    console.log(`Output: ${this.outputDir}`);
    console.log(`Organism: ${this.organism}`);
    console.log(`Data Type: ${this.dataType.toUpperCase()}`);
    console.log(`Mode: ${this.forceAutomation ? 'AUTOMATION (no agents)' : this.singleAgent ? `SINGLE-AGENT (${this.singleAgent})` : 'MULTI-AGENT'}`);
    console.log('');

    // Initialize
    this.initialize();

    // Stage 1
    const stage1Result = await this.runStage1();
    if (!stage1Result.proceed) {
      this.finalize(false, 'Stage 1 failed');
      return { success: false, failedAt: 1 };
    }

    // Stage 2
    const stage2Result = await this.runStage2();
    if (!stage2Result.proceed) {
      this.finalize(false, 'Stage 2 failed');
      return { success: false, failedAt: 2 };
    }

    // Branch based on data type
    if (this.dataType === 'scrna') {
      // scRNA-seq pipeline: Stage 3a → 3b → 3c → 4
      console.log('');
      console.log('🔬 scRNA-seq Pipeline: Running cell cycle correction stages...');

      // Stage 3a: Cell Cycle Scoring
      const stage3aResult = await this.runStage3a();
      if (!stage3aResult.proceed) {
        this.finalize(false, 'Stage 3a failed');
        return { success: false, failedAt: '3a' };
      }

      // Stage 3b: Cell Cycle Regression
      const stage3bResult = await this.runStage3b();
      if (!stage3bResult.proceed) {
        this.finalize(false, 'Stage 3b failed');
        return { success: false, failedAt: '3b' };
      }

      // Stage 3c: Clustering + UMAP
      const stage3cResult = await this.runStage3c();
      if (!stage3cResult.proceed) {
        this.finalize(false, 'Stage 3c failed');
        return { success: false, failedAt: '3c' };
      }
    } else {
      // Bulk RNA-seq pipeline: Stage 3 → 4
      const stage3Result = await this.runStage3();
      if (!stage3Result.proceed) {
        this.finalize(false, 'Stage 3 failed');
        return { success: false, failedAt: 3 };
      }
    }

    // Stage 4 (common for both bulk and scRNA-seq)
    const stage4Result = await this.runStage4();
    if (!stage4Result.proceed) {
      this.finalize(false, 'Stage 4 failed');
      return { success: false, failedAt: 4 };
    }

    // Success!
    this.finalize(true);
    return { success: true, stageResults: this.stageResults };
  }

  /**
   * Finalize the analysis
   */
  finalize(success, reason = '') {
    console.log('');
    console.log('#'.repeat(60));
    console.log(`  Pipeline ${success ? 'COMPLETE' : 'FAILED'}`);
    if (reason) console.log(`  Reason: ${reason}`);
    console.log('#'.repeat(60));
    console.log('');

    // Finalize logger
    if (this.logger) {
      this.logger.finalize({
        dataInfo: this.dataInfo,
        config: this.config,
        success,
        reason,
        stageResults: this.stageResults
      });
    }

    // Summary
    console.log('Stage Summary:');
    console.log(`  Stage 1 (FASTQ Validation): ${this.stageResults.stage1?.proceed ? 'PASS' : 'FAIL/SKIP'}`);
    console.log(`  Stage 2 (Alignment QC):     ${this.stageResults.stage2?.proceed ? 'PASS' : 'FAIL/SKIP'}`);

    if (this.dataType === 'scrna') {
      console.log(`  Stage 3a (Cell Cycle Score): ${this.stageResults.stage3a?.proceed ? 'PASS' : 'FAIL/SKIP'}`);
      console.log(`  Stage 3b (Cell Cycle Regr):  ${this.stageResults.stage3b?.proceed ? 'PASS' : 'FAIL/SKIP'}`);
      console.log(`  Stage 3c (Clustering/UMAP):  ${this.stageResults.stage3c?.proceed ? 'PASS' : 'FAIL/SKIP'}`);
    } else {
      console.log(`  Stage 3 (Quantification):   ${this.stageResults.stage3?.proceed ? 'PASS' : 'FAIL/SKIP'}`);
    }

    console.log(`  Stage 4 (DE Analysis):      ${this.stageResults.stage4?.proceed ? 'PASS' : 'FAIL/SKIP'}`);
    console.log('');

    // User input summary
    const userInputStages = [];
    if (this.stageResults.stage1?.userDecision) userInputStages.push('Stage 1');
    if (this.stageResults.stage2?.userDecision) userInputStages.push('Stage 2');
    if (this.stageResults.stage3?.userDecision) userInputStages.push('Stage 3');
    if (this.stageResults.stage4?.userDecision) userInputStages.push('Stage 4');

    if (userInputStages.length > 0) {
      console.log(`User input required at: ${userInputStages.join(', ')}`);
    } else {
      console.log('User input required at: None (full automation)');
    }
    console.log('');
  }
}

/**
 * CLI entry point for staged executor
 */
export async function runStagedPipeline(options) {
  const executor = new StagedExecutor(options);
  return executor.run();
}

export default StagedExecutor;
