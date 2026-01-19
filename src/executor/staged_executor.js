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

    // State tracking
    this.stageResults = {
      stage1: null,
      stage2: null,
      stage3: null,
      stage4: null
    };
    this.dataInfo = null;
    this.config = null;

    // Initialize coordinator
    this.coordinator = new Coordinator({
      verbose: this.verbose,
      singleAgent: this.singleAgent
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

    const files = fs.readdirSync(this.inputDir).filter(f => f.endsWith('.fastq.gz'));
    const samples = [];
    const pairedEnd = files.some(f => f.includes('_R2_'));

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
      comparison: this.comparison
    };

    this.log(`Detected ${samples.length} samples (${pairedEnd ? 'paired-end' : 'single-end'})`);
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
      forceAutomation: this.forceAutomation
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
      userDecision = await handleStage1UserDecision(reviewResult, stage1Output);
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
      userDecision = await handleStage2UserDecision(reviewResult, stage2Output);
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

    // 4. Skip agents if forceAutomation
    if (this.forceAutomation) {
      this.log('FORCE AUTOMATION - Skipping agent review');
      const proceed = stage3Output.overall_status !== 'ERROR';
      // Use defaults: simpleEdger if no batch effects, keep all samples
      const deMethod = stage3Output.batch_effect_detected ? 'batch_effect_edger' : 'simpleEdger';
      this.stageResults.stage3 = {
        output: stage3Output,
        proceed,
        automated: true,
        deMethod,
        batchSpecification: 'auto',
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
      userDecision = await handleStage3UserDecision(reviewResult, stage3Output);
      finalProceed = userDecision.proceed;
      finalDeMethod = userDecision.deMethod || finalDeMethod;
      finalBatchSpecification = userDecision.batchSpecification || finalBatchSpecification;
      finalOutlierAction = userDecision.outlierAction || finalOutlierAction;
      finalOutliersToRemove = userDecision.outliersToRemove || [];
      this.log(`User decided: ${finalProceed ? 'PROCEED' : 'ABORT'}`);
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

    // 9. Update dataInfo if outliers removed
    if (finalOutliersToRemove.length > 0) {
      this.dataInfo.samples = this.dataInfo.samples.filter(
        s => !finalOutliersToRemove.includes(s.name)
      );
      this.log(`Removed ${finalOutliersToRemove.length} outlier(s) from analysis`);
    }

    console.log('');
    console.log(`Stage 3 Result: ${finalProceed ? 'PASS' : 'FAIL'}`);
    console.log(`DE Method: ${finalDeMethod}`);
    if (finalOutliersToRemove.length > 0) {
      console.log(`Outliers removed: ${finalOutliersToRemove.join(', ')}`);
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
  async run() {
    console.log('');
    console.log('#'.repeat(60));
    console.log('  GeneExpert Staged Pipeline');
    console.log('#'.repeat(60));
    console.log('');
    console.log(`Input:  ${this.inputDir}`);
    console.log(`Output: ${this.outputDir}`);
    console.log(`Organism: ${this.organism}`);
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

    // Stage 3
    const stage3Result = await this.runStage3();
    if (!stage3Result.proceed) {
      this.finalize(false, 'Stage 3 failed');
      return { success: false, failedAt: 3 };
    }

    // Stage 4
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
    console.log(`  Stage 3 (Quantification):   ${this.stageResults.stage3?.proceed ? 'PASS' : 'FAIL/SKIP'}`);
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
