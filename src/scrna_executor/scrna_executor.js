/**
 * scRNA-seq Staged Executor
 *
 * Orchestrates 5-stage Seurat pipeline with 3 agent checkpoints:
 *   Stage 1: Load + QC (no agent)
 *   Stage 2: QC Filtering → Stats Agent
 *   Stage 3: Normalize + HVG (no agent)
 *   Stage 4: PCA → Stats Agent (PC selection)
 *   Stage 5: Clustering + Markers → Pipeline Agent
 */

import path from 'path';
import fs from 'fs';
import { execSync } from 'child_process';

// Stage modules
import {
  generateStage1Script,
  parseStage1Output,
  formatStage1ForDisplay
} from '../scrna_stages/stage1_load_qc.js';

import {
  generateStage2Script,
  parseStage2Output,
  formatStage2ForAgents
} from '../scrna_stages/stage2_filter_qc.js';

import {
  generateStage3Script,
  parseStage3Output,
  formatStage3ForDisplay
} from '../scrna_stages/stage3_normalize_hvg.js';

import {
  generateStage4Script,
  parseStage4Output,
  formatStage4ForAgents
} from '../scrna_stages/stage4_pca.js';

import {
  generateStage5Script,
  parseStage5Output,
  formatStage5ForAgents
} from '../scrna_stages/stage5_cluster_markers.js';

// Coordinator and utilities (REUSE from bulk RNA)
import { Coordinator } from '../coordinator/orchestrator.js';
import { Logger } from '../utils/logger.js';
import {
  handleScRNAStage2UserDecision,
  handleScRNAStage4UserDecision,
  handleScRNAStage5UserDecision
} from '../utils/user_input.js';

/**
 * scRNA-seq Staged Executor Class
 */
export class ScRNAExecutor {
  constructor(options = {}) {
    this.inputDir = options.inputDir;
    this.outputDir = options.outputDir;
    this.organism = options.organism || 'mouse';
    this.datasetName = options.datasetName || path.basename(options.inputDir);
    this.verbose = options.verbose !== false;
    this.dryRun = options.dryRun || false;
    this.singleAgent = options.singleAgent || null;
    this.forceAutomation = options.forceAutomation || false;
    this.sequentialChain = options.sequentialChain || false;
    this.roleAssignments = options.roleAssignments || {  // NEW: Role assignments for ablation study
      gptRole: 'stats',
      claudeRole: 'pipeline',
      geminiRole: 'biology'
    };

    // Stage results tracking
    this.stageResults = {
      stage1: null,
      stage2: null,
      stage3: null,
      stage4: null,
      stage5: null
    };

    this.dataInfo = null;
    this.config = null;

    // Initialize coordinator (REUSE from bulk RNA)
    this.coordinator = new Coordinator({
      verbose: this.verbose,
      singleAgent: this.singleAgent,
      sequentialChain: this.sequentialChain,
      roleAssignments: this.roleAssignments  // NEW: Pass role assignments for ablation study
    });

    this.logger = null;
  }

  log(message) {
    if (this.verbose) {
      console.log(`[scRNA Executor] ${message}`);
    }
  }

  /**
   * Detect scRNA data (supports .h5, 10x directory, or CSV formats)
   */
  detectDataInfo() {
    this.log('Detecting scRNA-seq data...');

    const files = fs.readdirSync(this.inputDir);

    // Option 1: Look for .h5 file
    const h5File = files.find(f => f.endsWith('.h5'));
    if (h5File) {
      this.dataInfo = {
        type: 'scrna',
        h5File: path.join(this.inputDir, h5File),
        organism: this.organism,
        datasetName: this.datasetName
      };
      this.log(`Found 10x HDF5 data: ${h5File}`);
      return this.dataInfo;
    }

    // Option 2: Look for 10x directory format (barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz)
    // Case 2a: Files directly in input directory
    const hasBarcodes = files.some(f => f === 'barcodes.tsv.gz' || f === 'barcodes.tsv');
    const hasFeatures = files.some(f => f === 'features.tsv.gz' || f === 'features.tsv' || f === 'genes.tsv.gz' || f === 'genes.tsv');
    const hasMatrix = files.some(f => f === 'matrix.mtx.gz' || f === 'matrix.mtx');

    if (hasBarcodes && hasFeatures && hasMatrix) {
      this.dataInfo = {
        type: 'scrna',
        inputDir: this.inputDir,  // Directory with 10x files
        organism: this.organism,
        datasetName: this.datasetName
      };
      this.log(`Found 10x directory format (files in root)`);
      return this.dataInfo;
    }

    // Case 2b: Files in filtered_feature_bc_matrix/ subdirectory
    const matrixDir = path.join(this.inputDir, 'filtered_feature_bc_matrix');
    if (fs.existsSync(matrixDir) && fs.statSync(matrixDir).isDirectory()) {
      this.dataInfo = {
        type: 'scrna',
        inputDir: matrixDir,  // Subdirectory with 10x files
        organism: this.organism,
        datasetName: this.datasetName
      };
      this.log(`Found 10x directory format (filtered_feature_bc_matrix/)`);
      return this.dataInfo;
    }

    // Option 3: Look for CSV/TSV files
    const csvFile = files.find(f => /\.(csv|txt|tsv)$/i.test(f));
    if (csvFile) {
      this.dataInfo = {
        type: 'scrna',
        csvFile: path.join(this.inputDir, csvFile),
        organism: this.organism,
        datasetName: this.datasetName
      };
      this.log(`Found CSV/TSV data: ${csvFile}`);
      return this.dataInfo;
    }

    // No valid format found
    throw new Error(
      'No scRNA-seq data found. Expected one of:\n' +
      '  - filtered_feature_bc_matrix.h5 (10x HDF5 format)\n' +
      '  - barcodes.tsv.gz + features.tsv.gz + matrix.mtx.gz (10x directory format)\n' +
      '  - filtered_feature_bc_matrix/ subdirectory (10x directory format)\n' +
      '  - .csv/.txt/.tsv file (genes as rows, cells as columns)'
    );
  }

  /**
   * Initialize executor
   */
  initialize() {
    const absOutputDir = path.resolve(this.outputDir);
    if (!fs.existsSync(absOutputDir)) {
      fs.mkdirSync(absOutputDir, { recursive: true });
    }

    this.config = {
      input: path.resolve(this.inputDir),
      output: absOutputDir,
      organism: this.organism,
      datasetName: this.datasetName,
      singleAgent: this.singleAgent,
      forceAutomation: this.forceAutomation
    };

    // Initialize logger
    this.logger = new Logger(absOutputDir, 'scrna_analysis', this.config);

    // Set dataset for coordinator
    this.coordinator.setDataset(this.datasetName);

    // Detect data
    this.detectDataInfo();

    return this;
  }

  /**
   * Execute a bash script
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
      return { success: false, exitCode: error.status, error: error.message };
    }
  }

  /**
   * Run Stage 1: Load + QC
   */
  async runStage1() {
    console.log('');
    console.log('='.repeat(60));
    console.log('  STAGE 1: Load + QC Metrics');
    console.log('='.repeat(60));
    console.log('');

    const absOutputDir = path.resolve(this.outputDir);
    const scriptPath = path.resolve(absOutputDir, `stage_1_scrna_load_qc.sh`);

    // Generate and execute script
    this.log('Generating Stage 1 script...');
    generateStage1Script(this.dataInfo, this.config, scriptPath);

    if (this.dryRun) {
      this.log('DRY RUN - Script generated but not executed');
      return { proceed: true, dryRun: true };
    }

    this.log('Executing Stage 1 script...');
    const execResult = this.executeScript(scriptPath, 'Load + QC');

    if (!execResult.success) {
      console.error(`Stage 1 execution failed: ${execResult.error}`);
      return { proceed: false, error: execResult.error };
    }

    // Parse output
    this.log('Parsing Stage 1 output...');
    const stage1Output = parseStage1Output(absOutputDir);

    // Display results (no agent for Stage 1)
    console.log('');
    console.log(formatStage1ForDisplay(stage1Output));

    // Log to JSONL (automated stage, no agents)
    this.logger.logStageDecision(
      1, 'Load + QC',
      this.dataInfo, stage1Output,
      {}, // No agent responses
      {
        decision: 'auto_proceed',
        confidence: 1.0,
        votes: { approve: 1, reject: 0, uncertain: 0 },
        reasoning: 'Automated stage - no agent review',
        recommendations: []
      },
      { proceed: true, automated: true }
    );

    this.stageResults.stage1 = { output: stage1Output, proceed: true };
    return { proceed: true, output: stage1Output };
  }

  /**
   * Run Stage 2: QC Filtering (AGENT CHECKPOINT 1)
   */
  async runStage2() {
    console.log('');
    console.log('='.repeat(60));
    console.log('  STAGE 2: QC Threshold Recommendation + Filtering');
    console.log('='.repeat(60));
    console.log('');

    const absOutputDir = path.resolve(this.outputDir);

    // STEP 1: Agent review of Stage 1 QC metrics → Recommend thresholds
    let thresholds;
    if (this.forceAutomation) {
      this.log('FORCE AUTOMATION - Using default thresholds');
      thresholds = {
        nFeature_min: 200,
        nFeature_max: 6000,
        percent_mt_max: 10
      };
    } else {
      this.log('Calling ALL 3 AGENTS to recommend QC thresholds...');
      const thresholdReview = await this.coordinator.reviewScRNAStage1ForThresholds(
        this.stageResults.stage1.output,
        this.dataInfo
      );

      if (thresholdReview.decision === 'INSUFFICIENT_DATA') {
        console.error('Agents determined dataset has insufficient data for analysis');
        return { proceed: false, error: 'Insufficient data (agent decision)' };
      }

      // Check if user decision needed (low consensus or disagreement on thresholds)
      let userDecision = null;
      const needsUserInput = thresholdReview.consensus.confidence < 0.7 ||
                             thresholdReview.consensus.decision.toLowerCase().includes('user_decision') ||
                             hasThresholdDisagreement(thresholdReview.agentResponses);

      if (needsUserInput) {
        this.log('Low consensus or threshold disagreement - escalating to user...');
        userDecision = await handleScRNAStage2UserDecision(
          thresholdReview,
          this.stageResults.stage1.output
        );

        if (!userDecision.proceed) {
          console.log('User aborted analysis');
          return { proceed: false, error: 'User aborted' };
        }

        thresholds = userDecision.thresholds;
      } else {
        thresholds = thresholdReview.thresholds;
      }

      // Log agent decision (Stage 2 checkpoint)
      if (thresholdReview.agentResponses && thresholdReview.consensus) {
        this.logger.logStageDecision(
          2, 'QC Filtering',
          this.dataInfo, this.stageResults.stage1.output,
          thresholdReview.agentResponses,
          thresholdReview.consensus,
          { proceed: true, thresholds, userDecision: userDecision || null }
        );
      }
      console.log('');
      console.log('Agent-Recommended Thresholds:');
      console.log(`  nFeature: ${thresholds.nFeature_min} - ${thresholds.nFeature_max}`);
      console.log(`  %MT: < ${thresholds.percent_mt_max}%`);
      console.log('');
    }

    // STEP 2: Generate and execute filtering script with agent thresholds
    const scriptPath = path.resolve(absOutputDir, `stage_2_scrna_filter_qc.sh`);

    this.log('Generating Stage 2 filtering script with agent thresholds...');
    generateStage2Script(this.stageResults.stage1.output, thresholds, this.config, scriptPath);

    if (this.dryRun) {
      return { proceed: true, dryRun: true };
    }

    this.log('Executing Stage 2 script...');
    const execResult = this.executeScript(scriptPath, 'QC Filtering');

    if (!execResult.success) {
      console.error(`Stage 2 execution failed: ${execResult.error}`);
      return { proceed: false, error: execResult.error };
    }

    // Parse output
    this.log('Parsing Stage 2 output...');
    const stage2Output = parseStage2Output(absOutputDir);

    console.log('');
    console.log(`Stage 2 Filtering Complete: ${stage2Output.cells_before} → ${stage2Output.cells_after} cells (${stage2Output.percent_removed}% removed)`);

    // Check if too many cells were removed (safety check)
    if (stage2Output.percent_removed > 80) {
      console.error('WARNING: >80% of cells removed! This may indicate problematic thresholds.');
    }

    const finalProceed = stage2Output.cells_after >= 100;  // Minimum 100 cells required

    this.stageResults.stage2 = { output: stage2Output, proceed: finalProceed, thresholds };

    console.log('');
    console.log(`Stage 2 Result: ${finalProceed ? 'PROCEED' : 'STOP (insufficient cells)'}`);
    return { proceed: finalProceed, output: stage2Output, thresholds };
  }

  /**
   * Run Stage 3: Normalize + HVG
   */
  async runStage3() {
    console.log('');
    console.log('='.repeat(60));
    console.log('  STAGE 3: Normalization + HVG');
    console.log('='.repeat(60));
    console.log('');

    const absOutputDir = path.resolve(this.outputDir);
    const scriptPath = path.resolve(absOutputDir, `stage_3_scrna_normalize_hvg.sh`);

    this.log('Generating Stage 3 script...');
    generateStage3Script(this.stageResults.stage2.output, this.config, scriptPath);

    if (this.dryRun) {
      return { proceed: true, dryRun: true };
    }

    this.log('Executing Stage 3 script...');
    const execResult = this.executeScript(scriptPath, 'Normalize + HVG');

    if (!execResult.success) {
      console.error(`Stage 3 execution failed: ${execResult.error}`);
      return { proceed: false, error: execResult.error };
    }

    const stage3Output = parseStage3Output(absOutputDir);
    console.log('');
    console.log(formatStage3ForDisplay(stage3Output));

    // Log to JSONL (automated stage, no agents)
    this.logger.logStageDecision(
      3, 'Normalize + HVG',
      this.dataInfo, stage3Output,
      {}, // No agent responses
      {
        decision: 'auto_proceed',
        confidence: 1.0,
        votes: { approve: 1, reject: 0, uncertain: 0 },
        reasoning: 'Automated stage - no agent review',
        recommendations: []
      },
      { proceed: true, automated: true }
    );

    this.stageResults.stage3 = { output: stage3Output, proceed: true };
    return { proceed: true, output: stage3Output };
  }

  /**
   * Run Stage 4: PCA (AGENT CHECKPOINT 2)
   */
  async runStage4() {
    console.log('');
    console.log('='.repeat(60));
    console.log('  STAGE 4: PCA + PC Selection');
    console.log('='.repeat(60));
    console.log('');

    const absOutputDir = path.resolve(this.outputDir);
    const scriptPath = path.resolve(absOutputDir, `stage_4_scrna_pca.sh`);

    this.log('Generating Stage 4 script...');
    generateStage4Script(this.stageResults.stage3.output, this.config, scriptPath);

    if (this.dryRun) {
      return { proceed: true, dryRun: true };
    }

    this.log('Executing Stage 4 script...');
    const execResult = this.executeScript(scriptPath, 'PCA');

    if (!execResult.success) {
      console.error(`Stage 4 execution failed: ${execResult.error}`);
      return { proceed: false, error: execResult.error };
    }

    const stage4Output = parseStage4Output(absOutputDir);

    // Agent review (Stats Agent for PC selection)
    if (this.forceAutomation) {
      this.log('FORCE AUTOMATION - Using default PCs 1-20');
      const pcSelection = { min_pc: 1, max_pc: 20 };
      this.stageResults.stage4 = { output: stage4Output, pcSelection, automated: true };
      return { proceed: true, pcSelection };
    }

    this.log('Calling ALL 3 AGENTS for PC selection...');
    const reviewResult = await this.coordinator.reviewScRNAStage4(stage4Output, this.dataInfo);

    // Check if user decision needed (agents disagree significantly on PC range)
    let userDecision = null;
    let pcSelection = { min_pc: 1, max_pc: 20 };  // Default

    const needsUserInput = reviewResult.consensus.confidence < 0.7 ||
                           reviewResult.consensus.decision.toLowerCase().includes('user_decision') ||
                           hasPCRangeDisagreement(reviewResult.agentResponses);

    if (needsUserInput) {
      this.log('PC range disagreement - escalating to user...');
      userDecision = await handleScRNAStage4UserDecision(reviewResult, stage4Output);

      if (!userDecision.proceed) {
        console.log('User aborted analysis');
        return { proceed: false, error: 'User aborted' };
      }

      pcSelection = userDecision.pcSelection;
    } else {
      // Use agent recommendation
      if (reviewResult.decision === 'SELECT_PC_RANGE') {
        pcSelection = { min_pc: reviewResult.min_pc, max_pc: reviewResult.max_pc };
      }
    }

    // Log agent decision (Stage 4 checkpoint)
    if (reviewResult.agentResponses && reviewResult.consensus) {
      this.logger.logStageDecision(
        4, 'PCA + PC Selection',
        this.dataInfo, stage4Output,
        reviewResult.agentResponses,
        reviewResult.consensus,
        { proceed: true, pcSelection, userDecision: userDecision || null }
      );
    }

    const finalProceed = reviewResult.decision !== 'STOP_AND_REVIEW';

    this.stageResults.stage4 = { output: stage4Output, review: reviewResult, pcSelection, proceed: finalProceed };

    console.log('');
    console.log(`Stage 4 Result: ${reviewResult.decision} (Using PCs ${pcSelection.min_pc}-${pcSelection.max_pc})`);
    return { proceed: finalProceed, output: stage4Output, pcSelection };
  }

  /**
   * Run Stage 5: Clustering + Markers (AGENT CHECKPOINT 3)
   */
  async runStage5() {
    console.log('');
    console.log('='.repeat(60));
    console.log('  STAGE 5: Clustering + Markers');
    console.log('='.repeat(60));
    console.log('');

    const absOutputDir = path.resolve(this.outputDir);
    const scriptPath = path.resolve(absOutputDir, `stage_5_scrna_cluster_markers.sh`);
    const pcSelection = this.stageResults.stage4.pcSelection;

    this.log('Generating Stage 5 script...');
    generateStage5Script(this.stageResults.stage4.output, pcSelection, this.config, scriptPath);

    if (this.dryRun) {
      return { proceed: true, dryRun: true };
    }

    this.log('Executing Stage 5 script...');
    const execResult = this.executeScript(scriptPath, 'Clustering + Markers', 7200000);  // 2 hour timeout

    if (!execResult.success) {
      console.error(`Stage 5 execution failed: ${execResult.error}`);
      return { proceed: false, error: execResult.error };
    }

    const stage5Output = parseStage5Output(absOutputDir);

    // Agent review (Pipeline Agent for clustering validation)
    if (this.forceAutomation) {
      this.log('FORCE AUTOMATION - Accepting clustering');
      this.stageResults.stage5 = { output: stage5Output, proceed: true, automated: true };
      return { proceed: true, output: stage5Output };
    }

    this.log('Calling ALL 3 AGENTS for clustering validation...');
    const reviewResult = await this.coordinator.reviewScRNAStage5(stage5Output, this.dataInfo);

    // Check if user decision needed (agents flag suspicious or disagree)
    let userDecision = null;
    let finalProceed = reviewResult.decision === 'ACCEPT_CLUSTERING';

    // FIX: Don't escalate if all 3 agents unanimously agreed to accept
    const unanimousAccept = reviewResult.consensus.votes.approve === 3;

    const needsUserInput = !unanimousAccept && (
                           reviewResult.consensus.confidence < 0.7 ||
                           reviewResult.consensus.decision.toLowerCase().includes('user_decision') ||
                           reviewResult.decision === 'FLAG_SUSPICIOUS' ||
                           reviewResult.decision === 'ADJUST_RESOLUTION' ||
                           reviewResult.consensus.votes.reject >= 2
                           );

    if (needsUserInput) {
      this.log('Clustering concerns detected - escalating to user...');
      userDecision = await handleScRNAStage5UserDecision(reviewResult, stage5Output);

      if (!userDecision.proceed) {
        console.log('User aborted or requested re-clustering');
        // Store action for potential cell cycle correction
        this.stageResults.stage5 = {
          output: stage5Output,
          review: reviewResult,
          proceed: false,
          userAction: userDecision.action || 'abort'
        };
        return { proceed: false, error: 'User aborted or requested re-analysis' };
      }

      finalProceed = userDecision.proceed;
    }

    // Log agent decision (Stage 5 checkpoint)
    if (reviewResult.agentResponses && reviewResult.consensus) {
      this.logger.logStageDecision(
        5, 'Clustering + Markers',
        this.dataInfo, stage5Output,
        reviewResult.agentResponses,
        reviewResult.consensus,
        {
          proceed: finalProceed,
          clusteringDecision: reviewResult.decision,
          userDecision: userDecision || null,
          userAction: userDecision?.action || null
        }
      );
    }

    this.stageResults.stage5 = { output: stage5Output, review: reviewResult, proceed: finalProceed };

    console.log('');
    console.log(`Stage 5 Result: ${reviewResult.decision}`);
    return { proceed: finalProceed, output: stage5Output, review: reviewResult };
  }

  /**
   * Main execution entry point
   */
  async run() {
    console.log('');
    console.log('#'.repeat(60));
    console.log('  scRNA-seq Seurat Pipeline (GeneExpert)');
    console.log('#'.repeat(60));
    console.log('');
    console.log(`Input:  ${this.inputDir}`);
    console.log(`Output: ${this.outputDir}`);
    console.log(`Organism: ${this.organism}`);
    console.log('');

    try {
      // Initialize
      this.initialize();

      // Stage 1: Load + QC (no agent)
      const stage1Result = await this.runStage1();
      if (!stage1Result.proceed) {
        throw new Error('Stage 1 failed');
      }

      // Stage 2: QC Filtering (Stats Agent)
      const stage2Result = await this.runStage2();
      if (!stage2Result.proceed) {
        throw new Error('Stage 2 rejected by agent');
      }

      // Stage 3: Normalize + HVG (no agent)
      const stage3Result = await this.runStage3();
      if (!stage3Result.proceed) {
        throw new Error('Stage 3 failed');
      }

      // Stage 4: PCA (Stats Agent for PC selection)
      const stage4Result = await this.runStage4();
      if (!stage4Result.proceed) {
        throw new Error('Stage 4 rejected by agent');
      }

      // Stage 5: Clustering + Markers (Pipeline Agent)
      const stage5Result = await this.runStage5();
      if (!stage5Result.proceed) {
        throw new Error('Stage 5 rejected by agent');
      }

      // Success!
      console.log('');
      console.log('#'.repeat(60));
      console.log('  Pipeline Complete!');
      console.log('#'.repeat(60));
      console.log('');
      console.log('Final Outputs:');
      console.log(`  - Seurat object: ${this.outputDir}/stage5_cluster_markers/seurat_stage5_clustered.rds`);
      console.log(`  - Markers: ${this.outputDir}/stage5_cluster_markers/markers_stage5.csv`);
      console.log('');

    } catch (error) {
      console.error('');
      console.error('#'.repeat(60));
      console.error('  Pipeline FAILED');
      console.error(`  Reason: ${error.message}`);
      console.error('#'.repeat(60));
      console.error('');
      throw error;
    }
  }
}

/**
 * ============================================================================
 * HELPER FUNCTIONS FOR USER DECISION ESCALATION
 * ============================================================================
 */

/**
 * Check if agents disagree significantly on QC thresholds
 */
function hasThresholdDisagreement(agentResponses) {
  if (!agentResponses.gpt5_2 || !agentResponses.claude || !agentResponses.gemini) {
    return false;
  }

  const extractThreshold = (text, pattern) => {
    const match = text.match(pattern);
    return match ? parseInt(match[1]) : null;
  };

  const gptMin = extractThreshold(agentResponses.gpt5_2.content, /nFeature[_\s]*min[:\s]*(\d+)/i);
  const gptMax = extractThreshold(agentResponses.gpt5_2.content, /nFeature[_\s]*max[:\s]*(\d+)/i);
  const claudeMin = extractThreshold(agentResponses.claude.content, /nFeature[_\s]*min[:\s]*(\d+)/i);
  const claudeMax = extractThreshold(agentResponses.claude.content, /nFeature[_\s]*max[:\s]*(\d+)/i);
  const geminiMin = extractThreshold(agentResponses.gemini.content, /nFeature[_\s]*min[:\s]*(\d+)/i);
  const geminiMax = extractThreshold(agentResponses.gemini.content, /nFeature[_\s]*max[:\s]*(\d+)/i);

  // Check if thresholds vary by >50% or >2000 absolute
  if (gptMin && claudeMin && geminiMin) {
    const minValues = [gptMin, claudeMin, geminiMin];
    const maxMin = Math.max(...minValues);
    const minMin = Math.min(...minValues);
    if (maxMin - minMin > 2000 || maxMin / minMin > 1.5) {
      return true;
    }
  }

  if (gptMax && claudeMax && geminiMax) {
    const maxValues = [gptMax, claudeMax, geminiMax];
    const maxMax = Math.max(...maxValues);
    const minMax = Math.min(...maxValues);
    if (maxMax - minMax > 5000 || maxMax / minMax > 1.5) {
      return true;
    }
  }

  return false;
}

/**
 * Check if agents disagree significantly on PC range
 */
function hasPCRangeDisagreement(agentResponses) {
  if (!agentResponses.gpt5_2 || !agentResponses.claude || !agentResponses.gemini) {
    return false;
  }

  const extractMaxPC = (text) => {
    const match = text.match(/max[_\s]*pc[:\s]*(\d+)/i);
    return match ? parseInt(match[1]) : null;
  };

  const gptMax = extractMaxPC(agentResponses.gpt5_2.content);
  const claudeMax = extractMaxPC(agentResponses.claude.content);
  const geminiMax = extractMaxPC(agentResponses.gemini.content);

  if (gptMax && claudeMax && geminiMax) {
    const pcValues = [gptMax, claudeMax, geminiMax];
    const maxPC = Math.max(...pcValues);
    const minPC = Math.min(...pcValues);

    // If range differs by >10 PCs, escalate
    return (maxPC - minPC) > 10;
  }

  return false;
}

export default ScRNAExecutor;
