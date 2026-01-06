/**
 * Pipeline Planner - Plans analysis steps based on input data type
 *
 * Creates a step-by-step execution plan that the Coordinator will orchestrate
 */

import dotenv from 'dotenv';
dotenv.config();

// Get scripts path from environment or use default
const SCRIPTS_PATH = process.env.SCRIPTS_PATH || '/users/ha00014/Halimas_projects/multi_llm_mcp/bio_informatics/scripts';

/**
 * Plan pipeline steps based on detected data
 * @param {Object} dataInfo - Output from detectInputData()
 * @param {Object} config - Analysis configuration
 * @returns {Array} - Ordered list of pipeline steps
 */
export function planPipeline(dataInfo, config) {
  console.log('[Pipeline Planner] Planning analysis steps...');

  const steps = [];

  // Plan based on input data type
  switch (dataInfo.type) {
    case 'fastq':
      steps.push(...planFromFASTQ(dataInfo, config));
      break;

    case 'bam':
      steps.push(...planFromBAM(dataInfo, config));
      break;

    case 'counts':
      steps.push(...planFromCounts(dataInfo, config));
      break;

    default:
      throw new Error(`Unknown data type: ${dataInfo.type}`);
  }

  console.log(`[Pipeline Planner] Planned ${steps.length} steps:`);
  console.log('[Pipeline Planner] ');
  console.log('[Pipeline Planner] 📋 STANDARD BULK RNA-SEQ PIPELINE (Lab Protocol):');
  steps.forEach((step, i) => {
    console.log(`   ${i + 1}. ${step.name} ${step.requiresDebate ? '🤔' : ''}`);
  });
  console.log('');

  return steps;
}

/**
 * Plan pipeline starting from FASTQ files
 * Agent will decide automation vs adaptation approach
 */
function planFromFASTQ(dataInfo, config) {

  return [
    // Step 1: QC on Raw Reads (FastQC)
    {
      name: 'FastQC',
      description: 'Quality control on raw FASTQ files',
      tool: 'fastqc',
      script: 'fastqc', // System command
      inputs: dataInfo.files,
      requiresDebate: false,
      agent: 'pipeline',
      outputs: ['*_fastqc.html', '*_fastqc.zip']
    },

    // Step 2: Alignment (fastq2bam - NO .sh extension!)
    {
      name: 'Alignment (fastq2bam)',
      description: 'Align FASTQ to reference genome',
      tool: 'fastq2bam',
      script: `${SCRIPTS_PATH}/fastq2bam`,
      inputs: dataInfo.files,
      requiresDebate: false,
      agent: 'pipeline',
      params: {
        genome: config.organism === 'mouse' ? 'mm10' : 'hg38',
        readType: dataInfo.pairedEnd ? 'paired' : 'single'
      },
      outputs: ['*.bam']
    },

    // Step 2: Feature Counts (featurecounts.R)
    {
      name: 'Feature Counts',
      description: 'Count reads per gene → output_name.count.txt',
      tool: 'featurecounts',
      script: `${SCRIPTS_PATH}/featurecounts.R`,
      inputs: 'bam_files',
      requiresDebate: false,
      agent: 'pipeline',
      params: {
        genome: config.organism === 'mouse' ? 'mm10' : 'hg38',
        outputName: config.comparison || 'analysis'
      },
      outputs: ['*.count.txt']
    },

    // Step 3: Filter Bad IDs (FilterIDs.R)
    {
      name: 'Filter Bad IDs',
      description: 'Remove bad genes → .count.txt.badIds.txt',
      tool: 'filterIDS',
      script: `${SCRIPTS_PATH}/filterIDS.R`,
      inputs: 'count_txt',
      requiresDebate: false,
      agent: 'pipeline',
      outputs: ['*.count.txt.badIds.txt']
    },

    // ========================================================================
    // IMPORTANT: Dual-path workflow from here:
    //
    // Path 1 (VISUALIZATION): badIds → RPKM → entrz → QC plots (PCA, MDS)
    // Path 2 (STATISTICS):    badIds → edgeR (uses raw counts + internal TMM)
    //
    // RPKM is ONLY for visualization. edgeR uses RAW filtered counts!
    // ========================================================================

    // Step 4: RPKM Normalization (RPKM.R) - For visualization only, NOT for edgeR
    {
      name: 'RPKM Normalization',
      description: 'Calculate RPKM → outRPKM.txt (for visualization)',
      tool: 'rpkm',
      script: `${SCRIPTS_PATH}/RPKM.R`,
      inputs: 'badIds_txt',  // Uses filtered counts from filterIDS.R
      requiresDebate: false,
      agent: 'pipeline',
      outputs: ['outRPKM.txt'],
      note: 'RPKM is for PCA/plots only. edgeR uses raw filtered counts.'
    },

    // Step 6: Add Gene Annotations (entrz.R)
    {
      name: 'Add Gene Symbols',
      description: 'Add Entrez IDs and gene symbols → outentrz.txt',
      tool: 'entrz',
      script: `${SCRIPTS_PATH}/entrz.R`,
      inputs: 'outRPKM.txt',
      requiresDebate: false,
      agent: 'pipeline',
      params: {
        genome: config.organism === 'mouse' ? 'mm10' : 'hg38'
      },
      outputs: ['outentrz.txt']
    },

    // QC Plots Step - Agent can read these and assess quality
    {
      name: 'Generate QC Plots',
      description: 'PCA, MDS, density plots for quality assessment',
      tool: 'qc_plots',
      script: `${SCRIPTS_PATH}/qc_plots.R`, // May need to create this
      inputs: 'outentrz.txt',
      requiresDebate: false,
      agent: 'pipeline',
      outputs: ['PCA.pdf', 'MDS.pdf', 'density.pdf'],
      note: 'Agent can read these plots via MCP to check for outliers/batch effects'
    },

    // NOTE: Agent decides if intervention needed after reviewing QC
    // Agent can write custom script if problems detected

    // Step 7: Differential Expression (simpleEdger3.R)
    {
      name: 'DE Analysis (edgeR)',
      description: 'edgeR differential expression analysis (uses RAW filtered counts)',
      tool: 'edger',
      script: `${SCRIPTS_PATH}/simpleEdger3.R`,
      inputs: 'badIds_txt',  // Uses RAW filtered counts from filterIDS.R, NOT RPKM!
      requiresDebate: false,
      agent: 'pipeline',
      params: {
        genome: config.organism === 'mouse' ? 'mm10' : 'hg38'
      },
      outputs: ['edger_results.txt'],
      note: 'edgeR expects raw counts and does its own TMM normalization internally'
    },

    // Step 8: Excel Conversion (edger3xl.R)
    {
      name: 'Convert to Excel',
      description: 'Create Excel file with FDR/logFC formulas',
      tool: 'edger_excel',
      script: `${SCRIPTS_PATH}/edger3xl.R`,
      inputs: 'badIds_txt',  // Same input as edgeR
      requiresDebate: false,
      agent: 'pipeline',
      outputs: ['*.xlsx']
    },

    // Step 9: Merge entrz with Excel (manual step - can add script later)
    {
      name: 'Merge Annotations',
      description: 'Merge outentrz.txt with Excel file',
      tool: 'merge',
      inputs: ['outentrz.txt', 'excel_file'],
      requiresDebate: false,
      agent: 'pipeline',
      note: 'Combines gene symbols with DE results'
    }
  ];
}

/**
 * Plan pipeline starting from BAM files (skip alignment)
 */
function planFromBAM(dataInfo, config) {
  return [
    // Skip FastQC and alignment, start from quantification
    {
      name: 'featureCounts',
      description: 'Count reads per gene from BAM files',
      tool: 'featurecounts',
      inputs: dataInfo.samples,
      requiresDebate: false,
      agent: 'pipeline',
      requiresGTF: true
    },

    {
      name: 'Filter Low Counts',
      description: 'Remove lowly expressed genes',
      tool: 'filterIDS',
      inputs: 'count_matrix',
      requiresDebate: false,
      agent: 'pipeline'
    },

    {
      name: 'RPKM Normalization',
      description: 'Normalize expression values',
      tool: 'rpkm',
      inputs: 'filtered_counts',
      requiresDebate: false,
      agent: 'pipeline'
    },

    // DECISION POINT: QC Review
    {
      name: 'QC Review',
      description: 'Multi-agent review of QC plots',
      tool: 'qc_plots',
      inputs: 'normalized_counts',
      requiresDebate: true,
      agent: 'multi',
      decisionType: 'qc_review'
    },

    // DECISION POINT: Threshold Selection
    {
      name: 'Threshold Selection',
      description: 'Multi-agent debate on thresholds',
      tool: null,
      inputs: { sampleSize: dataInfo.samples.length },
      requiresDebate: true,
      agent: 'multi',
      decisionType: 'threshold'
    },

    {
      name: `DE Analysis (${config.deTool})`,
      description: 'Differential expression analysis',
      tool: config.deTool,
      inputs: 'filtered_counts',
      requiresDebate: false,
      agent: 'pipeline',
      usesThresholds: true
    },

    {
      name: 'Generate Plots',
      description: 'Visualization',
      tool: 'visualization',
      inputs: 'de_results',
      requiresDebate: false,
      agent: 'pipeline'
    }
  ];
}

/**
 * Plan pipeline starting from count matrix (skip alignment and counting)
 */
function planFromCounts(dataInfo, config) {
  return [
    {
      name: 'Filter Low Counts',
      description: 'Remove lowly expressed genes',
      tool: 'filterIDS',
      inputs: dataInfo.files[0],
      requiresDebate: false,
      agent: 'pipeline'
    },

    {
      name: 'RPKM Normalization',
      description: 'Normalize expression values',
      tool: 'rpkm',
      inputs: 'filtered_counts',
      requiresDebate: false,
      agent: 'pipeline'
    },

    // DECISION POINT: QC Review
    {
      name: 'QC Review',
      description: 'Multi-agent review of QC plots',
      tool: 'qc_plots',
      inputs: 'normalized_counts',
      requiresDebate: true,
      agent: 'multi',
      decisionType: 'qc_review'
    },

    // DECISION POINT: Threshold Selection
    {
      name: 'Threshold Selection',
      description: 'Multi-agent debate on thresholds',
      tool: null,
      inputs: { sampleSize: dataInfo.samples.length },
      requiresDebate: true,
      agent: 'multi',
      decisionType: 'threshold'
    },

    {
      name: `DE Analysis (${config.deTool})`,
      description: 'Differential expression analysis',
      tool: config.deTool,
      inputs: 'filtered_counts',
      requiresDebate: false,
      agent: 'pipeline',
      usesThresholds: true
    },

    {
      name: 'Generate Plots',
      description: 'Visualization',
      tool: 'visualization',
      inputs: 'de_results',
      requiresDebate: false,
      agent: 'pipeline'
    }
  ];
}
