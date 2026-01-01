/**
 * Pipeline Planner - Plans analysis steps based on input data type
 *
 * Creates a step-by-step execution plan that the Coordinator will orchestrate
 */

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
  steps.forEach((step, i) => {
    console.log(`   ${i + 1}. ${step.name} ${step.requiresDebate ? '🤔' : ''}`);
  });
  console.log('');

  return steps;
}

/**
 * Plan pipeline starting from FASTQ files
 */
function planFromFASTQ(dataInfo, config) {
  return [
    // Step 1: Quality Control (Optional - can skip if user wants)
    {
      name: 'FastQC',
      description: 'Quality control on raw reads',
      tool: 'fastqc',
      script: 'fastqc', // System tool
      inputs: dataInfo.files,
      requiresDebate: false,
      agent: 'pipeline',
      optional: true
    },

    // Step 2: Alignment (using fastq2bam from /data/scripts)
    {
      name: 'Alignment (Subread)',
      description: 'Align reads to reference genome using Subread',
      tool: 'fastq2bam',
      script: '/data/scripts/fastq2bam', // User's script!
      inputs: dataInfo.samples,
      requiresDebate: false,
      agent: 'pipeline',
      params: {
        genome: config.organism === 'mouse' ? 'mm10' : 'hg38',
        readType: dataInfo.pairedEnd ? 'paired' : 'single'
      }
    },

    // Step 3: Quantification (using featurecounts.R from /data/scripts)
    {
      name: 'featureCounts',
      description: 'Count reads per gene using featureCounts',
      tool: 'featurecounts',
      script: '/data/scripts/featurecounts.R', // User's script!
      inputs: 'alignment_outputs', // BAM files from step 2
      requiresDebate: false,
      agent: 'pipeline',
      params: {
        annotation: config.organism === 'mouse' ? 'mm10' : 'hg38'
      }
    },

    // Step 4: Filtering (using filterIDS.R from /data/scripts)
    {
      name: 'Filter Low Counts',
      description: 'Remove lowly expressed genes',
      tool: 'filterIDS',
      script: '/data/scripts/filterIDS.R', // User's script!
      inputs: 'count_matrix',
      requiresDebate: false,
      agent: 'pipeline'
    },

    // Step 5: Normalization (using RPKM.R from /data/scripts)
    {
      name: 'RPKM Normalization',
      description: 'Normalize expression values',
      tool: 'rpkm',
      script: '/data/scripts/RPKM.R', // User's script!
      inputs: 'filtered_counts',
      requiresDebate: false,
      agent: 'pipeline'
    },

    // DECISION POINT: QC Review
    {
      name: 'QC Review',
      description: 'Multi-agent review of PCA/MDS plots for outliers and batch effects',
      tool: 'qc_plots',
      inputs: 'normalized_counts',
      requiresDebate: true, // Multi-agent debate
      agent: 'multi',
      decisionType: 'qc_review'
    },

    // DECISION POINT: Threshold Selection
    {
      name: 'Threshold Selection',
      description: 'Multi-agent debate on FDR and logFC thresholds',
      tool: null, // Just a decision
      inputs: { sampleSize: dataInfo.samples.length },
      requiresDebate: true, // Multi-agent debate
      agent: 'multi',
      decisionType: 'threshold'
    },

    // Step 6: Differential Expression (using simpleEdger3.R from /data/scripts)
    {
      name: 'DE Analysis (edgeR)',
      description: 'Identify differentially expressed genes using edgeR',
      tool: 'edger',
      script: '/data/scripts/simpleEdger3.R', // User's script!
      inputs: 'filtered_counts',
      requiresDebate: false,
      agent: 'pipeline',
      usesThresholds: true // Uses thresholds from previous decision
    },

    // Step 7: Visualization
    {
      name: 'Generate Plots',
      description: 'Volcano plots, MA plots, heatmaps',
      tool: 'visualization',
      inputs: 'de_results',
      requiresDebate: false,
      agent: 'pipeline'
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
