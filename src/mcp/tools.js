/**
 * GeneExpert MCP Tools
 * Wraps bioinformatics scripts from /data/scripts/ and /destiny/halima/dory/my_script/
 */

import { exec } from 'child_process';
import { promisify } from 'util';
import path from 'path';

const execAsync = promisify(exec);

// Script paths from environment or defaults
const SCRIPTS_PATH = process.env.SCRIPTS_PATH || '/data/scripts';
const CUSTOM_SCRIPTS_PATH = process.env.CUSTOM_SCRIPTS_PATH || '/destiny/halima/dory/my_script';

/**
 * Execute a script and return the result
 */
async function executeScript(scriptPath, args = [], options = {}) {
  try {
    const command = `${scriptPath} ${args.join(' ')}`;
    console.log(`Executing: ${command}`);

    const { stdout, stderr } = await execAsync(command, {
      maxBuffer: 10 * 1024 * 1024, // 10MB buffer for large outputs
      cwd: options.workingDir || process.cwd(),
      ...options
    });

    return {
      success: true,
      stdout: stdout.trim(),
      stderr: stderr.trim(),
      command
    };
  } catch (error) {
    return {
      success: false,
      error: error.message,
      stdout: error.stdout?.trim() || '',
      stderr: error.stderr?.trim() || '',
      command: scriptPath
    };
  }
}

/**
 * Tool definitions for MCP server
 */
export const tools = [
  // ============================================
  // PHASE 1: PRE-PROCESSING TOOLS
  // ============================================

  {
    name: 'run_fastqc',
    description: 'Run FastQC quality control on FASTQ files',
    inputSchema: {
      type: 'object',
      properties: {
        input_files: {
          type: 'array',
          items: { type: 'string' },
          description: 'Array of FASTQ file paths'
        },
        output_dir: {
          type: 'string',
          description: 'Output directory for QC reports'
        },
        threads: {
          type: 'number',
          description: 'Number of threads to use',
          default: 4
        }
      },
      required: ['input_files', 'output_dir']
    },
    handler: async (args) => {
      const files = args.input_files.join(' ');
      const threads = args.threads || 4;
      const result = await execAsync(
        `fastqc -t ${threads} -o ${args.output_dir} ${files}`
      );
      return {
        content: [{
          type: 'text',
          text: `FastQC completed for ${args.input_files.length} files\nOutput: ${args.output_dir}\n\n${result.stdout}`
        }]
      };
    }
  },

  {
    name: 'run_alignment',
    description: 'Align FASTQ files to reference genome using STAR/HISAT2',
    inputSchema: {
      type: 'object',
      properties: {
        fastq_r1: {
          type: 'string',
          description: 'Path to R1 FASTQ file (or single-end FASTQ)'
        },
        fastq_r2: {
          type: 'string',
          description: 'Path to R2 FASTQ file (paired-end only)'
        },
        genome: {
          type: 'string',
          description: 'Genome build (mm10, hg38, etc.)',
          default: 'mm10'
        },
        output_prefix: {
          type: 'string',
          description: 'Output BAM file prefix'
        },
        threads: {
          type: 'number',
          description: 'Number of threads',
          default: 8
        }
      },
      required: ['fastq_r1', 'genome', 'output_prefix']
    },
    handler: async (args) => {
      const script = path.join(SCRIPTS_PATH, 'fastq2bam');
      const scriptArgs = [
        args.fastq_r1,
        args.fastq_r2 || '',
        args.genome,
        args.output_prefix,
        args.threads || 8
      ].filter(Boolean);

      const result = await executeScript(script, scriptArgs);

      return {
        content: [{
          type: 'text',
          text: result.success
            ? `Alignment completed!\nBAM file: ${args.output_prefix}.bam\n\n${result.stdout}`
            : `Alignment failed: ${result.error}\n${result.stderr}`
        }]
      };
    }
  },

  {
    name: 'run_featurecounts',
    description: 'Quantify gene expression from BAM files using featureCounts',
    inputSchema: {
      type: 'object',
      properties: {
        bam_files: {
          type: 'array',
          items: { type: 'string' },
          description: 'Array of BAM file paths'
        },
        genome: {
          type: 'string',
          description: 'Genome annotation (mm10, hg38)',
          default: 'mm10'
        },
        output_file: {
          type: 'string',
          description: 'Output count matrix file'
        },
        paired: {
          type: 'boolean',
          description: 'Paired-end reads',
          default: true
        }
      },
      required: ['bam_files', 'genome', 'output_file']
    },
    handler: async (args) => {
      const script = path.join(CUSTOM_SCRIPTS_PATH, 'featurecounts_edited.R');
      const scriptArgs = [
        args.genome,
        args.output_file,
        args.paired ? 'paired' : 'single',
        ...args.bam_files
      ];

      const result = await executeScript(`Rscript ${script}`, scriptArgs);

      return {
        content: [{
          type: 'text',
          text: result.success
            ? `FeatureCounts completed!\nCount matrix: ${args.output_file}\nGenes quantified: ${result.stdout.match(/\d+/)?.[0] || 'N/A'}\n\n${result.stdout}`
            : `FeatureCounts failed: ${result.error}\n${result.stderr}`
        }]
      };
    }
  },

  // ============================================
  // PHASE 2: NORMALIZATION & QC TOOLS
  // ============================================

  {
    name: 'run_filter',
    description: 'Filter out lowly expressed genes (< 10 counts)',
    inputSchema: {
      type: 'object',
      properties: {
        count_matrix: {
          type: 'string',
          description: 'Input count matrix file'
        },
        output_file: {
          type: 'string',
          description: 'Output filtered matrix file'
        },
        min_count: {
          type: 'number',
          description: 'Minimum count threshold',
          default: 10
        }
      },
      required: ['count_matrix', 'output_file']
    },
    handler: async (args) => {
      const script = path.join(SCRIPTS_PATH, 'filterIDS.R');
      const scriptArgs = [
        args.count_matrix,
        args.output_file,
        args.min_count || 10
      ];

      const result = await executeScript(`Rscript ${script}`, scriptArgs);

      return {
        content: [{
          type: 'text',
          text: result.success
            ? `Filtering completed!\nFiltered matrix: ${args.output_file}\n\n${result.stdout}`
            : `Filtering failed: ${result.error}\n${result.stderr}`
        }]
      };
    }
  },

  {
    name: 'run_rpkm',
    description: 'Normalize counts to RPKM (Reads Per Kilobase Million)',
    inputSchema: {
      type: 'object',
      properties: {
        count_matrix: {
          type: 'string',
          description: 'Input count matrix file'
        },
        output_file: {
          type: 'string',
          description: 'Output RPKM matrix file'
        },
        genome: {
          type: 'string',
          description: 'Genome build (for gene lengths)',
          default: 'mm10'
        }
      },
      required: ['count_matrix', 'output_file']
    },
    handler: async (args) => {
      const script = path.join(SCRIPTS_PATH, 'RPKM.R');
      const scriptArgs = [
        args.count_matrix,
        args.output_file,
        args.genome || 'mm10'
      ];

      const result = await executeScript(`Rscript ${script}`, scriptArgs);

      return {
        content: [{
          type: 'text',
          text: result.success
            ? `RPKM normalization completed!\nNormalized matrix: ${args.output_file}\n\n${result.stdout}`
            : `RPKM normalization failed: ${result.error}\n${result.stderr}`
        }]
      };
    }
  },

  {
    name: 'run_annotation',
    description: 'Add gene symbols and descriptions to matrix',
    inputSchema: {
      type: 'object',
      properties: {
        matrix_file: {
          type: 'string',
          description: 'Input matrix file'
        },
        output_file: {
          type: 'string',
          description: 'Output annotated matrix file'
        },
        genome: {
          type: 'string',
          description: 'Genome (mm10, hg38)',
          default: 'mm10'
        }
      },
      required: ['matrix_file', 'output_file']
    },
    handler: async (args) => {
      const script = path.join(SCRIPTS_PATH, 'entrz.R');
      const scriptArgs = [
        args.matrix_file,
        args.output_file,
        args.genome || 'mm10'
      ];

      const result = await executeScript(`Rscript ${script}`, scriptArgs);

      return {
        content: [{
          type: 'text',
          text: result.success
            ? `Annotation completed!\nAnnotated file: ${args.output_file}\n\n${result.stdout}`
            : `Annotation failed: ${result.error}\n${result.stderr}`
        }]
      };
    }
  },

  {
    name: 'run_qc_plots',
    description: 'Generate PCA, MDS, and density plots for quality control (CRITICAL before DE analysis)',
    inputSchema: {
      type: 'object',
      properties: {
        rpkm_matrix: {
          type: 'string',
          description: 'RPKM-normalized matrix file'
        },
        output_dir: {
          type: 'string',
          description: 'Output directory for plots'
        },
        sample_groups: {
          type: 'array',
          items: { type: 'string' },
          description: 'Sample group labels (same order as matrix columns)'
        }
      },
      required: ['rpkm_matrix', 'output_dir']
    },
    handler: async (args) => {
      // PCA plot
      const pcaScript = path.join(SCRIPTS_PATH, 'PCA.R');
      const mdsScript = path.join(CUSTOM_SCRIPTS_PATH, 'MDS_mdfy.R');
      const densityScript = path.join(SCRIPTS_PATH, 'densityplot.R');

      const pcaResult = await executeScript(`Rscript ${pcaScript}`, [
        args.rpkm_matrix,
        path.join(args.output_dir, 'pca_plot.pdf')
      ]);

      const mdsResult = await executeScript(`Rscript ${mdsScript}`, [
        args.rpkm_matrix,
        path.join(args.output_dir, 'mds_plot.pdf')
      ]);

      const densityResult = await executeScript(`Rscript ${densityScript}`, [
        args.rpkm_matrix,
        path.join(args.output_dir, 'density_plot.pdf')
      ]);

      return {
        content: [{
          type: 'text',
          text: `QC plots generated in ${args.output_dir}/\n\n` +
                `PCA: ${pcaResult.success ? '✓' : '✗'}\n` +
                `MDS: ${mdsResult.success ? '✓' : '✗'}\n` +
                `Density: ${densityResult.success ? '✓' : '✗'}\n\n` +
                `⚠️ CRITICAL: Review these plots before proceeding to DE analysis!\n` +
                `Multi-agent review recommended for outlier detection.`
        }]
      };
    }
  },

  // ============================================
  // PHASE 3: DIFFERENTIAL EXPRESSION TOOLS
  // ============================================

  {
    name: 'run_edger',
    description: 'Run differential expression analysis using edgeR',
    inputSchema: {
      type: 'object',
      properties: {
        count_matrix: {
          type: 'string',
          description: 'Filtered count matrix file'
        },
        groups: {
          type: 'array',
          items: { type: 'string' },
          description: 'Sample group assignments (e.g., ["control", "control", "treatment", "treatment"])'
        },
        comparison: {
          type: 'string',
          description: 'Comparison to make (e.g., "treatment-control")'
        },
        output_file: {
          type: 'string',
          description: 'Output DEG results file'
        },
        fdr_cutoff: {
          type: 'number',
          description: 'FDR threshold',
          default: 0.05
        },
        logfc_cutoff: {
          type: 'number',
          description: 'Log fold-change threshold',
          default: 1.0
        }
      },
      required: ['count_matrix', 'groups', 'comparison', 'output_file']
    },
    handler: async (args) => {
      const script = path.join(SCRIPTS_PATH, 'simpleEdger3.R');
      const scriptArgs = [
        args.count_matrix,
        args.groups.join(','),
        args.comparison,
        args.output_file,
        args.fdr_cutoff || 0.05,
        args.logfc_cutoff || 1.0
      ];

      const result = await executeScript(`Rscript ${script}`, scriptArgs);

      return {
        content: [{
          type: 'text',
          text: result.success
            ? `Differential expression analysis completed!\n` +
              `Results: ${args.output_file}\n` +
              `FDR < ${args.fdr_cutoff}, |logFC| > ${args.logfc_cutoff}\n\n` +
              `${result.stdout}`
            : `edgeR analysis failed: ${result.error}\n${result.stderr}`
        }]
      };
    }
  },

  {
    name: 'export_to_excel',
    description: 'Export comprehensive results (DEGs + annotations + RPKM) to Excel',
    inputSchema: {
      type: 'object',
      properties: {
        deg_results: {
          type: 'string',
          description: 'DEG results file from edgeR'
        },
        rpkm_matrix: {
          type: 'string',
          description: 'RPKM normalized matrix'
        },
        annotations: {
          type: 'string',
          description: 'Gene annotations file'
        },
        output_excel: {
          type: 'string',
          description: 'Output Excel file path'
        }
      },
      required: ['deg_results', 'rpkm_matrix', 'output_excel']
    },
    handler: async (args) => {
      const script = path.join(SCRIPTS_PATH, 'edger3xl.R');
      const scriptArgs = [
        args.deg_results,
        args.rpkm_matrix,
        args.annotations || '',
        args.output_excel
      ].filter(Boolean);

      const result = await executeScript(`Rscript ${script}`, scriptArgs);

      return {
        content: [{
          type: 'text',
          text: result.success
            ? `Excel export completed!\nFile: ${args.output_excel}\n\n${result.stdout}`
            : `Excel export failed: ${result.error}\n${result.stderr}`
        }]
      };
    }
  },

  // ============================================
  // PHASE 4: VISUALIZATION TOOLS
  // ============================================

  {
    name: 'run_volcano_plot',
    description: 'Generate volcano plot of differential expression results',
    inputSchema: {
      type: 'object',
      properties: {
        deg_results: {
          type: 'string',
          description: 'DEG results file'
        },
        output_plot: {
          type: 'string',
          description: 'Output plot file path (.pdf or .png)'
        },
        fdr_line: {
          type: 'number',
          description: 'FDR significance line',
          default: 0.05
        },
        logfc_line: {
          type: 'number',
          description: 'LogFC threshold line',
          default: 1.0
        }
      },
      required: ['deg_results', 'output_plot']
    },
    handler: async (args) => {
      const script = path.join(CUSTOM_SCRIPTS_PATH, 'volcano.R');
      const scriptArgs = [
        args.deg_results,
        args.output_plot,
        args.fdr_line || 0.05,
        args.logfc_line || 1.0
      ];

      const result = await executeScript(`Rscript ${script}`, scriptArgs);

      return {
        content: [{
          type: 'text',
          text: result.success
            ? `Volcano plot created: ${args.output_plot}`
            : `Volcano plot failed: ${result.error}\n${result.stderr}`
        }]
      };
    }
  },

  {
    name: 'run_ma_plot',
    description: 'Generate MA plot (mean expression vs log fold-change)',
    inputSchema: {
      type: 'object',
      properties: {
        deg_results: {
          type: 'string',
          description: 'DEG results file'
        },
        output_plot: {
          type: 'string',
          description: 'Output plot file path'
        }
      },
      required: ['deg_results', 'output_plot']
    },
    handler: async (args) => {
      const script = path.join(CUSTOM_SCRIPTS_PATH, 'maplot.R');
      const scriptArgs = [args.deg_results, args.output_plot];

      const result = await executeScript(`Rscript ${script}`, scriptArgs);

      return {
        content: [{
          type: 'text',
          text: result.success
            ? `MA plot created: ${args.output_plot}`
            : `MA plot failed: ${result.error}\n${result.stderr}`
        }]
      };
    }
  },

  {
    name: 'run_venn',
    description: 'Create Venn diagram comparing DEG lists',
    inputSchema: {
      type: 'object',
      properties: {
        deg_lists: {
          type: 'array',
          items: { type: 'string' },
          description: 'Array of DEG result files to compare (2-3 files)'
        },
        list_names: {
          type: 'array',
          items: { type: 'string' },
          description: 'Names for each list'
        },
        output_plot: {
          type: 'string',
          description: 'Output Venn diagram file'
        }
      },
      required: ['deg_lists', 'list_names', 'output_plot']
    },
    handler: async (args) => {
      const script = path.join(SCRIPTS_PATH, 'venn_diagram.R');
      const scriptArgs = [
        ...args.deg_lists,
        ...args.list_names,
        args.output_plot
      ];

      const result = await executeScript(`Rscript ${script}`, scriptArgs);

      return {
        content: [{
          type: 'text',
          text: result.success
            ? `Venn diagram created: ${args.output_plot}`
            : `Venn diagram failed: ${result.error}\n${result.stderr}`
        }]
      };
    }
  },

  // ============================================
  // ATAC-seq SPECIFIC TOOLS
  // ============================================

  {
    name: 'run_atac_alignment',
    description: 'Run ATAC-seq alignment pipeline',
    inputSchema: {
      type: 'object',
      properties: {
        fastq_r1: {
          type: 'string',
          description: 'R1 FASTQ file'
        },
        fastq_r2: {
          type: 'string',
          description: 'R2 FASTQ file'
        },
        genome: {
          type: 'string',
          description: 'Genome build',
          default: 'mm10'
        },
        output_prefix: {
          type: 'string',
          description: 'Output prefix'
        }
      },
      required: ['fastq_r1', 'fastq_r2', 'output_prefix']
    },
    handler: async (args) => {
      const script = path.join(SCRIPTS_PATH, 'atac_seq_alignment');
      const scriptArgs = [
        args.fastq_r1,
        args.fastq_r2,
        args.genome || 'mm10',
        args.output_prefix
      ];

      const result = await executeScript(script, scriptArgs);

      return {
        content: [{
          type: 'text',
          text: result.success
            ? `ATAC-seq alignment completed: ${args.output_prefix}`
            : `ATAC-seq alignment failed: ${result.error}\n${result.stderr}`
        }]
      };
    }
  },

  {
    name: 'run_bedgraph',
    description: 'Convert BAM to bedGraph for genome browser visualization',
    inputSchema: {
      type: 'object',
      properties: {
        bam_file: {
          type: 'string',
          description: 'Input BAM file'
        },
        output_bedgraph: {
          type: 'string',
          description: 'Output bedGraph file'
        },
        genome: {
          type: 'string',
          description: 'Genome build',
          default: 'mm10'
        }
      },
      required: ['bam_file', 'output_bedgraph']
    },
    handler: async (args) => {
      const script = path.join(SCRIPTS_PATH, 'setgraph4ucsc');
      const scriptArgs = [
        args.bam_file,
        args.output_bedgraph,
        args.genome || 'mm10'
      ];

      const result = await executeScript(script, scriptArgs);

      return {
        content: [{
          type: 'text',
          text: result.success
            ? `bedGraph created: ${args.output_bedgraph}`
            : `bedGraph conversion failed: ${result.error}\n${result.stderr}`
        }]
      };
    }
  }
];

// Export helper function
export { executeScript };
