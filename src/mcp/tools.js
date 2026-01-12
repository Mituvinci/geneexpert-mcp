/**
 * GeneExpert MCP Tools
 * Wraps bioinformatics scripts (configured via SCRIPTS_PATH environment variable)
 */

import { exec } from 'child_process';
import { promisify } from 'util';
import path from 'path';
import fs from 'fs';

const execAsync = promisify(exec);

// Script paths from environment or defaults
const SCRIPTS_PATH = process.env.SCRIPTS_PATH || '/users/ha00014/Halimas_projects/multi_llm_mcp/bio_informatics/scripts';
const CUSTOM_SCRIPTS_PATH = process.env.CUSTOM_SCRIPTS_PATH || '/users/ha00014/Halimas_projects/multi_llm_mcp/bio_informatics/scripts';

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
  },

  // ============================================
  // ANALYSIS TOOLS - For agents to READ outputs
  // ============================================

  {
    name: 'read_bam_summary',
    description: 'Read BAM alignment statistics from log files. Returns mapping rates, aligned reads, etc.',
    inputSchema: {
      type: 'object',
      properties: {
        bam_directory: {
          type: 'string',
          description: 'Directory containing BAM files and .log files'
        }
      },
      required: ['bam_directory']
    },
    handler: async (args) => {
      try {
        const fs = await import('fs');
        const logFiles = fs.readdirSync(args.bam_directory).filter(f => f.endsWith('.log'));

        let summary = '=== BAM Alignment Summary ===\n\n';

        for (const logFile of logFiles) {
          const logPath = path.join(args.bam_directory, logFile);
          const logContent = fs.readFileSync(logPath, 'utf-8');

          // Extract key statistics from Subread log
          const totalReads = logContent.match(/Total reads\s*:\s*(\d+)/)?.[1];
          const mappedReads = logContent.match(/Mapped reads\s*:\s*(\d+)/)?.[1];
          const mappingRate = logContent.match(/Successfully assigned\s*:\s*([\d.]+)%/)?.[1];

          summary += `Sample: ${logFile.replace('.log', '')}\n`;
          summary += `  Total reads: ${totalReads || 'N/A'}\n`;
          summary += `  Mapped reads: ${mappedReads || 'N/A'}\n`;
          summary += `  Mapping rate: ${mappingRate || 'N/A'}%\n`;
          summary += `  Status: ${mappingRate && parseFloat(mappingRate) > 70 ? '✓ Good' : '⚠️ Low mapping rate!'}\n\n`;
        }

        return {
          content: [{
            type: 'text',
            text: summary
          }]
        };
      } catch (error) {
        return {
          content: [{
            type: 'text',
            text: `Failed to read BAM summary: ${error.message}`
          }]
        };
      }
    }
  },

  {
    name: 'read_count_summary',
    description: 'Read count matrix and provide summary statistics (total genes, counts distribution, etc.)',
    inputSchema: {
      type: 'object',
      properties: {
        count_file: {
          type: 'string',
          description: 'Path to count matrix file'
        }
      },
      required: ['count_file']
    },
    handler: async (args) => {
      try {
        const fs = await import('fs');
        const content = fs.readFileSync(args.count_file, 'utf-8');
        const lines = content.trim().split('\n');

        const header = lines[0].split(/\s+/);
        const samples = header.slice(2); // Skip ID and Length columns
        const nGenes = lines.length - 1;

        let summary = '=== Count Matrix Summary ===\n\n';
        summary += `File: ${path.basename(args.count_file)}\n`;
        summary += `Genes: ${nGenes}\n`;
        summary += `Samples: ${samples.length}\n`;
        summary += `Sample names: ${samples.join(', ')}\n\n`;

        // Calculate total counts per sample
        summary += 'Total counts per sample:\n';
        for (let i = 0; i < samples.length; i++) {
          let totalCounts = 0;
          for (let j = 1; j < lines.length; j++) {
            const cols = lines[j].split(/\s+/);
            totalCounts += parseInt(cols[i + 2]) || 0;
          }
          summary += `  ${samples[i]}: ${totalCounts.toLocaleString()}\n`;
        }

        return {
          content: [{
            type: 'text',
            text: summary
          }]
        };
      } catch (error) {
        return {
          content: [{
            type: 'text',
            text: `Failed to read count matrix: ${error.message}`
          }]
        };
      }
    }
  },

  {
    name: 'read_file',
    description: 'Read any text file content. Use this to inspect R script outputs, log files, or results.',
    inputSchema: {
      type: 'object',
      properties: {
        file_path: {
          type: 'string',
          description: 'Path to file to read'
        },
        max_lines: {
          type: 'number',
          description: 'Maximum number of lines to read (default: 100)',
          default: 100
        }
      },
      required: ['file_path']
    },
    handler: async (args) => {
      try {
        const fs = await import('fs');

        // Try multiple file paths (many scripts don't have .sh extension!)
        let filePath = args.file_path;
        let content;
        let actualPath;

        // Try original path first
        if (fs.existsSync(filePath)) {
          content = fs.readFileSync(filePath, 'utf-8');
          actualPath = filePath;
        }
        // If ends with .sh, try without it
        else if (filePath.endsWith('.sh')) {
          const pathWithoutSh = filePath.replace(/\.sh$/, '');
          if (fs.existsSync(pathWithoutSh)) {
            content = fs.readFileSync(pathWithoutSh, 'utf-8');
            actualPath = pathWithoutSh;
          }
        }
        // If doesn't end with .sh, try adding it
        else {
          const pathWithSh = filePath + '.sh';
          if (fs.existsSync(pathWithSh)) {
            content = fs.readFileSync(pathWithSh, 'utf-8');
            actualPath = pathWithSh;
          }
        }

        // If still not found, throw error
        if (!content) {
          throw new Error(`File not found: ${filePath} (also tried with/without .sh extension)`);
        }

        const lines = content.split('\n');
        const maxLines = args.max_lines || 100;

        const output = lines.slice(0, maxLines).join('\n');
        const truncated = lines.length > maxLines;

        return {
          content: [{
            type: 'text',
            text: `=== File: ${path.basename(actualPath)} ===\n${actualPath !== filePath ? `(Found as: ${path.basename(actualPath)})\n` : ''}\n${output}${truncated ? `\n\n... (${lines.length - maxLines} more lines truncated)` : ''}`
          }]
        };
      } catch (error) {
        return {
          content: [{
            type: 'text',
            text: `Failed to read file: ${error.message}`
          }]
        };
      }
    }
  },

  {
    name: 'write_custom_script',
    description: 'Write a custom R or bash script. Use this when standard pipeline needs adaptation (batch correction, special QC, etc.)',
    inputSchema: {
      type: 'object',
      properties: {
        script_path: {
          type: 'string',
          description: 'Where to save the script (e.g., /tmp/custom_batch_correction.R)'
        },
        script_content: {
          type: 'string',
          description: 'The complete script content'
        },
        description: {
          type: 'string',
          description: 'What this custom script does'
        }
      },
      required: ['script_path', 'script_content', 'description']
    },
    handler: async (args) => {
      try {
        const fs = await import('fs');
        fs.writeFileSync(args.script_path, args.script_content, 'utf-8');
        fs.chmodSync(args.script_path, '755'); // Make executable

        return {
          content: [{
            type: 'text',
            text: `✓ Custom script created: ${args.script_path}\n\nDescription: ${args.description}\n\nScript is ready to execute.`
          }]
        };
      } catch (error) {
        return {
          content: [{
            type: 'text',
            text: `Failed to write custom script: ${error.message}`
          }]
        };
      }
    }
  },

  {
    name: 'execute_custom_script',
    description: 'Execute a custom script that was written to handle edge cases',
    inputSchema: {
      type: 'object',
      properties: {
        script_path: {
          type: 'string',
          description: 'Path to custom script'
        },
        args: {
          type: 'array',
          items: { type: 'string' },
          description: 'Script arguments',
          default: []
        },
        working_dir: {
          type: 'string',
          description: 'Working directory',
          default: '/tmp'
        }
      },
      required: ['script_path']
    },
    handler: async (args) => {
      const scriptArgs = args.args || [];
      const result = await executeScript(args.script_path, scriptArgs, {
        workingDir: args.working_dir || '/tmp'
      });

      return {
        content: [{
          type: 'text',
          text: result.success
            ? `✓ Custom script executed successfully!\n\n${result.stdout}`
            : `✗ Custom script failed: ${result.error}\n${result.stderr}`
        }]
      };
    }
  },

  {
    name: 'list_available_scripts',
    description: 'List all available R and bash scripts in the bioinformatics scripts directory. Use this BEFORE generating pipeline scripts to validate which scripts exist.',
    inputSchema: {
      type: 'object',
      properties: {
        pattern: {
          type: 'string',
          description: 'Optional glob pattern to filter scripts (e.g., "*.R", "edger*", "qc*")',
          default: '*'
        }
      }
    },
    handler: async (args) => {
      try {
        const pattern = args.pattern || '*';
        const files = fs.readdirSync(SCRIPTS_PATH);

        // Filter by pattern
        const matchingFiles = files.filter(file => {
          if (pattern === '*') return true;

          // Simple glob matching
          const regex = new RegExp('^' + pattern.replace(/\*/g, '.*').replace(/\?/g, '.') + '$');
          return regex.test(file);
        });

        // Categorize scripts
        const rScripts = matchingFiles.filter(f => f.endsWith('.R'));
        const bashScripts = matchingFiles.filter(f => !f.endsWith('.R') && !f.endsWith('.sh'));
        const shellScripts = matchingFiles.filter(f => f.endsWith('.sh'));

        let output = `📂 Available scripts in ${SCRIPTS_PATH}\n\n`;

        output += `R Scripts (${rScripts.length}):\n`;
        rScripts.sort().forEach(script => {
          output += `  - ${script}\n`;
        });

        if (bashScripts.length > 0) {
          output += `\nBash Scripts (${bashScripts.length}):\n`;
          bashScripts.sort().forEach(script => {
            output += `  - ${script}\n`;
          });
        }

        if (shellScripts.length > 0) {
          output += `\nShell Scripts (${shellScripts.length}):\n`;
          shellScripts.sort().forEach(script => {
            output += `  - ${script}\n`;
          });
        }

        output += `\n✓ Total: ${matchingFiles.length} scripts found`;

        if (pattern !== '*') {
          output += ` (matching pattern: "${pattern}")`;
        }

        return {
          content: [{
            type: 'text',
            text: output
          }]
        };
      } catch (error) {
        return {
          content: [{
            type: 'text',
            text: `Failed to list scripts: ${error.message}`
          }]
        };
      }
    }
  },

  {
    name: 'validate_fastq',
    description: 'Validate FASTQ input files for quality, completeness, and paired-end consistency. Call this FIRST before any analysis to ensure data integrity.',
    inputSchema: {
      type: 'object',
      properties: {
        input_dir: {
          type: 'string',
          description: 'Directory containing FASTQ files'
        },
        output_report: {
          type: 'string',
          description: 'Path to save validation report TSV',
          default: '/tmp/fastq_validation.tsv'
        }
      },
      required: ['input_dir']
    },
    handler: async (args) => {
      const inputDir = args.input_dir;
      const reportPath = args.output_report || '/tmp/fastq_validation.tsv';

      try {
        // Check if validate_fastq.sh exists
        const validatorPath = `${SCRIPTS_PATH}/validate_fastq.sh`;
        if (!fs.existsSync(validatorPath)) {
          return {
            content: [{
              type: 'text',
              text: `❌ Validator script not found: ${validatorPath}\n\nPlease create this script or use manual validation.`
            }]
          };
        }

        // Run validation script
        const result = await executeScript(validatorPath, [inputDir, reportPath]);

        if (!result.success) {
          return {
            content: [{
              type: 'text',
              text: `❌ Validation script failed:\n${result.stderr}`
            }]
          };
        }

        // Read validation report
        const reportContent = fs.readFileSync(reportPath, 'utf-8');
        const lines = reportContent.trim().split('\n');

        // Parse TSV
        const header = lines[0].split('\t');
        const rows = lines.slice(1).map(line => {
          const values = line.split('\t');
          return Object.fromEntries(header.map((h, i) => [h, values[i]]));
        });

        // Analyze results
        const totalFiles = rows.length;
        const validFiles = rows.filter(r => r.status === 'VALID').length;
        const invalidFiles = rows.filter(r => r.status !== 'VALID');
        const mismatchFiles = rows.filter(r => r.status === 'MISMATCH');

        let summary = `📊 FASTQ Validation Report\n\n`;
        summary += `Input Directory: ${inputDir}\n`;
        summary += `Total Files: ${totalFiles}\n`;
        summary += `Valid Files: ${validFiles}\n`;
        summary += `Invalid Files: ${invalidFiles.length}\n`;
        summary += `Mismatched Pairs: ${mismatchFiles.length}\n\n`;

        if (invalidFiles.length > 0) {
          summary += `⚠️  ISSUES FOUND:\n`;
          invalidFiles.forEach(f => {
            summary += `  - ${f.file}: ${f.status} (reads: ${f.reads})\n`;
          });
          summary += `\n`;
        }

        summary += `Full Report:\n`;
        summary += reportContent;

        summary += `\n\n📋 DECISION:\n`;
        if (invalidFiles.length === 0) {
          summary += `✅ All files valid - PROCEED with analysis`;
        } else if (invalidFiles.length < totalFiles / 2) {
          summary += `⚠️  Some files invalid - CONSIDER excluding bad files and proceeding`;
        } else {
          summary += `❌ Too many invalid files - STOP and fix data`;
        }

        return {
          content: [{
            type: 'text',
            text: summary
          }]
        };

      } catch (error) {
        return {
          content: [{
            type: 'text',
            text: `Failed to validate FASTQ files: ${error.message}`
          }]
        };
      }
    }
  }
];

// Export helper function
export { executeScript };
