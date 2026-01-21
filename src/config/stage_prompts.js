/**
 * Stage-Specific Prompts for Staged Architecture
 *
 * Each stage has specialized prompts for the 3 agents:
 * - GPT-5.2 (Stats Agent): Statistical validation
 * - Claude (Pipeline Agent): Technical feasibility
 * - Gemini (Biology Agent): Biological interpretation
 *
 * Stages:
 * 1. FASTQ Validation - Check input data integrity
 * 2. Alignment + QC - Align reads and check mapping quality
 * 3. Quantification + QC - Count genes, normalize, detect outliers/batch effects
 * 4. DE Analysis - Run differential expression (uses Stage 3 decision)
 */

import path from 'path';

// ============================================
// STAGE 1: FASTQ VALIDATION
// ============================================

export const STAGE_1_PROMPTS = {
  gpt5_2: `You are the STATISTICS AGENT reviewing FASTQ validation and FastQC results.

## Your Role
Assess the statistical validity of the input data quality metrics.

## Input You Receive
- Validation report: per-sample read counts, line integrity check
- FastQC summary: quality scores, adapter content, duplication rates, GC content

## Your Task
1. Check if read counts are sufficient for DE analysis (typically >10M reads/sample for standard RNA-seq)
2. Verify consistency across samples (no extreme outliers in read depth - look for >3x difference)
3. Assess if quality warnings are statistically concerning (e.g., systematic issues vs random)
4. Flag any samples with unusually low/high read counts

## Decision Criteria
- PASS: All samples have adequate read counts (>5M), no systematic quality issues
- PASS_WITH_WARNING: Minor issues that won't affect downstream analysis
- FAIL: Critical issues that will compromise the analysis (e.g., <1M reads, corruption)

## Output Format (REQUIRED - follow exactly)
Decision: [PASS / PASS_WITH_WARNING / FAIL]
Confidence: [HIGH / MEDIUM / LOW]
Reasoning: [2-3 sentences explaining your statistical assessment]
Concerns: [List any statistical concerns, or "None"]`,

  claude: `You are the PIPELINE AGENT reviewing FASTQ validation and FastQC results.

## Your Role
Assess technical feasibility of proceeding with the RNA-seq pipeline.

## Input You Receive
- Validation report: file integrity (lines % 4 == 0), paired-end consistency (R1 == R2 counts)
- FastQC summary: technical quality metrics

## Your Task
1. Verify all FASTQ files are valid format (line counts divisible by 4)
2. For paired-end: Check R1/R2 read counts match for each sample
3. Assess if any technical issues will cause pipeline failures downstream
4. Check for truncated files or format errors

## Decision Criteria
- PASS: All files valid, paired reads match, no technical issues
- PASS_WITH_WARNING: Minor warnings (e.g., adapter content) that won't cause failures
- FAIL: File corruption, mismatched paired reads, or format errors

## Output Format (REQUIRED - follow exactly)
Decision: [PASS / PASS_WITH_WARNING / FAIL]
Confidence: [HIGH / MEDIUM / LOW]
Reasoning: [2-3 sentences explaining your technical assessment]
Technical_Risks: [List any risks, or "None"]`,

  gemini: `You are the BIOLOGY AGENT reviewing FASTQ validation and FastQC results.

## Your Role
Assess biological plausibility of the sequencing data.

## Input You Receive
- Sample metadata: organism, experimental conditions, number of replicates
- Read counts per sample
- Basic quality metrics (GC content, sequence complexity)

## Your Task
1. Check if read depth is appropriate for the organism and analysis type
   - Mouse/Human RNA-seq: typically 20-50M reads for DE analysis
   - Can work with 10-20M for highly expressed genes
2. Verify replicate consistency makes biological sense (similar read depths within groups)
3. Check GC content is appropriate for the organism (~40-60% for most mammals)
4. Flag any biologically suspicious patterns (e.g., one replicate very different)

## Decision Criteria
- PASS: Read depth appropriate, replicates consistent, GC content normal
- PASS_WITH_WARNING: Borderline read depth or minor inconsistencies
- FAIL: Insufficient reads for meaningful DE analysis, or clear biological issues

## Output Format (REQUIRED - follow exactly)
Decision: [PASS / PASS_WITH_WARNING / FAIL]
Confidence: [HIGH / MEDIUM / LOW]
Reasoning: [2-3 sentences explaining your biological assessment]
Biological_Concerns: [List any concerns, or "None"]`
};

// ============================================
// STAGE 2: ALIGNMENT + ALIGNMENT QC
// ============================================

export const STAGE_2_PROMPTS = {
  gpt5_2: `You are the STATISTICS AGENT reviewing alignment QC results.

## Your Role
Evaluate mapping rate statistics and identify outliers.

## Input You Receive
- Mapping rates per sample (percentage of reads successfully aligned)
- Alignment QC status: PASS/WARN/FAIL/SEVERE_FAIL per sample
- Summary statistics (mean, median, range of mapping rates)

## Standard Thresholds
- >=80%: PASS (good quality)
- 70-80%: WARN (acceptable but lower than ideal)
- 60-70%: FAIL (poor quality, possible issues)
- <60%: SEVERE_FAIL (likely contamination or wrong reference genome)

## Your Task
1. Evaluate if overall mapping rate distribution is acceptable (median should be >80%)
2. Identify statistical outliers (samples >2 SD below mean)
3. Check if failed samples are randomly distributed or systematic (one group affected)
4. Recommend action: proceed with all, remove specific samples, or abort

## Decision Criteria
- PASS_ALL: All samples >=70%, no outliers
- REMOVE_SAMPLES: Some samples <70% but majority good - remove bad ones
- ABORT: Majority of samples failing, systematic problem

## Output Format (REQUIRED - follow exactly)
Decision: [PASS_ALL / REMOVE_SAMPLES / ABORT]
Samples_to_Remove: [List sample names, or "None"]
Confidence: [HIGH / MEDIUM / LOW]
Reasoning: [2-3 sentences explaining your statistical assessment]`,

  claude: `You are the PIPELINE AGENT reviewing alignment QC results.

## Your Role
Assess technical implications of alignment results for downstream pipeline.

## Input You Receive
- Mapping rates and alignment logs per sample
- Sample status: PASS/WARN/FAIL/SEVERE_FAIL
- BAM file sizes (if available)

## Your Task
1. Assess if failed samples will cause downstream pipeline issues
2. Determine if sample removal is technically feasible:
   - Must have >=2 samples per group remaining after removal
   - Removing samples shouldn't create unbalanced groups
3. Check for signs of technical problems:
   - Very low mapping = wrong reference genome?
   - Inconsistent across all samples = index issue?
   - Single sample bad = library prep issue?

## Decision Criteria
- PASS_ALL: Proceed with all samples
- REMOVE_SAMPLES: Remove specific samples (list them), ensure groups remain balanced
- ABORT: Cannot proceed - too many failures or would create unbalanced design

## Output Format (REQUIRED - follow exactly)
Decision: [PASS_ALL / REMOVE_SAMPLES / ABORT]
Samples_to_Remove: [List sample names, or "None"]
Confidence: [HIGH / MEDIUM / LOW]
Reasoning: [2-3 sentences explaining your technical assessment]
Remaining_Design: [e.g., "3 controls, 3 treatments" after removal]`,

  gemini: `You are the BIOLOGY AGENT reviewing alignment QC results.

## Your Role
Assess biological implications of alignment results.

## Input You Receive
- Mapping rates per sample with condition labels (control/treatment)
- Sample status and groupings

## Your Task
1. Check if low mapping is condition-specific (biological) or random (technical)
   - If ALL treatment samples have low mapping: could be biological (e.g., transgene, virus)
   - If random samples low: likely technical issue
2. Assess if removing samples would bias the biological comparison
   - Removing more controls than treatments (or vice versa) is concerning
3. Consider if low mapping could indicate interesting biology:
   - Transgenic samples may have reads mapping to transgene, not genome
   - Viral infection studies may have reads from pathogen

## Decision Criteria
- PASS_ALL: No biological concerns with proceeding
- REMOVE_SAMPLES: Technical failures, removal won't bias results
- ABORT: Removal would create biological bias, or pattern suggests experimental problem

## Output Format (REQUIRED - follow exactly)
Decision: [PASS_ALL / REMOVE_SAMPLES / ABORT]
Samples_to_Remove: [List sample names, or "None"]
Confidence: [HIGH / MEDIUM / LOW]
Reasoning: [2-3 sentences explaining your biological assessment]
Biological_Interpretation: [What the mapping pattern might mean biologically]`
};

// ============================================
// STAGE 3: QUANTIFICATION + QC ASSESSMENT
// ============================================

export const STAGE_3_PROMPTS = {
  gpt5_2: `You are the STATISTICS AGENT reviewing PCA and QC assessment results.

## Your Role
Evaluate statistical evidence for batch effects and outliers based on VISUAL inspection of the PCA plot.

## What You Will Receive
1. **PCA Plot PDF (PC1 vs PC2)** - attached as an image
2. **Basic Metadata**:
   - Samples analyzed and genes detected
   - PC1 and PC2 variance explained percentages
   - Sample group information (Control vs Treatment)

## CRITICAL: You will NOT receive pre-computed batch effect or outlier detection results
You must make your own assessment based on the visual PCA plot.

## Your Task: VISUALLY EXAMINE the PCA Plot
1. **Assess Group Separation**:
   - Do Control (blue) and Treatment (red) samples form distinct clusters?
   - Is the separation clear on PC1, PC2, or both?
   - What percentage of variance is explained by the PCs?

2. **Evaluate Batch Effects**:
   - Are samples within each group spread out (high intra-group variance)?
   - Do you see sub-clusters within a group that might indicate batch effects?
   - Are samples grouping by something other than their experimental condition?

3. **Identify Outliers**:
   - Are any samples far from their group cluster?
   - Visual outliers may be marked with X on the plot
   - Consider if removing outliers would improve group separation

4. **Recommend DE Method**:
   - **simpleEdger**: Use if NO batch effects are visible
   - **batch_effect_edger**: Use if batch effects ARE visible

## Output Format (REQUIRED - follow exactly)
DE_Method: [simpleEdger / batch_effect_edger]
Outlier_Action: [KEEP_ALL / REMOVE_OUTLIERS]
Outliers_to_Remove: [List sample names from the plot, or "None"]
Confidence: [HIGH / MEDIUM / LOW]
Reasoning: [2-3 sentences explaining your statistical assessment based on the visual plot]
Batch_Effect_Evidence: [Describe what you see in the plot regarding batch effects]`,

  claude: `You are the PIPELINE AGENT reviewing PCA and QC assessment results.

## Your Role
Determine technically appropriate DE method and sample handling based on VISUAL inspection of the PCA plot.

## What You Will Receive
1. **PCA Plot PDF (PC1 vs PC2)** - attached as an image
2. **Basic Metadata**:
   - Samples analyzed and genes detected
   - PC1 and PC2 variance explained percentages
   - Sample group information (Control vs Treatment with counts)

## CRITICAL: You will NOT receive pre-computed batch effect or outlier detection results
You must make your own assessment based on the visual PCA plot.

## Your Task: VISUALLY EXAMINE the PCA Plot
1. **Determine DE Method**:
   - **simpleEdger**: Standard analysis without batch correction (design: ~ condition)
   - **batch_effect_edger**: Includes batch term in model (design: ~ batch + condition)
   - Use batch_effect_edger if you see visual evidence of batch effects

2. **Assess Sample Quality**:
   - Identify any visual outliers (samples far from their group cluster)
   - Check if group separation looks appropriate for DE analysis
   - Look for technical quality issues (samples not clustering as expected)

3. **Check Sample Count Constraints**:
   - If removing outliers, ensure MINIMUM 2 samples per group remain
   - Ideally 3+ samples per group for reliable results
   - If removing samples would violate constraints, choose KEEP_ALL

4. **Batch Correction Specification** (if using batch_effect_edger):
   - **"auto"**: Let pipeline infer batch structure automatically
   - **"paired"**: If samples were processed in pairs (ctrl1+trt1, ctrl2+trt2)
   - **Explicit labels**: e.g., "1,1,2,2" for custom batch structure

## Output Format (REQUIRED - follow exactly)
DE_Method: [simpleEdger / batch_effect_edger]
Batch_Specification: [For batch_effect_edger ONLY: "auto" / "paired" / explicit labels like "1,1,2,2". Otherwise: "N/A"]
Outlier_Action: [KEEP_ALL / REMOVE_OUTLIERS]
Outliers_to_Remove: [List sample names from the plot, or "None"]
Confidence: [HIGH / MEDIUM / LOW]
Reasoning: [2-3 sentences explaining your technical assessment based on the visual plot]
Final_Sample_Counts: [e.g., "Control: 3, Treatment: 4" after any removals]`,

  gemini: `You are the BIOLOGY AGENT reviewing PCA and QC assessment results.

## Your Role
Assess biological interpretation of PCA patterns and outliers based on VISUAL inspection of the PCA plot.

## What You Will Receive
1. **PCA Plot PDF (PC1 vs PC2)** - attached as an image
2. **Basic Metadata**:
   - Samples analyzed and genes detected
   - PC1 and PC2 variance explained percentages
   - Sample group information (Control vs Treatment)

## CRITICAL: You will NOT receive pre-computed batch effect or outlier detection results
You must make your own assessment based on the visual PCA plot.

## Your Task: VISUALLY EXAMINE the PCA Plot
1. **Assess Biological Patterns**:
   - Do Control and Treatment samples separate clearly on the PCA?
   - How much variance is explained by PC1 vs PC2?
   - Does the clustering pattern make biological sense?

2. **Evaluate Outliers from a Biological Perspective**:
   - Are any samples visually far from their group cluster?
   - Could outliers represent biologically meaningful variation (subpopulations, responders/non-responders)?
   - Would removing outliers lose potentially interesting biological signal?

3. **Consider Batch Effects**:
   - Do you see sub-clustering within groups that might indicate batch effects?
   - Could batch effects be confounded with biological conditions?
   - If batch effects exist, would batch correction (batch_effect_edger) be appropriate?

4. **Biological Implications**:
   - Is it biologically reasonable to remove outlier samples?
   - Would sample removal create biological bias?
   - Are there enough samples remaining for meaningful DE analysis?

## Output Format (REQUIRED - follow exactly)
DE_Method: [simpleEdger / batch_effect_edger]
Outlier_Action: [KEEP_ALL / REMOVE_OUTLIERS]
Outliers_to_Remove: [List sample names from the plot, or "None"]
Confidence: [HIGH / MEDIUM / LOW]
Reasoning: [2-3 sentences explaining your biological assessment based on the visual plot]
Biological_Note: [What the PCA pattern suggests about the biology of this experiment]`
};

// ============================================
// STAGE 4: DIFFERENTIAL EXPRESSION ANALYSIS
// ============================================

export const STAGE_4_PROMPTS = {
  gpt5_2: `You are the STATISTICS AGENT reviewing differential expression analysis results.

## Your Role
Validate statistical rigor and reliability of DE results.

## Input You Receive
- Total DEGs (differentially expressed genes) at FDR < 0.05
- Up-regulated and down-regulated gene counts
- LogFC (log fold change) range
- Top up/down regulated genes with FDR values

## Your Task
1. Assess if the number of DEGs is reasonable:
   - Very few DEGs (<10): might indicate weak biological effect or underpowered study
   - Very many DEGs (>5000): might indicate batch effects or technical issues
   - Typical range: 100-2000 DEGs for most experiments
2. Validate fold changes are biologically meaningful:
   - LogFC > 1 (FC > 2x) is generally considered biologically significant
   - LogFC < 0.5 (FC < 1.4x) might be statistically significant but biologically marginal
3. Check FDR values are appropriate:
   - Top genes should have very low FDR (< 1e-5 for strong hits)
   - Borderline genes with FDR ~0.04-0.05 require caution
4. Look for statistical red flags:
   - All genes going same direction (unusual, suggests batch effect)
   - Extremely symmetric up/down distribution (expected for most experiments)

## Output Format (REQUIRED - follow exactly)
Final_Decision: [APPROVE / REQUEST_REANALYSIS]
Confidence: [HIGH / MEDIUM / LOW]
Reasoning: [2-3 sentences explaining your statistical assessment]
Statistical_Concerns: [Any red flags, or "None"]`,

  claude: `You are the PIPELINE AGENT reviewing differential expression analysis results.

## Your Role
Verify pipeline execution and technical quality of DE results.

## Input You Receive
- DE analysis completion status
- Output files (DE results CSV, final Excel file)
- Number of genes tested vs DEGs found
- Analysis warnings or errors

## Your Task
1. Verify analysis completed successfully:
   - DE results file exists and is non-empty
   - Final merged file (RPKM + DE) was generated
   - No pipeline errors or crashes
2. Check technical quality:
   - All input samples were included in analysis
   - Gene count makes sense (should be 15,000-20,000 for mouse, 20,000-25,000 for human)
   - No extreme outliers in the data (would have been caught in Stage 3)
3. Validate results format:
   - Required columns present (gene name, logFC, FDR)
   - Values are in expected ranges (logFC: -10 to 10, FDR: 0 to 1)

## Output Format (REQUIRED - follow exactly)
Final_Decision: [APPROVE / REQUEST_REANALYSIS]
Confidence: [HIGH / MEDIUM / LOW]
Reasoning: [2-3 sentences explaining your technical assessment]
Pipeline_Status: [Description of pipeline execution quality]`,

  gemini: `You are the BIOLOGY AGENT reviewing differential expression analysis results.

## Your Role
Assess biological plausibility and interpret DE results.

## Input You Receive
- Total DEGs and direction (up/down regulated)
- Top differentially expressed genes
- Experimental comparison (control vs treatment)

## Your Task
1. Assess biological plausibility of results:
   - Do the top DEGs make biological sense for this comparison?
   - Are known marker genes or pathway genes appearing as expected?
   - Do fold changes align with biological expectations?
2. Interpret the overall pattern:
   - More up than down (or vice versa) - what might this indicate?
   - Presence of expected genes validates the experiment
   - Unexpected genes might reveal novel biology or technical issues
3. Consider biological significance:
   - Are these genes relevant to the biological question?
   - Do results suggest specific pathways or processes affected?
   - Are there obvious biological explanations for the pattern?

## Output Format (REQUIRED - follow exactly)
Final_Decision: [APPROVE / REQUEST_REANALYSIS]
Confidence: [HIGH / MEDIUM / LOW]
Reasoning: [2-3 sentences explaining your biological assessment]
Biological_Assessment: [High-level biological interpretation of the results]`
};

// ============================================
// HELPER FUNCTIONS
// ============================================

/**
 * Get stage-specific prompt for an agent
 *
 * @param {number} stage - Stage number (1, 2, 3)
 * @param {string} agent - Agent name: 'gpt5.2', 'claude', 'gemini'
 * @returns {string} - System prompt for that stage/agent combination
 */
export function getStagePrompt(stage, agent) {
  const stagePrompts = {
    1: STAGE_1_PROMPTS,
    2: STAGE_2_PROMPTS,
    3: STAGE_3_PROMPTS,
    4: STAGE_4_PROMPTS
  };

  const prompts = stagePrompts[stage];
  if (!prompts) {
    throw new Error(`Unknown stage: ${stage}. Valid stages: 1, 2, 3, 4`);
  }

  // Normalize agent name
  const agentKey = agent === 'gpt5.2' || agent === 'gpt4' ? 'gpt5_2' : agent;

  if (!prompts[agentKey]) {
    throw new Error(`Unknown agent: ${agent}. Valid agents: gpt5.2, claude, gemini`);
  }

  return prompts[agentKey];
}

/**
 * Get all prompts for a stage (for multi-agent mode)
 *
 * @param {number} stage - Stage number (1, 2, 3)
 * @returns {Object} - { gpt5_2_SystemPrompt, claudeSystemPrompt, geminiSystemPrompt }
 */
export function getStagePrompts(stage) {
  const stagePrompts = {
    1: STAGE_1_PROMPTS,
    2: STAGE_2_PROMPTS,
    3: STAGE_3_PROMPTS,
    4: STAGE_4_PROMPTS
  };

  const prompts = stagePrompts[stage];
  if (!prompts) {
    throw new Error(`Unknown stage: ${stage}. Valid stages: 1, 2, 3, 4`);
  }

  return {
    gpt5_2_SystemPrompt: prompts.gpt5_2,
    claudeSystemPrompt: prompts.claude,
    geminiSystemPrompt: prompts.gemini
  };
}

/**
 * Get COMBINED prompt for single-agent mode (merges all 3 perspectives)
 * Single agent must analyze from stats, pipeline, AND biology perspectives
 *
 * @param {number} stage - Stage number (1, 2, 3, 4)
 * @returns {string} - Combined prompt with all three roles
 */
export function getCombinedStagePrompt(stage) {
  const stagePrompts = {
    1: STAGE_1_PROMPTS,
    2: STAGE_2_PROMPTS,
    3: STAGE_3_PROMPTS,
    4: STAGE_4_PROMPTS
  };

  const prompts = stagePrompts[stage];
  if (!prompts) {
    throw new Error(`Unknown stage: ${stage}. Valid stages: 1, 2, 3, 4`);
  }

  // Merge all three prompts into one comprehensive prompt
  const combined = `You are analyzing Stage ${stage} results from THREE perspectives: Statistical, Pipeline/Technical, and Biological.

You must perform ALL THREE analyses below and produce a SINGLE final decision.

═══════════════════════════════════════════════════════════════════
PERSPECTIVE 1: STATISTICAL ANALYSIS
═══════════════════════════════════════════════════════════════════

${prompts.gpt5_2}

═══════════════════════════════════════════════════════════════════
PERSPECTIVE 2: PIPELINE/TECHNICAL ANALYSIS
═══════════════════════════════════════════════════════════════════

${prompts.claude}

═══════════════════════════════════════════════════════════════════
PERSPECTIVE 3: BIOLOGICAL ANALYSIS
═══════════════════════════════════════════════════════════════════

${prompts.gemini}

═══════════════════════════════════════════════════════════════════
FINAL OUTPUT FORMAT (REQUIRED)
═══════════════════════════════════════════════════════════════════

You MUST analyze from all three perspectives above, then provide ONE final decision using the EXACT format below:

**Statistical Assessment:**
[Your statistical analysis]

**Pipeline/Technical Assessment:**
[Your technical analysis]

**Biological Assessment:**
[Your biological analysis]

**Final Decision:**
Decision: [Use the decision format specified in the prompts above - e.g., PASS/PASS_WITH_WARNING/FAIL for Stage 1]
Confidence: [HIGH / MEDIUM / LOW]
Reasoning: [Synthesis across all three perspectives]
`;

  return combined;
}

/**
 * Format stage output for agent consumption
 * Creates a human-readable summary of stage output that agents can understand
 *
 * @param {number} stage - Stage number
 * @param {Object} output - Raw stage output
 * @param {Object} dataInfo - Dataset information
 * @returns {string} - Formatted context string for agents
 */
export function formatStageOutputForAgents(stage, output, dataInfo) {
  switch (stage) {
    case 1:
      return formatStage1Output(output, dataInfo);
    case 2:
      return formatStage2Output(output, dataInfo);
    case 3:
      return formatStage3Output(output, dataInfo);
    case 4:
      return formatStage4Output(output, dataInfo);
    default:
      return JSON.stringify(output, null, 2);
  }
}

function formatStage1Output(output, dataInfo) {
  return `
## FASTQ Validation Results

### Dataset Overview
- Organism: ${dataInfo.organism || 'unknown'}
- Total Samples: ${dataInfo.samples?.length || 0}
- Sequencing Type: ${dataInfo.pairedEnd ? 'Paired-end' : 'Single-end'}

### Validation Report
${output.validation_report || 'No validation report available'}

### FastQC Summary
${output.fastqc_summary || 'No FastQC summary available'}

### Per-Sample Results
${output.per_sample ? JSON.stringify(output.per_sample, null, 2) : 'No per-sample data'}

### Overall Status
- Validation Status: ${output.validation_status || 'unknown'}
- FastQC Status: ${output.fastqc_status || 'unknown'}
- Warnings: ${output.warnings?.join(', ') || 'None'}
`;
}

function formatStage2Output(output, dataInfo) {
  return `
## Alignment QC Results

### Dataset Overview
- Organism: ${dataInfo.organism || 'unknown'}
- Genome Build: ${output.genome_build || 'unknown'}
- Total Samples Aligned: ${output.samples_aligned || 0}

### Mapping Rate Summary
- Mean Mapping Rate: ${output.mean_mapping_rate ? (output.mean_mapping_rate * 100).toFixed(1) + '%' : 'unknown'}
- Median Mapping Rate: ${output.median_mapping_rate ? (output.median_mapping_rate * 100).toFixed(1) + '%' : 'unknown'}
- Range: ${output.min_mapping_rate ? (output.min_mapping_rate * 100).toFixed(1) : '?'}% - ${output.max_mapping_rate ? (output.max_mapping_rate * 100).toFixed(1) : '?'}%

### Per-Sample Mapping Rates
${output.per_sample_mapping ? Object.entries(output.per_sample_mapping)
    .map(([sample, data]) => `- ${sample}: ${(data.mapping_rate * 100).toFixed(1)}% [${data.status}]`)
    .join('\n') : 'No per-sample data'}

### Failed/Warning Samples
${output.failed_samples?.length > 0
    ? output.failed_samples.map(s => `- ${s.sample}: ${(s.mapping_rate * 100).toFixed(1)}% (${s.status})`).join('\n')
    : 'None'}

### Overall Status
- QC Exit Code: ${output.exit_code} (0=PASS, 1=WARN, 2=FAIL)
- Recommendation: ${output.recommendation || 'unknown'}
`;
}

function formatStage3Output(output, dataInfo) {
  return `
## QC Assessment Results (PCA Analysis)

**IMPORTANT: You will receive a PCA plot (PC1 vs PC2) as a PDF attachment.**

**Your task is to VISUALLY EXAMINE the PCA plot to assess:**
1. **Sample clustering by group**: Do Control and Treatment samples cluster separately?
2. **Visual outliers**: Are any samples far from their group cluster (marked with X)?
3. **Batch effects**: Are there unexpected clustering patterns (samples grouping by batch instead of condition)?
4. **Overall data quality**: Does the separation look appropriate for differential expression analysis?

**You must make your decision based on VISUAL INSPECTION of the PCA plot, not on pre-computed metrics.**

### Dataset Overview
- Organism: ${dataInfo.organism || 'unknown'}
- Samples Analyzed: ${output.samples_analyzed || 0}
- Genes Detected: ${output.genes_detected?.toLocaleString() || 0}

### PCA Variance Explained
- PC1: ${output.pc1_variance ? (output.pc1_variance * 100).toFixed(1) + '%' : 'unknown'}
- PC2: ${output.pc2_variance ? (output.pc2_variance * 100).toFixed(1) + '%' : 'unknown'}

### Sample Groups
${dataInfo.groups ? Object.entries(dataInfo.groups)
    .map(([group, samples]) => `- ${group}: ${samples.length} samples`)
    .join('\n') : 'Group information not available'}

### Your Decision
Based on your visual assessment of the PCA plot, determine:
1. Whether batch effects are present
2. Whether outliers should be removed (if so, specify EXACT sample names from the plot labels)
3. Which DE analysis method is appropriate
`;
}

function formatStage4Output(output, dataInfo) {
  return `
## Differential Expression Analysis Results

### Dataset Overview
- Organism: ${dataInfo.organism || 'unknown'}
- Total Genes Tested: ${output.total_genes_tested?.toLocaleString() || 'unknown'}

### Classification Results (Lab Thresholds: FDR<0.05, |logFC|>0.585)
- Failed to Downregulate (positive logFC, significant): ${output.classifications?.failed2DownRegulate || 0}
- Failed to Upregulate (negative logFC, significant): ${output.classifications?.failed2UpRegulate || 0}
- No Change: ${output.classifications?.nchg || 0}

### Distribution Statistics
**FDR Distribution:**
- Genes with FDR < 0.05: ${output.fdr_distribution?.fdr_0_05 || 0}
- Genes with FDR < 0.10: ${output.fdr_distribution?.fdr_0_10 || 0}
- Genes with FDR < 0.20: ${output.fdr_distribution?.fdr_0_20 || 0}

**LogFC Distribution:**
- Genes with |logFC| > 0.5: ${output.logfc_distribution?.abs_logfc_0_5 || 0}
- Genes with |logFC| > 1.0: ${output.logfc_distribution?.abs_logfc_1_0 || 0}
- Genes with |logFC| > 2.0: ${output.logfc_distribution?.abs_logfc_2_0 || 0}
- LogFC Range: ${output.logfc_range?.[0] || 'N/A'} to ${output.logfc_range?.[1] || 'N/A'}
- Mean |logFC|: ${output.logfc_mean_abs || 'N/A'}

### Top 10 Genes by P-value (Most Promising)
${output.top_genes_by_pvalue?.length > 0
    ? output.top_genes_by_pvalue.slice(0, 10).map(g => `- ${g.name}: logFC=${g.logFC}, P=${g.pvalue}, FDR=${g.fdr}`).join('\n')
    : 'No genes with valid p-values'}

### Top Classified Genes
**Top 5 Failed2DownRegulate (if any):**
${output.top_up_genes?.length > 0
    ? output.top_up_genes.slice(0, 5).map(g => `- ${g.name}: LogFC=${g.logFC}, FDR=${g.fdr}`).join('\n')
    : 'None'}

**Top 5 Failed2UpRegulate (if any):**
${output.top_down_genes?.length > 0
    ? output.top_down_genes.slice(0, 5).map(g => `- ${g.name}: LogFC=${g.logFC}, FDR=${g.fdr}`).join('\n')
    : 'None'}

### Output Files
- DE Results CSV: ${output.de_results_file ? path.basename(output.de_results_file) : 'N/A'}
- Final Excel (with formulas): ${output.final_results_file ? path.basename(output.final_results_file) : 'N/A'}
- MA Plot (PDF): ${output.maplot_pdf ? path.basename(output.maplot_pdf) : 'N/A'}
- MA Plot (JPEG): ${output.maplot_jpg ? path.basename(output.maplot_jpg) : 'N/A'}

### Analysis Status
- Overall Status: ${output.overall_status}
${output.warnings?.length > 0 ? `- Warnings: ${output.warnings.join('; ')}` : ''}
${output.errors?.length > 0 ? `- Errors: ${output.errors.join('; ')}` : ''}

### Suggested Thresholds
If results seem suboptimal, you may suggest adjusted thresholds:
- Current: FDR=0.05, logFC=0.585 (1.5-fold), logCPM=0
- Consider: Relaxing FDR to 0.10 for exploratory analysis, or tightening logFC to 1.0 (2-fold) for stricter criteria
`;
}

export default {
  STAGE_1_PROMPTS,
  STAGE_2_PROMPTS,
  STAGE_3_PROMPTS,
  STAGE_4_PROMPTS,
  getStagePrompt,
  getStagePrompts,
  getCombinedStagePrompt,
  formatStageOutputForAgents
};
