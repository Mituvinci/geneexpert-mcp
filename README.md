# GeneExpert: Multi-Agent RNA-seq Analysis System

**Staged Multi-Agent Pipeline** where GPT-5.2, Claude Sonnet 4.5, and Gemini Pro collaborate at decision checkpoints throughout RNA-seq analysis.

**Key Innovation:** 4-stage architecture with agent review checkpoints - agents validate each stage, detect issues early, and decide whether to proceed or adjust the approach.

**Research Goal:** Multi-agent collaboration with user-in-loop tracking reduces errors 40%+ through staged validation and consensus-based decision making.

---

## Quick Start

### 0. npm installation

```bash
Node: v18.20.8
npm: 10.8.2
OpenAI SDK: v4.104.0
Anthropic (Claude) SDK: v0.32.1
Google Gemini SDK: v0.21.0
Model Context Protocol (MCP) SDK: v1.25.1
```

### 1. Requirements

```bash
Python: 3.10.19
R: R version 4.3.3 (2024-02-29)
Subread-align: v2.1.1
FastQC: v0.12.1
SAMtools: samtools 1.22.1
```

### 2. Configure API Keys

```bash
cp .env.example .env
nano .env  # Add your OpenAI, Anthropic, Google API keys
```

### 3. Run Analysis (Staged Architecture)

```bash
node bin/geneexpert.js analyze data/your_input_folder \
  --staged \
  --organism mouse \
  --comparison "stroke_vs_control" \
  --control-keyword "cont" \
  --treatment-keyword "ips" \
  --output results/your_output_folder
```

**What happens:**
1. Stage 1: FASTQ Validation - Agents review file integrity and quality
2. Stage 2: Alignment + QC - Agents review mapping rates, detect contamination
3. Stage 3: Quantification + PCA QC - Agents review outliers, batch effects, choose DE method
4. Stage 4: DE Analysis - Agents review differential expression results
5. User input tracked when agents can't reach consensus (critical for research evaluation)

### Two Multi-Agent Modes

**Parallel Independent Voting (Default):**
```bash
# All agents respond simultaneously with independent perspectives
node bin/geneexpert.js analyze data/your_input \
  --staged --organism mouse --output results/parallel
```

**Sequential Chain:**
```bash
# Agents respond sequentially: GPT-5.2 → Gemini → Claude
# Each agent sees previous responses for informed synthesis
node bin/geneexpert.js analyze data/your_input \
  --staged --sequential-chain --organism mouse --output results/sequential
```

**Research Question:** Does sequential information passing improve decisions, or does parallel independent voting reduce errors better?


---

## What Makes This Novel

### Staged Architecture with Agent Checkpoints

**Traditional pipelines:**
```
Run all steps → Hope everything works → Debug if fails at the end
```

**GeneExpert Staged Multi-Agent System:**
```
Stage 1: FASTQ Validation
  Generate script → Execute → Parse results → Format for agents
  ↓
  3 Agents Review:
  [GPT-5.2 Stats]:    "Read depth sufficient, 15M reads/sample"
  [Claude Pipeline]:  "File integrity good, paired-end consistent"
  [Gemini Biology]:   "Quality scores acceptable, proceed to alignment"
  ↓
  Consensus: PASS (3/3 agents agree) → Proceed to Stage 2
  ↓
Stage 2: Alignment + Alignment QC
  Generate script → Execute → Parse alignment QC → Format for agents
  ↓
  3 Agents Review:
  [GPT-5.2 Stats]:    "Mapping rate 85% - excellent"
  [Claude Pipeline]:  "All samples >80%, no contamination"
  [Gemini Biology]:   "Genome match confirmed"
  ↓
  Consensus: PASS_ALL (keep all samples) → Proceed to Stage 3
  ↓
Stage 3: Quantification + PCA/QC Assessment
  Generate script → Execute → Parse PCA/outliers/batch → Format for agents
  ↓
  3 Agents Review:
  [GPT-5.2 Stats]:    "No statistical outliers, no batch effects"
  [Claude Pipeline]:  "PCA shows good separation by condition"
  [Gemini Biology]:   "Biological groups cluster appropriately"
  ↓
  Consensus: DE_Method=simpleEdger, KEEP_ALL → Proceed to Stage 4
  ↓
Stage 4: Differential Expression Analysis
  Generate script → Execute → Parse DE results → Format for agents
  ↓
  3 Agents Review:
  [GPT-5.2 Stats]:    "1,234 DEGs at FDR<0.05 - reasonable"
  [Claude Pipeline]:  "LogFC range appropriate, no red flags"
  [Gemini Biology]:   "Top genes biologically meaningful"
  ↓
  Consensus: APPROVE → Analysis Complete!
```

**This is TRUE staged validation:** Issues caught early, decisions made with full context, user input tracked when needed!

---

## System Architecture

### 4-Stage Pipeline

**Stage 1: FASTQ Validation**
```
Script Steps:
  1. validate_fastq.sh - File integrity check (lines % 4 == 0)
  2. fastqc - Quality control metrics

Agent Decision:
  - PASS: Sufficient read depth, no corruption
  - PASS_WITH_WARNING: Minor issues, user decides
  - FAIL: Critical issues, abort analysis

Outputs: stage1_validation/ with validation reports, FastQC summaries
```

**Stage 2: Alignment + Alignment QC**
```
Script Steps:
  1. subread-align - Align FASTQ to BAM (genome: mm10/hg38/rn6)
  2. alignment_qc_screening.R - Parse alignment logs for mapping rates

QC Thresholds:
  >=80%: PASS
  70-80%: WARN
  60-70%: FAIL
  <60%: SEVERE_FAIL (likely contamination or wrong genome)

Agent Decision:
  - PASS_ALL: All samples >=70%, proceed with all
  - REMOVE_SAMPLES: Some <70%, remove failed samples
  - ABORT: Too many failures or unbalanced design

Outputs: stage2_alignment/bam_files/, alignment QC summary
```

**Stage 3: Quantification + PCA/QC Assessment**
```
Script Steps:
  1. featureCounts - Count reads per gene
  2. filterIDS.R - Remove problematic gene IDs
  3. RPKM.R - Normalize for visualization
  4. entrz.R - Add gene symbols
  5. qc_assessment_pca.R - PCA analysis + outlier/batch detection

QC Detection:
  - Outliers: Robust distance-based (>3 MAD from group centroid)
  - Batch Effects: Median distance to centroid + MAD threshold

Agent Decision:
  - DE_Method: simpleEdger (no batch) vs batch_effect_edger (with batch correction)
  - Batch_Specification: auto / paired / explicit labels
  - Outlier_Action: KEEP_ALL vs REMOVE_OUTLIERS
  - Outliers_to_Remove: List of specific samples (unanimous vote required)

Outputs: stage3_quantification/ with count matrices, RPKM, PCA plots, QC metrics
```

**Stage 4: Differential Expression Analysis**
```
Script Steps:
  1. DE Analysis:
     - simpleEdger3.R (design: ~ condition) OR
     - batch_effect_edgeR_v3.R (design: ~ batch + condition)
  2. merge_results.R - Combine RPKM + DE results

Agent Decision:
  - APPROVE: Results validated, analysis complete
  - REQUEST_REANALYSIS: Issues detected, re-run with different parameters

Outputs: stage4_de_analysis/ with DE results CSV, final Excel file
```

### Multi-Agent Decision System

**3 Specialized Agents:**

- **GPT-5.2 (Stats Agent)**: Statistical validation, threshold selection
  - Validates mapping rates, DEG counts, statistical thresholds
  - Pure reasoning, no tool access

- **Claude Sonnet 4.5 (Pipeline Agent)**: Technical feasibility assessment
  - Validates pipeline execution, file formats, technical correctness
  - Pure reasoning, no tool access

- **Gemini Pro (Biology Agent)**: Biological interpretation
  - Validates biological assumptions, gene patterns, QC criteria
  - Pure reasoning, no tool access

**Two Multi-Agent Modes:**

**Mode 1: Parallel Independent Voting (Default)**
```
All 3 agents receive SAME input → Respond simultaneously → Majority vote
- Fast execution (1 round of API calls)
- True independent perspectives
- Reduces anchoring bias
- Default behavior (no flag needed)
```

**Mode 2: Sequential Chain**
```
GPT-5.2 (Stats) → Gemini (Biology) → Claude (Pipeline)
- Each agent sees previous responses
- Claude synthesizes with full context
- Slower (3x sequential calls)
- More informed decisions
- Use flag: --sequential-chain
```

**Consensus Voting Rules:**
- **Parallel mode**: Majority vote (2/3 required)
- **Sequential mode**: Final agent (Claude) makes synthesis decision
- **Sample removal**: Unanimous (3/3 required) in parallel mode
- **No consensus**: Escalate to user, track user input for evaluation

**User Input Tracking (Critical for Research):**
- Every decision tracked with `user_input_required` flag
- When agents disagree, user decides
- All user decisions logged in JSON with reasoning
- Metrics: % autonomy, user agreement rate, human-in-loop impact

---

## Key Features

### Implemented & Tested:

1. **4-Stage Architecture with Agent Checkpoints**
   - Each stage: Generate → Execute → Parse → Format → Agent Review → Decision
   - Early issue detection (catch bad samples in Stage 2, not Stage 4)
   - Progressive refinement (remove outliers before DE analysis)

2. **User Input Tracking (Research Critical)**
   - JSON logs track when user input required
   - User decisions recorded with reasoning
   - Enables research metrics: autonomy rate, agreement rate, impact analysis

3. **Agent Consensus at Each Stage**
   - Stage 1: PASS / PASS_WITH_WARNING / FAIL
   - Stage 2: PASS_ALL / REMOVE_SAMPLES / ABORT
   - Stage 3: DE method + outlier removal decisions
   - Stage 4: APPROVE / REQUEST_REANALYSIS

4. **Sample Propagation**
   - Stage 2 can remove failed samples (low mapping rate)
   - Stage 3 can remove outliers (PCA-based detection)
   - Stage 4 uses final refined sample set

5. **Batch Effect Handling**
   - Stage 3 agents detect batch effects from PCA
   - Choose simpleEdger (no batch) vs batch_effect_edger (with correction)
   - Automatic design matrix adjustment

6. **Complete Logging**
   - Stage-by-stage JSON logs
   - Agent conversations tracked
   - User decisions recorded
   - All metrics for research evaluation

7. **Flexible Execution**
   - Run full 4-stage pipeline: `--staged` flag
   - Test individual stages: `test_stage1.js`, `test_stage2.js`, etc.
   - Skip agent review: `--force-automation` flag
   - Single-agent baseline: `--single-agent gpt5.2|claude|gemini`
   - Sequential chain mode: `--sequential-chain` flag

8. **Hybrid Approach with Lab-Standard Formulas**
   - Excel file generation with dynamic formulas and thresholds
   - Expr formula checks if gene is expressed (max group average logRPKM > 2 AND logCPM > 0)
   - CG classification categorizes genes as failed2DownRegulate, failed2UpRegulate, or nchg
   - Formula-based auto-counting for each classification category
   - Lab thresholds: FDR=0.05, logCPM=0, logFC=0.585 (1.5-fold), logRPKM=2
   - MA plot visualization with customizable thresholds, PDF + JPEG formats
   - Enhanced agent summaries with distribution statistics, top genes by p-value, classification counts
   - Agent-driven threshold adjustments based on data characteristics

---

## Experimental Modes (Research Evaluation)

GeneExpert supports 5 experimental systems for ICML 2026 evaluation:

### System 1: No-Agent (Template Only)
```bash
node bin/geneexpert.js analyze <dataset> \
  --staged \
  --force-automation \
  --output results/no_agent
```
Baseline with no agent review, template-based decisions only.

### System 2: Single-Agent (GPT-5.2)
```bash
node bin/geneexpert.js analyze <dataset> \
  --staged \
  --single-agent gpt5.2 \
  --output results/single_gpt
```
GPT-5.2 performs ALL roles: statistics + pipeline + biology.

### System 3: Single-Agent (Claude)
```bash
node bin/geneexpert.js analyze <dataset> \
  --staged \
  --single-agent claude \
  --output results/single_claude
```
Claude performs ALL roles: statistics + pipeline + biology.

### System 4: Multi-Agent Parallel (Default)
```bash
node bin/geneexpert.js analyze <dataset> \
  --staged \
  --output results/multi_agent_parallel
```
3 agents review simultaneously with independent voting.

### System 5: Multi-Agent Sequential Chain
```bash
node bin/geneexpert.js analyze <dataset> \
  --staged \
  --sequential-chain \
  --output results/multi_agent_sequential
```
3 agents review sequentially: GPT-5.2 → Gemini → Claude (informed synthesis).

**Evaluation Metrics:**
- Decision accuracy (correct decisions / total)
- Error reduction (baseline errors - multi-agent errors) / baseline errors
- Success rate (successful analyses / total)
- Cost efficiency (total cost / successes)
- Inter-agent agreement (Cohen's kappa for multi-agent modes)
- User input frequency (% decisions requiring user input)
- Error propagation rate (sequential mode - measures if later agents correct or amplify first agent's errors)
- Reasoning quality (qualitative assessment)

**Evaluation Script:**
```bash
# Calculate all metrics
node bin/evaluate.js metrics

# Compare parallel vs sequential
node bin/evaluate.js compare

# Export results to CSV
node bin/evaluate.js metrics --output results.csv
```

---

## Project Structure

```
├── bin/geneexpert.js                    # CLI entry point (--staged flag)
├── src/
│   ├── executor/
│   │   └── staged_executor.js           # 4-stage orchestration
│   ├── stages/
│   │   ├── stage1_validation.js         # FASTQ validation
│   │   ├── stage2_alignment.js          # Alignment + QC
│   │   ├── stage3_quantification_qc.js  # Quantification + PCA QC
│   │   └── stage4_de_analysis.js        # Differential expression
│   ├── config/
│   │   └── stage_prompts.js             # Agent prompts for each stage
│   ├── coordinator/
│   │   └── orchestrator.js              # Multi-agent coordination
│   ├── utils/
│   │   ├── logger.js                    # JSON logging with user input tracking
│   │   ├── user_input.js                # User decision handling
│   │   └── llm_clients.js               # OpenAI, Anthropic, Google APIs
│   └── pipeline/                        # Old monolithic architecture (legacy)
├── test_stage1.js                       # Test Stage 1 independently
├── test_stage2.js                       # Test Stage 2 independently
├── test_stage3.js                       # Test Stage 3 independently
├── test_stage4.js                       # Test Stage 4 independently
├── experiments/
│   └── ground_truth.json                # Evaluation dataset labels
└── README.md                            # This file
```

---

## Cost Estimate

### Per-Stage Agent Review:
- 3 agents x 4 stages = 12 agent API calls per analysis
- ~$0.01 per agent call
- **Total agent cost: ~$0.12 per complete analysis**

### Script Execution:
- All bioinformatics tools run locally (subread-align, featureCounts, edgeR)
- **Execution cost: $0**

**Total: ~$0.12 per RNA-seq analysis with full multi-agent validation**



## Development Commands

```bash
# Run full staged pipeline
node bin/geneexpert.js analyze data/DA0036 \
  --staged \
  --organism mouse \
  --comparison "experiment_name" \
  --control-keyword "cont" \
  --treatment-keyword "ips" \
  --output results/DA0036

# Test individual stages
node test_stage1.js data/DA0036 results/test --organism=mouse
node test_stage2.js data/DA0036 results/test --organism=mouse
node test_stage3.js data/DA0036 results/test --organism=mouse
node test_stage4.js data/DA0036 results/test --organism=mouse

# Run with single agent (for  baseline)
node bin/geneexpert.js analyze data/DA0036 \
  --staged \
  --single-agent gpt5.2 \
  --output results/single_agent

# Run without agents (for  baseline)
node bin/geneexpert.js analyze data/DA0036 \
  --staged \
  --force-automation \
  --output results/no_agent

# Check logs
tail -f results/my_analysis/stage*_log.txt
```

---

## Research Contribution

**Hypothesis:** Staged multi-agent collaboration with user-in-loop tracking reduces errors 40%+ compared to single-agent or no-agent baseline systems.

**Novel Contributions:**

1. **Staged Multi-Agent Architecture**
   - 4 decision checkpoints throughout pipeline
   - Early issue detection (catch problems in Stage 2, not Stage 4)
   - Progressive refinement (remove outliers before DE analysis)

2. **Parallel vs Sequential Multi-Agent Comparison**
   - Parallel mode: Independent voting by all agents simultaneously
   - Sequential mode: Information passing (GPT-5.2 → Gemini → Claude)
   - Research question: Does sequential synthesis improve decisions or does parallel independence reduce errors?
   - Empirical comparison of two multi-agent architectures

3. **User-in-Loop Tracking**
   - JSON logs track when user input required
   - Enables empirical analysis of human-AI collaboration
   - Metrics: autonomy rate, agreement rate, impact on success

4. **Multi-Model Consensus**
   - Different foundation models (GPT-5.2, Claude Sonnet 4.5, Gemini Pro)
   - Voting system with confidence quantification
   - Disagreement signals uncertainty, triggers user input

5. **Condorcet's Jury Theorem Validation**
   - If each agent >50% accurate, ensemble >single agent
   - Empirical measurement of error reduction
   - Statistical significance testing

6. **Decision-Level Evaluation**
   - Not just final success/failure
   - Track correctness of intermediate decisions
   - Measure impact of each stage's agent review
   - Evaluate error propagation in sequential mode

---

## Evaluation Datasets

The system is evaluated on 9 RNA-seq datasets covering clean signals, batch effects, and contamination scenarios:

### Clean Biological Signal Datasets

**1. GSE52778 - Human Dexamethasone Response**
- Organism: Human cell lines (N61311, N052611, N080611, N061011)
- Design: Untreated vs dexamethasone-treated (4 vs 4 samples)
- Platform: Illumina paired-end RNA-seq
- Citation: Rüegger et al., Genome Biology, 2015
- URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778

**2. GSE114845 - Mouse Sleep Deprivation (Paired)**
- Organism: Mouse cortex
- Design: Control vs sleep-deprived, biologically paired by mouse (5 vs 5 samples)
- Platform: Illumina HiSeq 2500, single-end 101bp
- Citation: Bellesi et al., Journal of Neuroscience, 2018
- URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114845

**3. GSE113754 - Mouse Sleep Deprivation (Wildtype)**
- Organism: Mouse prefrontal cortex
- Design: Wildtype homecage control vs sleep deprivation (5 vs 5 samples)
- Platform: Illumina HiSeq 2500, paired-end
- Citation: Ingiosi et al., eLife, 2019 (PMID: 30973326)
- URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE113754

### Batch Effect Datasets

**4. GSE141496 - HeLa Technical Heterogeneity**
- Organism: Human HeLa cell line
- Design: Same biological condition across multiple lanes/runs (14 samples)
- Purpose: Stress-test for false-positive biological signals from technical variation
- Citation: Chen et al., Nature Biotechnology, 2019 (PMID: 30617342)
- URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE141496

**5. GSE47774 - SEQC Multi-Site Reference**
- Organism: Human reference RNA samples
- Design: Same RNA sequenced across multiple sites, flowcells, lanes (24 samples selected)
- Purpose: Multi-batch technical variation benchmark
- Citation: SEQC/MAQC-III Consortium, Nature Biotechnology, 2014
- URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE47774

### Contamination Source Datasets

**6. GSE193658 - Human Lab Experiment**
- Organism: Human cell line
- Design: Laboratory perturbation experiments
- Purpose: Source for synthetic contamination generation
- URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193658

**7. GSE48151 - E. coli Reference**
- Organism: Escherichia coli
- Design: Bacterial transcriptome under controlled conditions
- Purpose: Non-host organism for synthetic contamination (50% and 70% mixing)
- URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE48151

---

## License

MIT

---

## Links

- **GitHub:** https://github.com/Mituvinci/geneexpert-mcp
- **Status:** Staged architecture complete with Sequential Chain mode
- **Experiments:** 45 total (5 systems × 9 datasets) for ICML 2026
