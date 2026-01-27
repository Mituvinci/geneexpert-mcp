# GeneXpert-BioMCP : Poly-Foundational LLM Agent Orchestration for RNA-seq Upstream Intelligence

**Staged Multi-Agent Pipeline** where GPT-5.2, Claude Sonnet 4.5, and Gemini Pro collaborate at decision checkpoints throughout RNA-seq analysis.

**Key Innovation:** Dual-pipeline architecture (bulk RNA-seq 4-stage, scRNA-seq 7-stage) with agent review checkpoints - agents validate each stage, detect issues early, and decide whether to proceed or adjust the approach.

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
Seurat: v5.0.0 (for scRNA-seq analysis)
```

### 2. Configure API Keys

```bash
cp .env.example .env
nano .env  # Add your OpenAI, Anthropic, Google API keys
```

### 3. Run Analysis (Staged Architecture)

**Bulk RNA-seq Analysis:**
```bash
node bin/geneexpert.js analyze data/your_input_folder \
  --staged \
  --organism mouse \
  --comparison "stroke_vs_control" \
  --control-keyword "cont" \
  --treatment-keyword "ips" \
  --output results/your_output_folder
```

**Bulk RNA-seq Pipeline Flow:**
```
+---------------------------------------+
|  Input: FASTQ Files                   |
+---------------------------------------+
                  |
                  v
+---------------------------------------+
|  Stage 1: FASTQ Validation            |
|  - validate_fastq.sh                  |
|  - fastqc quality metrics             |
+---------------------------------------+
                  |
                  v
    +-------------------------------+
    |  AGENT CHECKPOINT #1          |
    |  Decision: PASS/WARN/FAIL     |
    +-------------------------------+
                  |
                  v
+---------------------------------------+
|  Stage 2: Alignment + QC              |
|  - subread-align (BAM files)          |
|  - alignment_qc_screening.R           |
+---------------------------------------+
                  |
                  v
    +--------------------------------------+
    |  AGENT CHECKPOINT #2                 |
    |  Decision: PASS_ALL/REMOVE/ABORT     |
    +--------------------------------------+
                  |
                  v
+---------------------------------------+
|  Stage 3: Quantification + PCA/QC     |
|  - featureCounts                      |
|  - RPKM normalization                 |
|  - PCA analysis                       |
+---------------------------------------+
                  |
                  v
    +------------------------------------------+
    |  AGENT CHECKPOINT #3                     |
    |  Decisions:                              |
    |  1. DE Method (simple/batch_effect)      |
    |  2. Outlier Action (KEEP/REMOVE)         |
    |  3. Batch Specification (auto/paired)    |
    +------------------------------------------+
                  |
                  v
+---------------------------------------+
|  Stage 4: Differential Expression     |
|  - edgeR analysis                     |
|  - merge_results.R                    |
+---------------------------------------+
                  |
                  v
    +--------------------------------------+
    |  AGENT CHECKPOINT #4                 |
    |  Decision: APPROVE/REANALYZE         |
    +--------------------------------------+
                  |
                  v
+---------------------------------------+
|  Output: DE Results + Excel Report    |
+---------------------------------------+
```
*User input tracked when agents can't reach consensus (critical for research evaluation)*

**Single-cell RNA-seq Analysis:**
```bash
node bin/scrna_geneexpert.js analyze data/scRNA_data/pbmc_healthy_human \
  --output results/scRNA_pbmc \
  --organism human \
  --verbose
```

**scRNA-seq Pipeline Flow (7 stages, 4 agent checkpoints):**
```
+---------------------------------------+
|  Input: 10x Genomics Data             |
|  (H5 / Directory / CSV)               |
+---------------------------------------+
                  |
                  v
+---------------------------------------+
|  Stage 1: Load + QC Metrics           |
|  - Load 10x data                      |
|  - Compute QC metrics                 |
|  - Auto-proceed (no agents)           |
+---------------------------------------+
                  |
                  v
+---------------------------------------+
|  Stage 2: QC Filtering                |
|  - Generate QC distribution plots     |
+---------------------------------------+
                  |
                  v
    +----------------------------------------+
    |  AGENT CHECKPOINT #1                   |
    |  Decision: SET_THRESHOLDS/USE_DEFAULT  |
    |  (nFeature_min/max, percent_mt_max)    |
    +----------------------------------------+
                  |
                  v
+---------------------------------------+
|  Stage 3A: Cell Cycle Scoring         |
|  - SCTransform normalization          |
|  - HVG selection                      |
|  - Cell cycle phase detection         |
+---------------------------------------+
                  |
                  v
    +----------------------------------------+
    |  AGENT CHECKPOINT #2                   |
    |  Decision: REMOVE_CELL_CYCLE/SKIP      |
    |  (Review PC-cycle correlation)         |
    +----------------------------------------+
                  |
                  v
+---------------------------------------+
|  Stage 3B: Cell Cycle Regression      |
|  - Conditional execution:             |
|    * REMOVE: Regress out S/G2M scores |
|    * SKIP: Scale without regression   |
|  - Auto-proceed based on 3A decision  |
+---------------------------------------+
                  |
                  v
+---------------------------------------+
|  Stage 4: PCA                         |
|  - Run PCA on HVGs                    |
|  - Generate elbow plot                |
+---------------------------------------+
                  |
                  v
    +----------------------------------------+
    |  AGENT CHECKPOINT #3                   |
    |  Decision: USE_DEFAULT/SELECT_PC_RANGE |
    |  (Select optimal # of PCs)             |
    +----------------------------------------+
                  |
                  v
+---------------------------------------+
|  Stage 5: Clustering + Markers        |
|  - Graph-based clustering             |
|  - UMAP visualization                 |
|  - Marker gene discovery              |
+---------------------------------------+
                  |
                  v
    +-------------------------------------------+
    |  AGENT CHECKPOINT #4                      |
    |  Decision: ACCEPT/ADJUST_RESOLUTION/FLAG  |
    |  (Validate cluster quality)               |
    +-------------------------------------------+
                  |
                  v
+---------------------------------------+
|  Output: Seurat Object + Clusters     |
|  + Markers + UMAP                     |
+---------------------------------------+
```
*Cell cycle correction (Stage 3A/3B) is CRITICAL for proper downstream clustering*

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
+---------------------------------------------------------------+
|  Stage 1: FASTQ Validation                                    |
|  Generate script → Execute → Parse results → Format for agents|
+---------------------------------------------------------------+
                            |
                            v
              +----------------------------------+
              |  3 Agents Review:                |
              |  [GPT-5.2 Stats]:                |
              |    "Read depth sufficient"       |
              |  [Claude Pipeline]:              |
              |    "File integrity good"         |
              |  [Gemini Biology]:               |
              |    "Quality acceptable"          |
              +----------------------------------+
                            |
                            v
              +----------------------------------+
              |  Consensus: PASS (3/3 agree)     |
              |  → Proceed to Stage 2            |
              +----------------------------------+
                            |
                            v
+---------------------------------------------------------------+
|  Stage 2: Alignment + QC                                      |
|  Generate script → Execute → Parse QC → Format for agents     |
+---------------------------------------------------------------+
                            |
                            v
              +----------------------------------+
              |  3 Agents Review:                |
              |  [GPT-5.2 Stats]:                |
              |    "Mapping rate 85% excellent"  |
              |  [Claude Pipeline]:              |
              |    "All samples >80%, no contam" |
              |  [Gemini Biology]:               |
              |    "Genome match confirmed"      |
              +----------------------------------+
                            |
                            v
              +----------------------------------+
              |  Consensus: PASS_ALL             |
              |  → Proceed to Stage 3            |
              +----------------------------------+
                            |
                            v
+---------------------------------------------------------------+
|  Stage 3: Quantification + PCA/QC                             |
|  Generate script → Execute → Parse PCA → Format for agents    |
+---------------------------------------------------------------+
                            |
                            v
              +----------------------------------+
              |  3 Agents Review:                |
              |  [GPT-5.2 Stats]:                |
              |    "No statistical outliers"     |
              |  [Claude Pipeline]:              |
              |    "PCA shows good separation"   |
              |  [Gemini Biology]:               |
              |    "Groups cluster properly"     |
              +----------------------------------+
                            |
                            v
              +----------------------------------+
              |  Consensus: simpleEdger, KEEP_ALL|
              |  → Proceed to Stage 4            |
              +----------------------------------+
                            |
                            v
+---------------------------------------------------------------+
|  Stage 4: Differential Expression                             |
|  Generate script → Execute → Parse DE → Format for agents     |
+---------------------------------------------------------------+
                            |
                            v
              +----------------------------------+
              |  3 Agents Review:                |
              |  [GPT-5.2 Stats]:                |
              |    "1,234 DEGs at FDR<0.05"      |
              |  [Claude Pipeline]:              |
              |    "LogFC range appropriate"     |
              |  [Gemini Biology]:               |
              |    "Top genes meaningful"        |
              +----------------------------------+
                            |
                            v
              +----------------------------------+
              |  Consensus: APPROVE              |
              |  → Analysis Complete!            |
              +----------------------------------+
```

**This is TRUE staged validation:** Issues caught early, decisions made with full context, user input tracked when needed!

---

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

### Auto-Resolution System (3-Tier Escalation)

When agents disagree, GeneXpert uses a tiered escalation system to resolve conflicts intelligently:

**Tier 1: Minor Disagreement (score < 0.3)**
- Auto-resolve using median/average of agent recommendations
- Example: PCs [25, 28, 30] → Use median 28
- Fast resolution, no user interruption

**Tier 2: Moderate Disagreement (0.3 ≤ score < 0.6)**
- Use highest confidence agent's recommendation
- Agents provide confidence scores: HIGH/MEDIUM/LOW
- Trust the most confident expert
- Example: Agent A (HIGH confidence): REMOVE, Agent B (LOW): SKIP → Use REMOVE

**Tier 3: Major Disagreement (score ≥ 0.6)**
- Escalate to user decision
- Fundamental disagreement signals uncertainty
- User input logged for research evaluation

**Auto-Resolution Modes:**
```bash
# Tiered escalation (default, recommended)
--auto-resolve auto

# Always use median (fastest, no escalation)
--auto-resolve median

# Always use highest confidence agent
--auto-resolve confidence

# Always escalate to user (most conservative)
--auto-resolve user
```

**Disagreement Scoring:**
- Numeric decisions: Coefficient of variation (std / mean)
- Binary decisions: Entropy-based scoring
- Categorical decisions: Simpson diversity index

**Implementation:**
- `/src/utils/auto_resolver.js` - Core escalation logic
- `/src/utils/consensus_helper.js` - Disagreement scoring, confidence extraction

---

### Role Swapping (Ablation Study Feature)

GeneXpert supports flexible role assignment to test which model excels in which role:

**Default Assignment:**
- GPT-5.2 → Stats Agent
- Claude Sonnet 4.5 → Pipeline Agent
- Gemini 2.0 Flash → Biology Agent

**Custom Role Assignment:**
```bash
# Swap roles for ablation study
node bin/geneexpert.js analyze <dataset> \
  --staged \
  --gpt-role biology \
  --claude-role stats \
  --gemini-role pipeline \
  --output results/role_swap

# Test all 6 permutations (3! = 6 possible assignments)
```

**How It Works:**
- Each agent checkpoint uses role-specific prompts
- Prompts dynamically assigned based on `--<model>-role` flag
- Agent responses formatted according to assigned role
- Enables research question: "Which model is best for statistical vs biological reasoning?"

**Supported Roles:**
- `stats`: Statistical validation, threshold selection
- `pipeline`: Technical feasibility, Seurat best practices
- `biology`: Biological interpretation, cell type knowledge

**Implementation:**
- `/src/coordinator/orchestrator.js` - `applyRoleSwapping()` function
- `/src/config/scrna_stage_prompts.js` - Role-specific prompt library

---

## Key Features

### Implemented & Tested:

1. **Dual-Pipeline Architecture with Agent Checkpoints**
   - Bulk RNA-seq: 4-stage pipeline, 4 agent checkpoints (Validation → Alignment → Quantification → DE)
   - scRNA-seq: 7-stage pipeline, 4 agent checkpoints (Load → Filter → Cell Cycle Scoring → Cell Cycle Regression → PCA → Cluster)
   - Each stage: Generate → Execute → Parse → Format → Agent Review → Decision
   - Early issue detection (catch bad samples early, not at the end)
   - Progressive refinement (remove outliers, correct cell cycle effects)
   - Conditional execution paths (Stage 3B branches based on 3A agent decision)

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

9. **3-Tier Auto-Resolution System**
   - Minor disagreements (score < 0.3): Auto-resolve using median/average
   - Moderate disagreements (0.3 ≤ score < 0.6): Use highest confidence agent
   - Major disagreements (score ≥ 0.6): Escalate to user
   - 4 modes: `auto` (tiered), `median` (fast), `confidence` (trust expert), `user` (conservative)
   - Disagreement scoring: CV for numeric, entropy for binary, diversity for categorical
   - Implementation: `/src/utils/auto_resolver.js`, `/src/utils/consensus_helper.js`

10. **Role Swapping for Ablation Studies**
   - Flexible role assignment: any model can play any role (stats/pipeline/biology)
   - Default: GPT-5.2 (Stats), Claude (Pipeline), Gemini (Biology)
   - Custom assignment via CLI: `--gpt-role`, `--claude-role`, `--gemini-role`
   - Role-specific prompt library (46KB of stage-specific prompts)
   - Enables research question: "Which model excels in which reasoning task?"
   - All 6 role permutations supported (3! = 6 combinations)

11. **Cell Cycle Analysis with Agent Decision**
   - Stage 3A: Agents review cell cycle correlation plots, decide REMOVE or SKIP
   - Stage 3B: Conditional execution based on 3A decision (regression vs skip)
   - Organism-specific marker genes (human: cc.genes.updated.2019, mouse: mm.cc.genes)
   - Correlation analysis between PC1 and S.Score/G2M.Score
   - Before/after plots for validation
   - Critical for preventing cell cycle-driven false clustering

---

## Experimental Modes (Research Evaluation)

GeneExpert supports 10 experimental systems for ICML 2026 evaluation:

### All Possible Role Combinations (6 permutations):

**System 1: GPT=stats, Claude=pipeline, Gemini=biology (DEFAULT)**
```bash
node bin/geneexpert.js analyze <dataset> --staged \
  --gpt-role stats --claude-role pipeline --gemini-role biology \
  --output results/perm1
```

**System 2: GPT=stats, Claude=biology, Gemini=pipeline**
```bash
node bin/geneexpert.js analyze <dataset> --staged \
  --gpt-role stats --claude-role biology --gemini-role pipeline \
  --output results/perm2
```

**System 3: GPT=pipeline, Claude=stats, Gemini=biology**
```bash
node bin/geneexpert.js analyze <dataset> --staged \
  --gpt-role pipeline --claude-role stats --gemini-role biology \
  --output results/perm3
```

**System 4: GPT=pipeline, Claude=biology, Gemini=stats**
```bash
node bin/geneexpert.js analyze <dataset> --staged \
  --gpt-role pipeline --claude-role biology --gemini-role stats \
  --output results/perm4
```

**System 5: GPT=biology, Claude=stats, Gemini=pipeline**
```bash
node bin/geneexpert.js analyze <dataset> --staged \
  --gpt-role biology --claude-role stats --gemini-role pipeline \
  --output results/perm5
```

**System 6: GPT=biology, Claude=pipeline, Gemini=stats**
```bash
node bin/geneexpert.js analyze <dataset> --staged \
  --gpt-role biology --claude-role pipeline --gemini-role stats \
  --output results/perm6
```

### Baseline Systems (4):

**System 7: No-Agent (Template Only)**
```bash
node bin/geneexpert.js analyze <dataset> --staged \
  --force-automation \
  --output results/no_agent
```
Baseline with no agent review, template-based decisions only.

**System 8: Single-LLM Multi-Agent (Claude Opus handles all roles)**
```bash
node bin/geneexpert.js analyze <dataset> --staged \
  --single-agent claude \
  --output results/single_claude
```
Claude called 3 times with different role prompts (stats, pipeline, biology).

**System 9: Single-LLM Multi-Agent (GPT-5.2 handles all roles)**
```bash
node bin/geneexpert.js analyze <dataset> --staged \
  --single-agent gpt5.2 \
  --output results/single_gpt
```
GPT-5.2 called 3 times with different role prompts (stats, pipeline, biology).

**System 10: Single-LLM Multi-Agent (Gemini Pro handles all roles)**
```bash
node bin/geneexpert.js analyze <dataset> --staged \
  --single-agent gemini \
  --output results/single_gemini
```
Gemini called 3 times with different role prompts (stats, pipeline, biology).

**Total Experiment Count:**
- 10 systems × 7 bulk RNA-seq datasets = 70 bulk analyses
- 10 systems × 6 scRNA-seq datasets = 60 scRNA analyses
- **Total: 130 analyses**

**Evaluation Metrics:**
- Decision accuracy (correct decisions / total)
- Error reduction (baseline errors - multi-agent errors) / baseline errors
- Success rate (successful analyses / total)
- Cost efficiency (total cost / successes)
- Inter-agent agreement (Cohen's kappa for multi-agent modes)
- User input frequency (% decisions requiring user input)
- Error propagation rate (sequential mode - measures if later agents correct or amplify first agent's errors)
- Role performance analysis (which model excels at which role: stats, pipeline, biology)
- Reasoning quality (qualitative assessment)

**Evaluation Scripts:**
```bash
# Calculate all metrics
node bin/evaluate.js metrics

# Compare parallel vs sequential
node bin/evaluate.js compare

# Export results to CSV
node bin/evaluate.js metrics --output results.csv

# Analyze costs from experiment logs
node bin/analyze_costs.js generate --results experiments/results --output costs.csv

# Filter cost analysis by system
node bin/analyze_costs.js generate --results experiments/results --system multi-agent-parallel
```

---

## Project Structure

```
├── bin/
│   ├── geneexpert.js                    # Bulk RNA-seq CLI entry point
│   ├── scrna_geneexpert.js              # scRNA-seq CLI entry point
│   ├── evaluate.js                      # Evaluation metrics calculator
│   ├── analyze_costs.js                 # Cost analysis from experiment logs
│   └── json_to_csv.js                   # JSON to CSV converter
├── src/
│   ├── executor/
│   │   └── staged_executor.js           # 4-stage bulk RNA-seq orchestration
│   ├── scrna_executor/
│   │   └── scrna_executor.js            # 7-stage scRNA-seq orchestration
│   ├── stages/
│   │   ├── stage1_validation.js         # FASTQ validation
│   │   ├── stage2_alignment.js          # Alignment + QC
│   │   ├── stage3_quantification_qc.js  # Quantification + PCA QC
│   │   └── stage4_de_analysis.js        # Differential expression
│   ├── scrna_stages/
│   │   ├── stage1_load_qc.js            # Load 10x data + QC metrics
│   │   ├── stage2_filter_qc.js          # QC filtering (agents recommend thresholds)
│   │   ├── stage3a_cell_cycle_scoring.js    # Cell cycle scoring (Agent Checkpoint #2)
│   │   ├── stage3b_cell_cycle_regression.js # Cell cycle regression or skip (conditional)
│   │   ├── stage4_pca.js                # PCA (agents select PC range)
│   │   └── stage5_cluster_markers.js    # Clustering + marker genes
│   ├── config/
│   │   ├── stage_prompts.js             # Bulk RNA-seq agent prompts
│   │   └── scrna_stage_prompts.js       # scRNA-seq agent prompts
│   ├── coordinator/
│   │   └── orchestrator.js              # Multi-agent coordination
│   ├── utils/
│   │   ├── logger.js                    # JSON logging with user input tracking
│   │   ├── user_input.js                # User decision handling
│   │   ├── llm_clients.js               # OpenAI, Anthropic, Google APIs
│   │   ├── auto_resolver.js             # 3-tier auto-resolution system
│   │   └── consensus_helper.js          # Disagreement scoring, confidence extraction
│   └── pipeline/                        # Old monolithic architecture (legacy)
└── README.md                            # This file
```

---

## Cost Estimate

### Bulk RNA-seq (4 Agent Checkpoints):
- 3 agents × 4 checkpoints = 12 agent API calls per analysis
- ~$0.01 per agent call
- **Total agent cost: ~$0.12 per complete bulk RNA-seq analysis**

### scRNA-seq (4 Agent Checkpoints):
- 3 agents × 4 checkpoints = 12 agent API calls per analysis
- ~$0.01 per agent call
- **Total agent cost: ~$0.12 per complete scRNA-seq analysis**

### Script Execution:
- All bioinformatics tools run locally (subread-align, featureCounts, edgeR, Seurat)
- **Execution cost: $0**

**Total: ~$0.12 per analysis with full multi-agent validation (both pipelines)**

---

## Research Contribution

**Hypothesis:** Staged multi-agent collaboration with user-in-loop tracking reduces errors 40%+ compared to single-agent or no-agent baseline systems.

**Novel Contributions:**

1. **Staged Multi-Agent Architecture**
   - 4 decision checkpoints throughout pipeline (bulk RNA-seq)
   - 4 decision checkpoints for scRNA-seq pipeline
   - Early issue detection (catch problems in Stage 2, not Stage 4)
   - Progressive refinement (remove outliers before DE analysis)

2. **Parallel vs Sequential Multi-Agent Comparison**
   - Parallel mode: Independent voting by all agents simultaneously
   - Sequential mode: Information passing (GPT-5.2 → Gemini → Claude)
   - Research question: Does sequential synthesis improve decisions or does parallel independence reduce errors?
   - Empirical comparison of two multi-agent architectures

3. **Comprehensive Role Swapping Analysis**
   - All 6 possible role permutations tested (3! = 6 combinations)
   - Empirical comparison: GPT-5.2 vs Claude vs Gemini in stats/pipeline/biology roles
   - Research questions:
     - Which model excels at statistical reasoning?
     - Which model excels at biological interpretation?
     - Which model excels at pipeline/technical decisions?
     - Is model-role alignment more important than raw model capability?
   - Enables identification of optimal multi-agent team composition

4. **User-in-Loop Tracking**
   - JSON logs track when user input required
   - Enables empirical analysis of human-AI collaboration
   - Metrics: autonomy rate, agreement rate, impact on success

5. **Multi-Model Consensus**
   - Different foundation models (GPT-5.2, Claude Sonnet 4.5, Gemini Pro)
   - Voting system with confidence quantification
   - Disagreement signals uncertainty, triggers user input

6. **Condorcet's Jury Theorem Validation**
   - If each agent >50% accurate, ensemble >single agent
   - Empirical measurement of error reduction
   - Statistical significance testing across 140 analyses

7. **Decision-Level Evaluation**
   - Not just final success/failure
   - Track correctness of intermediate decisions
   - Measure impact of each stage's agent review
   - Evaluate error propagation in sequential mode

8. **3-Tier Auto-Resolution System**
   - Intelligent disagreement resolution without always escalating to user
   - Tiered escalation based on disagreement severity (minor/moderate/major)
   - Confidence-weighted voting for moderate disagreements
   - Enables higher autonomy while maintaining safety
   - Research question: Can intelligent auto-resolution maintain accuracy while reducing user burden?

9. **Conditional Execution Architecture**
   - Stage 3B execution path determined by Stage 3A agent decision (scRNA)
   - Cell cycle regression vs skip based on agent consensus
   - Demonstrates adaptive pipeline with agent-driven branching
   - Research contribution: Agent decisions shape computational workflow, not just approve/reject

10. **Dual-Pipeline Validation**
   - Same multi-agent system validated on TWO distinct bioinformatics pipelines
   - Bulk RNA-seq: differential expression analysis
   - scRNA-seq: clustering and cell type identification
   - Tests generalizability of multi-agent approach across different computational biology domains

---

## Evaluation Datasets

The system is evaluated on 13 RNA-seq datasets (7 bulk + 6 single-cell) covering clean signals, batch effects, contamination, and cell cycle challenges:

### Bulk RNA-seq Datasets (7 datasets)

**1. GSE52778 - Human Dexamethasone Response**
- Organism: Human cell lines (N61311, N052611, N080611, N061011)
- Design: Untreated vs dexamethasone-treated (4 vs 4 samples)
- Platform: Illumina paired-end RNA-seq
- Category: Clean biological signal
- Citation: Rüegger et al., Genome Biology, 2015
- URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778

**2. GSE114845 - Mouse Sleep Deprivation (Paired)**
- Organism: Mouse cortex
- Design: Control vs sleep-deprived, biologically paired by mouse (5 vs 5 samples)
- Platform: Illumina HiSeq 2500, single-end 101bp
- Category: Clean biological signal
- Citation: Bellesi et al., Journal of Neuroscience, 2018
- URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114845

**3. GSE113754 - Mouse Sleep Deprivation (Wildtype)**
- Organism: Mouse prefrontal cortex
- Design: Wildtype homecage control vs sleep deprivation (5 vs 5 samples)
- Platform: Illumina HiSeq 2500, paired-end
- Category: Clean biological signal
- Citation: Ingiosi et al., eLife, 2019 (PMID: 30973326)
- URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE113754

**4. GSE141496 - HeLa Technical Heterogeneity**
- Organism: Human HeLa cell line
- Design: Same biological condition across multiple lanes/runs (14 samples)
- Platform: Illumina paired-end RNA-seq
- Category: Batch effect
- Purpose: Stress-test for false-positive biological signals from technical variation
- Citation: Chen et al., Nature Biotechnology, 2019 (PMID: 30617342)
- URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE141496

**5. GSE47774 - SEQC Multi-Site Reference**
- Organism: Human reference RNA samples
- Design: Same RNA sequenced across multiple sites, flowcells, lanes (24 samples selected)
- Platform: Illumina paired-end RNA-seq
- Category: Batch effect
- Purpose: Multi-batch technical variation benchmark
- Citation: SEQC/MAQC-III Consortium, Nature Biotechnology, 2014
- URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE47774

**6. GSE193658 - Human In-House Experiment**
- Organism: Human cell line
- Design: Laboratory perturbation experiments
- Platform: Illumina paired-end RNA-seq
- Category: In-house data
- Purpose: Real-world dataset generated in-house and published
- URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193658

**7. GSE114845_CONTAM70 - 70% E. coli Contamination**
- Base: GSE114845 single-end data + 70% synthetic E. coli contamination
- Contamination source: E. coli reads from GSE48151 (GSM1170025)
- Method: seqkit-based synthetic contamination (70% E. coli, 30% host)
- Category: Contamination
- Purpose: Test agent ability to detect severe cross-species contamination

### Single-cell RNA-seq Datasets (6 datasets)

**1. REH parental - Human Leukemia Cell Line**
- Organism: Human (leukemia cell line)
- Platform: 10x Chromium Multiome
- Category: Clean single-cell
- Accession: GSE293316
- URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE293316

**2. SUP-B15 parental - Human Leukemia Cell Line**
- Organism: Human (leukemia cell line)
- Platform: 10x Chromium Multiome
- Category: Clean single-cell
- Accession: GSE293316
- URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE293316

**3. PBMC Healthy Human - Peripheral Blood Mononuclear Cells**
- Organism: Human (healthy donor)
- Platform: 10x Chromium (10,194 cells)
- Category: Clean single-cell
- Cell cycle labels: 2,431 cells (1,963 G1, 127 S, 341 G2M) consensus labeled
- URL: https://www.10xgenomics.com/resources/datasets

**4. Mouse Brain Cells (Healthy) - 10k Brain Cells from E18 Mouse**
- Organism: Mouse (embryonic day 18)
- Platform: 10x Chromium v3 chemistry (11,843 cells)
- Category: Clean single-cell
- Cell cycle labels: 5,524 cells (3,830 G1, 1,116 S, 578 G2M) consensus labeled
- URL: https://www.10xgenomics.com/resources/datasets

**5. GSE75748 - Human Embryonic Stem Cells (hPSC/hESC)**
- Organism: Human embryonic stem cells
- Design: 1,776 cells across progenitor states and differentiation trajectory
- Category: Cell cycle challenge
- Purpose: Strong cell-cycle activity, developmental transitions
- Usage: Training/evaluating cell-cycle phase prediction
- URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75748

**6. GSE146773 - Human U-2 OS Cells with FUCCI Reporter**
- Organism: Human osteosarcoma cell line (U-2 OS)
- Design: FUCCI (Fluorescent Ubiquitination-based Cell Cycle Indicator) reporter
- Category: Cell cycle ground truth
- Purpose: Cell cycle phase ground truth validation
- URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146773

---

## License

MIT

---

## Links

- **GitHub:** https://github.com/Mituvinci/geneexpert-mcp
- **Status:** Dual-pipeline architecture complete (Bulk 4-stage + scRNA 7-stage) with 4 agent checkpoints each
- **New Features:** 3-tier auto-resolution, role swapping, conditional execution (Stage 3A/3B branching)
- **Experiments:** 45 bulk + scRNA experiments (5 systems × 9 datasets × 2 pipelines) for ICML 2026
