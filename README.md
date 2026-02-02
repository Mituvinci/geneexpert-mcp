# GeneExpert: Poly-Foundational LLM Multi-Agent System for RNA-seq Analysis

**Staged Multi-Agent Pipeline** where GPT-5.2, Claude Sonnet 4.5, and Gemini 2.0 Flash collaborate at decision checkpoints throughout RNA-seq analysis.

**Key Innovation:** Dual-pipeline architecture (bulk RNA-seq 4-stage, scRNA-seq 5-stage) with multi-agent consensus at each checkpoint. Agents validate each stage, detect issues early, and decide whether to proceed or adjust the approach.

**Research Goal:** Demonstrate that multi-agent collaboration reduces analysis errors by 40%+ compared to single-agent and no-agent baselines through staged validation and consensus-based decision making.

---

## Quick Start

### Installation

```bash
# Requirements
Node: v18.20.8
npm: 10.8.2
Python: 3.10.19
R: 4.3.3
Subread-align: v2.1.1
FastQC: v0.12.1
SAMtools: 1.22.1
Seurat: v5.0.0 (for scRNA-seq)

# Install dependencies
npm install

# Configure API keys
cp .env.example .env
nano .env  # Add OpenAI, Anthropic, Google API keys
```

### Run Analysis

**Bulk RNA-seq:**
```bash
node bin/geneexpert.js analyze data/your_dataset \
  --staged \
  --organism mouse \
  --comparison "condition1_vs_condition2" \
  --control-keyword "control" \
  --treatment-keyword "treatment" \
  --output results/output_folder
```

**Single-cell RNA-seq:**
```bash
node bin/scrna_geneexpert.js analyze data/scRNA_data/your_dataset \
  --output results/scRNA_output \
  --organism human \
  --verbose
```

---

## Evaluation Datasets

The system is evaluated on **14 RNA-seq datasets** (7 bulk + 7 single-cell) covering clean signals, batch effects, contamination, and cell cycle challenges.

### Bulk RNA-seq Datasets (7 datasets)

For reproducibility, we provide the exact GSM sample IDs downloaded from GEO and their renamed filenames used in our experiments.

#### 1. GSE52778 - Human Dexamethasone Response

- **Organism:** Human cell lines
- **Design:** Untreated vs dexamethasone-treated (4 vs 4 samples)
- **Platform:** Illumina paired-end RNA-seq
- **Category:** Clean biological signal
- **Citation:** Rüegger et al., Genome Biology, 2015
- **URL:** https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778

**Sample Mapping (GSM ID → Renamed Filename):**

Control group (untreated):
- GSM1275862 → N61311_untreated_1
- GSM1275866 → N052611_untreated_2
- GSM1275870 → N080611_untreated_3
- GSM1275874 → N061011_untreated_4

Treatment group (Dex):
- GSM1275863 → N61311_Dex_1
- GSM1275867 → N052611_Dex_2
- GSM1275871 → N080611_Dex_3
- GSM1275875 → N061011_Dex_4

---

#### 2. GSE114845 - Mouse Sleep Deprivation (Paired)

- **Organism:** Mouse cortex
- **Design:** Control vs sleep-deprived, biologically paired by mouse (5 vs 5 samples)
- **Platform:** Illumina HiSeq 2500, single-end 101bp
- **Category:** Clean biological signal
- **Citation:** Bellesi et al., Journal of Neuroscience, 2018
- **URL:** https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114845

**Sample Mapping (GSM ID → Renamed Filename):**

Control group:
- GSM3151467 → 05_Cortex_Control_R1_001
- GSM3151469 → 100_Cortex_Control_R1_001
- GSM3151471 → 101_Cortex_Control_R1_001
- GSM3151473 → 103_Cortex_Control_R1_001
- GSM3151475 → 29_Cortex_Control_R1_001

Sleep-deprived group:
- GSM3151468 → 05_Cortex_SleepDeprived_R1_001
- GSM3151470 → 100_Cortex_SleepDeprived_R1_001
- GSM3151472 → 101_Cortex_SleepDeprived_R1_001
- GSM3151474 → 103_Cortex_SleepDeprived_R1_001
- GSM3151476 → 29_Cortex_SleepDeprived_R1_001

---

#### 3. GSE113754 - Mouse Sleep Deprivation (Wildtype)

- **Organism:** Mouse prefrontal cortex
- **Design:** Wildtype homecage control vs sleep deprivation (5 vs 5 samples)
- **Platform:** Illumina HiSeq 2500, paired-end
- **Category:** Clean biological signal
- **Citation:** Ingiosi et al., eLife, 2019 (PMID: 30973326)
- **URL:** https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE113754

**Sample Mapping (GSM ID → Renamed Filename):**

Control group (homecage):
- GSM3118997 → WT_Homecage_Control_rep1
- GSM3118998 → WT_Homecage_Control_rep2
- GSM3118999 → WT_Homecage_Control_rep3
- GSM3119000 → WT_Homecage_Control_rep4
- GSM3119001 → WT_Homecage_Control_rep5

Sleep-deprived group:
- GSM3119007 → WT_SleepDeprived_rep1
- GSM3119008 → WT_SleepDeprived_rep2
- GSM3119009 → WT_SleepDeprived_rep3
- GSM3119010 → WT_SleepDeprived_rep4
- GSM3119011 → WT_SleepDeprived_rep5

---

#### 4. GSE141496 - HeLa Technical Heterogeneity

- **Organism:** Human HeLa cell line
- **Design:** Same biological condition across multiple lanes/runs (8 samples: 4 control + 4 treatment)
- **Platform:** Illumina paired-end RNA-seq
- **Category:** Batch effect
- **Purpose:** Stress-test for false-positive biological signals from technical variation
- **Citation:** Chen et al., Nature Biotechnology, 2019 (PMID: 30617342)
- **URL:** https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE141496

**Sample Mapping (GSM ID → Renamed Filename):**

Control group:
- GSM4204829 → Control_24h_post_HU_0.25m_rep1
- GSM4204830 → Control_24h_post_HU_0.25m_rep2
- GSM4204833 → Control_rep1
- GSM4204834 → Control_rep2

Treatment group (TAB182 KD):
- GSM4204839 → TAB182_KD_24h_post_HU_0.25m_rep1
- GSM4204840 → TAB182_KD_24h_post_HU_0.25m_rep2
- GSM4204841 → TAB182_KD_rep1
- GSM4204842 → TAB182_KD_rep2

---

#### 5. GSE47774 - SEQC Multi-Site Reference

- **Organism:** Human reference RNA samples
- **Design:** Same RNA sequenced across multiple sites, flowcells, lanes (22 samples selected: 11 AC + 11 BC)
- **Platform:** Illumina paired-end RNA-seq
- **Category:** Batch effect
- **Purpose:** Multi-batch technical variation benchmark
- **Citation:** SEQC/MAQC-III Consortium, Nature Biotechnology, 2014
- **URL:** https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE47774

**Sample Mapping (GSM ID → Renamed Filename):**

Sample group AC:
- GSM1156797 → Batch_BGI_A_L01_AC
- GSM1156799 → Batch_BGI_A_L02_AC
- GSM1156801 → Batch_BGI_A_L03_AC
- GSM1156803 → Batch_BGI_A_L04_AC
- GSM1156813 → Batch_BGI_A2_L01_AC
- GSM1156815 → Batch_BGI_A2_L02_AC
- GSM1156817 → Batch_BGI_A2_L03_AC
- GSM1156819 → Batch_BGI_A2_L04_AC
- GSM1156829 → Batch_BGI_A3_L01_AC
- GSM1156833 → Batch_BGI_A3_L03_AC
- GSM1156835 → Batch_BGI_A3_L04_AC

Sample group BC:
- GSM1156798 → Batch_BGI_A_L01_BC
- GSM1156800 → Batch_BGI_A_L02_BC
- GSM1156802 → Batch_BGI_A_L03_BC
- GSM1156804 → Batch_BGI_A_L04_BC
- GSM1156814 → Batch_BGI_A2_L01_BC
- GSM1156816 → Batch_BGI_A2_L02_BC
- GSM1156818 → Batch_BGI_A2_L03_BC
- GSM1156820 → Batch_BGI_A2_L04_BC
- GSM1156830 → Batch_BGI_A3_L01_BC
- GSM1156834 → Batch_BGI_A3_L03_BC
- GSM1156836 → Batch_BGI_A3_L04_BC

---

#### 6. GSE193658 - Human In-House Experiment

- **Organism:** Human MM1.S cell line
- **Design:** Laboratory perturbation experiments (3 mono + 3 TSW)
- **Platform:** Illumina paired-end RNA-seq
- **Category:** In-house data
- **Purpose:** Real-world dataset generated in-house and published
- **URL:** https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193658

**Sample Mapping (GSM ID → Renamed Filename):**

Mono group:
- GSM5815699 → RNA_MM1S_MONO1
- GSM5815700 → RNA_MM1S_MONO2
- GSM5815701 → RNA_MM1S_MONO3

TSW group:
- GSM5815702 → RNA_MM1S_TSW1
- GSM5815703 → RNA_MM1S_TSW2
- GSM5815704 → RNA_MM1S_TSW3

---

#### 7. GSE114845_CONTAM70 - 70% E. coli Contamination

- **Base Dataset:** GSE114845 single-end data + 70% synthetic E. coli contamination
- **Contamination Source:** E. coli reads from GSE48151 (GSM1170025)
- **Method:** seqkit-based synthetic contamination (70% E. coli, 30% host)
- **Category:** Contamination challenge
- **Purpose:** Test agent ability to detect severe cross-species contamination

**Sample Mapping (GSM ID → Renamed Filename - same as GSE114845):**

Control group:
- GSM3151467 → 05_Cortex_Control_R1_001
- GSM3151469 → 100_Cortex_Control_R1_001
- GSM3151471 → 101_Cortex_Control_R1_001
- GSM3151473 → 103_Cortex_Control_R1_001
- GSM3151475 → 29_Cortex_Control_R1_001

Sleep-deprived group:
- GSM3151468 → 05_Cortex_SleepDeprived_R1_001
- GSM3151470 → 100_Cortex_SleepDeprived_R1_001
- GSM3151472 → 101_Cortex_SleepDeprived_R1_001
- GSM3151474 → 103_Cortex_SleepDeprived_R1_001
- GSM3151476 → 29_Cortex_SleepDeprived_R1_001

---

### Single-cell RNA-seq Datasets (7 datasets)

#### 1. REH Parental - Human Leukemia Cell Line

- **Organism:** Human (leukemia cell line)
- **Platform:** 10x Chromium Multiome
- **Category:** Clean single-cell
- **Accession:** GSE293316
- **URL:** https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE293316

---

#### 2. SUP-B15 Parental - Human Leukemia Cell Line

- **Organism:** Human (leukemia cell line)
- **Platform:** 10x Chromium Multiome
- **Category:** Clean single-cell
- **Accession:** GSE293316
- **URL:** https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE293316

---

#### 3. PBMC Healthy Human - Peripheral Blood Mononuclear Cells

- **Organism:** Human (healthy donor)
- **Platform:** 10x Chromium (10,194 cells)
- **Category:** Clean single-cell
- **Cell cycle labels:** 2,431 cells (1,963 G1, 127 S, 341 G2M) consensus labeled
- **URL:** https://www.10xgenomics.com/resources/datasets

---

#### 4. Mouse Brain Cells (Healthy) - 10k Brain Cells from E18 Mouse

- **Organism:** Mouse (embryonic day 18)
- **Platform:** 10x Chromium v3 chemistry (11,843 cells)
- **Category:** Clean single-cell
- **Cell cycle labels:** 5,524 cells (3,830 G1, 1,116 S, 578 G2M) consensus labeled
- **URL:** https://www.10xgenomics.com/resources/datasets

---

#### 5. GSE75748 - Human Embryonic Stem Cells (hPSC/hESC)

- **Organism:** Human embryonic stem cells
- **Design:** 1,776 cells across progenitor states and differentiation trajectory
- **Category:** Cell cycle challenge
- **Purpose:** Strong cell-cycle activity, developmental transitions
- **Usage:** Training/evaluating cell-cycle phase prediction
- **URL:** https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75748

---

#### 6. GSE146773 - Human U-2 OS Cells with FUCCI Reporter

- **Organism:** Human osteosarcoma cell line (U-2 OS)
- **Design:** FUCCI (Fluorescent Ubiquitination-based Cell Cycle Indicator) reporter
- **Category:** Cell cycle ground truth
- **Purpose:** Cell cycle phase ground truth validation
- **URL:** https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146773

---

#### 7. GSE64016 - Human ESCs with FUCCI Reporter

- **Organism:** Human embryonic stem cells (H9 hESCs)
- **Design:** FUCCI reporter system for cell cycle phase identification
- **Category:** Cell cycle ground truth
- **Purpose:** Gold standard for cell cycle phase classification
- **URL:** https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64016

---

## Experimental Setup

GeneExpert supports **10 experimental systems** for comprehensive evaluation:

### Systems Evaluated

**6 Role Permutation Configurations (Multi-Agent):**
1. GPT-5.2 (Stats) + Claude (Pipeline) + Gemini (Biology) - **DEFAULT**
2. GPT-5.2 (Stats) + Claude (Biology) + Gemini (Pipeline)
3. GPT-5.2 (Pipeline) + Claude (Stats) + Gemini (Biology)
4. GPT-5.2 (Pipeline) + Claude (Biology) + Gemini (Stats)
5. GPT-5.2 (Biology) + Claude (Stats) + Gemini (Pipeline)
6. GPT-5.2 (Biology) + Claude (Pipeline) + Gemini (Stats)

**4 Baseline Systems:**
7. **No-Agent:** Template-only automation (no LLM review)
8. **Single-LLM (Claude):** Claude called 3x with different role prompts
9. **Single-LLM (GPT-5.2):** GPT called 3x with different role prompts
10. **Single-LLM (Gemini):** Gemini called 3x with different role prompts

### Running Experiments

**Run All 10 Systems on One Bulk RNA Dataset:**
```bash
bash bin/run_bulk_rna.sh data/your_dataset
```

**Run All 10 Systems on One scRNA Dataset:**
```bash
bash bin/run_sc_rna.sh data/scRNA_data/your_dataset
```

**Run Individual System (Example: System 1 - Default Roles):**
```bash
node bin/geneexpert.js analyze data/your_dataset \
  --staged \
  --gpt-role stats \
  --claude-role pipeline \
  --gemini-role biology \
  --output results/system1_output
```

**Run Baseline: No-Agent (System 7):**
```bash
node bin/geneexpert.js analyze data/your_dataset \
  --staged \
  --force-automation \
  --output results/no_agent_output
```

**Run Baseline: Single-Agent Claude (System 8):**
```bash
node bin/geneexpert.js analyze data/your_dataset \
  --staged \
  --single-agent claude \
  --output results/single_claude_output
```

### Total Experiments

- **Bulk RNA-seq:** 10 systems × 7 datasets = 70 analyses
- **scRNA-seq:** 10 systems × 7 datasets = 70 analyses
- **Total:** 140 analyses

---

## System Architecture

### Dual-Pipeline Design

**Bulk RNA-seq (4 stages, 4 agent checkpoints):**
```
Stage 1: FASTQ Validation → Agent Checkpoint
Stage 2: Alignment + QC → Agent Checkpoint (sample filtering)
Stage 3: Quantification + PCA → Agent Checkpoint (batch detection, outlier removal)
Stage 4: Differential Expression → Agent Checkpoint (result approval)
```

**Single-cell RNA-seq (5 stages, 4 agent checkpoints):**
```
Stage 1: Load 10x Data + QC → Auto-proceed (no agents)
Stage 2: QC Filtering → Agent Checkpoint (threshold selection)
Stage 3: Normalization + HVG → Auto-proceed
Stage 3A: Cell Cycle Scoring → Agent Checkpoint (decide regression)
Stage 3B: Cell Cycle Regression → Conditional execution based on 3A decision
Stage 4: PCA → Agent Checkpoint (PC selection)
Stage 5: Clustering + Markers → Agent Checkpoint (cluster validation)
```

### Multi-Agent Consensus

**3 Specialized Agents:**
- **Stats Agent (GPT-5.2):** Statistical validation, threshold selection
- **Pipeline Agent (Claude Sonnet 4.5):** Technical feasibility, best practices
- **Biology Agent (Gemini 2.0 Flash):** Biological interpretation, domain knowledge

**Two Collaboration Modes:**
1. **Parallel (Default):** All 3 agents respond independently → Majority vote
2. **Sequential Chain:** GPT → Gemini → Claude (each sees previous responses)

**Consensus Rules:**
- Majority (2/3): Proceed automatically
- Unanimous (3/3): Required for sample removal
- No consensus: Escalate to user

### Auto-Resolution System (3-Tier Escalation)

- **Minor disagreement (score < 0.3):** Auto-resolve using median/average
- **Moderate disagreement (0.3 ≤ score < 0.6):** Use highest confidence agent
- **Major disagreement (score ≥ 0.6):** Escalate to user decision

---

## Key Features

1. **Staged Validation:** Issues caught early (Stage 2) rather than at the end (Stage 4)
2. **Progressive Refinement:** Remove outliers and batch effects before final analysis
3. **Conditional Execution:** Stage 3B (cell cycle regression) branches based on Stage 3A decision
4. **User Input Tracking:** All user decisions logged with `user_input_required` flag for research evaluation
5. **Role Swapping:** Test all 6 permutations of agent role assignments
6. **Cost Tracking:** Per-decision API cost logging
7. **Comprehensive Logging:** JSON/JSONL logs with full agent conversations and decisions

---

## Evaluation Metrics

- **Decision Accuracy:** (correct decisions) / total → Target >90%
- **Error Reduction:** (baseline errors - multi-agent errors) / baseline → Target >40%
- **Success Rate:** (successful analyses) / total → Target >95%
- **Cost Efficiency:** Total cost / successes → Target <$0.10 per analysis
- **Inter-Agent Agreement:** Cohen's kappa → Target >0.7
- **User Input Frequency:** % decisions requiring human intervention
- **Role Performance:** Which model excels at which role (stats/pipeline/biology)

**Run Evaluation:**
```bash
# Calculate metrics for bulk RNA-seq
node bin/evaluate_bulk.js

# Calculate metrics for scRNA-seq
node bin/evaluate_scrna.js

# Statistical significance testing
python bin/run_ttest.py
```

---

## Project Structure

```
├── bin/
│   ├── geneexpert.js                    # Bulk RNA-seq CLI
│   ├── scrna_geneexpert.js              # scRNA-seq CLI
│   ├── evaluate_bulk.js                 # Bulk evaluation
│   ├── evaluate_scrna.js                # scRNA evaluation
│   └── [93 experiment/evaluation scripts]
│
├── src/
│   ├── executor/
│   │   └── staged_executor.js           # Bulk 4-stage orchestration
│   ├── scrna_executor/
│   │   └── scrna_executor.js            # scRNA 5-stage orchestration
│   ├── stages/                          # Bulk RNA-seq stages
│   ├── scrna_stages/                    # scRNA-seq stages
│   ├── config/
│   │   ├── stage_prompts.js             # Bulk prompts
│   │   └── scrna_stage_prompts.js       # scRNA prompts
│   ├── coordinator/
│   │   ├── orchestrator.js              # Multi-agent coordination
│   │   └── consensus.js                 # Voting logic
│   └── utils/
│       ├── llm_clients.js               # API wrappers
│       ├── logger.js                    # Comprehensive logging
│       ├── auto_resolver.js             # 3-tier escalation
│       └── metrics.js                   # Evaluation metrics
│
├── bio_informatics/
│   ├── scripts/                         # R/bash analysis scripts (30 files)
│   └── reference_data/                  # Genome annotations
│
└── experiments/
    ├── results/                         # Bulk RNA-seq outputs
    ├── scrna_results/                   # scRNA-seq outputs
    ├── bulk_rna_ground_truth.json       # Ground truth (7 datasets)
    └── scrna_ground_truth.json          # Ground truth (7 datasets)
```

---

## Cost Estimate

### Per Analysis:
- **Bulk RNA-seq:** 3 agents × 4 checkpoints = 12 API calls → ~$0.12
- **scRNA-seq:** 3 agents × 4 checkpoints = 12 API calls → ~$0.12
- **Script Execution:** Free (runs locally)

**Total: ~$0.12 per complete analysis with full multi-agent validation**

---

## Citation

If you use this system in your research, please cite:

```bibtex
@software{geneexpert2026,
  title={GeneExpert: Poly-Foundational LLM Multi-Agent System for RNA-seq Analysis},
  author={[Your Name]},
  year={2026},
  url={https://github.com/[your-username]/geneexpert-mcp}
}
```

---

## License

MIT

---

## Links

- **GitHub:** https://github.com/[your-username]/geneexpert-mcp
- **Status:** Production-ready (140 experiment pipeline complete)
- **Experiments:** 10 systems × 14 datasets = 140 total analyses
