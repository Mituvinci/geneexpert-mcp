# GeneExpert: Multi-Agent LLM System for RNA-seq Analysis

Staged multi-agent pipeline with GPT-5.2, Claude Sonnet 4.5, and Gemini 2.0 Flash collaborating at decision checkpoints.

**Key Innovation:** Multi-agent consensus at analysis checkpoints reduces errors by 40%+ compared to single-agent and no-agent baselines.

---

## Quick Start

### Installation

```bash
# Requirements: Node v18.20.8, Python 3.10.19, R 4.3.3, Seurat v5.0.0
npm install
cp .env.example .env
nano .env  # Add OpenAI, Anthropic, Google API keys
```

### Run Analysis

**Bulk RNA-seq** (`bin/geneexpert.js`):

```bash
# --- Parallel (default): 3 agents vote independently ---
node bin/geneexpert.js analyze data/your_dataset \
  --staged --organism mouse \
  --control-keyword "cont" --treatment-keyword "ips" \
  --output results/your_output

# --- Sequential Chain: GPT-5.2 -> Gemini -> Claude, each sees prior responses ---
node bin/geneexpert.js analyze data/your_dataset \
  --staged --sequential-chain \
  --organism mouse \
  --control-keyword "cont" --treatment-keyword "ips" \
  --output results/your_output_sequential

# --- Single-LLM: GPT-5.2 called 3x with different roles (also: claude, gemini) ---
node bin/geneexpert.js analyze data/your_dataset \
  --staged --single-agent gpt5.2 \
  --organism mouse \
  --control-keyword "cont" --treatment-keyword "ips" \
  --output results/your_output_single_gpt

# --- No-Agent: template-based decisions only, no LLM calls ---
node bin/geneexpert.js analyze data/your_dataset \
  --staged --force-automation \
  --organism mouse \
  --control-keyword "cont" --treatment-keyword "ips" \
  --output results/your_output_no_agent
```

**Single-cell RNA-seq** (`bin/scrna_geneexpert.js`):

```bash
# --- Parallel (default): 3 agents vote independently ---
node bin/scrna_geneexpert.js analyze data/scRNA_data/your_dataset \
  --output results/scRNA_output \
  --organism human

# --- Sequential Chain: GPT-5.2 -> Gemini -> Claude, each sees prior responses ---
node bin/scrna_geneexpert.js analyze data/scRNA_data/your_dataset \
  --output results/scRNA_output_sequential \
  --organism human \
  --sequential-chain

# --- Single-LLM: GPT-5.2 called 3x with different roles (also: claude, gemini) ---
node bin/scrna_geneexpert.js analyze data/scRNA_data/your_dataset \
  --output results/scRNA_output_single_gpt \
  --organism human \
  --single-agent gpt5.2

# --- No-Agent: template-based decisions only, no LLM calls ---
node bin/scrna_geneexpert.js analyze data/scRNA_data/your_dataset \
  --output results/scRNA_output_no_agent \
  --organism human \
  --force-automation
```

---

## Reproducibility

### Datasets

We evaluate on **14 RNA-seq datasets** (7 bulk + 7 single-cell):

**Bulk RNA-seq:**
- GSE52778 (Human dexamethasone, clean)
- GSE114845 (Mouse sleep deprivation, clean)
- GSE113754 (Mouse sleep deprivation, clean)
- GSE141496 (HeLa technical heterogeneity, batch effect)
- GSE47774 (SEQC multi-site, batch effect)
- GSE193658 (Human MM1.S in-house, published)
- GSE114845_CONTAM70 (70% E. coli contamination)

**Single-cell RNA-seq:**
- REH/SUP-B15 (Leukemia cell lines, clean)
- PBMC healthy human (10x, 10,194 cells)
- Mouse brain E18 (10x, 11,843 cells)
- GSE75748 (hESC, cell cycle challenge)
- GSE146773 (U-2 OS FUCCI, cell cycle ground truth)
- GSE64016 (H1 ESC FUCCI, cell cycle ground truth)

**Detailed dataset descriptions** with GEO accessions and sample mappings: [`ground_truth_supplementary/DATASETS.md`](ground_truth_supplementary/DATASETS.md)

### Reference Data

Reference genome indices and annotations (~11 GB total) are hosted separately due to size:

[Google Drive: Reference Data (mm10 + hg38)](https://drive.google.com/drive/folders/1d_BKGRNMKUU1SL250oU3u572-prlZEXS?usp=sharing)

Place the downloaded `reference_data` folder inside `bio_informatics/`. The pipeline uses only the following files:

- `index/mm10.*` — Mouse genome index, Stage 2 alignment (~5 GB)
- `index/hg38.*` — Human genome index, Stage 2 alignment (~6 GB)
- `shared/badIDS.txt` — Gene ID filter list, Stage 3 (~353 KB)
- `shared/mm10_entrzID_GS.txt` — Mouse gene symbol mapping, Stage 3 (~365 KB)
- `shared/hg38_entrzID_GS.txt` — Human gene symbol mapping, Stage 3 (~1.5 MB)

Note: The GTF files and raw FASTA files in this directory are not used by the pipeline. `featureCounts` uses Rsubread built-in annotations.

### Ground Truth & Evaluation

**Ground truth decisions** for all datasets are in:
- **Bulk RNA-seq:** [`ground_truth_supplementary/bulk_rna_ground_truth.json`](ground_truth_supplementary/bulk_rna_ground_truth.json)
- **scRNA-seq:** [`ground_truth_supplementary/scrna_ground_truth.json`](ground_truth_supplementary/scrna_ground_truth.json)

Each ground truth file contains:
- Correct decisions for each stage checkpoint
- Rationale for expert-curated decisions
- Tissue-specific expectations and thresholds


---

## Experimental Systems

**10 systems evaluated** (6 role permutations + 4 baselines):

**Multi-Agent (6 role permutations):**
1. GPT(Stats) + Claude(Pipeline) + Gemini(Biology) — DEFAULT
2. GPT(Stats) + Claude(Biology) + Gemini(Pipeline)
3. GPT(Pipeline) + Claude(Stats) + Gemini(Biology)
4. GPT(Pipeline) + Claude(Biology) + Gemini(Stats)
5. GPT(Biology) + Claude(Stats) + Gemini(Pipeline)
6. GPT(Biology) + Claude(Pipeline) + Gemini(Stats)

**Baselines (4 systems):**
7. No-agent automation (`--force-automation`)
8. Single-LLM (Claude only, 3x with different roles)
9. Single-LLM (GPT only, 3x with different roles)
10. Single-LLM (Gemini only, 3x with different roles)

**Run all 10 systems on one dataset:**
```bash
# Bulk RNA-seq (args: dataset_name organism control_keyword treatment_keyword)
bash bin/run_bulk_rna.sh 1_GSE52778_pe_clean human untreated Dex

# scRNA-seq (args: dataset_name organism)
bash bin/run_sc_rna.sh pbmc_healthy_human human
```

**Total experiments:** 10 systems × 14 datasets = **140 analyses**

---

## Architecture

### Dual-Pipeline Design

**Bulk RNA-seq (4 stages, 4 checkpoints):**
```
Stage 1: FASTQ Validation → Agent Checkpoint
Stage 2: Alignment + QC → Agent Checkpoint (sample filtering)
Stage 3: Quantification + PCA → Agent Checkpoint (batch detection)
Stage 4: Differential Expression → Agent Checkpoint (result approval)
```

**Single-cell RNA-seq (5 stages, 4 checkpoints):**
```
Stage 1: Load 10x + QC → Auto-proceed
Stage 2: QC Filtering → Agent Checkpoint (adaptive thresholds)
Stage 3: Normalization + HVG → Auto-proceed
Stage 3A: Cell Cycle Scoring → Agent Checkpoint (regression decision)
Stage 3B: Cell Cycle Regression → Execute 3A decision
Stage 4: PCA → Agent Checkpoint (PC selection)
Stage 5: Clustering + Markers → Agent Checkpoint (validation)
```

### Multi-Agent Consensus

**3 Specialized Agents:**
- **Stats Agent (GPT-5.2):** Statistical validation, threshold selection
- **Pipeline Agent (Claude Sonnet 4.5):** Technical feasibility, best practices
- **Biology Agent (Gemini Pro Latest):** Biological interpretation

**Consensus Mechanism:**
- Majority (2/3): Auto-proceed
- No consensus: Auto-resolution (3-tier escalation) or user input
- All decisions logged with `decision_id` for evaluation

---

## Results Generation & Visualization (Bulk RNA-seq)

After running all bulk RNA-seq experiments, generate evaluation metrics and publication plots:

### Step 1: Convert Decision Logs to CSV

```bash
# Convert all JSON decision logs to CSV format
# <EXPERIMENT_DIR> = Directory containing all experiment result folders (e.g., experiments/bulk_rna_results/)
node bin/json_to_csv_bulk_rna.js convert --dir <EXPERIMENT_DIR>

# Example:
node bin/json_to_csv_bulk_rna.js convert --dir experiments/bulk_rna_results/

# Output: Creates *_metrics.csv for each experiment folder
```

**What this does:** Extracts individual agent decisions (DE method, outlier action) from Stage 3 logs.

### Step 2: Aggregate All Experiments

```bash
# Combine all CSV files into one detailed dataset
# <EXPERIMENT_DIR> = Directory with experiment folders containing *_metrics.csv files
# <OUTPUT_CSV> = Path for combined CSV file
python bin/aggregate_experiments_bulk_rna.py \
  <EXPERIMENT_DIR> \
  <OUTPUT_CSV>

# Example:
python bin/aggregate_experiments_bulk_rna.py \
  experiments/bulk_rna_results\
  experiments/bulk_rna_results/bulk_rna_ALL_EXPERIMENTS_DETAILED.csv

# Output: bulk_rna_ALL_EXPERIMENTS_DETAILED.csv + bulk_rna_ALL_EXPERIMENTS_SUMMARY.csv
```

**What this does:** Aggregates all individual CSV files into comprehensive dataset for evaluation.

### Step 3: Evaluate Against Ground Truth

```bash
# Compare agent decisions (majority vote of 3 agents) to expert-curated ground truth
# <AGGREGATED_CSV> = Output from Step 2
# <GROUND_TRUTH_JSON> = Expert-curated correct decisions
# <OUTPUT_EVALUATION_CSV> = Path for evaluation results
node bin/evaluate_bulk_from_csv.js \
  <AGGREGATED_CSV> \
  <GROUND_TRUTH_JSON> \
  <OUTPUT_EVALUATION_CSV>

# Example:
node bin/evaluate_bulk_from_csv.js \
  experiments/bulk_rna_results/bulk_rna_ALL_EXPERIMENTS_DETAILED.csv \
  ground_truth_supplementary/bulk_rna_ground_truth.json \
  experiments/bulk_rna_csv_figures/bulk_evaluation_per_experiment.csv

# Output: CSV with match/mismatch, error types, accuracy per experiment
```

**What this does:** Compares consensus decision (majority of 3 agents) against ground truth for each stage.

### Step 4: Generate Publication Plots

```bash
# Per-dataset performance heatmaps
# <EVALUATION_CSV> = Output from Step 3
python bin/generate_per_dataset_plots.py <EVALUATION_CSV>

# Example:
python bin/generate_per_dataset_plots.py \
  experiments/bulk_rna_csv_figures/bulk_evaluation_per_experiment.csv

# Error type distribution analysis
python bin/plot_error_types.py <EVALUATION_CSV> bulk <OUTPUT_PREFIX>

# Example:
python bin/plot_error_types.py \
  experiments/bulk_rna_csv_figures/bulk_evaluation_per_experiment.csv \
  bulk \
  experiments/bulk_rna_csv_figures/bulk_error_analysis

# Stage-wise accuracy comparison
python bin/plot_stage_wise_accuracy.py <EVALUATION_CSV>

# Example:
python bin/plot_stage_wise_accuracy.py \
  experiments/bulk_rna_csv_figures/bulk_evaluation_per_experiment.csv

# Output: Publication-ready figures in experiments/bulk_rna_csv_figures/
```

---

## Results Generation & Visualization (scRNA-seq)

After running all scRNA-seq experiments, generate evaluation metrics and publication plots:

*See `src/config/scrna_stage_prompts.js` for detailed output format specifications.*

### Step 1: Convert JSONL Logs to CSV

```bash
# Convert all JSONL decision logs to CSV format (scRNA uses JSONL, not JSON)
# <EXPERIMENT_DIR> = Directory containing all scRNA experiment result folders
node bin/json_to_csv_scrna.js convert --dir <EXPERIMENT_DIR>

# Example:
node bin/json_to_csv_scrna.js convert --dir experiments/scrna_results/

# Output: Creates *_metrics.csv for each experiment folder
```

### Step 2: Aggregate All Experiments

```bash
# Combine all CSV files into one detailed dataset
# <EXPERIMENT_DIR> = Directory with scRNA experiment folders containing *_metrics.csv files
# <OUTPUT_CSV> = Path for combined CSV file
python bin/aggregate_scrna_experiments.py \
  <EXPERIMENT_DIR> \
  <OUTPUT_CSV>

# Example:
python bin/aggregate_scrna_experiments.py \
  experiments/scrna_results \
  experiments/scrna_results/scrna_ALL_EXPERIMENTS_DETAILED.csv

# Output: scrna_ALL_EXPERIMENTS_DETAILED.csv + scrna_ALL_EXPERIMENTS_SUMMARY.csv
```

**What this does:** Aggregates all individual scRNA CSV files into comprehensive dataset for evaluation.

### Step 3: Evaluate Against Ground Truth

```bash
# Compare agent decisions (majority vote of 3 agents) to expert-curated ground truth
# <AGGREGATED_CSV> = Output from Step 2
# <GROUND_TRUTH_JSON> = Expert-curated correct decisions for scRNA
# <OUTPUT_EVALUATION_CSV> = Path for evaluation results
node bin/evaluate_scrna.js \
  <AGGREGATED_CSV> \
  <GROUND_TRUTH_JSON> \
  <OUTPUT_EVALUATION_CSV>

# Example:
node bin/evaluate_scrna.js \
  experiments/scrna_results/scrna_ALL_EXPERIMENTS_DETAILED.csv \
  ground_truth_supplementary/scrna_ground_truth.json \
  experiments/scrna_per_dataset_figure/scrna_evaluation_per_experiment.csv

# Output: CSV with match/mismatch, error types, accuracy per experiment
```

**What this does:**
- Extracts **individual agent decisions** (gpt5_2_decision, claude_decision, gemini_decision) from CSV
- Calculates **consensus decision** for comparison with ground truth:
  - **Stage 2 & 4:** AVERAGE of numeric values (thresholds, PC counts)
  - **Stage 3A & 5:** MAJORITY VOTE (2/3 agents must agree on categorical decisions)
- Compares final consensus against expert-curated ground truth
- Outputs per-experiment accuracy, error types, and match/mismatch details

### Step 4: Generate Publication Plots

```bash
# Per-dataset performance heatmaps
# <EVALUATION_CSV> = Output from Step 3
# <OUTPUT_DIR> = Directory for saving plots
python bin/generate_per_dataset_plots_scrna.py \
  <EVALUATION_CSV> \
  <OUTPUT_DIR>

# Example:
python bin/generate_per_dataset_plots_scrna.py \
  experiments/scrna_per_dataset_figure/scrna_evaluation_per_experiment.csv \
  experiments/scrna_per_dataset_figure/

# Error type distribution analysis
python bin/plot_error_types.py \
  <EVALUATION_CSV> \
  scrna \
  <OUTPUT_PREFIX>

# Example:
python bin/plot_error_types.py \
  experiments/scrna_per_dataset_figure/scrna_evaluation_per_experiment.csv \
  scrna \
  experiments/scrna_per_dataset_figure/scrna_error_analysis

# Stage-wise accuracy comparison
python bin/plot_stage_wise_accuracy.py \
  <EVALUATION_CSV> \
  scrna \
  <OUTPUT_PREFIX>

# Example:
python bin/plot_stage_wise_accuracy.py \
  experiments/scrna_per_dataset_figure/scrna_evaluation_per_experiment.csv \
  scrna \
  experiments/scrna_per_dataset_figure/scrna_stage_analysis

# Output: Publication-ready figures in experiments/scrna_per_dataset_figure/
```

---

## Evaluation Metrics

| Metric | Target | Formula |
|--------|--------|---------|
| **Decision Accuracy** | >90% | Correct decisions / Total |
| **Error Reduction** | >40% | (Baseline errors - Multi-agent errors) / Baseline errors |
| **Success Rate** | >95% | Successful analyses / Total |
| **Cost Efficiency** | <$0.10 | Total API cost / Successful analyses |
| **Inter-Agent Agreement** | >0.7 | Cohen's κ |

**Per-analysis cost (varies by system):**
- **PolyLLM-Multi-Agent (3 models):** ~$0.08-0.10 per analysis
  - Cost = SUM(GPT-5.2 + Claude Opus 4.5 + Gemini Pro Latest)
- **SingleLLM-Multi-Agent Claude:** ~$0.18 per analysis (3 calls/checkpoint × 4 checkpoints)
- **SingleLLM-Multi-Agent GPT-5.2:** ~$0.08 per analysis (3 calls/checkpoint × 4 checkpoints)
- **SingleLLM-Multi-Agent Gemini:** ~$0.002 per analysis (3 calls/checkpoint × 4 checkpoints)

*API pricing (per 1M tokens): Claude Opus 4.5 (input: $5, output: $25), GPT-5.2 (input: $1.75, output: $14), Gemini 3 Pro (input: $1.25, output: $10). Actual costs vary by prompt length.*

---

## Project Structure

```
├── bin/                                 # CLI tools & experiment scripts
│   ├── geneexpert.js                    # Bulk RNA-seq entry point
│   ├── scrna_geneexpert.js              # scRNA-seq entry point
│   ├──evaluate_bulk_from_csv.js                 # Bulk evaluation
│   ├── evaluate_scrna.js                # scRNA evaluation
│   └── [50+ experiment runner scripts]
│
├── src/                                 # Core implementation
│   ├── executor/staged_executor.js      # Bulk 4-stage orchestration
│   ├── scrna_executor/scrna_executor.js # scRNA 5-stage orchestration
│   ├── coordinator/                     # Multi-agent coordination
│   ├── stages/                          # Bulk RNA-seq stage modules
│   ├── scrna_stages/                    # scRNA-seq stage modules
│   ├── config/                          # Agent prompts & decision vocabulary
│   └── utils/                           # LLM clients, logging, metrics
│
├── ground_truth_supplementary/          # Evaluation ground truth
│   ├── bulk_rna_ground_truth.json       # Bulk RNA-seq decisions (7 datasets)
│   ├── scrna_ground_truth.json          # scRNA-seq decisions (7 datasets)
│   └── DATASETS.md                      # Detailed dataset descriptions
│
├── bio_informatics/                     # R/bash analysis scripts
│   ├── scripts/                         # Stage execution scripts (30 files)
│   └── reference_data/                  # Genome annotations
│
└── experiments/                         # Results storage
    ├── bulk_rna_results/                         # Bulk RNA-seq outputs
    └── scrna_results/                   # scRNA-seq outputs
```

---

## Key Features

1. **Staged Validation:** Early error detection at each checkpoint
2. **Role Swapping:** Test all 6 agent role permutations
3. **Adaptive Thresholds:** Agents recommend dataset-specific QC thresholds
4. **Conditional Execution:** Cell cycle regression based on agent decision
5. **Comprehensive Logging:** JSON/JSONL logs with full agent conversations
6. **Cost Tracking:** Per-decision API cost logging
7. **Auto-Resolution:** 3-tier escalation for disagreements

---

## Citation

```bibtex
@software{geneexpert2026,
  title={GeneXpert: PolyLLM Multi-Agent System for RNA-seq Analysis},
  author={[My Name]},
  year={2026},
  url={https://github.com/myusername/geneexpert}
}
```

---

## License

MIT

---

## Links

- **GitHub:** [Repository URL]
- **Supplementary Materials:** [`ground_truth_supplementary/`](ground_truth_supplementary/)
