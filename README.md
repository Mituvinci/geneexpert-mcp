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

**Bulk RNA-seq:**
```bash
node bin/geneexpert.js analyze data/your_dataset \
  --staged --organism mouse \
  --control-keyword "control" --treatment-keyword "treatment" \
  --output results/output_folder
```

**Single-cell RNA-seq:**
```bash
node bin/scrna_geneexpert.js analyze data/scRNA_data/your_dataset \
  --output results/scRNA_output --organism human
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

### Ground Truth & Evaluation

**Ground truth decisions** for all datasets are in:
- **Bulk RNA-seq:** [`ground_truth_supplementary/bulk_rna_ground_truth.json`](ground_truth_supplementary/bulk_rna_ground_truth.json)
- **scRNA-seq:** [`ground_truth_supplementary/scrna_ground_truth.json`](ground_truth_supplementary/scrna_ground_truth.json)

Each ground truth file contains:
- Correct decisions for each stage checkpoint
- Rationale for expert-curated decisions
- Tissue-specific expectations and thresholds

**Run evaluation:**
```bash
# Bulk RNA-seq evaluation
node bin/evaluate_bulk.js

# scRNA-seq evaluation
node bin/evaluate_scrna.js

# Statistical significance testing
python bin/run_ttest.py
```

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

**Run all systems:**
```bash
# Bulk RNA-seq
bash bin/run_bulk_rna.sh data/your_dataset

# scRNA-seq
bash bin/run_sc_rna.sh data/scRNA_data/your_dataset
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
- **Biology Agent (Gemini 2.0 Flash):** Biological interpretation

**Consensus Mechanism:**
- Majority (2/3): Auto-proceed
- No consensus: Auto-resolution (3-tier escalation) or user input
- All decisions logged with `decision_id` for evaluation

---

## Evaluation Metrics

| Metric | Target | Formula |
|--------|--------|---------|
| **Decision Accuracy** | >90% | Correct decisions / Total |
| **Error Reduction** | >40% | (Baseline errors - Multi-agent errors) / Baseline errors |
| **Success Rate** | >95% | Successful analyses / Total |
| **Cost Efficiency** | <$0.10 | Total API cost / Successful analyses |
| **Inter-Agent Agreement** | >0.7 | Cohen's κ |

**Per-analysis cost:** ~$0.12 (3 agents × 4 checkpoints × $0.01/call)

---

## Project Structure

```
├── bin/                                 # CLI tools & experiment scripts
│   ├── geneexpert.js                    # Bulk RNA-seq entry point
│   ├── scrna_geneexpert.js              # scRNA-seq entry point
│   ├── evaluate_bulk.js                 # Bulk evaluation
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
    ├── results/                         # Bulk RNA-seq outputs
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
  title={GeneExpert: Multi-Agent LLM System for RNA-seq Analysis},
  author={[Your Name]},
  year={2026},
  url={https://github.com/yourusername/geneexpert}
}
```

---

## License

MIT

---

## Links

- **GitHub:** [Repository URL]
- **Supplementary Materials:** [`ground_truth_supplementary/`](ground_truth_supplementary/)
