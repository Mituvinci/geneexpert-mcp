# Supplementary Material: Dataset Descriptions

**GeneExpert: Multi-Agent LLM System for RNA-seq Analysis**

This document provides detailed descriptions, GEO accessions, and sample mappings for all 14 RNA-seq datasets used in the GeneExpert evaluation.

For a concise overview, see the main [README.md](../README.md).

---

### Data Preparation (Bulk RNA-seq)

Each dataset in `data/download_new_bulk_RNA_Data/` has its own download and rename scripts. Run them in order from anywhere in the repo:

```bash
# Step 1: Download FASTQs from NCBI SRA
bash data/download_new_bulk_RNA_Data/1_GSE52778_pe_clean/1_GSE52778_download.sh

# Step 2: Rename to pipeline format (label-only, no GSM ID)
bash data/download_new_bulk_RNA_Data/1_GSE52778_pe_clean/1_GSE52778_rename.sh
```

After rename, the dataset folder contains:

```
1_GSE52778_pe_clean/
├── gsm_map_GSE52778.tsv
├── 1_GSE52778_download.sh
├── 1_GSE52778_rename.sh
│
├── N61311_untreated_R1_001.fastq.gz      # Control (4 samples x 2 reads = 8 files)
├── N61311_untreated_R2_001.fastq.gz
├── N052611_untreated_R1_001.fastq.gz
├── N052611_untreated_R2_001.fastq.gz
├── N080611_untreated_R1_001.fastq.gz
├── N080611_untreated_R2_001.fastq.gz
├── N061011_untreated_R1_001.fastq.gz
├── N061011_untreated_R2_001.fastq.gz
│
├── N61311_Dex_R1_001.fastq.gz            # Treatment (4 samples x 2 reads = 8 files)
├── N61311_Dex_R2_001.fastq.gz
├── N052611_Dex_R1_001.fastq.gz
├── N052611_Dex_R2_001.fastq.gz
├── N080611_Dex_R1_001.fastq.gz
├── N080611_Dex_R2_001.fastq.gz
├── N061011_Dex_R1_001.fastq.gz
└── N061011_Dex_R2_001.fastq.gz
```

The pipeline identifies groups by keyword: `--control-keyword "untreated"` matches control files, `--treatment-keyword "Dex"` matches treatment files. All datasets follow the same two-step download-then-rename pattern:

| Dataset folder | Type | Organism | Control keyword | Treatment keyword |
|----------------|------|----------|-----------------|-------------------|
| `1_GSE52778_pe_clean` | Paired-end | human | `untreated` | `Dex` |
| `2_GSE114845_se_clean` | Single-end | mouse | `Control` | `SleepDeprived` |
| `3_GSE113754_pe_clean` | Paired-end | mouse | `Control` | `SleepDeprived` |
| `4_GSE141496_batch_effect` | Paired-end | human | `Control` | `TAB182_KD` |
| `5_GSE47774_batch_effect` | Paired-end | human | `AC` | `BC` |
| `6_GSE193658_Lab_data` | Paired-end | human | `MONO` | `TSW` |
| `E.coli_GSE48151` | Single-end | — | — | — |

Note: `E.coli_GSE48151` is the contamination source used to generate the `7_GSE114845_CONTAM70` dataset — it is not run as a standalone analysis.


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
- **URL:** https://www.10xgenomics.com/datasets/pbm-cs-from-a-healthy-donor-whole-transcriptome-analysis-3-1-standard-4-0-0

---

#### 4. Mouse Brain Cells (Healthy) - 10k Brain Cells from E18 Mouse

- **Organism:** Mouse (embryonic day 18)
- **Platform:** 10x Chromium v3 chemistry (11,843 cells)
- **Category:** Clean single-cell
- **Cell cycle labels:** 5,524 cells (3,830 G1, 1,116 S, 578 G2M) consensus labeled
- **URL:** https://www.10xgenomics.com/datasets/10-k-brain-cells-from-an-e-18-mouse-v-3-chemistry-3-standard-3-0-0

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
