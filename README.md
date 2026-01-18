# GeneExpert: Multi-Agent RNA-seq Analysis System

**AUTOMATION + ADAPTATION Hybrid System** where GPT-4, Claude (MCP-enabled), and Gemini collaborate to analyze RNA-seq data intelligently.

**Key Innovation:** Agents DECIDE the approach dynamically - fast template-based automation for clean data, intelligent MCP-powered adaptation for edge cases.

**Research Goal:** ICML 2026 paper - Multi-agent collaboration reduces errors 40%+ through adaptive decision-making

---

## Quick Start

### 1. Requirement

```bash
Node: v18.20.8
npm: 10.8.2
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

### 3. Run Analysis

```bash
node bin/geneexpert.js analyze data/DA0036 \
  --organism mouse \
  --comparison "stroke_vs_control" \
  --control-keyword "cont" \
  --treatment-keyword "ips" \
  --output results/my_analysis
```

**What happens:**
1.  System detects data type (FASTQ/BAM/counts), paired-end status, sample groups
2.  **3 agents analyze your data** (GPT-4, Claude, Gemini)
3.  **Agents vote: AUTOMATION or ADAPTATION?**
   - Small n, batch effects, outliers → ADAPTATION (custom script)
   - Clean data, standard design → AUTOMATION (fast template)
4.  **Beautiful plan displayed** - you confirm with Y/N
5.  **Script executes** with real-time output streaming
6.  Results saved with full logs

---

##  What Makes This Novel

### AUTOMATION vs ADAPTATION: Agents Decide Dynamically

**Traditional pipelines:**
```
Run fixed steps → Hope it works → Manual debugging if fails
```

**GeneExpert Multi-Agent System:**
```
Step 1: Detect data characteristics
  ↓
Step 2: Multi-agent analysis
  [GPT-4]: "n=2 per group is risky for edgeR"
  [Claude]: "Low replicates need custom approach"
  [Gemini]: "Statistical power too low for standard pipeline"
  ↓
Step 3: Consensus voting
  Vote: 3/3 ADAPTATION (100% confidence)
  ↓
Step 4: Decision implementation

  IF AUTOMATION (template-based):
    - Generate bash script from template
    - Standard 10-step RNA-seq pipeline
    - Fast, $0 cost, runs in minutes

  IF ADAPTATION (MCP-powered):
    - MCP Claude agent READS actual R/bash scripts
    - Discovers parameters, understands workflow
    - WRITES custom script addressing concerns
    - Intelligent, ~$0.50 cost, solves edge cases
  ↓
Step 5: User confirms → Execute → Results!
```

**This is TRUE intelligence**: Agents see your data, debate the approach, and write solutions - not just templates!

---

##  Architecture Overview

### Two Execution Paths

#### Path 1: AUTOMATION (Fast & Free)
```
When: Clean data, standard design (n≥3), no edge cases
How: Template-based bash script
Cost: $0 (no API calls during execution)
Speed: Fast (no agent overhead)

Pipeline: FASTQ Validation → FastQC → Alignment → featureCounts → filterIDS
          → RPKM → entrz → QC Plots → edgeR → merge_results → Excel
```

#### Path 2: ADAPTATION (Intelligent - Agent Writes Code!)
```
When: Small n, batch effects, outliers, edge cases
How: MCP-enabled Claude READS existing scripts & WRITES custom solution
Cost: ~$0.05 (MCP tool calls)
Speed: Slower (agent analysis + custom script generation)

Process:
  Step 1: Claude MCP agent READS existing R/bash scripts
          - Calls read_file("/bio_informatics/scripts/featurecounts.R")
          - Discovers parameters: need annotation, output, control, treatment, *.bam
          - Understands workflow requirements

  Step 2: Claude ANALYZES agent concerns from debate
          - GPT-4: "Small n=2, use exact test"
          - Gemini: "Need extra QC for low replicates"

  Step 3: Claude WRITES custom bash script from scratch
          - Not modifying existing files!
          - Creates NEW orchestration script: geneexpert_v1.sh
          - Adjusts parameter order, adds validations, modifies logic

  Step 4: System EXECUTES Claude's generated script
          - Bash executes the custom script
          - Claude monitors output via stdout/stderr
          - If fails → Feedback loop (Claude writes v2)

IMPORTANT: Claude doesn't modify your R scripts - it writes NEW bash
orchestration scripts that CALL your existing R scripts intelligently!
```

### Multi-Agent Decision System

**3 Specialized Agents:**
-  **GPT-4 (Stats Agent)**: Statistical validation, threshold selection
  - **No tool access** - Pure reasoning only
  - Validates thresholds, sample sizes, statistical methods

- 🔧 **Claude Sonnet 4.5 (Pipeline Agent)**: Technical execution, MCP tool calling
  - **HAS MCP tool access** - Can read files and write code
  - **In ADAPTATION mode:** Reads existing R/bash scripts via MCP tools
  - **Writes custom bash scripts** that orchestrate pipeline with modifications
  - **Executes generated scripts** and monitors output
  - **Can modify workflow logic**: skip steps, add validations, adjust parameters

-  **Gemini Pro (Biology Agent)**: Biological interpretation, pathway analysis
  - **No tool access** - Pure reasoning only
  - Validates biological assumptions and QC criteria

**How Agents Communicate:**
1. **Coordinator** sends same question to all 3 agents **in parallel**
2. Each agent analyzes independently (no cross-talk)
3. **Coordinator** collects all responses
4. **Consensus mechanism** applies voting rules
5. **Synthesized decision** returned to user

**Consensus Voting Rules:**
- **Standard execution**: Pipeline Agent only (no vote needed)
- **Threshold selection**: Majority vote (2/3 required)
- **Sample removal**: Unanimous (3/3 required)
- **Disagreement**: Escalate to user for final decision

---

##  Which Agent Writes Code?

**ONLY Claude Sonnet 4.5 (Pipeline Agent) can write and execute code!**

### What Each Agent Does:

| Agent | Can Read Files? | Can Write Code? | Can Execute Code? | Can View Images? | Role |
|-------|----------------|-----------------|-------------------|------------------|------|
| **GPT-4 (Stats)** |  No |  No |  No |  YES | Pure reasoning - validates statistics |
| **Claude (Pipeline)** |  YES (MCP) |  YES (writes bash) |  YES (executes scripts) |  YES | Code generation & execution |
| **Gemini (Biology)** |  No |  No |  No |  YES | Pure reasoning - validates biology |

**NEW! All 3 agents have VISION capability** - they can view PCA/MDS plots and vote on outliers/batch effects!

### What Claude Writes in ADAPTATION Mode:

**Claude WRITES:** New bash orchestration scripts
```bash
# Example: geneexpert_stroke_vs_control_v1.sh
#!/bin/bash
set -e

# Claude wrote this entire script!
# It orchestrates existing R scripts with custom logic

# Step 1: FastQC
fastqc data/*.fastq.gz -o results/fastqc/

# Step 2: Alignment (Claude adjusted parameters based on agent debate)
cd data/
for i in *_R1_001.fastq.gz; do
  fname=$(basename "$i" _R1_001.fastq.gz)
  subread-align -t 0 -i /genome/mm10 \
    -r "${fname}_R1_001.fastq.gz" \
    -R "${fname}_R2_001.fastq.gz" \
    -T 8 -o "${fname}.bam"
done

# Step 3: Feature Counts (Claude got parameters from reading featurecounts.R)
Rscript /scripts/featurecounts.R mm10 stroke_vs_control cont ips *.bam

# ... continues with all 9 steps
```

**Claude DOES NOT WRITE:** New R scripts or Python code
- Your existing R scripts (`featurecounts.R`, `RPKM.R`, etc.) are NEVER modified
- Claude just orchestrates them with correct parameters

### How Code Generation Works:

```
AUTOMATION Mode:
  Coordinator → Template engine → Fills in config variables
  Result: Standard bash script (no agent writing needed)

ADAPTATION Mode:
  Coordinator → Sends debate context to Claude MCP agent
  Claude → Reads existing scripts via read_file tool
  Claude → Understands parameters and requirements
  Claude → WRITES custom bash script addressing concerns
  System → EXECUTES Claude's script
  If fails → Claude WRITES v2 with fixes
```

**Key Point:** Claude is the ONLY agent that can read your files, write new scripts, and trigger execution. GPT-4 and Gemini provide analysis and recommendations, but can't touch code.

---

## Key Features (All Working! )

###  Implemented & Tested:

1. **Agent-Driven Decision Making**
   - Agents analyze data and vote AUTOMATION vs ADAPTATION
   - No hardcoded rules - truly dynamic
   - Shows vote breakdown and confidence (0-100%)

2. **User Confirmation Flow**
   - Beautiful plan display before execution
   - Shows: experiment details, agent decision, pipeline steps
   - Y/N confirmation prevents unwanted execution

3. **Script Execution with Live Output**
   - Real-time stdout/stderr streaming (using `spawn()`)
   - Exit code detection (success/failure)
   - All output captured to log files

4. **Smart File Handling**
   - Handles scripts with/without `.sh` extension automatically
   - MCP `read_file` tool tries multiple paths
   - No renaming of existing lab scripts needed

5. **Auto-Detection**
   - Data type: FASTQ, BAM, or count matrix
   - Paired-end: R1/R2 or _1/_2 patterns
   - Sample groups: User-specified keywords or common patterns

6. **Genome Build Translation**
   - User-friendly: `--organism mouse`
   - Scripts use technical: `mm10`
   - Prevents parameter errors

7. **Complete Logging**
   - All terminal output saved
   - Agent conversations logged separately
   - Session history tracked

8. **Error-Propagating Feedback Loop** (NEW! Jan 11, 2026)
   - Script fails → Agents analyze error automatically
   - Detect completed vs failed steps
   - Generate v2 script with fixes (skips completed work!)
   - Retry up to 3 attempts (v1 → v2 → v3)
   - Loop ends: Success OR max retries reached

9. **Robust Consensus Voting**
   - Fixed confidence calculation bug (Jan 6, 2026)
   - Correctly handles all vote types
   - Shows accurate reasoning ("3/3 agents recommend...")

10. **QC Plots with Agent Vision** (NEW! Jan 15, 2026)
    - `qc_plots.R` generates PCA, MDS, density plots as PNG
    - ALL 3 agents can VIEW images (GPT-4, Claude, Gemini have vision)
    - Agents VOTE on outlier detection (z-score > 2.5 SD from group centroid)
    - Agents VOTE on batch effects (separation score, PC1 correlation)
    - Unanimous vote required for sample removal

11. **Batch Effect Correction** (NEW! Jan 15, 2026)
    - `batch_effect_edgeR_v3.R` - same input format as standard edgeR
    - Uses design matrix `~batch + condition` when batch effects detected
    - Agents decide: standard edgeR OR batch-corrected edgeR

12. **Decision-Level Evaluation Framework** (NEW! Jan 15, 2026)
    - Every agent decision gets unique `decision_id` for ICML evaluation
    - Format: `{dataset}_{step}_{type}` (e.g., `DA0036_step0_mode`)
    - Runtime logging: decisions can be joined with ground truth offline
    - Decision-type aware disagreement scoring
    - `TOTAL_AGENTS` constant for scalability (currently 3, extensible to 4+)

---

##  QC Plots: Agent Vision for Outlier/Batch Detection

**The Innovation:** Agents don't just read numbers - they VIEW PCA plots and VOTE on sample quality!

### How It Works:

```
Step 7: QC Plots (after normalization)
    ↓
qc_plots.R generates:
  - PCA plot (PNG) - sample clustering
  - MDS plot (PNG) - distance visualization
  - Density plot (PNG) - library size distribution
  - qc_summary.json - structured metrics
    ↓
ALL 3 AGENTS VIEW THE IMAGES:
  [GPT-4]:   "Sample_3 is 3.2 SD from group centroid - outlier"
  [Claude]:  "Separation score 1.2 suggests batch effects"
  [Gemini]:  "Groups overlap significantly - check batch"
    ↓
CONSENSUS VOTING:
  - Outlier detection: Unanimous required (3/3) for removal
  - Batch effects: Majority vote (2/3) to use corrected edgeR
    ↓
IF batch effects detected:
  → Use batch_effect_edgeR_v3.R (design: ~batch + condition)
ELSE:
  → Use simpleEdger3.R (standard design: ~condition)
```

### Outlier Detection Criteria:
- **Z-score > 2.5** from group centroid in PCA space
- Flagged samples shown to agents with reasoning
- User always has final say (unanimous vote OR escalation)

### Batch Effect Detection Criteria:
- **Separation score < 1.5** (samples spread within groups)
- **PC1 correlation with group < 0.8**
- Either triggers batch correction recommendation

---

##  Error-Propagating Feedback Loop

**The Innovation:** When scripts fail, agents debug and generate fixed versions automatically.

### How It Works:

```
ATTEMPT 1: Execute generated script v1.sh
    ↓
  [FAILS: Exit code 1]
    ↓
FEEDBACK LOOP ACTIVATED:
    ↓
Step 1: Detect Progress
  - Completed: Steps 1-3 
  - Failed: Step 4 
    ↓
Step 2: Capture Error Context
  - Exit code: 1
  - Error message: "cannot open file test_run.count.txt"
  - STDOUT: Last 50 lines
  - STDERR: Full error output
    ↓
Step 3: Multi-Agent Debugging
  Coordinator asks all 3 agents:

  "PREVIOUS SCRIPT FAILED:
   - Completed: Steps 1-3
   - Failed: Step 4 (Filter Bad IDs)
   - Error: cannot open file test_run.count.txt

   YOUR TASK:
   1. SKIP completed steps (1-3)
   2. RESUME from Step 4
   3. Fix the error (likely wrong file extension)
   4. Continue with remaining steps"
    ↓
  GPT-4: "File extension mismatch. Use .csv not .txt"
  Claude: [calls list_available_scripts] "Confirmed: output is .csv"
  Gemini: "Agree with file extension fix"
    ↓
Step 4: Generate v2 Script
  - Skips Steps 1-3 (already done!)
  - Starts from Step 4 with fix
  - Continues to completion
    ↓
ATTEMPT 2: Execute v2.sh
    ↓
  [SUCCESS: Exit code 0] 
    ↓
ANALYSIS COMPLETE!
```

### Loop Stopping Criteria:

**The loop ENDS when:**
1.  **Script succeeds** (exit code 0) → Analysis complete!
2.  **Max retries reached** (3 attempts) → Report failure to user
3.  **User cancels** during confirmation → Stop execution

### Maximum Retry Attempts: **3**

```javascript
const MAX_RETRIES = 3;

Attempt 1: v1.sh (generated from agent decision)
Attempt 2: v2.sh (fixed based on v1 errors)
Attempt 3: v3.sh (fixed based on v2 errors)

If v3 fails → STOP and report to user
```

### Why This Matters:

**Traditional Pipeline:**
```
Script fails at Step 7 → User must manually debug → Hours wasted
```

**GeneExpert Feedback Loop:**
```
Script fails at Step 7 → Agents debug automatically → v2 generated
  ↓
Skips Steps 1-6 (saves 20-40 minutes!)
  ↓
Resumes from Step 7 with fix → Success!
```

### Example Error Correction:

**v1 Error:**
```bash
# v1.sh (WRONG):
Rscript filterIDS.R test_run.count.txt

# Error: File not found
```

**Agents Analyze:**
- GPT-4: "File extension mismatch detected"
- Claude: Calls `list_available_scripts` → sees .csv output
- Gemini: "All R scripts use CSV format since Jan 7"

**v2 Fix:**
```bash
# v2.sh (FIXED):
Rscript filterIDS.R test_run.count.csv

# Success! 
```

---

## ️ MCP Tools

### Implemented Tools:

**File Reading (for ADAPTATION):**
-  `read_file` - Read any text file (R scripts, logs, results)
  - Smart extension handling (tries with/without .sh)
  - Auto-discovery of script parameters

**Pipeline Execution (for both):**
-  `run_fastqc` - Quality control
-  `run_alignment` - STAR/HISAT2 alignment
-  `run_featurecounts` - Gene quantification
-  `run_filter` - Filter low counts
-  `run_rpkm` - RPKM normalization
-  `run_annotation` - Add gene symbols
-  `run_edger` - Differential expression
-  `export_to_excel` - Format results
-  Visualization: volcano, MA plots, Venn diagrams

### In Development:

**Analysis Tools (for future ADAPTATION enhancements):**
-  `read_bam_summary` - Parse alignment statistics
-  `read_count_summary` - Analyze count matrices
-  `write_custom_script` - Save custom R scripts
-  `execute_custom_script` - Run custom solutions

---

##  Project Structure

```
├── bin/geneexpert.js                    # CLI entry point
├── src/
│   ├── pipeline/
│   │   ├── executor.js                  # Main orchestration
│   │   ├── planner.js                   # Pipeline planning (10 steps)
│   │   ├── script_generator.js          # AUTOMATION & ADAPTATION scripts
│   │   ├── modification_engine.js       # Script modification for ADAPTATION
│   │   └── data_detector.js             # Auto-detect data type
│   ├── agents/
│   │   └── mcp_claude_agent.js          # MCP-enabled Claude (ADAPTATION)
│   ├── coordinator/
│   │   ├── orchestrator.js              # Multi-agent coordination + decision_id
│   │   └── consensus.js                 # Voting, disagreement scoring, TOTAL_AGENTS
│   ├── mcp/
│   │   ├── server.js                    # MCP server
│   │   └── tools.js                     # 14+ MCP tools
│   └── utils/
│       ├── llm_clients.js               # OpenAI, Anthropic, Google APIs
│       └── logger.js                    # Logging system
├── bio_informatics/
│   ├── scripts/                         # Lab R/bash scripts (existing)
│   │   ├── qc_plots.R                   # QC: PCA, MDS, outlier detection (NEW!)
│   │   ├── batch_effect_edgeR_v3.R      # Batch-corrected DE analysis (NEW!)
│   │   └── ...                          # Other existing scripts
│   └── myprog/                          # Reference data
├── experiments/
│   └── ground_truth.json                # Template for ICML evaluation
├── IMPLEMENTATION_STATUS.md              # Detailed progress tracking
├── HYBRID_ARCHITECTURE.md                # System design doc
└── README.md                             # This file
```

---

##  Example: Small Sample Size (n=2)

### Agent Analysis:

```bash
$ node bin/geneexpert.js analyze data/DA0036 \
    --organism mouse \
    --comparison "DA0036_stroke_vs_control" \
    --control-keyword "cont" \
    --treatment-keyword "ips" \
    --output results/DA0036_test
```

**What happens:**

```
[Data Detector] Detected: 4 samples, paired-end, 2 groups (cont vs ips)

[Coordinator] Consulting agents to decide approach...

[GPT-4 Stats Agent] ✓ Response received
────────────────────────────────────────────────────────────
I recommend ADAPTATION. The primary reason is the small sample
size (n=2 per group), which is particularly challenging for
differential expression studies. Small sample sizes increase the
risk of both Type I and Type II errors...
────────────────────────────────────────────────────────────

[Claude Pipeline Agent] ✓ Response received
────────────────────────────────────────────────────────────
I would recommend using ADAPTATION. While the proposed pipeline
covers standard steps, there are issues that warrant a customized
approach:
1. Small sample size: Only 2 replicates per group
2. Potential batch effects: Need careful assessment
────────────────────────────────────────────────────────────

[Gemini Biology Agent] ✓ Response received
────────────────────────────────────────────────────────────
The dataset's primary issue is the very low number of replicates
(n=2 per group). This severely limits statistical power.
Therefore, ADAPTATION is recommended to address these limitations.
────────────────────────────────────────────────────────────

[Coordinator]  Synthesizing consensus...
[Coordinator] Decision: ADAPTATION
[Coordinator] Confidence: 100%
[Coordinator] Reasoning: 3/3 agents recommend ADAPTATION

[Script Generator] Generating ADAPTATION script...
[MCP Claude Agent] Reading actual scripts to discover parameters...
[MCP Claude Agent] Tool call: read_file
  Reading: /path/to/bio_informatics/scripts/fastq2bam
  (Found as: fastq2bam - no .sh extension)

[Script Generator] ✓ Custom script saved to:
  results/DA0036_test/geneexpert_DA0036_stroke_vs_control_v1.sh

─────────────────────────────────────────────────────────────
 ANALYSIS PLAN
─────────────────────────────────────────────────────────────
Experiment: DA0036_stroke_vs_control
Organism: mouse (mm10)
Samples: 4 (paired-end)
Groups: cont (n=2) vs ips (n=2)

 Agent Decision: ADAPTATION (100% confidence)
   Reason: 3/3 agents recommend custom approach for small n

 Generated Script:
   results/DA0036_test/geneexpert_DA0036_stroke_vs_control_v1.sh

Pipeline Steps:
  0. FASTQ Validation
  1. FastQC
  2. Alignment (Subread-align)
  3. Feature Counts
  4. Filter Bad IDs
  5. RPKM Normalization (for visualization)
  6. Add Gene Symbols
  7. QC Plots (PCA, MDS) → Agents VIEW & VOTE
  8. DE Analysis (edgeR with raw counts)
  9. Merge RPKM + DE Results into Excel
  10. Complete

  This will execute the generated script!

Proceed? (Y/N):
```

**User types Y:**

```
[Executor] Running script...

Step 1/10: Running FastQC...
[live output streaming...]

Step 2/10: Alignment (FASTQ to BAM)...
[live output streaming...]

...

[Executor]  Analysis complete!
 Results saved to: results/DA0036_test/
```

**This is true intelligence!** Agents saw n=2, debated the risk, and chose ADAPTATION.

---

##  Development Status (Jan 15, 2026)

###  SYSTEM OPERATIONAL - QC & Evaluation Framework Complete!

**Completed:**
-  Multi-agent decision making (AUTOMATION vs ADAPTATION)
-  Consensus voting with accurate confidence calculation
-  AUTOMATION script generation (template-based, 10 steps)
-  ADAPTATION script generation (MCP-powered, intelligent)
-  User confirmation flow with beautiful plan display
-  Script execution with real-time output streaming
-  Smart file handling (with/without .sh extension)
-  Data detection and auto-grouping
-  Genome build translation (mouse → mm10)
-  Complete logging system
-  CLI interface

**Bugs Fixed (Jan 6, 2026):**
-  Confidence calculation bug (0% → 100% for unanimous votes)
-  File extension handling (fastq2bam.sh vs fastq2bam)
-  MCP agent file reading errors

**Major Achievement (Jan 11, 2026):**
-  First successful end-to-end run (DA0036 dataset)
-  All 10 steps completed from FASTQ → Excel
-  Feedback loops working (error detection + v2 generation)
-  Intelligent step detection (skip completed steps on retry)

**Major Achievement (Jan 15, 2026):**
-  QC Plots with Agent Vision (`qc_plots.R`)
  - PCA, MDS, density plots as PNG for agents to VIEW
  - Outlier detection (z-score > 2.5 SD from group centroid)
  - Batch effect detection (separation score, PC1 correlation)
-  Batch Effect Correction (`batch_effect_edgeR_v3.R`)
  - Same input format as `simpleEdger3.R`
  - Uses `~batch + condition` design matrix
-  Decision-Level Evaluation Framework
  - Unique `decision_id` for every agent decision
  - Runtime logging for ICML paper evaluation
  - Decision-type aware disagreement scoring
  - `TOTAL_AGENTS` constant for scalability (3 → 4+ agents)

**Next Steps:**
- Week 3 (Jan 15-17): Test with different datasets (clean, problematic, edge cases)
- Week 3 (Jan 15-17): Collect baseline metrics (no-agent, single-agent, multi-agent)
- Week 4 (Jan 18-22): Measure error reduction, consensus quality, cost
- Week 5 (Jan 23-28): Write ICML paper

**Days to ICML Deadline:** 13 days

---

##  Cost Estimate

### AUTOMATION Path:
- Agent decision-making: ~$0.03 (one-time, 3 agents analyze data)
- Script generation: $0 (template-based)
- Execution: $0 (no API calls)
- **Total: ~$0.03 per analysis**

### ADAPTATION Path:
- Agent decision-making: ~$0.03
- MCP tool calls (read scripts): ~$0.20 (Claude reads 5-10 scripts)
- Custom script generation: ~$0.30 (Claude writes custom bash)
- Execution: $0
- **Total: ~$0.53 per analysis**

**Both very affordable!** Even ADAPTATION is <$1 per complete RNA-seq analysis.

---

##  Documentation

- **IMPLEMENTATION_STATUS.md** - Complete progress tracking (updated Jan 6, 2026)
- **HYBRID_ARCHITECTURE.md** - Detailed system design
- **PROJECT_PURPOSE.md** - Research goals & ICML contribution
- **Guide.md** - User guide

---

##  Development Commands

```bash
# Test LLM APIs
npm run test:apis

# Run analysis (full pipeline)
node bin/geneexpert.js analyze <input> \
  --organism mouse|human|rat \
  --comparison "experiment_name" \
  --control-keyword "ctrl" \
  --treatment-keyword "treat" \
  --output <output_dir>

# Check logs
tail -f results/my_analysis/analysis_*.log
```

---

##  Research Contribution (ICML 2026)

**Hypothesis:** Multi-foundation model collaboration with dynamic approach selection reduces errors 40%+ compared to single-agent or static pipeline systems.

**Novel Contributions:**

1. **AUTOMATION + ADAPTATION Hybrid System**
   - Agents dynamically choose execution strategy
   - Template-based for speed, MCP-powered for intelligence

2. **MCP-Enabled Intelligence**
   - Agents READ actual scripts (not just documentation)
   - Discover parameters through file inspection
   - Write truly custom solutions

3. **Multi-Model Consensus**
   - Different foundation models (GPT-4, Claude, Gemini)
   - Voting system with confidence quantification
   - Disagreement signals uncertainty

4. **Adaptive Execution**
   - Custom script generation for edge cases
   - Addresses specific data characteristics
   - No one-size-fits-all assumptions

5. **User-in-the-Loop**
   - Transparent decision display
   - Confirmation before execution
   - Trust through visibility

**Target Conference:** ICML 2026 (International Conference on Machine Learning)
**Submission Deadline:** January 28, 2026
**Status:** System complete, ready for benchmark testing

---

##  License

MIT

---

##  Links

- **GitHub:** https://github.com/Mituvinci/geneexpert-mcp
- **Status:**  System operational - First successful end-to-end run complete!
- **Last Updated:** January 15, 2026

---

**Ready to revolutionize bioinformatics analysis with multi-agent intelligence!** 

See `IMPLEMENTATION_STATUS.md` for detailed progress tracking.
