# GeneExpert: HYBRID Multi-Agent RNA-seq Analysis

Multi-agent system where Claude (with MCP), GPT-4, and Gemini collaborate to analyze RNA-seq data. **HYBRID approach:** fast automation for standard steps, true agentic behavior at decision points.

**Research Goal:** ICML 2026 paper - Multi-agent collaboration reduces errors 40%+ through adaptive decision-making

---

## 🚀 Quick Start

### 1. Install Dependencies

```bash
npm install
```

### 2. Configure API Keys

```bash
cp .env.example .env
nano .env  # Add your OpenAI, Anthropic, Google API keys
```

### 3. Run Analysis

```bash
node bin/geneexpert.js analyze /path/to/fastqs/ \
  --output results/my_analysis \
  --organism mouse \
  --comparison "treatment_vs_control"
```

**That's it!** The system:
- Detects your data automatically
- Runs full pipeline (alignment, counting, DE analysis)
- Agents debate at decision points
- Saves results with full logs

---

## 🎯 What Makes This Novel

### Standard Automation vs HYBRID Agentic System

**Traditional automation:**
```
Run pipeline → Fixed steps → Done
(Can't adapt to batch effects, outliers, or edge cases)
```

**GeneExpert HYBRID:**
```
Run standard steps (fast automation)
    ↓
DECISION POINT: QC Review
    ↓
MCP Claude Agent:
- read_bam_summary → Sees actual mapping rates
- read_count_summary → Sees actual counts
- DETECTS: "Sample 2 has 45% mapping - investigate!"
    ↓
Multi-agent debate based on REAL data
    ↓
If edge case: Write custom script to solve
    ↓
Continue with corrected data
```

**Key Innovation:** Agents READ actual data via MCP and ADAPT to problems, not just run predefined steps!

---

## 📊 Architecture

### Standard Steps (Automation - $0 cost)
- FastQC → fastq2bam (Subread) → featureCounts → filterIDS → RPKM
- JavaScript executes `/data/scripts/` via bash
- Fast, cheap, reliable

### Decision Points (Agentic - ~$0.03 per decision)
1. **QC Review:**
   - MCP Claude reads alignment stats & count summaries
   - Detects problems (low mapping, batch effects, outliers)
   - Multi-agent debate: GPT-4 + Claude + Gemini
   - If edge case: Claude writes custom R script

2. **Threshold Selection:**
   - Multi-agent debate on FDR/logFC based on sample size
   - Small N → More stringent thresholds

**Total cost:** ~$0.055 per analysis

---

## 🛠️ MCP Tools (19 tools)

### Execution Tools (14):
- RNA-seq pipeline: fastq2bam, featureCounts, filterIDS, RPKM, edgeR
- QC plots: PCA, MDS, density
- Visualization: volcano, MA plots, Venn diagrams

### Analysis Tools (5 - NEW! For agentic behavior):
- `read_bam_summary` - Read alignment statistics
- `read_count_summary` - Analyze count matrix
- `read_file` - Read any output
- `write_custom_script` - Write custom R/bash scripts
- `execute_custom_script` - Execute custom solutions

**Agents use analysis tools to SEE real data and ADAPT!**

---

## 📁 Project Structure

```
├── bin/geneexpert.js           # CLI entry point
├── src/
│   ├── pipeline/               # Orchestration & planning
│   ├── agents/
│   │   ├── pipeline_agent.js   # Bash executor
│   │   └── mcp_claude_agent.js # MCP-enabled Claude (agentic!)
│   ├── coordinator/            # Multi-agent coordination
│   ├── mcp/                    # MCP server & 19 tools
│   └── utils/                  # LLM APIs, logging
├── AGENTS.md                    # Full architecture
├── HYBRID_ARCHITECTURE.md       # HYBRID system details
├── PROJECT_PURPOSE.md           # Research goals
└── README.md                    # This file
```

---

## 🧪 Example: Handling Batch Effect (Edge Case)

### Standard System:
```
Run pipeline → Generate PCA → Run DE → Done
(Batch effect goes undetected - ERRORS!)
```

### GeneExpert HYBRID:
```
Run pipeline → QC Review Decision Point
    ↓
[MCP Claude Agent] read_count_summary
[MCP Claude Agent] DETECTED: Samples cluster by batch, not treatment!
    ↓
Multi-agent debate:
  Stats: "Need batch correction"
  Biology: "Agree, use ComBat or batch term in model"
  Claude: "I'll write custom edgeR script"
    ↓
Claude writes:
  design <- model.matrix(~batch + treatment)
    ↓
execute_custom_script → Problem solved! ✅
```

**This is TRUE agentic behavior!**

---

## 📈 Development Status

### ✅ Complete (Week 1):
- HYBRID architecture implemented
- MCP server with 19 tools
- Pipeline Agent (bash executor)
- MCP Claude Agent (tool calling)
- Multi-agent coordinator
- Logging system
- CLI interface

### 🔄 Current (Week 2):
- Testing on real data (DA0036 - stroke vs control)
- Validating MCP tool calling
- Verifying edge case handling

### 📋 Upcoming:
- Week 3: Run 3 benchmark datasets
- Week 4: Measure error reduction (target: 40%+)
- Week 5: Write ICML paper (due Jan 28, 2026)

---

## 💡 Key Features

✅ **ONE command** - Full analysis from FASTQ to DEGs
✅ **Auto-detection** - FASTQ, BAM, or count matrix
✅ **Fast automation** - Standard steps run without API calls
✅ **MCP-enabled** - Agents READ actual data at decision points
✅ **Multi-agent debate** - GPT-4 + Claude + Gemini collaborate
✅ **Adaptive** - Writes custom scripts for edge cases
✅ **Full logging** - All agent conversations saved
✅ **Cheap** - ~$0.055 per analysis

---

## 📊 API Requirements

- **OpenAI API** (GPT-4): Stats validation
- **Anthropic API** (Claude): Pipeline + MCP tool calling
- **Google AI API** (Gemini): Biology interpretation

**Cost per analysis:** ~$0.055
**Gemini free tier:** 20 requests/day (10 analyses/day)

---

## 🔬 Research Contribution

**Hypothesis:** Multi-foundation model collaboration reduces errors 40%+ through adaptive decision-making

**Novel aspects:**
1. HYBRID approach (automation + agentic decision points)
2. MCP tool calling (agents see real data, not summaries)
3. Multi-model consensus (different LLMs debate)
4. Adaptive execution (custom scripts for edge cases)
5. Uncertainty quantification (disagreement signals ambiguity)

**Target:** ICML 2026

---

## 📚 Documentation

- **AGENTS.md** - Full architecture & agent roles
- **HYBRID_ARCHITECTURE.md** - Detailed HYBRID system explanation
- **PROJECT_PURPOSE.md** - Research goals & contribution

---

## 🛠️ Development

```bash
# Test LLM APIs
npm run test:apis

# Test coordinator (simulated)
npm run test:coordinator

# Run analysis
node bin/geneexpert.js analyze <input> --output <dir>
```

---

## 📝 License

MIT

---

## 🔗 Links

- **GitHub:** https://github.com/Mituvinci/geneexpert-mcp
- **Status:** HYBRID system complete, testing phase

---

**Ready to test!** See `AGENTS.md` for full details.
