# GeneExpert: Multi-Agent Bioinformatics System

## Project Overview
Multi-agent MCP server where Claude, GPT-4, and Gemini collaborate as specialized teammates to solve bioinformatics problems. Each agent has a specific role; they debate decisions and reach consensus.

**Goal:** Build a practical tool people can install + ICML 2026 paper (deadline: Jan 28, 2026)  
**Research Focus:** Multi-agent collaboration reduces errors by 40%+ vs single-agent

---

## System Environment

**Server:** OpenSUSE Linux at `/data/`  
**Project:** `/data/halimaakhter/multi_llm_mcp/`  
**Lab Scripts:** `/data/scripts/` (all bioinformatics tools pre-installed)  
**Custom Scripts:** `/data/halimaakhter/my_script/`

**Pre-installed Tools:**
- Alignment: STAR, HISAT2, Bowtie2
- QC: FastQC, MultiQC
- Quantification: featureCounts, RSEM
- DE Analysis: edgeR, DESeq2, limma
- Utilities: SAMtools, BEDtools, deepTools
- R libraries: All loaded
- Conda environments: Available

**Requirements:**
- No GPU needed (CPU-based server)
- API keys: OpenAI, Anthropic, Google AI
- Node.js + npm for MCP server
- Python 3.8+ for orchestration

---

## Architecture

### Multi-Agent Team
```
Coordinator Agent (orchestrates)
    ↓
├─ Pipeline Agent (Claude + MCP) - Executes tools
├─ Stats Agent (GPT-4) - Validates statistics
├─ Biology Agent (Gemini) - Interprets biology
└─ QC Agent (Claude) - Monitors quality
    ↓
Agents debate → Vote → Consensus → Result
```

### Agent Roles

**Coordinator:** Routes tasks, synthesizes consensus, presents results  
**Pipeline Agent:** Runs bioinformatics tools via MCP, manages workflows  
**Stats Agent:** Validates thresholds (FDR, logFC), checks assumptions, detects outliers  
**Biology Agent:** Interprets pathways, validates coherence, suggests follow-ups  
**QC Agent:** Monitors PCA/MDS plots, flags technical issues, recommends sample removal

### Collaboration Rules

| Decision Type | Consensus Needed |
|---------------|------------------|
| Run standard tool | Pipeline Agent only |
| Select threshold | Majority vote (3/4) |
| Remove sample | Unanimous (4/4) |
| Change methods | Majority + user approval |

**Disagreement handling:**
- Minor (1 dissent): Proceed, note dissent
- Major (2-2 split): Pause, present options to user
- Critical safety: Immediate stop, require confirmation

---

## RNA-seq Pipeline (Biologically Correct)

### Phase 1: Pre-processing
1. FastQC (quality control)
2. STAR/HISAT2 (alignment → BAM)
3. featureCounts (quantification → count matrix)

### Phase 2: Normalization & QC
4. filterIDS.R (remove low counts)
5. RPKM.R (normalize)
6. entrz.R (add gene annotations)
7. **QC plots** (PCA, MDS, density) ← MUST DO BEFORE DE

### Phase 3: Differential Expression
8. simpleEdger3.R (DE analysis)
9. Export to Excel (combine DEGs + annotations + RPKM)

### Phase 4: Visualization
10. Volcano plots, MA plots, heatmaps
11. venn_diagram.R (compare DEG lists)

### Phase 5: Downstream (optional)
12. GO enrichment, pathway analysis
13. setgraph4ucsc (bedGraph for genome browser)

**Critical:** Always run QC plots (step 7) before DE analysis (step 8). Multi-agent review required at QC stage.

---

## Project Structure

```
/data/halimaakhter/multi_llm_mcp/
├── AGENTS.md                    # This file (Claude Code reads)
├── RESEARCH_CONTEXT.md          # Paper details (reference when needed)
├── package.json
├── .env                         # API keys (never commit)
│
├── src/
│   ├── coordinator/
│   │   ├── orchestrator.js      # Main coordinator
│   │   └── consensus.js         # Voting logic
│   ├── agents/
│   │   ├── pipeline_agent.js    # Claude + MCP
│   │   ├── stats_agent.js       # GPT-4
│   │   ├── biology_agent.js     # Gemini
│   │   └── qc_agent.js          # Claude
│   ├── subagents/
│   │   ├── pathway_analyzer.js
│   │   ├── literature_search.js
│   │   └── qc_diagnostics.js
│   ├── mcp/
│   │   ├── server.js            # MCP server
│   │   └── tools.js             # Wrap /data/scripts
│   ├── collaboration/
│   │   ├── debate.js
│   │   ├── voting.js
│   │   └── escalation.js
│   └── utils/
│       ├── llm_clients.js       # API wrappers
│       └── logging.js
│
├── config/
│   ├── agent_roles.yaml
│   ├── consensus_rules.yaml
│   └── escalation_policy.yaml
│
├── data/                        # Experimental datasets
├── experiments/                 # Baselines & tests
├── results/                     # Pipeline outputs
└── paper/                       # ICML manuscript
```

---

## MCP Tools to Implement

### Core Bioinformatics (wrap existing scripts)
1. `run_fastqc` → FastQC
2. `run_alignment` → fastq2bam
3. `run_featurecounts` → `/data/halimaakhter/my_script/featurecounts_edited.R`
4. `run_filter` → `filterIDS.R`
5. `run_rpkm` → `RPKM.R`
6. `run_annotation` → `entrz.R`
7. `run_qc_plots` → PCA, MDS, density
8. `run_edger` → `simpleEdger3.R`
9. `export_to_excel` → Format DEG results with annotations

### Visualization
10. `run_volcano_plot`
11. `run_ma_plot`
12. `run_heatmap`
13. `run_venn` → `venn_diagram.R`

### Multi-LLM Review
14. `review_with_gpt4` → Send results to GPT-4 API
15. `review_with_claude` → Send to Claude API
16. `review_with_gemini` → Send to Gemini API
17. `compare_reviews` → Analyze agreement/disagreement
18. `consensus_report` → Generate decision summary

---

## Tech Stack

**MCP Server:** Node.js with `@modelcontextprotocol/sdk`  
**Orchestration:** Python or JavaScript  
**LLM APIs:**
- OpenAI SDK (`openai`)
- Anthropic SDK (`anthropic`)
- Google Generative AI (`google-generativeai`)

**Bioinformatics:** R scripts (already in `/data/scripts/`)

---

## API Configuration

Store in `.env` (never commit):
```bash
OPENAI_API_KEY=sk-...
ANTHROPIC_API_KEY=sk-ant-...
GOOGLE_API_KEY=...
```

**Cost estimate:** ~$10-30 per full RNA-seq analysis

---

## Development Priorities

### Week 1 (Dec 26 - Jan 1): Infrastructure
1. Build MCP server wrapping `/data/scripts/` tools
2. Implement Coordinator Agent
3. Test basic tool execution

### Week 2 (Jan 2 - Jan 8): Agents
1. Pipeline Agent (Claude + MCP)
2. Stats Agent (GPT-4 API)
3. Biology Agent (Gemini API)
4. QC Agent (Claude API)

### Week 3 (Jan 9 - Jan 15): Collaboration
1. Debate/voting system
2. Consensus mechanisms
3. Test on sample dataset

### Week 4 (Jan 16 - Jan 22): Experiments
1. Run 3 datasets (clean, corrupted, ambiguous)
2. Collect metrics

### Week 5 (Jan 23 - Jan 28): Paper
1. Write manuscript
2. Submit to ICML 2026

---

## Key Constraints

**Timeline:** 33 days until Jan 28, 2026  
**Token efficiency:** Minimize API calls where possible  
**Biological correctness:** Always run QC before DE analysis  
**User approval:** Required for sample removal, method changes

---

## Installation Target

Final product should install like:
```bash
npm install -g geneexpert-mcp
geneexpert config --setup
geneexpert analyze ./my_rnaseq_data/
```

---

## Research Contribution

**Main hypothesis:** Multi-agent collaboration reduces errors by 40%+ vs single-agent

**What makes this novel:**
- First multi-LLM collaborative system for science (not just coding)
- Specialized agent roles (each LLM handles what it's best at)
- Consensus mechanisms for critical decisions
- Disagreement = uncertainty quantification

**NOT novel:** Using MCP (just infrastructure), using LLMs for bioinformatics (already exists)

---

## Important Notes

### For Claude Code Development

**Focus on:**
- Building working MCP server first
- Implementing agent communication
- Getting one complete analysis working end-to-end
- Code quality > premature optimization

**Don't worry about:**
- Perfect accuracy (baseline first, improve later)
- Paper writing (separate phase)
- UI polish (terminal-based is fine)

### Biological Best Practices

**Must follow:**
- QC plots BEFORE differential expression
- Check PCA for outliers before proceeding
- Export results with full annotations (gene symbols, RPKM, etc.)
- Never skip quality control steps

### Multi-Agent Principles

**Agents should:**
- Actively collaborate (not just review)
- Debate decisions in real-time
- Reach consensus before proceeding
- Escalate to user when split

**Avoid:**
- Silent failures
- Single agent making critical decisions
- Skipping consensus on important choices

---

## Quick Reference

**Main lab scripts:** `/data/scripts/`  
**Custom scripts:** `/data/halimaakhter/my_script/`  
**Example output:** `liver_KOECvsKOCT_WTECvsWTCT_edger.xlsx`

**Key R scripts:**
- `featurecounts_edited.R` - Gene quantification
- `simpleEdger3.R` - Differential expression
- `filterIDS.R` - Filter low counts
- `RPKM.R` - Normalization
- `entrz.R` - Add gene annotations
- `venn_diagram.R` - Compare DEG lists

**ATAC-seq scripts:**
- `atac_seq_alignment`
- `atac_fastq_to_bam`
- `setgraph4ucsc` - Bedgraph for visualization

---

## For More Details

See `RESEARCH_CONTEXT.md` for:
- Full experimental design
- ICML paper outline
- Detailed workflow examples
- Literature review
- Comparison to related work

---

**Last Updated:** December 26, 2024  
**Status:** Ready to start building  
**Next Step:** Implement MCP server wrapping `/data/scripts/` tools
