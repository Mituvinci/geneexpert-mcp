# GeneExpert MCP Server

Multi-agent bioinformatics system where Claude, GPT-4, and Gemini collaborate to analyze RNA-seq and ATAC-seq data.

## Quick Start

### 1. Configuration

Copy the example environment file and add your API keys:

```bash
cp .env.example .env
nano .env  # Add your API keys
```

### 2. Test the MCP Server

```bash
npm test
```

You should see all 14 tools loaded successfully.

### 3. Start the MCP Server

```bash
npm start
```

The server runs on stdio and communicates via the MCP protocol.

## Available Tools

### RNA-seq Pipeline (14 tools)

**Phase 1: Pre-processing**
- `run_fastqc` - Quality control on FASTQ files
- `run_alignment` - Align reads to reference genome (STAR/HISAT2)
- `run_featurecounts` - Count reads per gene

**Phase 2: Normalization & QC**
- `run_filter` - Remove lowly expressed genes
- `run_rpkm` - RPKM normalization
- `run_annotation` - Add gene symbols and descriptions
- `run_qc_plots` - **CRITICAL**: PCA, MDS, density plots

**Phase 3: Differential Expression**
- `run_edger` - Statistical analysis (edgeR)
- `export_to_excel` - Comprehensive results file

**Phase 4: Visualization**
- `run_volcano_plot` - DE visualization
- `run_ma_plot` - Mean vs fold-change
- `run_venn` - Compare DEG lists

**ATAC-seq**
- `run_atac_alignment` - ATAC-seq pipeline
- `run_bedgraph` - Genome browser tracks

## Project Structure

```
/data/halimaakhter/multi_llm_mcp/
├── src/
│   ├── mcp/
│   │   ├── server.js       # MCP server (✓ DONE)
│   │   ├── tools.js        # Tool wrappers (✓ DONE)
│   │   └── test.js         # Test script
│   ├── coordinator/        # TODO: Multi-agent orchestrator
│   ├── agents/             # TODO: Claude, GPT-4, Gemini agents
│   ├── collaboration/      # TODO: Debate & consensus
│   └── utils/              # TODO: LLM API clients
├── config/                 # TODO: Agent configs
├── data/                   # Experimental datasets
├── results/                # Pipeline outputs
└── reviews/                # Multi-LLM reviews
```

## Next Steps

### Week 1 (Current)
- [x] Build MCP server wrapping tools
- [ ] Implement Coordinator Agent
- [ ] Test basic tool execution on sample data

### Week 2
- [ ] Pipeline Agent (Claude + MCP)
- [ ] Stats Agent (GPT-4)
- [ ] Biology Agent (Gemini)
- [ ] QC Agent (Claude)

### Week 3
- [ ] Debate/voting system
- [ ] Consensus mechanisms
- [ ] End-to-end testing

## Wrapped Scripts

**Lab scripts** (`/data/scripts/`):
- FastQC, STAR, HISAT2, Bowtie2
- featureCounts, SAMtools
- edgeR, DESeq2 wrappers
- Visualization tools

**Custom scripts** (`/destiny/halima/dory/my_script/`):
- `featurecounts_edited.R`
- `MDS_mdfy.R`
- `maplot.R`
- `volcano.R`

## API Keys Required

- **OpenAI** (GPT-4): Stats validation
- **Anthropic** (Claude): Pipeline execution & QC
- **Google AI** (Gemini): Biology interpretation

Store in `.env` file (never commit!)

## Research Goal

**ICML 2026 Paper**: Demonstrate that multi-agent collaboration reduces errors by 40%+ vs single-agent approaches in genomic analysis.

**Key Innovation**: Specialized agent roles + consensus-based decision making + disagreement as uncertainty signal.

## Development

```bash
# Test configuration
npm test

# Start MCP server
npm start

# Start with auto-reload (development)
npm run dev
```

## License

MIT

---

**Status**: Week 1 - MCP Server Complete ✓
**Next**: Build Coordinator Agent & test pipeline execution
