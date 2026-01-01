# DA0036 Analysis - Stroke (ips) vs Control (cont)

## Dataset Information
- **Study:** DA0036
- **Comparison:** Stroke side (ips) vs Non-stroke side (cont)
- **Design:** 2 vs 2 samples (paired-end RNA-seq)
- **Samples:**
  - **Stroke (ips):** GD700, GD701
  - **Control (cont):** GD702, GD703

## Data Location
- **Input (READ-ONLY):** `/destiny/fastqs/DA0036/`
- **Output (YOUR FOLDER):** `/data/halimaakhter/multi_llm_mcp/test_results/DA0036/`

⚠️ **IMPORTANT:** Never write outputs to `/destiny/` - always to `/data/halimaakhter/`

---

## Analysis Pipeline

### Step 1: Quality Control (FastQC)
```bash
# Run the script
./01_run_fastqc.sh

# Save output to file
./01_run_fastqc.sh > 01_fastqc_output.txt 2>&1
```

**Output:** `fastqc/*.html` reports

---

### Step 2: Consult Multi-Agent System (After FastQC)
After reviewing the FastQC HTML reports, ask the agents for advice:

```bash
cd /data/halimaakhter/multi_llm_mcp
npm run test:coordinator > test_results/DA0036/02_agent_consultation.txt 2>&1
```

---

### Step 3: Alignment (TBD)
- STAR or HISAT2?
- Reference genome?

### Step 4: Quantification (TBD)
- featureCounts or RSEM?

### Step 5: Differential Expression (TBD)
- edgeR or DESeq2?

---

## Notes
- All commands run via `.sh` scripts for history tracking
- npm outputs saved to `.txt` files
- Small sample size (2 vs 2) - perfect test case for multi-agent consensus!
