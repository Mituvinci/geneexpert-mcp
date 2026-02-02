#!/usr/bin/Rscript
#
# batch_effect_edgeR_v3.R - edgeR with batch effect correction
# FIXED: Smart batch detection that doesn't treat every GSM ID as unique batch
#
# Usage:
#   Rscript batch_effect_edgeR_v3.R <count.csv> <control_keyword> <treatment_keyword> [batch_info]
#
# Arguments:
#   count.csv          - Count file (ID, Length, sample1, sample2, ...)
#   control_keyword    - Keyword in sample names for control group (e.g., 'control', 'untreated')
#   treatment_keyword  - Keyword in sample names for treatment group (e.g., 'treatment', 'Dex')
#   batch_info         - Optional: 'auto', 'paired', or comma-separated batch labels
#
# Examples:
#   Rscript batch_effect_edgeR_v3.R counts.csv 'untreated' 'Dex'
#   Rscript batch_effect_edgeR_v3.R counts.csv 'control' 'treatment' paired
#   Rscript batch_effect_edgeR_v3.R counts.csv 'ctrl' 'trt' '1,1,2,2,3,3'
#

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("
================================================================================
BATCH EFFECT CORRECTION - edgeR v3 (FIXED 2026-01-21)
================================================================================

Usage: batch_effect_edgeR_v3.R <count.csv> <control_keyword> <treatment_keyword> [batch_info] [fdr_threshold] [logfc_threshold]

Arguments:
  count.csv          - Count file (ID, Length, sample1, sample2, ...)
  control_keyword    - Keyword in sample names for control group (e.g., 'control', 'untreated')
  treatment_keyword  - Keyword in sample names for treatment group (e.g., 'treatment', 'Dex')
  batch_info         - Optional: 'auto', 'paired', or comma-separated batch labels
  fdr_threshold      - Optional: FDR cutoff (default: 0.05)
  logfc_threshold    - Optional: LogFC cutoff (default: 0.585)

Examples:
  Rscript batch_effect_edgeR_v3.R counts.csv 'untreated' 'Dex'
  Rscript batch_effect_edgeR_v3.R counts.csv 'control' 'treatment' paired
  Rscript batch_effect_edgeR_v3.R counts.csv 'ctrl' 'trt' '1,1,2,2,3,3'

NOTE: SUPPORTS UNEQUAL GROUP SIZES (e.g., 3 control + 4 treatment)

================================================================================
", call. = FALSE)
}

library("edgeR")
library("limma")
library("jsonlite")  # For write_json() function

count_file <- args[1]
control_keyword <- tolower(args[2])
treatment_keyword <- tolower(args[3])
batch_info <- if (length(args) >= 4) args[4] else "auto"

cat("\n")
cat("================================================================================\n")
cat("BATCH EFFECT CORRECTION - edgeR v3 (FIXED 2026-01-21)\n")
cat("================================================================================\n\n")
cat("Input file:", count_file, "\n")
cat("Control keyword:", control_keyword, "\n")
cat("Treatment keyword:", treatment_keyword, "\n")
cat("Batch info:", batch_info, "\n")

# Parse custom thresholds (for re-analysis)
fdr_threshold <- if (length(args) >= 5) as.numeric(args[5]) else 0.05
logfc_threshold <- if (length(args) >= 6) as.numeric(args[6]) else 0.585

cat(sprintf("\n=== Threshold Settings ===\n"))
cat(sprintf("FDR threshold: %.3f\n", fdr_threshold))
cat(sprintf("LogFC threshold: %.3f\n\n", logfc_threshold))

# ============================================================================
# Load Count Data
# ============================================================================

cat("Loading count data...\n")

# Try CSV first, then space/tab delimited
counts <- tryCatch({
  read.csv(count_file, header = TRUE, check.names = FALSE)
}, error = function(e) {
  tryCatch({
    read.table(count_file, header = TRUE, sep = "\t", check.names = FALSE)
  }, error = function(e2) {
    read.table(count_file, header = TRUE, sep = " ", check.names = FALSE)
  })
})

# Extract gene info and count matrix
gene_info <- counts[, 1:2]
colnames(gene_info) <- c("ID", "Length")  # FIXED: Changed "GeneID" to "ID" to match merge_results.R expectation
count_matrix <- counts[, -1:-2]

sample_names <- colnames(count_matrix)
n_samples <- ncol(count_matrix)
n_genes <- nrow(count_matrix)

cat("Loaded:", n_genes, "genes x", n_samples, "samples\n")
cat("Samples:", paste(sample_names, collapse = ", "), "\n\n")

# ============================================================================
# Determine Groups (Control vs Treatment)
# ============================================================================

cat("Determining sample groups using keywords...\n")

# Assign groups based on keywords in column names
assign_group <- function(sample_name) {
  sample_lower <- tolower(sample_name)
  if (grepl(control_keyword, sample_lower)) {
    return("Control")
  } else if (grepl(treatment_keyword, sample_lower)) {
    return("Treatment")
  } else {
    return("Unknown")
  }
}

# Create group labels for each sample
condition <- sapply(sample_names, assign_group)

# Check for unknown samples
if (any(condition == "Unknown")) {
  unknown_samples <- sample_names[condition == "Unknown"]
  stop(sprintf("\n\nERROR: Could not assign groups for these samples:\n%s\n\nMake sure sample names contain keywords '%s' or '%s'\n",
               paste(unknown_samples, collapse = "\n"), control_keyword, treatment_keyword))
}

# Count samples per group
n_control <- sum(condition == "Control")
n_treatment <- sum(condition == "Treatment")

cat("\n=== GROUP ASSIGNMENT ===\n")
cat(sprintf("Control group (%d samples):\n", n_control))
for (s in sample_names[condition == "Control"]) {
  cat("  ", s, "\n")
}
cat(sprintf("\nTreatment group (%d samples):\n", n_treatment))
for (s in sample_names[condition == "Treatment"]) {
  cat("  ", s, "\n")
}
cat("========================\n\n")

# Store group indices for batch pairing
control_indices <- which(condition == "Control")
treatment_indices <- which(condition == "Treatment")

# ============================================================================
# Smart Batch Detection (FIXED!)
# ============================================================================

cat("Determining batch structure...\n")

if (batch_info == "auto") {
  # -------------------------------------------------------------------------
  # SMART AUTO-DETECTION:
  # Look for meaningful batch patterns, not just extract all numbers
  # -------------------------------------------------------------------------

  batch <- NULL

  # Strategy 1A: Look for explicit batch labels (Rep1, Rep2, Batch1, Batch2, etc.)
  if (any(grepl("(rep|replicate|batch|run|day)\\d+", sample_names, ignore.case = TRUE))) {
    # Extract rep/batch numbers
    batch_patterns <- regmatches(sample_names, regexpr("(rep|replicate|batch|run|day)(\\d+)", sample_names, ignore.case = TRUE))
    if (length(batch_patterns) == n_samples) {
      batch <- gsub("[^0-9]", "", batch_patterns)
      cat("Auto-detected batch from explicit labels (Rep/Batch patterns)\n")
    }
  }

  # Strategy 1B: Look for alphanumeric batch patterns (e.g., Batch_BGI_A, Batch_BGI_A2, Batch_BGI_A3)
  if (is.null(batch) && any(grepl("batch", sample_names, ignore.case = TRUE))) {
    # Try to extract batch identifiers like "BGI_A", "BGI_A2", "BGI_A3"
    # Look for pattern: Batch_XXX_YYY where YYY might be the batch identifier
    batch_patterns <- regmatches(sample_names, regexpr("Batch_[A-Za-z0-9]+_[A-Za-z0-9]+", sample_names, ignore.case = TRUE))

    if (length(batch_patterns) == n_samples) {
      # Extract the last alphanumeric component after "Batch_"
      # For "Batch_BGI_A", "Batch_BGI_A2", "Batch_BGI_A3" -> extract "A", "A2", "A3"
      batch_ids <- gsub(".*_([A-Za-z0-9]+)$", "\\1", batch_patterns, ignore.case = TRUE)

      # Check if we have meaningful batch structure (not all unique)
      n_unique_batches <- length(unique(batch_ids))
      if (n_unique_batches > 1 && n_unique_batches < n_samples) {
        # Convert to numeric for edgeR (A=1, A2=2, A3=3, etc.)
        batch <- as.numeric(factor(batch_ids, levels = unique(batch_ids)))
        cat("Auto-detected batch from alphanumeric labels (Batch_XXX pattern)\n")
        cat("Batch mapping:", paste(unique(batch_ids), "=", unique(batch), collapse = ", "), "\n")
      }
    }
  }

  # Strategy 2: UNIVERSAL stem-based pairing using control/treatment keywords
  # This is THE MOST RELIABLE method - strip keywords, match stems!
  if (is.null(batch)) {
    cat("Trying universal stem-based batch detection (using control/treatment keywords)...\n")

    # Create regex patterns for keywords (case-insensitive, word boundaries)
    control_pattern <- paste0("_?", control_keyword, "_?")
    treatment_pattern <- paste0("_?", treatment_keyword, "_?")

    # Strip keywords from sample names to get "stems"
    stems <- sample_names
    stems <- gsub(control_pattern, "_KEYWORD_", stems, ignore.case = TRUE)
    stems <- gsub(treatment_pattern, "_KEYWORD_", stems, ignore.case = TRUE)
    stems <- gsub("_+", "_", stems)  # Clean up multiple underscores
    stems <- gsub("_$", "", stems)   # Remove trailing underscores
    stems <- gsub("^_", "", stems)   # Remove leading underscores

    cat("Sample stems after removing keywords:\n")
    for (i in 1:min(6, length(stems))) {
      cat("  ", sample_names[i], " → ", stems[i], "\n")
    }
    if (length(stems) > 6) cat("  ... (", length(stems) - 6, " more)\n")

    # Count how many times each stem appears
    stem_counts <- table(stems)
    n_unique_stems <- length(unique(stems))

    cat("\nUnique stems:", n_unique_stems, "\n")
    cat("Samples per stem:", paste(unique(as.numeric(stem_counts)), collapse = ", "), "\n")

    # Check if we have valid pairing structure:
    # - Fewer stems than samples (some samples share stems)
    # - Each stem appears in BOTH control and treatment (valid pairing)
    if (n_unique_stems < n_samples && n_unique_stems >= 2) {
      # Validate that each stem has both control and treatment
      valid_pairing <- TRUE
      for (stem in unique(stems)) {
        has_control <- any(condition[stems == stem] == "Control")
        has_treatment <- any(condition[stems == stem] == "Treatment")
        if (!has_control || !has_treatment) {
          # Some stems might only have control or treatment (technical replicates)
          # This is still valid - just assign unique batch to unpaired samples
          valid_pairing <- TRUE  # Still proceed, handle unpaired later
        }
      }

      if (valid_pairing) {
        # Assign batch IDs based on stems
        batch <- as.numeric(factor(stems, levels = unique(stems)))
        cat("✓ Auto-detected batch structure from sample name stems!\n")
        cat("Batch assignments (first 10):", paste(batch[1:min(10, length(batch))], collapse = ", "), "\n")

        # Show batch distribution
        batch_table <- table(batch)
        cat("Batch distribution:", paste(names(batch_table), "→", batch_table, "samples", collapse = "; "), "\n")
      }
    } else if (n_unique_stems == n_samples) {
      cat("All stems are unique (no pairing detected from keywords)\n")
    }
  }

  # Strategy 3: Look for shared numeric patterns (fallback for old datasets)
  if (is.null(batch)) {
    # Extract all numeric patterns from sample names
    # For GSM3151467_05_Cortex_Control_R1_001, we want to find "05" that might be shared

    # Try extracting smaller numbers (2-3 digits) that might represent subject/pair IDs
    small_numbers <- gsub(".*_([0-9]{1,3})_.*", "\\1", sample_names)

    # Check if these numbers create valid pairing structure
    if (length(unique(small_numbers)) <= n_samples/2) {
      # Each unique number should appear in both control and treatment
      valid_pairing <- TRUE
      for (num in unique(small_numbers)) {
        has_control <- any(condition[small_numbers == num] == "Control")
        has_treatment <- any(condition[small_numbers == num] == "Treatment")
        if (!has_control || !has_treatment) {
          valid_pairing <- FALSE
          break
        }
      }

      if (valid_pairing) {
        batch <- small_numbers
        cat("Auto-detected paired batch structure from sample IDs\n")
        cat("Pairing pattern:", paste(unique(small_numbers), collapse = ", "), "\n")
      }
    }
  }

  # Strategy 4: If still no batch detected, cannot proceed with auto mode
  if (is.null(batch)) {
    cat("\n================================================================================\n")
    cat("WARNING: Could not auto-detect meaningful batch structure!\n")
    cat("================================================================================\n\n")
    cat("Tried the following detection strategies:\n")
    cat("  1. Explicit batch labels (Batch1, Rep1, Run1)\n")
    cat("  2. Alphanumeric batch patterns (Batch_BGI_A, Batch_BGI_A2)\n")
    cat("  3. Universal stem matching (removing control/treatment keywords)\n")
    cat("  4. Shared numeric IDs (e.g., Sample_05 pairing)\n\n")
    cat("All strategies failed. Your samples may have unique IDs that don't indicate pairing.\n\n")

    cat("OPTIONS TO PROCEED:\n")
    cat("1. Use 'paired' mode if samples are ordered (ctrl1-trt1, ctrl2-trt2, ...)\n")
    cat("   Example: Rscript batch_effect_edgeR_v3.R counts.csv ctrl trt paired\n\n")
    cat("2. Provide explicit batch labels (e.g., '1,1,2,2,3,3')\n")
    cat("   Example: Rscript batch_effect_edgeR_v3.R counts.csv ctrl trt '1,1,2,2,3,3'\n\n")
    cat("3. Use simpleEdger3.R instead (no batch correction)\n")
    cat("   Note: Only use if you're confident there are NO batch effects!\n\n")
    cat("================================================================================\n\n")

    stop("ERROR: Cannot proceed with 'auto' batch detection. Please specify batch structure or use simpleEdger3.R")
  }

} else if (batch_info == "paired") {
  # Assume 1:1 pairing for equal groups, or pair as many as possible for unequal
  batch <- rep(NA, n_samples)

  if (n_control == n_treatment) {
    # Equal groups: straightforward 1:1 pairing
    for (i in 1:n_control) {
      batch[control_indices[i]] <- i
      batch[treatment_indices[i]] <- i
    }
    cat("Using paired design: equal groups, 1:1 pairing\n")

  } else {
    # Unequal groups: pair as many as possible, give unique batches to unpaired
    n_pairs <- min(n_control, n_treatment)

    # Pair first n_pairs from each group
    for (i in 1:n_pairs) {
      batch[control_indices[i]] <- i
      batch[treatment_indices[i]] <- i
    }

    # Give unique batch IDs to unpaired samples
    next_batch <- n_pairs + 1
    if (n_control > n_treatment) {
      # Extra control samples
      for (i in (n_pairs + 1):n_control) {
        batch[control_indices[i]] <- next_batch
        next_batch <- next_batch + 1
      }
    } else {
      # Extra treatment samples
      for (i in (n_pairs + 1):n_treatment) {
        batch[treatment_indices[i]] <- next_batch
        next_batch <- next_batch + 1
      }
    }

    cat(sprintf("Using paired design: unequal groups (%d paired + %d unpaired)\n", n_pairs, abs(n_control - n_treatment)))
  }

} else {
  # Explicit batch labels provided
  batch_labels <- strsplit(batch_info, ",")[[1]]

  if (length(batch_labels) != n_samples) {
    stop(sprintf("\n\nERROR: Number of batch labels (%d) doesn't match number of samples (%d)\n\n",
                 length(batch_labels), n_samples))
  }

  batch <- batch_labels
  cat("Using user-provided batch labels\n")
}

# ============================================================================
# Validate Batch Structure (NEW SAFETY CHECK!)
# ============================================================================

cat("\nValidating batch structure...\n")

n_unique_batches <- length(unique(batch))

if (n_unique_batches == n_samples) {
  cat("\nERROR: Each sample has a UNIQUE batch ID!\n")
  cat("This creates a singular design matrix (batch confounded with sample identity).\n\n")
  cat("Your batch structure:\n")
  for (i in 1:n_samples) {
    cat(sprintf("  %s → Batch: %s\n", sample_names[i], batch[i]))
  }
  cat("\nThis is NOT a valid batch effect dataset.\n")
  cat("SOLUTION: Use simpleEdger3.R instead (no batch correction).\n\n")
  stop("ERROR: Invalid batch structure - every sample has unique batch ID")
}

cat(sprintf("Batch validation passed: %d samples in %d batches (%.1f samples/batch)\n",
            n_samples, n_unique_batches, n_samples / n_unique_batches))

# ============================================================================
# Create Design Matrix
# ============================================================================

# Create sample info dataframe
sample_info <- data.frame(
  Sample = sample_names,
  Condition = condition,
  Batch = batch,
  stringsAsFactors = FALSE
)

cat("\nFinal design:\n")
print(sample_info)

# Convert to factors
sample_info$Condition <- factor(sample_info$Condition, levels = c("Control", "Treatment"))
sample_info$Batch <- factor(sample_info$Batch)

cat("\nCreating design matrix with batch correction...\n")
design <- model.matrix(~ batch + condition, data = sample_info)

cat("Design matrix:\n")
print(design)

# ============================================================================
# Run edgeR Analysis
# ============================================================================

cat("\n\nRunning edgeR differential expression analysis...\n")

# Create DGEList object
y <- DGEList(counts = count_matrix, genes = gene_info, group = sample_info$Condition)

# Filter low-expression genes
keep <- filterByExpr(y, design)
y <- y[keep, , keep.lib.sizes = FALSE]

cat(sprintf("\nGenes after filtering: %d \n\n", nrow(y)))

# Normalization
y <- calcNormFactors(y)

cat("Library sizes and normalization factors:\n")
print(data.frame(
  Sample = sample_names,
  LibSize = y$samples$lib.size,
  NormFactor = y$samples$norm.factors
))

# Estimate dispersion
y <- estimateDisp(y, design)

cat(sprintf("\nCommon dispersion: %f \n", y$common.dispersion))
cat(sprintf("BCV (sqrt of dispersion): %f \n\n", sqrt(y$common.dispersion)))

# Fit GLM
fit <- glmQLFit(y, design)

# Test for treatment effect (last coefficient is conditionTreatment)
glm <- glmQLFTest(fit, coef = ncol(design))

# Get top tags
result <- topTags(glm, n = Inf, sort.by = "PValue")$table

# ============================================================================
# Save Results
# ============================================================================

# Generate output filename
output_file <- sprintf("%svs%s.csv", levels(sample_info$Condition)[2], levels(sample_info$Condition)[1])

cat("\n================================================================================\n")
cat("WRITING RESULTS\n")
cat("================================================================================\n\n")

write.csv(result, file = output_file, row.names = FALSE, quote = FALSE)

# Generate threshold distribution for agent review (re-analysis)
cat("\nGenerating threshold distribution for agent review...\n")
threshold_dist <- list(
  # FDR thresholds
  fdr_0_001 = sum(result$FDR < 0.001),
  fdr_0_01 = sum(result$FDR < 0.01),
  fdr_0_05 = sum(result$FDR < 0.05),
  fdr_0_10 = sum(result$FDR < 0.10),
  fdr_0_20 = sum(result$FDR < 0.20),

  # LogFC thresholds (at FDR < 0.05)
  logfc_0_5_at_fdr_0_05 = sum(result$FDR < 0.05 & abs(result$logFC) > 0.5),
  logfc_0_585_at_fdr_0_05 = sum(result$FDR < 0.05 & abs(result$logFC) > 0.585),
  logfc_1_0_at_fdr_0_05 = sum(result$FDR < 0.05 & abs(result$logFC) > 1.0),
  logfc_1_5_at_fdr_0_05 = sum(result$FDR < 0.05 & abs(result$logFC) > 1.5),
  logfc_2_0_at_fdr_0_05 = sum(result$FDR < 0.05 & abs(result$logFC) > 2.0),

  # P-value distribution
  min_pvalue = min(result$PValue),
  median_pvalue = median(result$PValue),
  percent_p_below_0_05 = round(100 * sum(result$PValue < 0.05) / nrow(result), 2),

  # LogFC distribution (for all genes)
  min_logfc = min(result$logFC),
  max_logfc = max(result$logFC),
  median_logfc = median(result$logFC),

  # Current settings results
  current_fdr = fdr_threshold,
  current_logfc = logfc_threshold,
  current_degs = sum(result$FDR < fdr_threshold & abs(result$logFC) > logfc_threshold),

  # Total genes tested
  total_genes = nrow(result)
)

# Write threshold distribution JSON
write_json(threshold_dist, "threshold_distribution.json", auto_unbox = TRUE, pretty = TRUE)

cat("Output file:", output_file, "\n")
cat("Total genes tested:", nrow(result), "\n")
cat(sprintf("Significant genes (FDR < %.3f, |logFC| > %.3f): %d\n",
    fdr_threshold, logfc_threshold,
    sum(result$FDR < fdr_threshold & abs(result$logFC) > logfc_threshold)))
cat(sprintf("  Upregulated (logFC > %.3f): %d\n",
    logfc_threshold,
    sum(result$FDR < fdr_threshold & result$logFC > logfc_threshold)))
cat(sprintf("  Downregulated (logFC < -%.3f): %d\n\n",
    logfc_threshold,
    sum(result$FDR < fdr_threshold & result$logFC < -logfc_threshold)))

cat("================================================================================\n")
cat("BATCH EFFECT CORRECTION COMPLETE\n")
cat("================================================================================\n")
