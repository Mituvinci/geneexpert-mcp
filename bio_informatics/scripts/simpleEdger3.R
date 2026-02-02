#!/usr/bin/Rscript
# EDITED by Claude Code (2026-01-07):
#   - Read CSV input instead of space-delimited text
#   - Write CSV output instead of space-delimited text
# EDITED by Claude Code (2026-01-20):
#   - Support UNEQUAL group sizes (e.g., 3 control + 4 treatment)
#   - Automatically detect groups from column names using keywords
#
# read in .count.filtered.csv file
# generates edger3 for each possible group comparison

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("\n\nUsage: simpleEdger3.R input.count.filtered.csv control_keyword treatment_keyword [fdr_threshold] [logfc_threshold]\n\nThis script compares 2 groups and now SUPPORTS UNEQUAL GROUP SIZES.\n\nExample: simpleEdger3.R counts.csv \"control\" \"treatment\"\nExample: simpleEdger3.R counts.csv \"untreated\" \"Dex\" 0.05 0.585\n\nOptional arguments:\n  fdr_threshold: FDR cutoff (default: 0.05)\n  logfc_threshold: LogFC cutoff (default: 0.585)\n\nThe script will:\n1. Read column names from CSV\n2. Assign samples to groups based on keywords (case-insensitive)\n3. Handle any combination (2+3, 3+4, 4+5, etc.)\n\nPlease make sure the counts are formatted:\nID Length sample1 sample2 sample3...\n\n", call. = FALSE)
}

# Parse custom thresholds (for re-analysis)
fdr_threshold <- if (length(args) >= 4) as.numeric(args[4]) else 0.05
logfc_threshold <- if (length(args) >= 5) as.numeric(args[5]) else 0.585

cat(sprintf("\n=== Threshold Settings ===\n"))
cat(sprintf("FDR threshold: %.3f\n", fdr_threshold))
cat(sprintf("LogFC threshold: %.3f\n\n", logfc_threshold))

library("Rsubread")
library("edgeR")
library("limma")
library("jsonlite")  # For write_json() function

# EDIT: Read CSV instead of space-delimited
counts <- read.csv(args[1], header = TRUE)
control_keyword <- tolower(args[2])
treatment_keyword <- tolower(args[3])

# Get sample column names (skip ID and Length columns)
sample_cols <- colnames(counts[, -1:-2])

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
group_labels <- sapply(sample_cols, assign_group)

# Check for unknown samples
if (any(group_labels == "Unknown")) {
  unknown_samples <- sample_cols[group_labels == "Unknown"]
  stop(sprintf("\n\nERROR: Could not assign groups for these samples:\n%s\n\nMake sure sample names contain keywords '%s' or '%s'\n",
               paste(unknown_samples, collapse = "\n"), control_keyword, treatment_keyword))
}

# Count samples per group
n_control <- sum(group_labels == "Control")
n_treatment <- sum(group_labels == "Treatment")

cat("\n=== GROUP ASSIGNMENT ===\n")
cat(sprintf("Control group (%d samples): %s\n", n_control,
            paste(sample_cols[group_labels == "Control"], collapse = ", ")))
cat(sprintf("Treatment group (%d samples): %s\n", n_treatment,
            paste(sample_cols[group_labels == "Treatment"], collapse = ", ")))
cat("========================\n\n")

# Create factor for EdgeR (edgeR handles unequal groups automatically!)
grp <- factor(group_labels, levels = c("Control", "Treatment"))

y <- DGEList(counts = counts[, -1:-2], genes = counts[, 1:2], group = grp)

design <- model.matrix(~0 + grp)
colnames(design) <- levels(grp)

cat("Estimating dispersions...\n")
y <- estimateDisp(y, design)

cat("Fitting GLM model...\n")
fit <- glmQLFit(y, design)

# Save the names of the groups
groups <- levels(grp)  # Should be c("Control", "Treatment")

# Generate filename
fname <- paste(groups[1], "vs", groups[2], ".csv", sep = "")

cat("\nGenerating DE results:", fname, "\n")

# Generate contrast vector: Treatment vs Control (Treatment - Control)
# This means positive logFC = upregulated in Treatment, negative = downregulated in Treatment
my_contrast <- makeContrasts(
  TreatmentVsControl = Treatment - Control,
  levels = design
)

# Generate DE table
cat("Running QL F-test...\n")
glm <- glmQLFTest(fit, contrast = my_contrast)
table <- topTags(glm, n = length(y$counts[, 1]))$table

# Write CSV output
write.csv(table, file = fname, row.names = FALSE, quote = FALSE)

# Generate threshold distribution for agent review (re-analysis)
cat("\nGenerating threshold distribution for agent review...\n")
threshold_dist <- list(
  # FDR thresholds
  fdr_0_001 = sum(table$FDR < 0.001),
  fdr_0_01 = sum(table$FDR < 0.01),
  fdr_0_05 = sum(table$FDR < 0.05),
  fdr_0_10 = sum(table$FDR < 0.10),
  fdr_0_20 = sum(table$FDR < 0.20),

  # LogFC thresholds (at FDR < 0.05)
  logfc_0_5_at_fdr_0_05 = sum(table$FDR < 0.05 & abs(table$logFC) > 0.5),
  logfc_0_585_at_fdr_0_05 = sum(table$FDR < 0.05 & abs(table$logFC) > 0.585),
  logfc_1_0_at_fdr_0_05 = sum(table$FDR < 0.05 & abs(table$logFC) > 1.0),
  logfc_1_5_at_fdr_0_05 = sum(table$FDR < 0.05 & abs(table$logFC) > 1.5),
  logfc_2_0_at_fdr_0_05 = sum(table$FDR < 0.05 & abs(table$logFC) > 2.0),

  # P-value distribution
  min_pvalue = min(table$PValue),
  median_pvalue = median(table$PValue),
  percent_p_below_0_05 = round(100 * sum(table$PValue < 0.05) / nrow(table), 2),

  # LogFC distribution (for all genes)
  min_logfc = min(table$logFC),
  max_logfc = max(table$logFC),
  median_logfc = median(table$logFC),

  # Current settings results
  current_fdr = fdr_threshold,
  current_logfc = logfc_threshold,
  current_degs = sum(table$FDR < fdr_threshold & abs(table$logFC) > logfc_threshold),

  # Total genes tested
  total_genes = nrow(table)
)

# Write threshold distribution JSON
write_json(threshold_dist, "threshold_distribution.json", auto_unbox = TRUE, pretty = TRUE)

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Output written to:", fname, "\n")
cat("Total genes tested:", nrow(table), "\n")
cat(sprintf("Significant genes (FDR < %.3f, |logFC| > %.3f): %d\n",
    fdr_threshold, logfc_threshold,
    sum(table$FDR < fdr_threshold & abs(table$logFC) > logfc_threshold)))
cat(sprintf("  Upregulated (logFC > %.3f): %d\n",
    logfc_threshold,
    sum(table$FDR < fdr_threshold & table$logFC > logfc_threshold)))
cat(sprintf("  Downregulated (logFC < -%.3f): %d\n",
    logfc_threshold,
    sum(table$FDR < fdr_threshold & table$logFC < -logfc_threshold)))
cat("========================\n")
