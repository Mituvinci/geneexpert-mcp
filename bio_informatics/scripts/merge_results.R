#!/usr/bin/Rscript
# UPDATED by Claude Code (2026-01-19)
# Merges RPKM+gene symbols with edgeR results into final Excel file
# Includes exact formulas matching lab workflow:
#   - Expr: checks if gene is expressed (max group average logRPKM > 2 AND logCPM > 0)
#   - CG: classifies as failed2UpRegulate/failed2DownRegulate/nchg
#
# Usage: merge_results.R rpkm.entrz.csv edgeR_results.csv output_prefix control_keyword treatment_keyword

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
  stop("\nUsage: merge_results.R rpkm.entrz.csv edgeR_results.csv output_prefix control_keyword treatment_keyword\n\nExample:\n  merge_results.R DA0036.rpkm.entrz.csv DA0036_DE.csv final_results cont ips\n")
}

library("openxlsx")

# Parse arguments
rpkm_file <- args[1]
edger_file <- args[2]
output_prefix <- args[3]
control_keyword <- args[4]
treatment_keyword <- args[5]

# Read input files
cat("Reading RPKM data:", rpkm_file, "\n")
rpkm <- read.csv(rpkm_file, header = TRUE)

cat("Reading edgeR results:", edger_file, "\n")
edger <- read.csv(edger_file, header = TRUE)

# Merge on ID column (Entrzid in rpkm, ID in edger)
cat("Merging datasets on ID column...\n")
merged <- merge(rpkm, edger, by.x = "Entrzid", by.y = "ID", all = FALSE)

# Identify sample columns (all columns from rpkm except Entrzid and gene_symbol)
rpkm_sample_cols <- setdiff(colnames(rpkm), c("Entrzid", "gene_symbol"))

# Split samples into control and treatment groups based on keywords
control_samples <- grep(control_keyword, rpkm_sample_cols, value = TRUE, ignore.case = TRUE)
treatment_samples <- grep(treatment_keyword, rpkm_sample_cols, value = TRUE, ignore.case = TRUE)

cat("Identified", length(control_samples), "control samples (keyword:", control_keyword, ")\n")
cat("Identified", length(treatment_samples), "treatment samples (keyword:", treatment_keyword, ")\n")

# Build final column order: ID, Length, gene_symbol, RPKM samples, edgeR stats
edger_stat_cols <- c("logFC", "logCPM", "F", "PValue", "FDR")
final_cols <- c("Entrzid", "Length", "gene_symbol", rpkm_sample_cols, edger_stat_cols)
merged_ordered <- merged[, final_cols]

cat("Merged data dimensions:", nrow(merged_ordered), "genes x", ncol(merged_ordered), "columns\n")

# Create Excel workbook
cat("Creating Excel workbook with formulas...\n")
wb <- createWorkbook()
addWorksheet(wb, "Results")

# Write merged data
writeData(wb, "Results", merged_ordered, startRow = 1, startCol = 1)

# Calculate column positions
n_samples <- length(rpkm_sample_cols)
first_sample_col <- 4  # After ID, Length, gene_symbol
logfc_col <- 3 + n_samples + 1  # After ID, Length, gene_symbol, samples
logcpm_col <- logfc_col + 1
f_col <- logcpm_col + 1
pval_col <- f_col + 1
fdr_col <- pval_col + 1
expr_col <- fdr_col + 1
cg_col <- expr_col + 1

# Add Expr and CG column headers
writeData(wb, "Results", "Expr", startCol = expr_col, startRow = 1)
writeData(wb, "Results", "CG", startCol = cg_col, startRow = 1)

# Find column indices for control and treatment samples
control_col_indices <- which(colnames(merged_ordered) %in% control_samples)
treatment_col_indices <- which(colnames(merged_ordered) %in% treatment_samples)

# FIXED 2026-01-27: Calculate threshold column position FIRST (before formulas)
# This ensures formulas reference the correct dynamic columns
threshold_col_x <- max(24, cg_col + 4)  # Start at column X (24) or 4 columns after CG, whichever is further right
threshold_col_y <- threshold_col_x + 1  # Column Y (FDR value)
threshold_col_aa <- threshold_col_x + 3  # Column AA (logRPKM threshold)

# Build Expr formula: =AND(MAX(AVERAGE(control), AVERAGE(treatment)) > $AA$2, logCPM > $Y$2)
# AA2 = logRPKM threshold (2) - dynamically positioned
# Y2 = logCPM threshold (0) - dynamically positioned
n_rows <- nrow(merged_ordered)
expr_formula <- sapply(2:(n_rows + 1), function(r) {
  # Build AVERAGE ranges for control and treatment
  control_range <- paste0(int2col(control_col_indices[1]), r, ":", int2col(control_col_indices[length(control_col_indices)]), r)
  treatment_range <- paste0(int2col(treatment_col_indices[1]), r, ":", int2col(treatment_col_indices[length(treatment_col_indices)]), r)

  paste0("=AND(MAX(AVERAGE(", control_range, "),AVERAGE(", treatment_range, "))>$", int2col(threshold_col_aa), "$2,",
         int2col(logcpm_col), r, ">$", int2col(threshold_col_y), "$2)")
})
writeFormula(wb, "Results", expr_formula, startCol = expr_col, startRow = 2)

# Build CG formula: =IF(AND(logFC > $Y$3, FDR < $Y$1, Expr), "failed2DownRegulate", IF(AND(logFC < -$Y$3, FDR < $Y$1, Expr), "failed2UpRegulate", "nchg"))
# Y1 = FDR threshold (0.05) - dynamically positioned
# Y3 = logFC threshold (0.585) - dynamically positioned
cg_formula <- sapply(2:(n_rows + 1), function(r) {
  paste0("=IF(AND(", int2col(logfc_col), r, ">$", int2col(threshold_col_y), "$3,",
         int2col(fdr_col), r, "<$", int2col(threshold_col_y), "$1,",
         int2col(expr_col), r, "),\"failed2DownRegulate\",IF(AND(",
         int2col(logfc_col), r, "<-$", int2col(threshold_col_y), "$3,",
         int2col(fdr_col), r, "<$", int2col(threshold_col_y), "$1,",
         int2col(expr_col), r, "),\"failed2UpRegulate\",\"nchg\"))")
})
writeFormula(wb, "Results", cg_formula, startCol = cg_col, startRow = 2)

# Add threshold section (columns X-AA, rows 1-3)
# X1 = "FDR",      Y1 = 0.05
# X2 = "logCPM",   Y2 = 0
# X3 = "logFC",    Y3 = 0.585
# Z2 = "logRPKM",  AA2 = 2
# NOTE: threshold_col_x already calculated above with dynamic positioning

writeData(wb, "Results", data.frame(
  Label = c("FDR", "logCPM", "logFC"),
  Value = c(0.05, 0, 0.585)
), startCol = threshold_col_x, startRow = 1, colNames = FALSE)

# Add logRPKM threshold (Z2 = "logRPKM", AA2 = 2)
writeData(wb, "Results", "logRPKM", startCol = threshold_col_x + 2, startRow = 2)  # Z2
writeData(wb, "Results", 2, startCol = threshold_col_x + 3, startRow = 2)  # AA2

# Add summary count labels and formulas (rows 5-6)
writeData(wb, "Results", "failed2DownRegulate", startCol = threshold_col_x, startRow = 5)
writeData(wb, "Results", "failed2UpRegulate", startCol = threshold_col_x, startRow = 6)

down_formula <- paste0("=COUNTIFS(", int2col(cg_col), ":", int2col(cg_col), ",\"failed2DownRegulate\")")
up_formula <- paste0("=COUNTIFS(", int2col(cg_col), ":", int2col(cg_col), ",\"failed2UpRegulate\")")
writeFormula(wb, "Results", c(down_formula, up_formula), startCol = threshold_col_x + 1, startRow = 5)

# Save workbook
output_file <- paste0(output_prefix, "_final.xlsx")
saveWorkbook(wb, output_file, overwrite = TRUE)

cat("\n✓ Merge complete!\n")
cat("  Output file:", output_file, "\n")
cat("  Total genes:", nrow(merged_ordered), "\n")
cat("  Columns:\n")
cat("    - Entrzid, Length, gene_symbol\n")
cat("    - RPKM samples:", n_samples, "\n")
cat("    - edgeR stats: logFC, logCPM, F, PValue, FDR\n")
cat("    - Formulas: Expr, CG\n")
cat("    - Thresholds: FDR=0.05, logCPM=0, logFC=0.585, logRPKM=2\n")
cat("\nOpen Excel to see formulas and counts!\n")

# Helper function to convert column number to Excel column letter
int2col <- function(n) {
  result <- ""
  while (n > 0) {
    remainder <- (n - 1) %% 26
    result <- paste0(LETTERS[remainder + 1], result)
    n <- (n - remainder - 1) %/% 26
  }
  return(result)
}
