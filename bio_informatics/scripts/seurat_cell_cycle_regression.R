#!/usr/bin/Rscript
#
# seurat_cell_cycle_regression.R - Stage 3b: Cell Cycle Regression + Verification
#
# Usage:
#   Rscript seurat_cell_cycle_regression.R <seurat_scored.rds> <output_prefix>
#
# Arguments:
#   seurat_scored.rds - Seurat object from Stage 3a (with cell cycle scores)
#   output_prefix     - Prefix for output files
#
# Outputs:
#   <prefix>_seurat_regressed.rds        - Seurat object with regressed data
#   <prefix>_cell_cycle_after.pdf        - DimPlot (PCA) after regression
#   <prefix>_cell_cycle_after.jpg        - DimPlot for agents
#   <prefix>_cell_cycle_comparison.pdf   - Side-by-side before/after
#   <prefix>_cell_cycle_comparison.jpg   - Comparison for agents
#   <prefix>_regression_summary.txt      - Regression metrics
#

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("\n================================================================================
STAGE 3b: CELL CYCLE REGRESSION
================================================================================

Usage: seurat_cell_cycle_regression.R <seurat_scored.rds> <output_prefix>

Arguments:
  seurat_scored.rds - Seurat object from Stage 3a (with cell cycle scores)
  output_prefix     - Prefix for output files

Example:
  Rscript seurat_cell_cycle_regression.R my_experiment_seurat_scored.rds my_experiment

================================================================================
", call. = FALSE)
}

library(Seurat)
library(ggplot2)
library(patchwork)

rds_file <- args[1]
output_prefix <- args[2]

cat("\n================================================================================\n")
cat("STAGE 3b: CELL CYCLE REGRESSION\n")
cat("================================================================================\n\n")
cat("Input RDS:", rds_file, "\n")
cat("Output prefix:", output_prefix, "\n\n")

# ============================================================================
# Load Seurat Object from Stage 3a
# ============================================================================

cat("Loading Seurat object from Stage 3a...\n")
seurat_obj <- readRDS(rds_file)
cat("Loaded:", ncol(seurat_obj), "cells\n\n")

# Store "before" PCA embeddings and correlations
pc_before <- Embeddings(seurat_obj, reduction = "pca")
cor_s_pc1_before <- cor(seurat_obj$S.Score, pc_before[, 1])
cor_g2m_pc1_before <- cor(seurat_obj$G2M.Score, pc_before[, 1])
cor_s_pc2_before <- cor(seurat_obj$S.Score, pc_before[, 2])
cor_g2m_pc2_before <- cor(seurat_obj$G2M.Score, pc_before[, 2])

cat("Cell Cycle Correlations (BEFORE regression):\n")
cat(sprintf("  PC1 vs S.Score: r=%.3f\n", cor_s_pc1_before))
cat(sprintf("  PC1 vs G2M.Score: r=%.3f\n", cor_g2m_pc1_before))
cat(sprintf("  PC2 vs S.Score: r=%.3f\n", cor_s_pc2_before))
cat(sprintf("  PC2 vs G2M.Score: r=%.3f\n", cor_g2m_pc2_before))

# Store "before" plot
p_before <- DimPlot(seurat_obj, reduction = "pca", group.by = "Phase", pt.size = 1.2) +
  ggtitle("Before Regression") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))

# ============================================================================
# Scale Data with Cell Cycle Regression
# ============================================================================

cat("\nPerforming cell cycle regression...\n")
cat("Regressing out: S.Score, G2M.Score\n")

seurat_obj <- ScaleData(
  seurat_obj,
  vars.to.regress = c("S.Score", "G2M.Score"),
  features = rownames(seurat_obj)
)

cat("Regression complete.\n")

# ============================================================================
# Re-run PCA (after regression)
# ============================================================================

cat("\nRe-running PCA on regressed data...\n")
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))

# Get variance explained
pca_var <- (seurat_obj@reductions$pca@stdev)^2
pca_var_pct <- 100 * pca_var / sum(pca_var)

cat("\nPCA Variance Explained (AFTER regression):\n")
for (i in 1:5) {
  cat(sprintf("  PC%d: %.2f%%\n", i, pca_var_pct[i]))
}

# Calculate correlation after regression
pc_after <- Embeddings(seurat_obj, reduction = "pca")
cor_s_pc1_after <- cor(seurat_obj$S.Score, pc_after[, 1])
cor_g2m_pc1_after <- cor(seurat_obj$G2M.Score, pc_after[, 1])
cor_s_pc2_after <- cor(seurat_obj$S.Score, pc_after[, 2])
cor_g2m_pc2_after <- cor(seurat_obj$G2M.Score, pc_after[, 2])

cat("\nCell Cycle Correlations (AFTER regression):\n")
cat(sprintf("  PC1 vs S.Score: r=%.3f\n", cor_s_pc1_after))
cat(sprintf("  PC1 vs G2M.Score: r=%.3f\n", cor_g2m_pc1_after))
cat(sprintf("  PC2 vs S.Score: r=%.3f\n", cor_s_pc2_after))
cat(sprintf("  PC2 vs G2M.Score: r=%.3f\n", cor_g2m_pc2_after))

# ============================================================================
# Calculate Regression Success Metrics
# ============================================================================

cat("\n================================================================================\n")
cat("REGRESSION SUCCESS METRICS\n")
cat("================================================================================\n")

# Absolute correlation reduction
reduction_s_pc1 <- abs(cor_s_pc1_before) - abs(cor_s_pc1_after)
reduction_g2m_pc1 <- abs(cor_g2m_pc1_before) - abs(cor_g2m_pc1_after)
reduction_s_pc2 <- abs(cor_s_pc2_before) - abs(cor_s_pc2_after)
reduction_g2m_pc2 <- abs(cor_g2m_pc2_before) - abs(cor_g2m_pc2_after)

cat("\nCorrelation Reduction:\n")
cat(sprintf("  PC1 vs S.Score: %.3f → %.3f (reduced by %.3f)\n", cor_s_pc1_before, cor_s_pc1_after, reduction_s_pc1))
cat(sprintf("  PC1 vs G2M.Score: %.3f → %.3f (reduced by %.3f)\n", cor_g2m_pc1_before, cor_g2m_pc1_after, reduction_g2m_pc1))
cat(sprintf("  PC2 vs S.Score: %.3f → %.3f (reduced by %.3f)\n", cor_s_pc2_before, cor_s_pc2_after, reduction_s_pc2))
cat(sprintf("  PC2 vs G2M.Score: %.3f → %.3f (reduced by %.3f)\n", cor_g2m_pc2_before, cor_g2m_pc2_after, reduction_g2m_pc2))

# Assess success
success_threshold <- 0.3  # Correlation should be reduced to <0.3
max_cor_after <- max(abs(cor_s_pc1_after), abs(cor_g2m_pc1_after), abs(cor_s_pc2_after), abs(cor_g2m_pc2_after))

if (max_cor_after < success_threshold) {
  cat("\n✅ REGRESSION SUCCESS: Cell cycle effect successfully removed\n")
  cat(sprintf("   (Max correlation after regression: %.3f < %.1f threshold)\n", max_cor_after, success_threshold))
  regression_status <- "SUCCESS"
} else {
  cat("\n⚠ REGRESSION PARTIAL: Cell cycle effect reduced but still present\n")
  cat(sprintf("   (Max correlation after regression: %.3f ≥ %.1f threshold)\n", max_cor_after, success_threshold))
  regression_status <- "PARTIAL"
}

cat("================================================================================\n\n")

# ============================================================================
# Generate DimPlots (After Regression)
# ============================================================================

cat("Generating DimPlots (after regression)...\n")

# Plot after regression
p_after <- DimPlot(seurat_obj, reduction = "pca", group.by = "Phase", pt.size = 1.2) +
  ggtitle("After Regression") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))

# Save individual "after" plot
pdf_after <- paste0(output_prefix, "_cell_cycle_after.pdf")
jpg_after <- paste0(output_prefix, "_cell_cycle_after.jpg")
ggsave(pdf_after, p_after, width = 8, height = 6)
ggsave(jpg_after, p_after, width = 8, height = 6, dpi = 150)

cat("  PDF (after):", pdf_after, "\n")
cat("  JPG (after):", jpg_after, "\n")

# Side-by-side comparison
p_comparison <- p_before | p_after
p_comparison <- p_comparison +
  plot_annotation(
    title = "Cell Cycle Regression: Before vs After",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  )

pdf_comparison <- paste0(output_prefix, "_cell_cycle_comparison.pdf")
jpg_comparison <- paste0(output_prefix, "_cell_cycle_comparison.jpg")
ggsave(pdf_comparison, p_comparison, width = 14, height = 6)
ggsave(jpg_comparison, p_comparison, width = 14, height = 6, dpi = 150)

cat("  PDF (comparison):", pdf_comparison, "\n")
cat("  JPG (comparison):", jpg_comparison, "\n")

# ============================================================================
# Save Regressed Seurat Object + Summary
# ============================================================================

cat("\nSaving regressed Seurat object...\n")
rds_regressed <- paste0(output_prefix, "_seurat_regressed.rds")
saveRDS(seurat_obj, file = rds_regressed)
cat("  RDS:", rds_regressed, "\n")

# Write summary file
summary_file <- paste0(output_prefix, "_regression_summary.txt")
sink(summary_file)
cat("================================================================================\n")
cat("STAGE 3b: CELL CYCLE REGRESSION SUMMARY\n")
cat("================================================================================\n\n")
cat("Input RDS:", rds_file, "\n")
cat("Output prefix:", output_prefix, "\n")
cat("Cells processed:", ncol(seurat_obj), "\n\n")

cat("================================================================================\n")
cat("REGRESSION RESULTS\n")
cat("================================================================================\n\n")

cat("BEFORE Regression:\n")
cat(sprintf("  PC1 vs S.Score: r=%.3f\n", cor_s_pc1_before))
cat(sprintf("  PC1 vs G2M.Score: r=%.3f\n", cor_g2m_pc1_before))
cat(sprintf("  PC2 vs S.Score: r=%.3f\n", cor_s_pc2_before))
cat(sprintf("  PC2 vs G2M.Score: r=%.3f\n", cor_g2m_pc2_before))

cat("\nAFTER Regression:\n")
cat(sprintf("  PC1 vs S.Score: r=%.3f (reduced by %.3f)\n", cor_s_pc1_after, reduction_s_pc1))
cat(sprintf("  PC1 vs G2M.Score: r=%.3f (reduced by %.3f)\n", cor_g2m_pc1_after, reduction_g2m_pc1))
cat(sprintf("  PC2 vs S.Score: r=%.3f (reduced by %.3f)\n", cor_s_pc2_after, reduction_s_pc2))
cat(sprintf("  PC2 vs G2M.Score: r=%.3f (reduced by %.3f)\n", cor_g2m_pc2_after, reduction_g2m_pc2))

cat("\nRegression Status:", regression_status, "\n")
cat("Max correlation after regression:", round(max_cor_after, 3), "\n")

cat("\nOutput Files:\n")
cat("  Regressed Seurat object:", rds_regressed, "\n")
cat("  DimPlot After (PDF):", pdf_after, "\n")
cat("  DimPlot After (JPG):", jpg_after, "\n")
cat("  Comparison (PDF):", pdf_comparison, "\n")
cat("  Comparison (JPG):", jpg_comparison, "\n")
cat("  Summary:", summary_file, "\n")

cat("\n================================================================================\n")
cat("STAGE 3b COMPLETE - Proceeding to Stage 3c (Clustering + UMAP)\n")
cat("================================================================================\n")
sink()

cat("\nSummary written to:", summary_file, "\n")

cat("\n================================================================================\n")
cat("STAGE 3b COMPLETE\n")
cat("================================================================================\n")
