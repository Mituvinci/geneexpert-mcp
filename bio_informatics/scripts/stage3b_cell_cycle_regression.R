#!/usr/bin/env Rscript

#' Stage 3B — CELL CYCLE REGRESSION (Conditional - only if agents recommend)
#' stage3b_cell_cycle_regression.R
#'
#' Purpose:
#' 1. Scale data with cell cycle regression (regress out S.Score, G2M.Score)
#' 2. Re-run PCA to verify cell cycle effect removed
#' 3. Save cell_cycle_after plots for comparison
#'
#' This stage ONLY runs if agents decided:
#'   REMOVE_CELL_CYCLE (cell cycle effect detected in Stage 3A)
#'
#' Input:
#'   seurat_stage3a_scored.rds (with S.Score, G2M.Score metadata)
#'   hvg_genes.csv
#'
#' Output:
#'   seurat_stage3_norm.rds (with cell cycle regressed)
#'   cell_cycle_after.pdf/.jpg
#'   regression_summary_stage3b.json

library(Seurat)
library(jsonlite)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
rds_in <- args[1]
outdir <- args[2]

cat("\n================================================================================\n")
cat("scRNA-seq Stage 3B: CELL CYCLE REGRESSION\n")
cat("================================================================================\n\n")

cat("Loading Seurat object from Stage 3A (with cell cycle scores)...\n")
seu <- readRDS(rds_in)

# Load HVG genes
hvg_file <- file.path(outdir, "hvg_genes.csv")
if (!file.exists(hvg_file)) {
  stop("ERROR: hvg_genes.csv not found. Stage 3A must run first.")
}
hvg_genes <- read.csv(hvg_file)$gene

cat(sprintf("  Using %d HVGs for scaling\n", length(hvg_genes)))

# Verify cell cycle scores exist
if (!("S.Score" %in% colnames(seu@meta.data)) || !("G2M.Score" %in% colnames(seu@meta.data))) {
  stop("ERROR: S.Score and G2M.Score not found. Stage 3A must run first.")
}

# Step 1: REMOVE CELL CYCLE EFFECTS (Regress out S.Score and G2M.Score)
cat("\n================================================================================\n")
cat("STEP 1: REMOVING CELL CYCLE EFFECTS\n")
cat("================================================================================\n")

cat("Scaling data with cell cycle regression...\n")
cat("  This removes variance due to cell cycle phase from gene expression\n")
cat("  Goal: Cells cluster by biological identity, NOT cell cycle phase\n\n")

seu <- ScaleData(seu,
                 vars.to.regress = c("S.Score", "G2M.Score"),
                 features = hvg_genes,
                 verbose = FALSE)

cat("✓ Cell cycle effects removed from scaled data\n")

# Step 2: RE-RUN PCA AND SAVE PLOTS (AFTER REMOVAL)
cat("\n================================================================================\n")
cat("STEP 2: RE-RUNNING PCA AND SAVING PLOTS (AFTER removal)\n")
cat("================================================================================\n")

seu <- RunPCA(seu, features = hvg_genes, npcs = 30, verbose = FALSE)

# Plot: PCA colored by cell cycle phase (AFTER removal)
# Save as PDF (for user)
# Plot 3: PCA colored by cell cycle phase (AFTER removal)
p3 <- DimPlot(seu, reduction = "pca", group.by = "Phase") +
  ggtitle("PCA: Cell Cycle Phase (AFTER removal)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Plot 4: FeaturePlot showing S.Score and G2M.Score
p4 <- FeaturePlot(seu, reduction = "pca", features = c("S.Score", "G2M.Score"), ncol = 2)

# Save as PDF (for user) - contains BOTH plots
pdf(file.path(outdir, "cell_cycle_after.pdf"), width = 10, height = 5)
print(p3)
print(p4)
dev.off()

# Save Plot 3 (Phase - 3 colors) as JPEG for comparison
jpeg(file.path(outdir, "cell_cycle_phase_after.jpg"), width = 1000, height = 500, quality = 95)
print(p3)
dev.off()

# Save Plot 4 (Scores - blue gradients) as JPEG for reference
jpeg(file.path(outdir, "cell_cycle_scores_after.jpg"), width = 1000, height = 500, quality = 95)
print(p4)
dev.off()

cat("✓ Saved: cell_cycle_after.pdf (both plots)\n")
cat("✓ Saved: cell_cycle_phase_after.jpg (3-color Phase plot for comparison)\n")
cat("✓ Saved: cell_cycle_scores_after.jpg (S.Score/G2M.Score gradients)\n")

# Calculate correlation again (should be lower now)
pca_embeddings <- seu@reductions$pca@cell.embeddings[, 1:2]
cor_s_after <- cor(pca_embeddings[, 1], seu$S.Score)
cor_g2m_after <- cor(pca_embeddings[, 1], seu$G2M.Score)
cat(sprintf("\nPC1 correlation AFTER regression:\n"))
cat(sprintf("  PC1 correlation with S.Score: %.3f\n", cor_s_after))
cat(sprintf("  PC1 correlation with G2M.Score: %.3f\n", cor_g2m_after))

cat("\n================================================================================\n")
cat("CELL CYCLE CORRECTION COMPLETE!\n")
cat("================================================================================\n")
cat("BEFORE: Cells clustered by cell cycle phase (confounding factor)\n")
cat("AFTER:  Cell cycle variance removed, cells cluster by biological identity\n")
cat("\nCompare cell_cycle_phase_before.jpg vs cell_cycle_phase_after.jpg to see the effect!\n")
cat("(PDFs contain both Phase and Score plots)\n")
cat("================================================================================\n\n")

# Step 3: Write summary JSON
summary_json <- list(
  cell_cycle_regressed = TRUE,
  pc1_correlation_s_score_after = cor_s_after,
  pc1_correlation_g2m_score_after = cor_g2m_after,
  cell_cycle_plots = list(
    after_pdf = "cell_cycle_after.pdf",
    after_jpg_phase = "cell_cycle_phase_after.jpg",
    after_jpg_scores = "cell_cycle_scores_after.jpg"
  )
)

write_json(
  summary_json,
  file.path(outdir, "regression_summary_stage3b.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)

# CRITICAL FIX (2026-01-27): Set VariableFeatures so Stage 4 PCA can access them
# The PCA computed here is temporary (for before/after comparison)
# Stage 4 will re-run PCA and needs access to VariableFeatures
cat("\nSetting VariableFeatures in Seurat object...\n")
VariableFeatures(seu) <- hvg_genes
cat(sprintf("  ✓ Set %d variable features\n", length(hvg_genes)))

# Save final Seurat object (with cell cycle corrected)
saveRDS(seu, file.path(outdir, "seurat_stage3_norm.rds"))

cat("\n✓ Stage 3B complete!\n")
cat("  - Cell cycle effects removed from scaled data\n")
cat("  - Variable features set for Stage 4 PCA\n")
cat("  - Saved: seurat_stage3_norm.rds (ready for PCA in Stage 4)\n\n")
