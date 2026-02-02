#!/usr/bin/env Rscript

#' Stage 3B (Skip) — SCALE WITHOUT CELL CYCLE REGRESSION
#' stage3b_skip_cell_cycle.R
#'
#' Purpose:
#' Scale data WITHOUT cell cycle regression (agents determined no effect)
#'
#' This stage runs if agents decided:
#'   SKIP_CELL_CYCLE (no visible cell cycle effect in Stage 3A)
#'
#' Input:
#'   seurat_stage3a_scored.rds (with S.Score, G2M.Score metadata)
#'   hvg_genes.csv
#'
#' Output:
#'   seurat_stage3_norm.rds (scaled but NOT regressed)
#'   skip_summary_stage3b.json

library(Seurat)
library(jsonlite)

args <- commandArgs(trailingOnly = TRUE)
rds_in <- args[1]
outdir <- args[2]

cat("\n================================================================================\n")
cat("scRNA-seq Stage 3B (Skip): SCALING WITHOUT CELL CYCLE REGRESSION\n")
cat("================================================================================\n\n")

cat("Loading Seurat object from Stage 3A...\n")
seu <- readRDS(rds_in)

# Load HVG genes
hvg_file <- file.path(outdir, "hvg_genes.csv")
if (!file.exists(hvg_file)) {
  stop("ERROR: hvg_genes.csv not found. Stage 3A must run first.")
}
hvg_genes <- read.csv(hvg_file)$gene

cat(sprintf("  Using %d HVGs for scaling\n", length(hvg_genes)))

cat("\nAgents determined: NO significant cell cycle effect detected\n")
cat("Proceeding with standard scaling (no regression)\n\n")

# Scale data WITHOUT cell cycle regression
cat("Scaling data (without cell cycle regression)...\n")
seu <- ScaleData(seu, features = hvg_genes, verbose = FALSE)

cat("✓ Data scaled successfully\n")

# Step 2: Write summary JSON
summary_json <- list(
  cell_cycle_regressed = FALSE,
  reason = "Agents determined no significant cell cycle effect",
  action = "Scaled data without regression"
)

write_json(
  summary_json,
  file.path(outdir, "skip_summary_stage3b.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)

# CRITICAL FIX (2026-01-27): Set VariableFeatures so Stage 4 PCA can access them
# Without this, RunPCA() in Stage 4 fails with "No variable features" error
cat("\nSetting VariableFeatures in Seurat object...\n")
VariableFeatures(seu) <- hvg_genes
cat(sprintf("  ✓ Set %d variable features\n", length(hvg_genes)))

# Save Seurat object (scaled but NOT regressed)
saveRDS(seu, file.path(outdir, "seurat_stage3_norm.rds"))

cat("\n✓ Stage 3B (skip) complete!\n")
cat("  - Data scaled without cell cycle regression\n")
cat("  - Variable features set for Stage 4 PCA\n")
cat("  - Saved: seurat_stage3_norm.rds (ready for PCA in Stage 4)\n\n")
