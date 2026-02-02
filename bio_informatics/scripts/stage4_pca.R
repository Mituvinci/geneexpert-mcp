#!/usr/bin/env Rscript

#' Stage 4 — PCA + Elbow (Agent Gate)
#' stage4_pca.R
#'
#' Purpose:
#' Compute PCA and summarize variance.
#'
#' Agent decision:
#'   Select number of PCs (e.g. 1:20)
#'
#' Input:
#'   seurat_stage3_norm.rds
#'
#' Output:
#'   seurat_stage4_pca.rds
#'   pca_variance.json
#'   elbow_plot.pdf/.jpg (PCA variance elbow plot)

library(Seurat)
library(jsonlite)

args <- commandArgs(trailingOnly = TRUE)
rds_in <- args[1]
outdir <- args[2]

cat("\n================================================================================\n")
cat("scRNA-seq Stage 4: PCA + Elbow Plot\n")
cat("================================================================================\n\n")

cat("Loading Seurat object from Stage 3B...\n")
seu <- readRDS(rds_in)

# SAFETY CHECK (2026-01-27): Ensure VariableFeatures are set
hvg_count <- length(VariableFeatures(seu))
cat(sprintf("  Variable features in object: %d\n", hvg_count))

if (hvg_count == 0) {
  cat("\n⚠ WARNING: No variable features set in Seurat object!\n")
  cat("  Attempting to load from hvg_genes.csv (fallback)...\n")

  # Try to load HVG genes from Stage 3 output
  stage3_dir <- dirname(rds_in)
  hvg_file <- file.path(stage3_dir, "hvg_genes.csv")

  if (file.exists(hvg_file)) {
    hvg_genes <- read.csv(hvg_file)$gene
    VariableFeatures(seu) <- hvg_genes
    cat(sprintf("  ✓ Loaded %d variable features from CSV\n", length(hvg_genes)))
  } else {
    stop("ERROR: No variable features in object and hvg_genes.csv not found. Cannot run PCA.")
  }
}

cat("\nRunning PCA (50 components)...\n")
seu <- RunPCA(seu, npcs = 50, verbose = FALSE)
cat("  ✓ PCA complete\n")

cat("\nGenerating elbow plot for agent review...\n")
pdf(file.path(outdir, "elbow_plot.pdf"), width = 8, height = 6)
print(ElbowPlot(seu, ndims = 50))
dev.off()

jpeg(file.path(outdir, "elbow_plot.jpg"), width = 800, height = 600, quality = 95)
print(ElbowPlot(seu, ndims = 50))
dev.off()
cat("  ✓ Saved elbow_plot.pdf and elbow_plot.jpg\n")

stdev <- seu[["pca"]]@stdev
variance <- (stdev^2) / sum(stdev^2) * 100

cat("\nWriting PCA variance summary...\n")
write_json(
  list(
    n_cells = ncol(seu),
    pc_variance_percent = variance
  ),
  file.path(outdir, "pca_variance.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)
cat("  ✓ Saved pca_variance.json\n")

cat("\nSaving Seurat object with PCA...\n")
saveRDS(seu, file.path(outdir, "seurat_stage4_pca.rds"))
cat("  ✓ Saved seurat_stage4_pca.rds\n")

cat("\n✓ Stage 4 complete!\n")
cat(sprintf("  - %d cells\n", ncol(seu)))
cat(sprintf("  - %d variable features\n", length(VariableFeatures(seu))))
cat("  - 50 PCs computed\n")
cat("  - Elbow plot ready for agent review\n\n")
