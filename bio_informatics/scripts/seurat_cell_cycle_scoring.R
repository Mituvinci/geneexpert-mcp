#!/usr/bin/Rscript
#
# seurat_cell_cycle_scoring.R - Stage 3a: Cell Cycle Scoring + Visualization
#
# Usage:
#   Rscript seurat_cell_cycle_scoring.R <count_matrix.csv> <organism> <output_prefix>
#
# Arguments:
#   count_matrix.csv - Gene expression count matrix (genes x cells)
#   organism         - 'human' or 'mouse'
#   output_prefix    - Prefix for output files
#
# Outputs:
#   <prefix>_seurat_scored.rds           - Seurat object with cell cycle scores
#   <prefix>_cell_cycle_before.pdf       - DimPlot (PCA) colored by Phase
#   <prefix>_cell_cycle_before.jpg       - DimPlot for agents
#   <prefix>_cell_cycle_summary.txt      - Summary statistics
#

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("\n================================================================================
STAGE 3a: CELL CYCLE SCORING
================================================================================

Usage: seurat_cell_cycle_scoring.R <count_matrix.csv> <organism> <output_prefix>

Arguments:
  count_matrix.csv - Gene expression count matrix (genes x cells)
  organism         - 'human' or 'mouse'
  output_prefix    - Prefix for output files

Example:
  Rscript seurat_cell_cycle_scoring.R counts.csv human my_experiment

================================================================================
", call. = FALSE)
}

library(Seurat)
library(ggplot2)

count_file <- args[1]
organism <- tolower(args[2])
output_prefix <- args[3]

cat("\n================================================================================\n")
cat("STAGE 3a: CELL CYCLE SCORING\n")
cat("================================================================================\n\n")
cat("Input file:", count_file, "\n")
cat("Organism:", organism, "\n")
cat("Output prefix:", output_prefix, "\n\n")

# ============================================================================
# Load Data
# ============================================================================

cat("Loading count matrix...\n")
counts <- read.csv(count_file, header = TRUE, row.names = 1, check.names = FALSE)
cat("Loaded:", nrow(counts), "genes x", ncol(counts), "cells\n\n")

# ============================================================================
# Create Seurat Object + Basic QC
# ============================================================================

cat("Creating Seurat object...\n")
seurat_obj <- CreateSeuratObject(counts = counts, project = output_prefix, min.cells = 3, min.features = 200)

# Calculate mitochondrial percentage
if (organism == "human") {
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
} else if (organism == "mouse") {
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
} else {
  stop("ERROR: Organism must be 'human' or 'mouse'")
}

cat("\nBasic QC metrics:\n")
cat("  Median nFeature_RNA:", median(seurat_obj$nFeature_RNA), "\n")
cat("  Median nCount_RNA:", median(seurat_obj$nCount_RNA), "\n")
cat("  Median percent.mt:", round(median(seurat_obj$percent.mt), 2), "%\n\n")

# Filter cells (basic thresholds)
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)
cat("Cells after QC filtering:", ncol(seurat_obj), "\n\n")

# ============================================================================
# Normalization + Variable Features
# ============================================================================

cat("Normalizing data...\n")
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)

cat("Finding variable features...\n")
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

# ============================================================================
# Cell Cycle Scoring
# ============================================================================

cat("\nPerforming cell cycle scoring...\n")

# Load cell cycle gene lists (built-in Seurat lists)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# Convert to mouse genes if needed
if (organism == "mouse") {
  s.genes <- tolower(s.genes)
  g2m.genes <- tolower(g2m.genes)

  # Capitalize first letter (mouse gene format)
  s.genes <- paste0(toupper(substring(s.genes, 1, 1)), substring(s.genes, 2))
  g2m.genes <- paste0(toupper(substring(g2m.genes, 1, 1)), substring(g2m.genes, 2))
}

# Score cell cycle phases
seurat_obj <- CellCycleScoring(seurat_obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# Cell cycle distribution
phase_table <- table(seurat_obj$Phase)
cat("\nCell Cycle Distribution:\n")
for (phase in names(phase_table)) {
  cat(sprintf("  %s: %d cells (%.1f%%)\n", phase, phase_table[phase], 100 * phase_table[phase] / sum(phase_table)))
}

# ============================================================================
# PCA (before regression)
# ============================================================================

cat("\nRunning PCA...\n")
seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj))
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))

# Get variance explained
pca_var <- (seurat_obj@reductions$pca@stdev)^2
pca_var_pct <- 100 * pca_var / sum(pca_var)

cat("\nPCA Variance Explained:\n")
for (i in 1:5) {
  cat(sprintf("  PC%d: %.2f%%\n", i, pca_var_pct[i]))
}

# Calculate correlation between PCs and cell cycle scores
pc_embeddings <- Embeddings(seurat_obj, reduction = "pca")
cor_s_pc1 <- cor(seurat_obj$S.Score, pc_embeddings[, 1])
cor_g2m_pc1 <- cor(seurat_obj$G2M.Score, pc_embeddings[, 1])
cor_s_pc2 <- cor(seurat_obj$S.Score, pc_embeddings[, 2])
cor_g2m_pc2 <- cor(seurat_obj$G2M.Score, pc_embeddings[, 2])

cat("\nCorrelation with Cell Cycle Scores:\n")
cat(sprintf("  PC1 vs S.Score: r=%.3f\n", cor_s_pc1))
cat(sprintf("  PC1 vs G2M.Score: r=%.3f\n", cor_g2m_pc1))
cat(sprintf("  PC2 vs S.Score: r=%.3f\n", cor_s_pc2))
cat(sprintf("  PC2 vs G2M.Score: r=%.3f\n", cor_g2m_pc2))

# ============================================================================
# Generate DimPlot (PCA colored by Phase)
# ============================================================================

cat("\nGenerating DimPlot (before regression)...\n")

pdf_file <- paste0(output_prefix, "_cell_cycle_before.pdf")
jpg_file <- paste0(output_prefix, "_cell_cycle_before.jpg")

p <- DimPlot(seurat_obj, reduction = "pca", group.by = "Phase", pt.size = 1.2) +
  ggtitle("PCA: Cell Cycle Effect (Before Regression)") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "right"
  )

ggsave(pdf_file, p, width = 8, height = 6)
ggsave(jpg_file, p, width = 8, height = 6, dpi = 150)

cat("  PDF:", pdf_file, "\n")
cat("  JPG:", jpg_file, "\n")

# ============================================================================
# Save Seurat Object + Summary
# ============================================================================

cat("\nSaving Seurat object...\n")
rds_file <- paste0(output_prefix, "_seurat_scored.rds")
saveRDS(seurat_obj, file = rds_file)
cat("  RDS:", rds_file, "\n")

# Write summary file
summary_file <- paste0(output_prefix, "_cell_cycle_summary.txt")
sink(summary_file)
cat("================================================================================\n")
cat("STAGE 3a: CELL CYCLE SCORING SUMMARY\n")
cat("================================================================================\n\n")
cat("Input file:", count_file, "\n")
cat("Organism:", organism, "\n")
cat("Output prefix:", output_prefix, "\n\n")

cat("Cell Cycle Distribution:\n")
for (phase in names(phase_table)) {
  cat(sprintf("  %s: %d cells (%.1f%%)\n", phase, phase_table[phase], 100 * phase_table[phase] / sum(phase_table)))
}

cat("\nPCA Variance Explained:\n")
for (i in 1:5) {
  cat(sprintf("  PC%d: %.2f%%\n", i, pca_var_pct[i]))
}

cat("\nCorrelation with Cell Cycle Scores:\n")
cat(sprintf("  PC1 vs S.Score: r=%.3f\n", cor_s_pc1))
cat(sprintf("  PC1 vs G2M.Score: r=%.3f\n", cor_g2m_pc1))
cat(sprintf("  PC2 vs S.Score: r=%.3f\n", cor_s_pc2))
cat(sprintf("  PC2 vs G2M.Score: r=%.3f\n", cor_g2m_pc2))

cat("\nOutput Files:\n")
cat("  Seurat object:", rds_file, "\n")
cat("  DimPlot (PDF):", pdf_file, "\n")
cat("  DimPlot (JPG):", jpg_file, "\n")
cat("  Summary:", summary_file, "\n")

cat("\n================================================================================\n")
cat("STAGE 3a COMPLETE - Proceeding to Stage 3b (Cell Cycle Regression)\n")
cat("================================================================================\n")
sink()

cat("\nSummary written to:", summary_file, "\n")

cat("\n================================================================================\n")
cat("STAGE 3a COMPLETE\n")
cat("================================================================================\n")
