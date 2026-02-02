#!/usr/bin/env Rscript

#' Stage 5 — Clustering + UMAP + Markers (Agent Gate)
#' stage5_cluster_markers.R
#'
#' Purpose:
#' Final clustering and marker discovery.
#'
#' Agent decision:
#'   Accept clustering
#'   Adjust resolution
#'   Flag suspicious clusters
#'
#' Input:
#'   seurat_stage4_pca.rds
#'   PC range (e.g. 1:20)
#'
#' Output:
#'   markers_stage5.csv
#'   cluster_summary_stage5.json
#'   seurat_stage5_clustered.rds
#'   umap_plot.pdf/.jpg (UMAP cluster visualization)

library(Seurat)
library(jsonlite)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
rds_in <- args[1]
outdir <- args[2]
min_pc <- as.numeric(args[3])
max_pc <- as.numeric(args[4])
resolution <- if (length(args) >= 5) as.numeric(args[5]) else 0.5  # Optional resolution parameter
pcs <- min_pc:max_pc  # Create PC range (e.g., 1:21)

# Create output directory if it doesn't exist
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("Using clustering resolution: %.2f\n", resolution))

seu <- readRDS(rds_in)

seu <- FindNeighbors(seu, dims = pcs)
seu <- FindClusters(seu, resolution = resolution)
seu <- RunUMAP(seu, dims = pcs)

# Generate UMAP cluster plot for agent review
pdf(file.path(outdir, "umap_plot.pdf"), width = 10, height = 8)
print(DimPlot(seu, reduction = "umap", label = TRUE, pt.size = 0.5))
dev.off()

jpeg(file.path(outdir, "umap_plot.jpg"), width = 1000, height = 800, quality = 95)
print(DimPlot(seu, reduction = "umap", label = TRUE, pt.size = 0.5))
dev.off()

# CRITICAL FIX (2026-01-27): Seurat v5 requires joining data layers before FindAllMarkers
# Without this, FindAllMarkers fails with "data layers are not joined" warning
cat("\nJoining data layers for Seurat v5 compatibility...\n")
seu <- JoinLayers(seu)
cat("  ✓ Data layers joined\n")

cat("\nFinding marker genes for all clusters...\n")
markers <- FindAllMarkers(
  seu,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  verbose = FALSE
)

cat(sprintf("  ✓ Found %d marker genes across all clusters\n", nrow(markers)))

write.csv(markers, file.path(outdir, "markers_stage5.csv"), row.names = FALSE)

cluster_sizes <- table(seu$seurat_clusters)

# Get top 3 markers per cluster for summary
# Handle case where no markers were found
if (nrow(markers) > 0) {
  top_markers_per_cluster <- markers %>%
    group_by(cluster) %>%
    arrange(desc(avg_log2FC), p_val_adj) %>%
    slice_head(n = 3) %>%
    summarise(top_markers = paste(gene, collapse = ", "), .groups = 'drop')

  top_markers_summary <- paste(
    sprintf("Cluster %s: %s", top_markers_per_cluster$cluster, top_markers_per_cluster$top_markers),
    collapse = "\n  "
  )
} else {
  cat("\n⚠ WARNING: No marker genes identified!\n")
  cat("  This can happen with:\n")
  cat("  - Very homogeneous datasets (e.g., clonal cell lines)\n")
  cat("  - Low gene expression variance between clusters\n")
  cat("  - Overly stringent filtering thresholds\n\n")

  top_markers_summary <- "No marker genes identified (homogeneous dataset)"
}

summary_json <- list(
  n_cells = ncol(seu),
  n_clusters = length(cluster_sizes),
  total_markers = nrow(markers),
  median_markers_per_cluster = if (nrow(markers) > 0) median(table(markers$cluster)) else 0,
  resolution = resolution,  # Use actual resolution parameter
  cluster_sizes = as.list(cluster_sizes),
  top_markers_summary = top_markers_summary
)

write_json(summary_json,
           file.path(outdir, "cluster_summary_stage5.json"),
           pretty = TRUE,
           auto_unbox = TRUE)

saveRDS(seu, file.path(outdir, "seurat_stage5_clustered.rds"))
