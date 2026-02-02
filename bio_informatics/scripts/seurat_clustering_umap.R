#!/usr/bin/Rscript
#
# seurat_clustering_umap.R - Stage 3c: Clustering + UMAP Visualization
#
# Usage:
#   Rscript seurat_clustering_umap.R <seurat_regressed.rds> <resolution> <output_prefix>
#
# Arguments:
#   seurat_regressed.rds - Seurat object from Stage 3b (with regressed data)
#   resolution           - Clustering resolution (e.g., 0.5, 0.8, 1.0)
#   output_prefix        - Prefix for output files
#
# Outputs:
#   <prefix>_seurat_clustered.rds        - Seurat object with clusters + UMAP
#   <prefix>_umap_clusters.pdf           - UMAP colored by cluster
#   <prefix>_umap_clusters.jpg           - UMAP for agents
#   <prefix>_umap_phase.pdf              - UMAP colored by cell cycle phase
#   <prefix>_umap_phase.jpg              - UMAP Phase for agents
#   <prefix>_clustering_summary.txt      - Cluster sizes and metrics
#

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("\n================================================================================
STAGE 3c: CLUSTERING + UMAP
================================================================================

Usage: seurat_clustering_umap.R <seurat_regressed.rds> <resolution> <output_prefix>

Arguments:
  seurat_regressed.rds - Seurat object from Stage 3b (with regressed data)
  resolution           - Clustering resolution (e.g., 0.5, 0.8, 1.0)
  output_prefix        - Prefix for output files

Example:
  Rscript seurat_clustering_umap.R my_experiment_seurat_regressed.rds 0.8 my_experiment

================================================================================
", call. = FALSE)
}

library(Seurat)
library(ggplot2)

rds_file <- args[1]
resolution <- as.numeric(args[2])
output_prefix <- args[3]

cat("\n================================================================================\n")
cat("STAGE 3c: CLUSTERING + UMAP\n")
cat("================================================================================\n\n")
cat("Input RDS:", rds_file, "\n")
cat("Resolution:", resolution, "\n")
cat("Output prefix:", output_prefix, "\n\n")

# ============================================================================
# Load Seurat Object from Stage 3b
# ============================================================================

cat("Loading Seurat object from Stage 3b...\n")
seurat_obj <- readRDS(rds_file)
cat("Loaded:", ncol(seurat_obj), "cells\n\n")

# ============================================================================
# Find Neighbors + Clusters
# ============================================================================

cat("Finding neighbors...\n")
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)

cat("Finding clusters (resolution =", resolution, ")...\n")
seurat_obj <- FindClusters(seurat_obj, resolution = resolution)

# Cluster sizes
cluster_table <- table(Idents(seurat_obj))
n_clusters <- length(cluster_table)

cat("\nClustering Results:\n")
cat("  Number of clusters:", n_clusters, "\n")
cat("\nCluster Sizes:\n")
for (cluster in names(cluster_table)) {
  cat(sprintf("  Cluster %s: %d cells (%.1f%%)\n", cluster, cluster_table[cluster], 100 * cluster_table[cluster] / sum(cluster_table)))
}

# ============================================================================
# Run UMAP
# ============================================================================

cat("\nRunning UMAP...\n")
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)
cat("UMAP complete.\n")

# ============================================================================
# Generate UMAP Plots
# ============================================================================

cat("\nGenerating UMAP plots...\n")

# UMAP colored by cluster
p_clusters <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = 1.0) +
  ggtitle("UMAP: Clusters") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "right"
  )

pdf_clusters <- paste0(output_prefix, "_umap_clusters.pdf")
jpg_clusters <- paste0(output_prefix, "_umap_clusters.jpg")
ggsave(pdf_clusters, p_clusters, width = 8, height = 6)
ggsave(jpg_clusters, p_clusters, width = 8, height = 6, dpi = 150)

cat("  Clusters PDF:", pdf_clusters, "\n")
cat("  Clusters JPG:", jpg_clusters, "\n")

# UMAP colored by cell cycle phase (verify regression)
p_phase <- DimPlot(seurat_obj, reduction = "umap", group.by = "Phase", pt.size = 1.0) +
  ggtitle("UMAP: Cell Cycle Phase (Post-Regression)") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "right"
  )

pdf_phase <- paste0(output_prefix, "_umap_phase.pdf")
jpg_phase <- paste0(output_prefix, "_umap_phase.jpg")
ggsave(pdf_phase, p_phase, width = 8, height = 6)
ggsave(jpg_phase, p_phase, width = 8, height = 6, dpi = 150)

cat("  Phase PDF:", pdf_phase, "\n")
cat("  Phase JPG:", jpg_phase, "\n")

# ============================================================================
# Find Marker Genes (Top 5 per cluster)
# ============================================================================

cat("\nFinding marker genes for clusters...\n")
markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, max.cells.per.ident = 200)

if (nrow(markers) > 0) {
  cat("\nTop marker genes per cluster:\n")
  top_markers <- markers %>%
    group_by(cluster) %>%
    top_n(n = 5, wt = avg_log2FC)

  for (cluster_id in unique(top_markers$cluster)) {
    cluster_genes <- top_markers$gene[top_markers$cluster == cluster_id]
    cat(sprintf("  Cluster %s: %s\n", cluster_id, paste(cluster_genes, collapse = ", ")))
  }

  # Save full marker table
  marker_file <- paste0(output_prefix, "_cluster_markers.csv")
  write.csv(markers, file = marker_file, row.names = FALSE)
  cat("\nFull marker table saved to:", marker_file, "\n")
} else {
  cat("\nNo significant marker genes found.\n")
  marker_file <- NA
}

# ============================================================================
# Cell Cycle Phase Distribution per Cluster
# ============================================================================

cat("\nCell cycle phase distribution per cluster:\n")
phase_by_cluster <- table(Idents(seurat_obj), seurat_obj$Phase)
for (cluster_id in rownames(phase_by_cluster)) {
  cat(sprintf("  Cluster %s: G1=%d, S=%d, G2M=%d\n",
              cluster_id,
              phase_by_cluster[cluster_id, "G1"],
              phase_by_cluster[cluster_id, "S"],
              phase_by_cluster[cluster_id, "G2M"]))
}

# ============================================================================
# Save Clustered Seurat Object + Summary
# ============================================================================

cat("\nSaving clustered Seurat object...\n")
rds_clustered <- paste0(output_prefix, "_seurat_clustered.rds")
saveRDS(seurat_obj, file = rds_clustered)
cat("  RDS:", rds_clustered, "\n")

# Write summary file
summary_file <- paste0(output_prefix, "_clustering_summary.txt")
sink(summary_file)
cat("================================================================================\n")
cat("STAGE 3c: CLUSTERING + UMAP SUMMARY\n")
cat("================================================================================\n\n")
cat("Input RDS:", rds_file, "\n")
cat("Resolution:", resolution, "\n")
cat("Output prefix:", output_prefix, "\n")
cat("Cells processed:", ncol(seurat_obj), "\n\n")

cat("================================================================================\n")
cat("CLUSTERING RESULTS\n")
cat("================================================================================\n\n")
cat("Number of clusters:", n_clusters, "\n\n")

cat("Cluster Sizes:\n")
for (cluster in names(cluster_table)) {
  cat(sprintf("  Cluster %s: %d cells (%.1f%%)\n", cluster, cluster_table[cluster], 100 * cluster_table[cluster] / sum(cluster_table)))
}

cat("\nCell Cycle Phase Distribution per Cluster:\n")
for (cluster_id in rownames(phase_by_cluster)) {
  cat(sprintf("  Cluster %s: G1=%d, S=%d, G2M=%d\n",
              cluster_id,
              phase_by_cluster[cluster_id, "G1"],
              phase_by_cluster[cluster_id, "S"],
              phase_by_cluster[cluster_id, "G2M"]))
}

if (nrow(markers) > 0) {
  cat("\nTop Marker Genes per Cluster (Top 5):\n")
  for (cluster_id in unique(top_markers$cluster)) {
    cluster_genes <- top_markers$gene[top_markers$cluster == cluster_id]
    cat(sprintf("  Cluster %s: %s\n", cluster_id, paste(cluster_genes, collapse = ", ")))
  }
}

cat("\nOutput Files:\n")
cat("  Clustered Seurat object:", rds_clustered, "\n")
cat("  UMAP Clusters (PDF):", pdf_clusters, "\n")
cat("  UMAP Clusters (JPG):", jpg_clusters, "\n")
cat("  UMAP Phase (PDF):", pdf_phase, "\n")
cat("  UMAP Phase (JPG):", jpg_phase, "\n")
if (!is.na(marker_file)) {
  cat("  Cluster markers:", marker_file, "\n")
}
cat("  Summary:", summary_file, "\n")

cat("\n================================================================================\n")
cat("STAGE 3c COMPLETE - Ready for Stage 4 (Differential Expression)\n")
cat("================================================================================\n")
sink()

cat("\nSummary written to:", summary_file, "\n")

cat("\n================================================================================\n")
cat("STAGE 3c COMPLETE\n")
cat("================================================================================\n")
