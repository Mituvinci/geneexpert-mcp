#!/usr/bin/env Rscript

#' Stage 3 — Normalization + HVGs + CELL CYCLE ANALYSIS
#' stage3_normalize_hvg.R
#'
#' Purpose:
#' Standard Seurat preprocessing WITH CELL CYCLE CORRECTION
#' 1. Normalize data
#' 2. Find HVGs
#' 3. CELL CYCLE PHASE DETECTION
#' 4. SAVE CELL CYCLE PLOTS (before removal)
#' 5. REMOVE CELL CYCLE EFFECTS (regress out S.Score, G2M.Score)
#' 6. SCALE DATA
#' 7. RE-PLOT AND SAVE (after correction)
#'
#' No agent decision. Always runs (CRITICAL for downstream clustering).
#'
#' Input:
#'   seurat_stage2_filtered.rds
#'   organism (human or mouse)
#'
#' Output:
#'   seurat_stage3_norm.rds
#'   hvg_genes.csv
#'   cell_cycle_before.pdf (if cell cycle detected)
#'   cell_cycle_after.pdf (if cell cycle detected)
#'   hvg_summary_stage3.json

library(Seurat)
library(jsonlite)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
rds_in <- args[1]
outdir <- args[2]
organism <- if (length(args) >= 3) tolower(args[3]) else "human"

cat("\n================================================================================\n")
cat("scRNA-seq Stage 3: Normalization + HVGs + CELL CYCLE ANALYSIS\n")
cat("================================================================================\n\n")

cat("Loading Seurat object from Stage 2...\n")
seu <- readRDS(rds_in)

# Step 1: Normalize data
cat("Step 1: Normalizing data (LogNormalize)...\n")
seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)

# Step 2: Find highly variable genes
cat("Step 2: Finding highly variable genes (top 2000)...\n")
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
hvg_genes <- VariableFeatures(seu)

write.csv(
  data.frame(gene = hvg_genes),
  file.path(outdir, "hvg_genes.csv"),
  row.names = FALSE
)

cat(sprintf("  → Identified %d HVGs\n", length(hvg_genes)))

# Step 3: CELL CYCLE PHASE DETECTION
cat("\n================================================================================\n")
cat("STEP 3: CELL CYCLE PHASE DETECTION (CRITICAL!)\n")
cat("================================================================================\n")

# Load cell cycle markers (built into Seurat - human genes)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# Convert to mouse genes if needed
if (organism == "mouse") {
  cat("Converting cell cycle markers to mouse gene format...\n")
  # Convert to lowercase then capitalize first letter
  s.genes <- tolower(s.genes)
  g2m.genes <- tolower(g2m.genes)
  s.genes <- paste0(toupper(substring(s.genes, 1, 1)), substring(s.genes, 2))
  g2m.genes <- paste0(toupper(substring(g2m.genes, 1, 1)), substring(g2m.genes, 2))
}

cat(sprintf("Using %d S phase markers and %d G2M phase markers (%s)\n",
            length(s.genes), length(g2m.genes), organism))

# Check how many cell cycle genes are in the dataset
s_genes_present <- sum(s.genes %in% rownames(seu))
g2m_genes_present <- sum(g2m.genes %in% rownames(seu))

cat(sprintf("  S phase markers detected in data: %d/%d (%.1f%%)\n",
            s_genes_present, length(s.genes), 100 * s_genes_present / length(s.genes)))
cat(sprintf("  G2M phase markers detected in data: %d/%d (%.1f%%)\n",
            g2m_genes_present, length(g2m.genes), 100 * g2m_genes_present / length(g2m.genes)))

# Only proceed with cell cycle if we have sufficient markers
min_markers_threshold <- 5  # Need at least 5 markers per phase

if (s_genes_present < min_markers_threshold || g2m_genes_present < min_markers_threshold) {
  cat("\n⚠ WARNING: Insufficient cell cycle markers detected!\n")
  cat("  This can happen with:\n")
  cat("  - Small/shallow datasets (low read depth)\n")
  cat("  - Species-specific marker mismatches\n")
  cat("  - Highly filtered datasets\n")
  cat("\n  Proceeding WITHOUT cell cycle correction...\n\n")

  # Just scale without cell cycle correction
  cat("Scaling data (without cell cycle regression)...\n")
  seu <- ScaleData(seu)

  # Create summary
  summary_json <- list(
    n_hvgs = length(hvg_genes),
    cell_cycle_detected = FALSE,
    reason = sprintf("Insufficient markers: %d S-phase, %d G2M-phase (need %d minimum)",
                     s_genes_present, g2m_genes_present, min_markers_threshold)
  )

  write_json(
    summary_json,
    file.path(outdir, "hvg_summary_stage3.json"),
    pretty = TRUE,
    auto_unbox = TRUE
  )

  saveRDS(seu, file.path(outdir, "seurat_stage3_norm.rds"))

  cat("\n✓ Stage 3 complete (without cell cycle correction)\n")
  quit(save = "no", status = 0)
}

# Score cells for cell cycle phases
cat("\nScoring cells for S phase and G2M phase...\n")
seu <- CellCycleScoring(seu,
                        s.features = s.genes,
                        g2m.features = g2m.genes,
                        set.ident = FALSE)

# Check if scoring worked
if (!("S.Score" %in% colnames(seu@meta.data)) || !("G2M.Score" %in% colnames(seu@meta.data))) {
  cat("WARNING: Cell cycle scoring failed unexpectedly\n")
  cat("Proceeding without cell cycle correction...\n")

  seu <- ScaleData(seu)

  summary_json <- list(
    n_hvgs = length(hvg_genes),
    cell_cycle_detected = FALSE,
    reason = "Cell cycle scoring failed (CellCycleScoring returned no scores)"
  )

  write_json(
    summary_json,
    file.path(outdir, "hvg_summary_stage3.json"),
    pretty = TRUE,
    auto_unbox = TRUE
  )

  saveRDS(seu, file.path(outdir, "seurat_stage3_norm.rds"))
  quit(save = "no", status = 0)
}

# Get phase distribution
phase_table <- table(seu$Phase)
cat("\nCell cycle phase distribution:\n")
print(phase_table)

# Step 4: SAVE CELL CYCLE PLOTS (BEFORE REMOVAL)
cat("\n================================================================================\n")
cat("STEP 4: SAVING CELL CYCLE PLOTS (BEFORE removal)\n")
cat("================================================================================\n")

# Need to run PCA temporarily to visualize cell cycle effects
cat("Running temporary PCA for visualization...\n")
seu_temp <- ScaleData(seu, features = hvg_genes, verbose = FALSE)
seu_temp <- RunPCA(seu_temp, features = hvg_genes, npcs = 30, verbose = FALSE)

# Plot 1: PCA colored by cell cycle phase (BEFORE removal)
# Save as PDF (for user)
pdf(file.path(outdir, "cell_cycle_before.pdf"), width = 10, height = 5)

p1 <- DimPlot(seu_temp, reduction = "pca", group.by = "Phase") +
  ggtitle("PCA: Cell Cycle Phase (BEFORE removal)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

p2 <- FeaturePlot(seu_temp, reduction = "pca", features = c("S.Score", "G2M.Score"), ncol = 2)

print(p1)
print(p2)

dev.off()

# Save as JPEG (for potential agent review / reports)
jpeg(file.path(outdir, "cell_cycle_before.jpg"), width = 1000, height = 500, quality = 95)
print(p1)
print(p2)
dev.off()

cat("✓ Saved: cell_cycle_before.pdf + cell_cycle_before.jpg\n")

# Step 5: REMOVE CELL CYCLE EFFECTS (Regress out S.Score and G2M.Score)
cat("\n================================================================================\n")
cat("STEP 5: REMOVING CELL CYCLE EFFECTS (Regressing out S.Score, G2M.Score)\n")
cat("================================================================================\n")

cat("Scaling data with cell cycle regression...\n")
cat("  This removes variance due to cell cycle phase from gene expression\n")
cat("  Goal: Cells cluster by biological identity, NOT cell cycle phase\n\n")

seu <- ScaleData(seu,
                 vars.to.regress = c("S.Score", "G2M.Score"),
                 features = hvg_genes,
                 verbose = FALSE)

cat("✓ Cell cycle effects removed from scaled data\n")

# Step 6: RE-RUN PCA AND SAVE PLOTS (AFTER REMOVAL)
cat("\n================================================================================\n")
cat("STEP 6: RE-RUNNING PCA AND SAVING PLOTS (AFTER removal)\n")
cat("================================================================================\n")

seu <- RunPCA(seu, features = hvg_genes, npcs = 30, verbose = FALSE)

# Plot 2: PCA colored by cell cycle phase (AFTER removal)
# Save as PDF (for user)
pdf(file.path(outdir, "cell_cycle_after.pdf"), width = 10, height = 5)

p3 <- DimPlot(seu, reduction = "pca", group.by = "Phase") +
  ggtitle("PCA: Cell Cycle Phase (AFTER removal)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

p4 <- FeaturePlot(seu, reduction = "pca", features = c("S.Score", "G2M.Score"), ncol = 2)

print(p3)
print(p4)

dev.off()

# Save as JPEG (for potential agent review / reports)
jpeg(file.path(outdir, "cell_cycle_after.jpg"), width = 1000, height = 500, quality = 95)
print(p3)
print(p4)
dev.off()

cat("✓ Saved: cell_cycle_after.pdf + cell_cycle_after.jpg\n")

cat("\n================================================================================\n")
cat("CELL CYCLE CORRECTION COMPLETE!\n")
cat("================================================================================\n")
cat("BEFORE: Cells clustered by cell cycle phase (confounding factor)\n")
cat("AFTER:  Cell cycle variance removed, cells cluster by biological identity\n")
cat("\nCompare cell_cycle_before.pdf vs cell_cycle_after.pdf to see the effect!\n")
cat("================================================================================\n\n")

# Step 7: Write summary JSON
summary_json <- list(
  n_hvgs = length(hvg_genes),
  cell_cycle_detected = TRUE,
  phase_distribution = as.list(phase_table),
  s_score_range = c(min(seu$S.Score), max(seu$S.Score)),
  g2m_score_range = c(min(seu$G2M.Score), max(seu$G2M.Score)),
  cell_cycle_plots = list(
    before_pdf = "cell_cycle_before.pdf",
    before_jpg = "cell_cycle_before.jpg",
    after_pdf = "cell_cycle_after.pdf",
    after_jpg = "cell_cycle_after.jpg"
  )
)

write_json(
  summary_json,
  file.path(outdir, "hvg_summary_stage3.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)

# Save final Seurat object (with cell cycle corrected)
saveRDS(seu, file.path(outdir, "seurat_stage3_norm.rds"))

cat("\n✓ Stage 3 complete!\n")
cat(sprintf("  - HVGs: %d\n", length(hvg_genes)))
cat(sprintf("  - Cell cycle phases detected: %s\n", paste(names(phase_table), collapse = ", ")))
cat("  - Cell cycle effects removed from scaled data\n")
cat("  - Saved: seurat_stage3_norm.rds (ready for PCA in Stage 4)\n\n")
