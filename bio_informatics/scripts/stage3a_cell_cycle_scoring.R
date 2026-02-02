#!/usr/bin/env Rscript

#' Stage 3A — Normalization + HVGs + CELL CYCLE SCORING
#' stage3a_cell_cycle_scoring.R
#'
#' Purpose:
#' 1. Normalize data
#' 2. Find HVGs
#' 3. CELL CYCLE PHASE DETECTION (score cells)
#' 4. SAVE CELL CYCLE PLOTS (before removal) for agent review
#'
#' Agent decision AFTER this stage:
#'   REMOVE_CELL_CYCLE - Strong cell cycle effect visible
#'   SKIP_CELL_CYCLE - No visible cell cycle effect
#'   UNCERTAIN - Escalate to user
#'
#' Input:
#'   seurat_stage2_filtered.rds
#'   organism (human or mouse)
#'
#' Output:
#'   seurat_stage3a_scored.rds (with S.Score, G2M.Score, Phase metadata)
#'   hvg_genes.csv
#'   cell_cycle_before.pdf/.jpg (for agent review)
#'   cell_cycle_summary_stage3a.json

library(Seurat)
library(jsonlite)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
rds_in <- args[1]
outdir <- args[2]
organism <- if (length(args) >= 3) tolower(args[3]) else "human"

cat("\n================================================================================\n")
cat("scRNA-seq Stage 3A: Normalization + HVGs + CELL CYCLE SCORING\n")
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

# Check if we found enough HVGs (minimum threshold: 100)
# Homogeneous cell lines (e.g., cancer lines) may have very few HVGs
MIN_HVG_THRESHOLD <- 100

if (length(hvg_genes) < MIN_HVG_THRESHOLD) {
  cat(sprintf("\n⚠ WARNING: Only %d HVGs found (threshold: %d)\n", length(hvg_genes), MIN_HVG_THRESHOLD))
  cat("  This can happen with:\n")
  cat("  - Homogeneous cell lines (cancer lines, clonal populations)\n")
  cat("  - Small datasets with low cell diversity\n")
  cat("  - Technical issues (low read depth, over-filtered)\n")
  cat("\n  SOLUTION: Using all genes instead of HVGs for downstream analysis\n\n")

  # Use all genes instead
  hvg_genes <- rownames(seu)
  cat(sprintf("  → Using all %d genes for analysis\n", length(hvg_genes)))
}

write.csv(
  data.frame(gene = hvg_genes),
  file.path(outdir, "hvg_genes.csv"),
  row.names = FALSE
)

cat(sprintf("  → Final gene set: %d genes\n", length(hvg_genes)))

# Step 3: CELL CYCLE PHASE DETECTION
cat("\n================================================================================\n")
cat("STEP 3: CELL CYCLE PHASE DETECTION\n")
cat("================================================================================\n")

# Load cell cycle markers (built into Seurat - human genes)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# Convert to mouse genes if needed
if (organism == "mouse") {
  cat("Converting cell cycle markers to mouse gene format...\n")
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
min_markers_threshold <- 15  # Need at least 15 markers per phase (Seurat CellCycleScoring requires more than minimal)

if (s_genes_present < min_markers_threshold || g2m_genes_present < min_markers_threshold) {
  cat("\n⚠ WARNING: Insufficient cell cycle markers detected!\n")
  cat("  This can happen with:\n")
  cat("  - Small/shallow datasets (low read depth)\n")
  cat("  - Species-specific marker mismatches\n")
  cat("  - Highly filtered datasets\n")
  cat("  - Gene ID format mismatch (Ensembl IDs vs gene symbols)\n")
  cat("\n  RECOMMENDATION: Skip cell cycle correction\n\n")

  # Create summary
  summary_json <- list(
    n_hvgs = length(hvg_genes),
    cell_cycle_detected = FALSE,
    sufficient_markers = FALSE,
    s_markers_present = s_genes_present,
    g2m_markers_present = g2m_genes_present,
    recommendation = "SKIP_CELL_CYCLE",
    reason = sprintf("Insufficient markers: %d S-phase, %d G2M-phase (need %d minimum)",
                     s_genes_present, g2m_genes_present, min_markers_threshold),
    organism = organism
  )

  write_json(
    summary_json,
    file.path(outdir, "cell_cycle_summary_stage3a.json"),
    pretty = TRUE,
    auto_unbox = TRUE
  )

  saveRDS(seu, file.path(outdir, "seurat_stage3a_scored.rds"))

  cat("\n✓ Stage 3A complete (insufficient markers for cell cycle detection)\n")
  quit(save = "no", status = 0)
}

# Score cells for cell cycle phases (with error handling)
cat("\nScoring cells for S phase and G2M phase...\n")

# Use tryCatch to handle both errors and warnings without script exit
cell_cycle_success <- FALSE
seu_scored <- tryCatch(
  {
    # Suppress warnings but capture them
    withCallingHandlers(
      {
        CellCycleScoring(seu,
                        s.features = s.genes,
                        g2m.features = g2m.genes,
                        set.ident = FALSE)
      },
      warning = function(w) {
        cat(sprintf("  ⚠ Warning: %s\n", w$message))
        invokeRestart("muffleWarning")
      }
    )
  },
  error = function(e) {
    cat(sprintf("\n⚠ CellCycleScoring failed: %s\n", e$message))
    cat("  This is common for:\n")
    cat("  - Homogeneous cell lines (clonal populations like cancer lines)\n")
    cat("  - Gene name format mismatches (case sensitivity)\n")
    cat("  - Sparse expression of cell cycle genes\n")
    cat("  → Proceeding without cell cycle scoring\n\n")
    NULL  # Return NULL on error
  }
)

# Check if scoring succeeded
if (!is.null(seu_scored)) {
  seu <- seu_scored
  cell_cycle_success <- TRUE
  cat("✓ Cell cycle scoring completed successfully\n")
}

# Check if scoring succeeded by verifying metadata columns were added
has_s_score <- "S.Score" %in% colnames(seu@meta.data)
has_g2m_score <- "G2M.Score" %in% colnames(seu@meta.data)
has_phase <- "Phase" %in% colnames(seu@meta.data)

cat(sprintf("\nCell cycle scoring check:\n"))
cat(sprintf("  - S.Score column added: %s\n", if(has_s_score) "YES" else "NO"))
cat(sprintf("  - G2M.Score column added: %s\n", if(has_g2m_score) "YES" else "NO"))
cat(sprintf("  - Phase column added: %s\n", if(has_phase) "YES" else "NO"))

if (!cell_cycle_success || !has_s_score || !has_g2m_score || !has_phase) {
  cat("\n⚠ WARNING: Cell cycle scoring failed!\n")
  cat("  Likely causes:\n")
  cat("  - Gene name format mismatch (e.g., Ensembl IDs vs gene symbols)\n")
  cat("  - Too few cell cycle genes detected in dataset\n")
  cat("  - Species mismatch\n")
  cat(sprintf("  - Markers present: S=%d/%d, G2M=%d/%d\n",
              s_genes_present, length(s.genes), g2m_genes_present, length(g2m.genes)))
  cat("\n  RECOMMENDATION: Skip cell cycle correction\n\n")

  summary_json <- list(
    n_hvgs = length(hvg_genes),
    cell_cycle_detected = FALSE,
    sufficient_markers = s_genes_present >= min_markers_threshold && g2m_genes_present >= min_markers_threshold,
    s_markers_present = s_genes_present,
    g2m_markers_present = g2m_genes_present,
    recommendation = "SKIP_CELL_CYCLE",
    reason = "Cell cycle scoring failed (CellCycleScoring could not find enough features despite markers present)",
    organism = organism
  )

  write_json(
    summary_json,
    file.path(outdir, "cell_cycle_summary_stage3a.json"),
    pretty = TRUE,
    auto_unbox = TRUE
  )

  saveRDS(seu, file.path(outdir, "seurat_stage3a_scored.rds"))

  cat("\n✓ Stage 3A complete (cell cycle scoring failed - skipping correction)\n")
  quit(save = "no", status = 0)
}

# Get phase distribution
phase_table <- table(seu$Phase)
cat("\nCell cycle phase distribution:\n")
print(phase_table)

# Step 4: SAVE CELL CYCLE PLOTS (BEFORE REMOVAL) - FOR AGENT REVIEW
cat("\n================================================================================\n")
cat("STEP 4: SAVING CELL CYCLE PLOTS (BEFORE removal) - FOR AGENT REVIEW\n")
cat("================================================================================\n")

# Need to run PCA temporarily to visualize cell cycle effects
cat("Running temporary PCA for visualization...\n")
seu_temp <- ScaleData(seu, features = hvg_genes, verbose = FALSE)
seu_temp <- RunPCA(seu_temp, features = hvg_genes, npcs = 30, verbose = FALSE)

# Plot 1: PCA colored by cell cycle phase (BEFORE removal)
p1 <- DimPlot(seu_temp, reduction = "pca", group.by = "Phase") +
  ggtitle("PCA: Cell Cycle Phase (BEFORE removal)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Plot 2: FeaturePlot showing S.Score and G2M.Score
p2 <- FeaturePlot(seu_temp, reduction = "pca", features = c("S.Score", "G2M.Score"), ncol = 2)

# Save as PDF (for user) - contains BOTH plots
pdf(file.path(outdir, "cell_cycle_before.pdf"), width = 10, height = 5)
print(p1)
print(p2)
dev.off()

# Save Plot 1 (Phase - 3 colors) as JPEG for agent review ← THIS IS WHAT AGENTS SEE!
jpeg(file.path(outdir, "cell_cycle_phase_before.jpg"), width = 1000, height = 500, quality = 95)
print(p1)
dev.off()

# Save Plot 2 (Scores - blue gradients) as JPEG for reference
jpeg(file.path(outdir, "cell_cycle_scores_before.jpg"), width = 1000, height = 500, quality = 95)
print(p2)
dev.off()

cat("✓ Saved: cell_cycle_before.pdf (both plots)\n")
cat("✓ Saved: cell_cycle_phase_before.jpg (3-color Phase plot for agents)\n")
cat("✓ Saved: cell_cycle_scores_before.jpg (S.Score/G2M.Score gradients)\n")

# Calculate variance explained by cell cycle scores
cat("\nCalculating variance explained by cell cycle...\n")
pca_embeddings <- seu_temp@reductions$pca@cell.embeddings[, 1:2]
cor_s <- cor(pca_embeddings[, 1], seu_temp$S.Score)
cor_g2m <- cor(pca_embeddings[, 1], seu_temp$G2M.Score)
cat(sprintf("  PC1 correlation with S.Score: %.3f\n", cor_s))
cat(sprintf("  PC1 correlation with G2M.Score: %.3f\n", cor_g2m))

# Step 5: Write summary JSON
summary_json <- list(
  n_hvgs = length(hvg_genes),
  cell_cycle_detected = TRUE,
  sufficient_markers = TRUE,
  s_markers_present = s_genes_present,
  g2m_markers_present = g2m_genes_present,
  phase_distribution = as.list(phase_table),
  s_score_range = c(min(seu$S.Score), max(seu$S.Score)),
  g2m_score_range = c(min(seu$G2M.Score), max(seu$G2M.Score)),
  pc1_correlation_s_score = cor_s,
  pc1_correlation_g2m_score = cor_g2m,
  cell_cycle_plot_phase = "cell_cycle_phase_before.jpg",
  cell_cycle_plot_scores = "cell_cycle_scores_before.jpg",
  recommendation = if (abs(cor_s) > 0.3 || abs(cor_g2m) > 0.3) "REMOVE_CELL_CYCLE" else "SKIP_CELL_CYCLE",
  organism = organism
)

write_json(
  summary_json,
  file.path(outdir, "cell_cycle_summary_stage3a.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)

# Save Seurat object with cell cycle scores (but NOT regressed yet)
saveRDS(seu, file.path(outdir, "seurat_stage3a_scored.rds"))

cat("\n✓ Stage 3A complete!\n")
cat(sprintf("  - HVGs: %d\n", length(hvg_genes)))
cat(sprintf("  - Cell cycle phases detected: %s\n", paste(names(phase_table), collapse = ", ")))
cat("  - Cell cycle scores calculated (S.Score, G2M.Score)\n")
cat("  - Saved: cell_cycle_phase_before.jpg for agent review (3-color Phase plot)\n")
cat("\n→ Next: Agents will review cell_cycle_phase_before.jpg (3-color Phase plot) and decide:\n")
cat("    REMOVE_CELL_CYCLE (if strong effect) → Stage 3B\n")
cat("    SKIP_CELL_CYCLE (if no effect) → Stage 4\n\n")
