#!/usr/bin/env Rscript

#' Stage 2 — Filter + QC
#' stage2_filter_qc.R
#'
#' Purpose:
#' Apply QC thresholds and summarize effect.
#'
#' Agent decision:
#'   PROCEED
#'   PROCEED_WITH_WARNING
#'   STOP / ADJUST_THRESHOLDS
#'
#' Input:
#'   seurat_stage1_raw.rds
#'
#' Output:
#'   seurat_stage2_filtered.rds
#'   qc_summary_stage2.json
#'   qc_distributions_before.pdf/.jpg (violin plots before filtering)
#'   qc_distributions_after.pdf/.jpg (violin plots after filtering)


library(Seurat)
library(jsonlite)

args <- commandArgs(trailingOnly = TRUE)
rds_in <- args[1]
outdir <- args[2]

# Agent-recommended thresholds (dynamic)
nFeature_min <- as.numeric(args[3])
nFeature_max <- as.numeric(args[4])
percent_mt_max <- as.numeric(args[5])

cat(sprintf("Applying agent-recommended QC thresholds:\n"))
cat(sprintf("  nFeature_RNA: %d - %d genes\n", nFeature_min, nFeature_max))
cat(sprintf("  percent.mt: < %d%%\n", percent_mt_max))
cat("\n")

seu <- readRDS(rds_in)
before <- ncol(seu)

# Generate QC distribution plots BEFORE filtering
pdf(file.path(outdir, "qc_distributions_before.pdf"), width = 12, height = 4)
print(VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1))
dev.off()

jpeg(file.path(outdir, "qc_distributions_before.jpg"), width = 1200, height = 400, quality = 95)
print(VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1))
dev.off()

# Apply dynamic thresholds
seu <- subset(
  seu,
  subset = nFeature_RNA > nFeature_min &
           nFeature_RNA < nFeature_max &
           percent.mt < percent_mt_max
)

after <- ncol(seu)

# Generate QC distribution plots AFTER filtering
pdf(file.path(outdir, "qc_distributions_after.pdf"), width = 12, height = 4)
print(VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1))
dev.off()

jpeg(file.path(outdir, "qc_distributions_after.jpg"), width = 1200, height = 400, quality = 95)
print(VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1))
dev.off()

cat(sprintf("Filtered: %d → %d cells (%.2f%% removed)\n",
            before, after, 100 * (before - after) / before))

summary_json <- list(
  cells_before = before,
  cells_after = after,
  percent_removed = round(100 * (before - after) / before, 2),
  thresholds = list(
    min_features = nFeature_min,
    max_features = nFeature_max,
    max_percent_mt = percent_mt_max
  )
)

write_json(summary_json,
           file.path(outdir, "qc_summary_stage2.json"),
           pretty = TRUE,
           auto_unbox = TRUE)

saveRDS(seu, file.path(outdir, "seurat_stage2_filtered.rds"))
