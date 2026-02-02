#!/usr/bin/env Rscript

#' Stage 1 — Load + QC Metrics
#' stage1_load_qc.R
#'
#' Purpose:
#' Load 10x data, compute QC metrics.
#' No filtering. No agent decision.
#'
#' Input:
#'   filtered_feature_bc_matrix.h5
#'
#' Output:
#'   seurat_stage1_raw.rds
#'   qc_metrics_stage1.csv
#'   qc_summary_stage1.json


library(Seurat)
library(jsonlite)
library(Seurat)
library(sctransform)
library(tidyr)
library(dplyr)


args <- commandArgs(trailingOnly = TRUE)
input_path <- args[1]
outdir <- args[2]

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Auto-detect input format
if (dir.exists(input_path)) {
  # Directory with barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz
  cat("Loading from 10x directory format...\n")
  counts <- Read10X(input_path)
} else if (file.exists(input_path) && grepl("\\.h5$", input_path)) {
  # HDF5 file (requires hdf5r package and C library)
  cat("Loading from HDF5 format...\n")
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("hdf5r package not available. Please use directory format or CSV instead.")
  }
  counts <- Read10X_h5(input_path)
} else if (file.exists(input_path) && grepl("\\.(csv|txt|tsv)$", input_path, ignore.case = TRUE)) {
  # CSV/TSV file (genes as rows, cells as columns)
  cat("Loading from CSV/TSV format...\n")
  cat("  Expected format: genes as rows, cells as columns, first column = gene names\n")

  # Detect separator
  sep <- if (grepl("\\.csv$", input_path, ignore.case = TRUE)) "," else "\t"

  # Read the file
  counts <- read.table(input_path,
                      sep = sep,
                      header = TRUE,
                      row.names = 1,
                      check.names = FALSE,
                      stringsAsFactors = FALSE)

  # Convert to matrix
  counts <- as.matrix(counts)

  cat(sprintf("  Loaded %d genes x %d cells\n", nrow(counts), ncol(counts)))
} else {
  stop("Input must be:\n  - Directory (10x format)\n  - .h5 file (10x HDF5)\n  - .csv/.txt/.tsv file (genes=rows, cells=columns)")
}

seu <- CreateSeuratObject(
  counts = counts,
  project = "scRNA",
  min.cells = 3,
  min.features = 200
)

seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

qc <- seu@meta.data[, c("nFeature_RNA", "nCount_RNA", "percent.mt")]

write.csv(qc, file.path(outdir, "qc_metrics_stage1.csv"))

summary_json <- list(
  cells_total = ncol(seu),
  genes_total = nrow(seu),
  nFeature_median = median(qc$nFeature_RNA),
  nCount_median = median(qc$nCount_RNA),
  percent_mt_median = median(qc$percent.mt)
)

write_json(summary_json,
           file.path(outdir, "qc_summary_stage1.json"),
           pretty = TRUE,
           auto_unbox = TRUE)

saveRDS(seu, file.path(outdir, "seurat_stage1_raw.rds"))
