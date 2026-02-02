#!/usr/bin/env Rscript
#
# qc_plots.R - Generate QC plots (PCA, MDS, Density) with outlier/batch effect detection
#
# This script generates plots for agent interpretation and outputs a summary
# that helps detect outliers and batch effects.
#
# Usage:
#   Rscript qc_plots.R -i <rpkm_file> -g <groups> -o <output_dir> [-c <control_keyword>] [-t <treatment_keyword>]
#
# Example:
#   Rscript qc_plots.R -i outRPKM.txt -g "auto" -o ./qc_results -c "cont" -t "ips"
#
# Outputs:
#   - PCA_plot.png         : PCA plot (PC1 vs PC2)
#   - MDS_plot.png         : MDS plot
#   - density_plot.png     : Expression density plot
#   - qc_summary.json      : Summary for agent interpretation (outliers, batch effect)
#   - qc_summary.txt       : Human-readable summary
#

library(ggplot2)
library(optparse)
library(jsonlite)

# Suppress warnings for cleaner output
options(warn = -1)

# ============================================================================
# Command Line Arguments
# ============================================================================

option_list <- list(
  make_option(c("-i", "--infile"),
    type = "character", default = NULL,
    help = "Input RPKM/normalized expression file", metavar = "FILE"
  ),
  make_option(c("-g", "--groups"),
    type = "character", default = "auto",
    help = "Group labels: 'auto' or comma-separated (e.g., 'ctrl,ctrl,trt,trt')"
  ),
  make_option(c("-o", "--outdir"),
    type = "character", default = "./qc_results",
    help = "Output directory for plots and summary"
  ),
  make_option(c("-c", "--control_keyword"),
    type = "character", default = NULL,
    help = "Keyword to identify control samples (e.g., 'cont', 'ctrl', 'wt')"
  ),
  make_option(c("-t", "--treatment_keyword"),
    type = "character", default = NULL,
    help = "Keyword to identify treatment samples (e.g., 'ips', 'treated', 'ko')"
  ),
  make_option(c("-e", "--expression_threshold"),
    type = "double", default = 1,
    help = "Minimum expression value to count as 'expressed' [default: 1]"
  ),
  make_option(c("-n", "--num_expressed"),
    type = "integer", default = 2,
    help = "Minimum samples a gene needs to be expressed in [default: 2]"
  ),
  make_option(c("--outlier_threshold"),
    type = "double", default = 2.5,
    help = "SD threshold for outlier detection [default: 2.5]"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

# Validate required arguments
if (is.null(opt$infile)) {
  stop("Error: Input file (-i) is required!\n\nUsage: Rscript qc_plots.R -i <rpkm_file> -g <groups> -o <output_dir>")
}

# Create output directory
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

cat("\n========================================\n")
cat("QC Plots Generator - Outlier & Batch Effect Detection\n")
cat("========================================\n\n")
cat("Input file:", opt$infile, "\n")
cat("Output dir:", opt$outdir, "\n\n")

# ============================================================================
# Load and Prepare Data
# ============================================================================

cat("Loading expression data...\n")

# Try different separators
data <- tryCatch({
  read.table(opt$infile, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
}, error = function(e) {
  tryCatch({
    read.table(opt$infile, header = TRUE, sep = " ", row.names = 1, check.names = FALSE)
  }, error = function(e2) {
    read.csv(opt$infile, header = TRUE, row.names = 1, check.names = FALSE)
  })
})

# Remove non-numeric columns (like gene symbols)
numeric_cols <- sapply(data, is.numeric)
if (sum(numeric_cols) < ncol(data)) {
  cat("Removing non-numeric columns...\n")
  # Keep gene symbol column name if exists
  gene_col <- which(!numeric_cols)[1]
  if (!is.na(gene_col)) {
    gene_symbols <- data[, gene_col]
  }
  data <- data[, numeric_cols, drop = FALSE]
}

sample_names <- colnames(data)
n_samples <- ncol(data)
n_genes <- nrow(data)

cat("Loaded:", n_genes, "genes x", n_samples, "samples\n")
cat("Samples:", paste(sample_names, collapse = ", "), "\n\n")

# ============================================================================
# Determine Groups
# ============================================================================

cat("Determining sample groups...\n")

if (opt$groups == "auto") {
  # Try to auto-detect groups from sample names
  if (!is.null(opt$control_keyword) && !is.null(opt$treatment_keyword)) {
    # Use keywords to assign groups
    groups <- ifelse(
      grepl(opt$control_keyword, sample_names, ignore.case = TRUE),
      "Control",
      ifelse(
        grepl(opt$treatment_keyword, sample_names, ignore.case = TRUE),
        "Treatment",
        "Unknown"
      )
    )
  } else {
    # Default: use first character of sample names
    groups <- substr(sample_names, 1, 3)
  }
} else {
  # Use provided groups
  groups <- strsplit(opt$groups, ",")[[1]]
  if (length(groups) != n_samples) {
    stop("Number of groups (", length(groups), ") doesn't match samples (", n_samples, ")")
  }
}

cat("Groups assigned:", paste(unique(groups), collapse = ", "), "\n")
cat("Group mapping:\n")
for (i in 1:n_samples) {
  cat("  ", sample_names[i], "->", groups[i], "\n")
}
cat("\n")

# ============================================================================
# Filter Low Expression Genes
# ============================================================================

cat("Filtering lowly expressed genes...\n")
keep <- rowSums(data > opt$expression_threshold) >= opt$num_expressed
data_filtered <- data[keep, ]
cat("Kept", nrow(data_filtered), "of", n_genes, "genes after filtering\n\n")

# ============================================================================
# PCA Analysis
# ============================================================================

cat("Performing PCA...\n")

# Handle zero variance genes
gene_vars <- apply(data_filtered, 1, var)
data_pca <- data_filtered[gene_vars > 0, ]

pca_result <- prcomp(t(data_pca), scale. = TRUE, center = TRUE)
pca_var <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)

pca_df <- data.frame(
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  PC3 = if(ncol(pca_result$x) >= 3) pca_result$x[, 3] else 0,
  Sample = sample_names,
  Group = groups
)

cat("PC1 explains", pca_var[1], "% of variance\n")
cat("PC2 explains", pca_var[2], "% of variance\n\n")

# ============================================================================
# Outlier Detection
# ============================================================================

cat("Detecting outliers...\n")

outlier_info <- list()
outlier_samples <- c()

# Calculate centroid for each group
for (grp in unique(groups)) {
  grp_idx <- which(groups == grp)
  grp_samples <- sample_names[grp_idx]

  if (length(grp_idx) > 1) {
    # Calculate group centroid
    centroid_pc1 <- mean(pca_df$PC1[grp_idx])
    centroid_pc2 <- mean(pca_df$PC2[grp_idx])

    # Calculate distances from centroid
    distances <- sqrt((pca_df$PC1[grp_idx] - centroid_pc1)^2 +
                      (pca_df$PC2[grp_idx] - centroid_pc2)^2)

    # Calculate SD of distances
    dist_mean <- mean(distances)
    dist_sd <- sd(distances)

    # Flag outliers (> threshold SDs from centroid)
    for (j in 1:length(grp_idx)) {
      idx <- grp_idx[j]
      z_score <- if(dist_sd > 0) (distances[j] - dist_mean) / dist_sd else 0

      is_outlier <- abs(z_score) > opt$outlier_threshold

      outlier_info[[sample_names[idx]]] <- list(
        sample = sample_names[idx],
        group = grp,
        distance_from_centroid = round(distances[j], 3),
        z_score = round(z_score, 3),
        is_outlier = is_outlier
      )

      if (is_outlier) {
        outlier_samples <- c(outlier_samples, sample_names[idx])
        cat("  WARNING: Potential outlier detected:", sample_names[idx],
            "(z-score:", round(z_score, 2), ")\n")
      }
    }
  }
}

if (length(outlier_samples) == 0) {
  cat("  No outliers detected (threshold:", opt$outlier_threshold, "SD)\n")
}
cat("\n")

# ============================================================================
# Batch Effect Detection
# ============================================================================

cat("Checking for batch effects...\n")

# Method 1: Check if samples cluster by group
# Calculate silhouette-like score (within-group vs between-group distances)

within_group_dist <- c()
between_group_dist <- c()

for (i in 1:n_samples) {
  for (j in 1:n_samples) {
    if (i < j) {
      dist_ij <- sqrt((pca_df$PC1[i] - pca_df$PC1[j])^2 +
                      (pca_df$PC2[i] - pca_df$PC2[j])^2)
      if (groups[i] == groups[j]) {
        within_group_dist <- c(within_group_dist, dist_ij)
      } else {
        between_group_dist <- c(between_group_dist, dist_ij)
      }
    }
  }
}

avg_within <- if(length(within_group_dist) > 0) mean(within_group_dist) else 0
avg_between <- if(length(between_group_dist) > 0) mean(between_group_dist) else 0

# Separation score: higher = better separation
separation_score <- if(avg_within > 0) avg_between / avg_within else Inf

# Check for batch effect indicators
# Batch effect suspected if:
# 1. Within-group distances are high (samples spread out)
# 2. PC1 doesn't correlate well with groups

# Calculate correlation between PC1 and group (as numeric)
group_numeric <- as.numeric(as.factor(groups))
pc1_group_cor <- abs(cor(pca_df$PC1, group_numeric))

batch_effect_suspected <- FALSE
batch_effect_reason <- ""

if (separation_score < 1.5) {
  batch_effect_suspected <- TRUE
  batch_effect_reason <- "Low separation between groups (samples spread within groups)"
} else if (pc1_group_cor < 0.5 && pca_var[1] > 30) {
  batch_effect_suspected <- TRUE
  batch_effect_reason <- "PC1 explains high variance but doesn't correlate with groups"
}

cat("  Avg within-group distance:", round(avg_within, 3), "\n")
cat("  Avg between-group distance:", round(avg_between, 3), "\n")
cat("  Separation score:", round(separation_score, 3), "(higher = better)\n")
cat("  PC1-Group correlation:", round(pc1_group_cor, 3), "\n")

if (batch_effect_suspected) {
  cat("\n  WARNING: Potential BATCH EFFECT detected!\n")
  cat("  Reason:", batch_effect_reason, "\n")
  cat("  Recommendation: Use batch_effect_script_edgeR_v2.R for DE analysis\n")
} else {
  cat("  No batch effect detected.\n")
}
cat("\n")

# ============================================================================
# Generate PCA Plot
# ============================================================================

cat("Generating PCA plot...\n")

# Add outlier indicator to plot
pca_df$Outlier <- ifelse(pca_df$Sample %in% outlier_samples, "Outlier", "Normal")

pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, shape = Outlier)) +
  geom_point(size = 4) +
  geom_text(aes(label = Sample), vjust = -1, size = 3) +
  scale_shape_manual(values = c("Normal" = 16, "Outlier" = 17)) +
  labs(
    title = "PCA Plot - Sample Quality Assessment",
    subtitle = ifelse(batch_effect_suspected,
                      "WARNING: Potential batch effect detected!",
                      ifelse(length(outlier_samples) > 0,
                             paste("WARNING: Outliers detected:", paste(outlier_samples, collapse = ", ")),
                             "No outliers or batch effects detected")),
    x = paste0("PC1 (", pca_var[1], "%)"),
    y = paste0("PC2 (", pca_var[2], "%)")
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10,
                                 color = ifelse(batch_effect_suspected | length(outlier_samples) > 0,
                                               "red", "darkgreen"))
  )

ggsave(file.path(opt$outdir, "PCA_plot.png"), pca_plot, width = 10, height = 8, dpi = 150)
cat("  Saved:", file.path(opt$outdir, "PCA_plot.png"), "\n")

# ============================================================================
# Generate MDS Plot
# ============================================================================

cat("Generating MDS plot...\n")

# Classical MDS
dist_matrix <- dist(t(data_pca))
mds_result <- cmdscale(dist_matrix, k = 2)

mds_df <- data.frame(
  Dim1 = mds_result[, 1],
  Dim2 = mds_result[, 2],
  Sample = sample_names,
  Group = groups,
  Outlier = ifelse(sample_names %in% outlier_samples, "Outlier", "Normal")
)

mds_plot <- ggplot(mds_df, aes(x = Dim1, y = Dim2, color = Group, shape = Outlier)) +
  geom_point(size = 4) +
  geom_text(aes(label = Sample), vjust = -1, size = 3) +
  scale_shape_manual(values = c("Normal" = 16, "Outlier" = 17)) +
  labs(
    title = "MDS Plot - Sample Similarity",
    x = "Dimension 1",
    y = "Dimension 2"
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

ggsave(file.path(opt$outdir, "MDS_plot.png"), mds_plot, width = 10, height = 8, dpi = 150)
cat("  Saved:", file.path(opt$outdir, "MDS_plot.png"), "\n")

# ============================================================================
# Generate Density Plot
# ============================================================================

cat("Generating density plot...\n")

# Reshape for density plot
density_data <- stack(as.data.frame(data_filtered))
colnames(density_data) <- c("Expression", "Sample")
density_data$Group <- rep(groups, each = nrow(data_filtered))

density_plot <- ggplot(density_data, aes(x = Expression, color = Sample)) +
  geom_density(linewidth = 1) +
  labs(
    title = "Expression Density Plot",
    subtitle = "Distribution of expression values per sample",
    x = "Log2 Expression (RPKM)",
    y = "Density"
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
  xlim(0, quantile(density_data$Expression, 0.99))

ggsave(file.path(opt$outdir, "density_plot.png"), density_plot, width = 10, height = 8, dpi = 150)
cat("  Saved:", file.path(opt$outdir, "density_plot.png"), "\n\n")

# ============================================================================
# Generate Summary for Agents
# ============================================================================

cat("Generating summary for agent interpretation...\n")

summary_data <- list(
  input_file = opt$infile,
  n_samples = n_samples,
  n_genes = n_genes,
  n_genes_after_filter = nrow(data_filtered),

  samples = lapply(1:n_samples, function(i) {
    list(
      name = sample_names[i],
      group = groups[i],
      pc1 = round(pca_df$PC1[i], 3),
      pc2 = round(pca_df$PC2[i], 3)
    )
  }),

  pca_variance = list(
    PC1 = pca_var[1],
    PC2 = pca_var[2],
    PC3 = if(length(pca_var) >= 3) pca_var[3] else NA
  ),

  outlier_detection = list(
    threshold_sd = opt$outlier_threshold,
    outliers_detected = length(outlier_samples) > 0,
    outlier_samples = outlier_samples,
    outlier_details = outlier_info
  ),

  batch_effect_detection = list(
    suspected = batch_effect_suspected,
    reason = batch_effect_reason,
    metrics = list(
      avg_within_group_distance = round(avg_within, 3),
      avg_between_group_distance = round(avg_between, 3),
      separation_score = round(separation_score, 3),
      pc1_group_correlation = round(pc1_group_cor, 3)
    ),
    recommendation = ifelse(batch_effect_suspected,
                           "Use batch_effect_script_edgeR_v2.R",
                           "Use simpleEdger3.R")
  ),

  plots_generated = list(
    pca = file.path(opt$outdir, "PCA_plot.png"),
    mds = file.path(opt$outdir, "MDS_plot.png"),
    density = file.path(opt$outdir, "density_plot.png")
  ),

  agent_decision_required = list(
    outlier_removal = length(outlier_samples) > 0,
    batch_correction = batch_effect_suspected
  )
)

# Write JSON summary (for agents)
json_file <- file.path(opt$outdir, "qc_summary.json")
write_json(summary_data, json_file, pretty = TRUE, auto_unbox = TRUE)
cat("  Saved:", json_file, "\n")

# Write human-readable summary
txt_file <- file.path(opt$outdir, "qc_summary.txt")
sink(txt_file)
cat("========================================\n")
cat("QC SUMMARY REPORT\n")
cat("========================================\n\n")
cat("Input:", opt$infile, "\n")
cat("Samples:", n_samples, "\n")
cat("Genes (filtered):", nrow(data_filtered), "/", n_genes, "\n\n")

cat("--- SAMPLE GROUPS ---\n")
for (i in 1:n_samples) {
  cat(sample_names[i], "->", groups[i], "\n")
}

cat("\n--- PCA RESULTS ---\n")
cat("PC1:", pca_var[1], "% variance\n")
cat("PC2:", pca_var[2], "% variance\n")

cat("\n--- OUTLIER DETECTION ---\n")
if (length(outlier_samples) > 0) {
  cat("OUTLIERS DETECTED:\n")
  for (s in outlier_samples) {
    info <- outlier_info[[s]]
    cat("  -", s, "(z-score:", info$z_score, ")\n")
  }
  cat("\nACTION REQUIRED: Agents must vote (UNANIMOUS) to remove outliers.\n")
} else {
  cat("No outliers detected.\n")
}

cat("\n--- BATCH EFFECT DETECTION ---\n")
cat("Separation score:", round(separation_score, 3), "\n")
cat("PC1-Group correlation:", round(pc1_group_cor, 3), "\n")
if (batch_effect_suspected) {
  cat("\nBATCH EFFECT SUSPECTED!\n")
  cat("Reason:", batch_effect_reason, "\n")
  cat("\nACTION REQUIRED: Use batch_effect_script_edgeR_v2.R instead of simpleEdger3.R\n")
} else {
  cat("No batch effect detected.\n")
}

cat("\n--- RECOMMENDATION ---\n")
if (length(outlier_samples) > 0) {
  cat("1. Review PCA plot to confirm outliers\n")
  cat("2. If confirmed, remove:", paste(outlier_samples, collapse = ", "), "\n")
  cat("3. Re-run pipeline from featureCounts step\n")
} else if (batch_effect_suspected) {
  cat("1. Review PCA plot to confirm batch effect\n")
  cat("2. Use batch_effect_script_edgeR_v2.R for DE analysis\n")
} else {
  cat("Proceed with standard pipeline (simpleEdger3.R)\n")
}

cat("\n========================================\n")
cat("END OF REPORT\n")
cat("========================================\n")
sink()

cat("  Saved:", txt_file, "\n\n")

# ============================================================================
# Final Summary
# ============================================================================

cat("========================================\n")
cat("QC ANALYSIS COMPLETE\n")
cat("========================================\n\n")

cat("Output files:\n")
cat("  - PCA_plot.png\n")
cat("  - MDS_plot.png\n")
cat("  - density_plot.png\n")
cat("  - qc_summary.json (for agent parsing)\n")
cat("  - qc_summary.txt (human readable)\n\n")

if (length(outlier_samples) > 0 || batch_effect_suspected) {
  cat("ATTENTION REQUIRED:\n")
  if (length(outlier_samples) > 0) {
    cat("  - Outliers detected:", paste(outlier_samples, collapse = ", "), "\n")
  }
  if (batch_effect_suspected) {
    cat("  - Batch effect suspected\n")
  }
  cat("\nAgents should review the PCA plot and vote on action.\n")
} else {
  cat("All samples look good! Proceed with DE analysis.\n")
}

cat("\n")
