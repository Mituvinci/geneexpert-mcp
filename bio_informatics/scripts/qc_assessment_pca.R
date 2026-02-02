#!/usr/bin/env Rscript

#' Comprehensive QC Assessment using PCA
#'
#' Detects:
#' - Within-group spread (hidden batch effects / technical heterogeneity)
#' - Outlier samples using robust median + MAD distance
#'
#' Usage:
#'   Rscript qc_assessment_pca.R <input_matrix_or_pca> <output_dir> <control_keyword> <treatment_keyword> [n_pcs]
#'
#' Arguments:
#'   input_matrix_or_pca: Path to normalized expression matrix (.csv) OR PCA scores (.csv)
#'   output_dir: Directory to save outputs
#'   control_keyword: Keyword in sample names for control group (e.g., "cont", "ctrl")
#'   treatment_keyword: Keyword in sample names for treatment group (e.g., "treat", "ips")
#'   n_pcs: Number of PCs to use (default: 10, will cap at available PCs)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  cat("Usage: Rscript qc_assessment_pca.R <input_matrix_or_pca> <output_dir> <control_keyword> <treatment_keyword> [n_pcs]\n")
  quit(status = 1)
}

input_file <- args[1]
comparison_name <- args[2]  # This is the comparison name (prefix for output files)
control_kw <- tolower(args[3])
treatment_kw <- tolower(args[4])
n_pcs <- if (length(args) >= 5) as.numeric(args[5]) else 10

# Note: Files will be saved in current directory with comparison_name prefix
# No subdirectory created - bash script already cd'd into stage3_quantification/

cat("=== QC Assessment using PCA ===\n")
cat("Input:", input_file, "\n")
cat("Comparison name:", comparison_name, "\n")
cat("Control keyword:", control_kw, "\n")
cat("Treatment keyword:", treatment_kw, "\n")
cat("Number of PCs:", n_pcs, "\n\n")

# Read input
data <- read.csv(input_file, row.names = 1, check.names = FALSE)

# Check if first column is "Length" (from featureCounts output) and skip it
if(ncol(data) > 0 && colnames(data)[1] == "Length") {
  cat("Detected 'Length' column - skipping for PCA analysis\n")
  data <- data[, -1, drop = FALSE]  # Remove first column (Length)
  cat("Using", ncol(data), "sample columns for analysis\n")
}

# Determine if input is expression matrix or PCA scores
is_pca_input <- any(grepl("^PC[0-9]+$", colnames(data)))

if (is_pca_input) {
  cat("Input detected as PCA scores CSV\n")
  pca_scores <- data
  sample_names <- rownames(pca_scores)
} else {
  cat("Input detected as expression matrix - computing PCA\n")

  # Transpose (genes as rows, samples as columns)
  # Assume input is samples x genes, transpose to genes x samples
  if (nrow(data) < ncol(data)) {
    cat("Transposing matrix (assuming samples as rows)\n")
    expr_matrix <- t(data)
  } else {
    expr_matrix <- data
  }

  # Log2 transform if not already (assume RPKM or similar)
  cat("Applying log2(x+1) transformation\n")
  expr_log <- log2(expr_matrix + 1)

  # Remove zero-variance genes
  gene_var <- apply(expr_log, 1, var)
  expr_log <- expr_log[gene_var > 0, ]

  cat("Running PCA on", nrow(expr_log), "genes and", ncol(expr_log), "samples\n")

  # Run PCA
  pca_result <- prcomp(t(expr_log), center = TRUE, scale. = FALSE)
  pca_scores <- pca_result$x
  sample_names <- rownames(pca_scores)

  # Save variance explained
  var_explained <- summary(pca_result)$importance[2, ]
  write.csv(var_explained, paste0(comparison_name, "_pca_variance_explained.csv"))
}

# Assign group labels based on keywords
assign_group <- function(sample_name, control_kw, treatment_kw) {
  sample_lower <- tolower(sample_name)
  if (grepl(control_kw, sample_lower)) {
    return("Control")
  } else if (grepl(treatment_kw, sample_lower)) {
    return("Treatment")
  } else {
    return("Unknown")
  }
}

group_labels <- sapply(sample_names, assign_group, control_kw, treatment_kw)

# ============================================================================
# CREATE SAMPLE NAME DICTIONARY (A1-An for control, B1-Bn for treatment)
# ============================================================================

cat("\n=== CREATING SAMPLE NAME DICTIONARY ===\n")

# Separate samples by group
control_samples <- sample_names[group_labels == "Control"]
treatment_samples <- sample_names[group_labels == "Treatment"]
unknown_samples <- sample_names[group_labels == "Unknown"]

# Create short labels
short_labels <- character(length(sample_names))
sample_dictionary <- data.frame(
  short_label = character(),
  full_name = character(),
  group = character(),
  stringsAsFactors = FALSE
)

# Assign A1, A2, A3... to control samples
for (i in seq_along(control_samples)) {
  idx <- which(sample_names == control_samples[i])
  short_label <- paste0("A", i)
  short_labels[idx] <- short_label
  sample_dictionary <- rbind(sample_dictionary, data.frame(
    short_label = short_label,
    full_name = control_samples[i],
    group = "Control",
    stringsAsFactors = FALSE
  ))
}

# Assign B1, B2, B3... to treatment samples
for (i in seq_along(treatment_samples)) {
  idx <- which(sample_names == treatment_samples[i])
  short_label <- paste0("B", i)
  short_labels[idx] <- short_label
  sample_dictionary <- rbind(sample_dictionary, data.frame(
    short_label = short_label,
    full_name = treatment_samples[i],
    group = "Treatment",
    stringsAsFactors = FALSE
  ))
}

# Handle unknown samples (if any) - assign C1, C2...
for (i in seq_along(unknown_samples)) {
  idx <- which(sample_names == unknown_samples[i])
  short_label <- paste0("C", i)
  short_labels[idx] <- short_label
  sample_dictionary <- rbind(sample_dictionary, data.frame(
    short_label = short_label,
    full_name = unknown_samples[i],
    group = "Unknown",
    stringsAsFactors = FALSE
  ))
}

# Save dictionary to file
dict_file <- paste0(comparison_name, "_sample_dictionary.txt")
sink(dict_file)
cat("======================================\n")
cat("SAMPLE NAME DICTIONARY\n")
cat("======================================\n\n")
cat("Use this dictionary to identify samples by their plot labels:\n\n")

if (length(control_samples) > 0) {
  cat(sprintf("Control Group (%s, n=%d):\n", control_kw, length(control_samples)))
  for (i in seq_along(control_samples)) {
    cat(sprintf("  A%d = %s\n", i, control_samples[i]))
  }
  cat("\n")
}

if (length(treatment_samples) > 0) {
  cat(sprintf("Treatment Group (%s, n=%d):\n", treatment_kw, length(treatment_samples)))
  for (i in seq_along(treatment_samples)) {
    cat(sprintf("  B%d = %s\n", i, treatment_samples[i]))
  }
  cat("\n")
}

if (length(unknown_samples) > 0) {
  cat(sprintf("Unknown Group (n=%d):\n", length(unknown_samples)))
  for (i in seq_along(unknown_samples)) {
    cat(sprintf("  C%d = %s\n", i, unknown_samples[i]))
  }
  cat("\n")
}
sink()
cat("  Saved:", dict_file, "\n")

# Also save as CSV for easy parsing
write.csv(sample_dictionary, paste0(comparison_name, "_sample_dictionary.csv"), row.names = FALSE)
cat("  Saved:", paste0(comparison_name, "_sample_dictionary.csv\n"))

# Cap n_pcs at available PCs
n_pcs <- min(n_pcs, ncol(pca_scores))
cat("Using top", n_pcs, "PCs for analysis\n\n")

pca_subset <- pca_scores[, 1:n_pcs, drop = FALSE]

# ============================================================================
# A) WITHIN-GROUP SPREAD ANALYSIS (Hidden Batch Effects)
# ============================================================================

cat("=== A) WITHIN-GROUP SPREAD ANALYSIS ===\n")

euclidean_dist <- function(x, y) {
  sqrt(sum((x - y)^2))
}

# Compute group centroids
groups <- unique(group_labels[group_labels != "Unknown"])
group_stats <- list()

for (g in groups) {
  samples_in_group <- group_labels == g
  group_pcs <- pca_subset[samples_in_group, , drop = FALSE]

  # Centroid
  centroid <- colMeans(group_pcs)

  # Distance from each sample to centroid
  dists <- apply(group_pcs, 1, function(x) euclidean_dist(x, centroid))

  # Group dispersion score (median distance to centroid)
  S_g <- median(dists)

  group_stats[[g]] <- list(
    centroid = centroid,
    dispersion = S_g,
    n_samples = sum(samples_in_group),
    dists = dists
  )

  cat(sprintf("  %s group: n=%d, dispersion=%.2f\n", g, sum(samples_in_group), S_g))
}

# Overall reference (all samples)
all_centroid <- colMeans(pca_subset)
all_dists <- apply(pca_subset, 1, function(x) euclidean_dist(x, all_centroid))
M <- median(all_dists)
MAD <- median(abs(all_dists - M))

cat(sprintf("\nOverall reference: median=%.2f, MAD=%.2f\n", M, MAD))

# Flag groups with excessive spread
spread_threshold <- M + 2 * MAD
spread_flags <- list()
for (g in groups) {
  is_spread <- group_stats[[g]]$dispersion > spread_threshold
  spread_flags[[g]] <- is_spread
  if (is_spread) {
    cat(sprintf("  WARNING: %s group is TOO SPREAD (%.2f > %.2f)\n",
                g, group_stats[[g]]$dispersion, spread_threshold))
  }
}

# Spread imbalance ratio
if (length(groups) == 2) {
  dispersions <- sapply(group_stats, function(x) x$dispersion)
  R <- max(dispersions) / min(dispersions)
  cat(sprintf("\nSpread imbalance ratio: %.2f", R))
  if (R > 1.5) {
    cat(" (IMBALANCED - possible batch effect)\n")
  } else {
    cat(" (balanced)\n")
  }
} else {
  R <- NA
}

# ============================================================================
# B) OUTLIER DETECTION (Robust distance within group)
# ============================================================================

cat("\n=== B) OUTLIER DETECTION ===\n")

outlier_results <- data.frame(
  sample_id = sample_names,
  group_label = group_labels,
  dist_to_group_centroid = NA,
  mean_dist_to_own_group = NA,
  outlier_flag = FALSE,
  stringsAsFactors = FALSE
)

for (g in groups) {
  samples_in_group <- which(group_labels == g)
  n_g <- length(samples_in_group)

  cat(sprintf("\nGroup: %s (n=%d)\n", g, n_g))

  if (n_g < 3) {
    cat("  Too few replicates (<3) for robust outlier detection\n")
    outlier_results$dist_to_group_centroid[samples_in_group] <- group_stats[[g]]$dists
    next
  }

  group_pcs <- pca_subset[samples_in_group, , drop = FALSE]

  # For each sample in group, compute mean distance to all other samples in group
  d_i_g <- numeric(n_g)
  for (i in 1:n_g) {
    sample_i <- group_pcs[i, ]
    dists_to_others <- apply(group_pcs[-i, , drop = FALSE], 1,
                              function(x) euclidean_dist(sample_i, x))
    d_i_g[i] <- mean(dists_to_others)
  }

  # Robust outlier threshold: median + 3*MAD
  m_g <- median(d_i_g)
  mad_g <- median(abs(d_i_g - m_g))
  threshold <- m_g + 3 * mad_g

  outliers <- d_i_g > threshold

  # Store results
  outlier_results$dist_to_group_centroid[samples_in_group] <- group_stats[[g]]$dists
  outlier_results$mean_dist_to_own_group[samples_in_group] <- d_i_g
  outlier_results$outlier_flag[samples_in_group] <- outliers

  cat(sprintf("  Median distance: %.2f, MAD: %.2f, Threshold: %.2f\n",
              m_g, mad_g, threshold))

  if (any(outliers)) {
    outlier_samples <- sample_names[samples_in_group[outliers]]
    cat(sprintf("  OUTLIERS DETECTED: %s\n", paste(outlier_samples, collapse=", ")))
  } else {
    cat("  No outliers detected\n")
  }
}

# ============================================================================
# OUTPUTS
# ============================================================================

cat("\n=== SAVING OUTPUTS ===\n")

# 1) PCA plot (PC1 vs PC2)
# Function to generate the plot (to avoid code duplication)
plot_pca <- function() {
  # Use SHORT labels (A1, B1, etc.) - dictionary maps to full names
  plot_labels <- short_labels

  # Calculate data ranges for PC1 and PC2
  x_range <- range(pca_scores[, 1])
  y_range <- range(pca_scores[, 2])

  # Add 50% buffer on each side to prevent label clipping at boundaries
  x_buffer <- diff(x_range) * 0.5
  y_buffer <- diff(y_range) * 0.5

  # Expanded axis limits
  xlim <- c(x_range[1] - x_buffer, x_range[2] + x_buffer)
  ylim <- c(y_range[1] - y_buffer, y_range[2] + y_buffer)

  # Create plot with larger bottom margin for legend
  par(mar = c(9, 4, 4, 2) + 0.1)  # bottom, left, top, right

  plot(pca_scores[, 1], pca_scores[, 2],
       col = ifelse(group_labels == "Control", "blue", "red"),  # Only blue and red
       pch = 19, cex = 2,
       xlab = "PC1", ylab = "PC2",
       main = "PCA Plot - QC Assessment",
       xlim = xlim,  # Expanded X-axis limits (50% buffer)
       ylim = ylim)  # Expanded Y-axis limits (50% buffer)

  # Add SHORT sample labels LEFT of each point (A1, B1, etc.)
  text(pca_scores[, 1], pca_scores[, 2],
       labels = plot_labels,
       pos = 2,  # Position: 1=below, 2=left, 3=above, 4=right
       cex = 1.6,  # Larger font for readability
       col = ifelse(group_labels == "Control", "darkblue", "darkred"))  # Only blue and red

  # No X marks for outliers - agents will decide based on plot

  # Count samples per group
  control_count <- sum(group_labels == "Control")
  treatment_count <- sum(group_labels == "Treatment")

  # Legend outside plot box, horizontal, below X-axis with better spacing
  legend("bottom",
         legend = c(
           paste0(control_kw, " (n=", control_count, ")"),
           paste0(treatment_kw, " (n=", treatment_count, ")")
         ),
         col = c("blue", "red"),
         pch = c(19, 19),
         horiz = TRUE,  # Horizontal layout
         xpd = TRUE,  # Allow legend outside plot region
         inset = c(0, -0.22),  # Move further below X-axis to avoid overlap
         cex = 0.9,
         bty = "n")  # No box around legend
}

# Save as PDF (for user - no _qc suffix) - standard size
pdf(paste0(comparison_name, "_pca_plot.pdf"), width = 8, height = 6)
plot_pca()
dev.off()
cat("  Saved:", paste0(comparison_name, "_pca_plot.pdf\n"))

# Save as JPEG (for LLM agents - GPT-5.2 and Claude require JPEG/PNG, with _qc suffix)
# Standard size since we use short labels now
jpeg(paste0(comparison_name, "_pca_plot_qc.jpg"), width = 800, height = 600, quality = 95)
plot_pca()
dev.off()
cat("  Saved:", paste0(comparison_name, "_pca_plot_qc.jpg\n"))

# 2) CSV report
write.csv(outlier_results, paste0(comparison_name, "_qc_assessment_report.csv"),
          row.names = FALSE)
cat("  Saved:", paste0(comparison_name, "_qc_assessment_report.csv\n"))

# 3) QC Summary text
summary_file <- paste0(comparison_name, "_qc_summary.txt")
sink(summary_file)

cat("======================================\n")
cat("QC ASSESSMENT SUMMARY\n")
cat("======================================\n\n")

cat("DATASET OVERVIEW:\n")
cat(sprintf("  Total samples: %d\n", nrow(pca_subset)))
for (g in groups) {
  cat(sprintf("  %s: n=%d\n", g, group_stats[[g]]$n_samples))
}

cat("\nWITHIN-GROUP SPREAD (Batch Effect Check):\n")
for (g in groups) {
  status <- if (spread_flags[[g]]) "TOO SPREAD (possible batch effect)" else "Normal"
  cat(sprintf("  %s: dispersion=%.2f [%s]\n", g, group_stats[[g]]$dispersion, status))
}

if (!is.na(R)) {
  status <- if (R > 1.5) "IMBALANCED (concern)" else "Balanced"
  cat(sprintf("  Spread ratio: %.2f [%s]\n", R, status))
}

cat("\nOUTLIER DETECTION:\n")
n_outliers <- sum(outlier_results$outlier_flag)
if (n_outliers > 0) {
  cat(sprintf("  %d outlier(s) detected:\n", n_outliers))
  outlier_samples <- outlier_results$sample_id[outlier_results$outlier_flag]
  for (s in outlier_samples) {
    cat(sprintf("    - %s (%s)\n", s,
                outlier_results$group_label[outlier_results$sample_id == s]))
  }
} else {
  cat("  No outliers detected\n")
}

cat("\nRECOMMENDATIONS:\n")
if (any(unlist(spread_flags)) || n_outliers > 0) {
  cat("  ⚠ QUALITY CONCERNS DETECTED\n")
  if (any(unlist(spread_flags))) {
    cat("    - High within-group spread suggests hidden batch effects\n")
    cat("    - Consider using batch-aware DE analysis (batch_effect_edger.R)\n")
  }
  if (n_outliers > 0) {
    cat("    - Outlier samples detected\n")
    cat("    - Consult with analysts: remove or keep these samples?\n")
  }
} else {
  cat("  ✓ No major quality concerns\n")
  cat("    - Standard DE analysis can proceed (simpleEdger3.R)\n")
}

sink()
cat("  Saved: qc_summary.txt\n")

# Print summary to console
cat("\n")
cat(readLines(summary_file), sep = "\n")

cat("\n=== QC ASSESSMENT COMPLETE ===\n")

# Exit with status indicating QC issues
if (any(unlist(spread_flags)) || n_outliers > 0) {
  quit(status = 2)  # QC issues detected
} else {
  quit(status = 0)  # All good
}
