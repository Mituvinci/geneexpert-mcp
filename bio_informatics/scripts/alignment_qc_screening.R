#!/usr/bin/env Rscript

#' Alignment-Based Contamination/Quality Screening
#'
#' Parses Subread alignment log files to detect possible contamination,
#' wrong reference genome, or poor library quality.
#'
#' Usage:
#'   Rscript alignment_qc_screening.R <bam_dir> <output_dir>
#'
#' Arguments:
#'   bam_dir: Directory containing .log files from subread-align
#'   output_dir: Directory to save screening reports
#'
#' Outputs:
#'   - alignment_qc_summary.csv: All samples with metrics and status
#'   - alignment_qc_report.txt: Detailed report of WARN/FAIL samples
#'
#' Exit codes:
#'   0 = All samples PASS
#'   1 = Some samples WARN
#'   2 = Some samples FAIL or SEVERE_FAIL

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  cat("Usage: Rscript alignment_qc_screening.R <bam_dir> <output_dir>\n")
  quit(status = 1)
}

bam_dir <- args[1]
output_dir <- args[2]

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== Alignment QC Screening ===\n")
cat("BAM directory:", bam_dir, "\n")
cat("Output directory:", output_dir, "\n\n")

# ============================================================================
# Helper: Parse Subread log file
# ============================================================================

parse_subread_log <- function(log_file) {
  # Read file
  lines <- readLines(log_file, warn = FALSE)

  # Find Summary block
  summary_start <- grep("^//=+\\s+Summary\\s+=+", lines)
  if (length(summary_start) == 0) {
    return(list(
      sample = basename(log_file),
      total_fragments = NA,
      mapped = NA,
      unmapped = NA,
      properly_paired = NA,
      not_properly_paired = NA,
      singleton = NA,
      chimeric = NA,
      unexpected_fraglen = NA,
      is_paired_end = NA,
      status = "INSUFFICIENT_DATA"
    ))
  }

  # Extract ~25 lines after Summary header
  summary_lines <- lines[summary_start:(summary_start + 25)]
  summary_text <- paste(summary_lines, collapse = "\n")

  # Parse key metrics (handle variable whitespace)
  extract_metric <- function(pattern, text) {
    match <- regmatches(text, regexpr(pattern, text, perl = TRUE))
    if (length(match) == 0) return(NA)
    # Extract number (first capturing group or first number)
    num <- as.numeric(gsub("[^0-9]", "", strsplit(match, ":")[[1]][2]))
    return(num)
  }

  # FIX: For single-end, Subread uses "Total reads" not "Total fragments"
  total_fragments <- extract_metric("Total fragments\\s*:\\s*([0-9,]+)", summary_text)
  if (is.na(total_fragments)) {
    total_fragments <- extract_metric("Total reads\\s*:\\s*([0-9,]+)", summary_text)
  }

  mapped <- extract_metric("Mapped\\s*:\\s*([0-9,]+)", summary_text)
  unmapped <- extract_metric("Unmapped\\s*:\\s*([0-9,]+)", summary_text)

  # FIX: If still NA, calculate from mapped + unmapped
  if (is.na(total_fragments) && !is.na(mapped) && !is.na(unmapped)) {
    total_fragments <- mapped + unmapped
  }
  properly_paired <- extract_metric("Properly paired\\s*:\\s*([0-9,]+)", summary_text)
  not_properly_paired <- extract_metric("Not properly paired\\s*:\\s*([0-9,]+)", summary_text)
  singleton <- extract_metric("Singleton\\s*:\\s*([0-9,]+)", summary_text)
  chimeric <- extract_metric("Chimeric\\s*:\\s*([0-9,]+)", summary_text)
  unexpected_fraglen <- extract_metric("Unexpected.*length\\s*:\\s*([0-9,]+)", summary_text)

  # Auto-detect paired-end (check if Input file 2 exists in settings)
  is_paired_end <- any(grepl("Input file 2\\s*:", lines))

  return(list(
    sample = gsub("\\.log$", "", basename(log_file)),
    total_fragments = total_fragments,
    mapped = mapped,
    unmapped = unmapped,
    properly_paired = properly_paired,
    not_properly_paired = not_properly_paired,
    singleton = singleton,
    chimeric = chimeric,
    unexpected_fraglen = unexpected_fraglen,
    is_paired_end = is_paired_end,
    status = NA  # To be computed
  ))
}

# ============================================================================
# Helper: Compute screening status
# ============================================================================

compute_status <- function(metrics) {
  total <- metrics$total_fragments

  # Check for insufficient data
  if (is.na(total) || total < 100) {
    return(list(
      mapped_pct = NA,
      properly_paired_pct = NA,
      singleton_pct = NA,
      chimeric_pct = NA,
      unexpected_fraglen_pct = NA,
      status = "INSUFFICIENT_DATA",
      reasons = "Total fragments < 100 or missing"
    ))
  }

  # Compute percentages
  mapped_pct <- (metrics$mapped / total) * 100
  properly_paired_pct <- if (!is.na(metrics$properly_paired)) (metrics$properly_paired / total) * 100 else NA
  singleton_pct <- if (!is.na(metrics$singleton)) (metrics$singleton / total) * 100 else NA
  chimeric_pct <- if (!is.na(metrics$chimeric)) (metrics$chimeric / total) * 100 else NA
  unexpected_fraglen_pct <- if (!is.na(metrics$unexpected_fraglen)) (metrics$unexpected_fraglen / total) * 100 else NA

  # Apply thresholds
  reasons <- c()

  # Primary: mapping rate
  if (is.na(mapped_pct)) {
    status <- "INSUFFICIENT_DATA"
    reasons <- c(reasons, "Mapping rate missing")
  } else if (mapped_pct < 60) {
    status <- "SEVERE_FAIL"
    reasons <- c(reasons, sprintf("Mapping rate critically low: %.1f%% < 60%%", mapped_pct))
  } else if (mapped_pct < 70) {
    status <- "FAIL"
    reasons <- c(reasons, sprintf("Mapping rate low: %.1f%% < 70%%", mapped_pct))
  } else if (mapped_pct < 80) {
    status <- "WARN"
    reasons <- c(reasons, sprintf("Mapping rate borderline: %.1f%% < 80%%", mapped_pct))
  } else {
    status <- "PASS"
  }

  # Secondary checks (paired-end only)
  if (metrics$is_paired_end && status != "INSUFFICIENT_DATA") {
    # Properly paired
    if (!is.na(properly_paired_pct)) {
      if (properly_paired_pct < 40) {
        status <- ifelse(status %in% c("SEVERE_FAIL", "FAIL"), status, "FAIL")
        reasons <- c(reasons, sprintf("Properly paired critically low: %.1f%% < 40%%", properly_paired_pct))
      } else if (properly_paired_pct < 60) {
        status <- ifelse(status %in% c("SEVERE_FAIL", "FAIL"), status, "WARN")
        reasons <- c(reasons, sprintf("Properly paired low: %.1f%% < 60%%", properly_paired_pct))
      }
    }

    # Singleton
    if (!is.na(singleton_pct)) {
      if (singleton_pct > 20) {
        status <- ifelse(status %in% c("SEVERE_FAIL", "FAIL"), status, "FAIL")
        reasons <- c(reasons, sprintf("Singleton rate high: %.1f%% > 20%%", singleton_pct))
      } else if (singleton_pct > 10) {
        status <- ifelse(status %in% c("SEVERE_FAIL", "FAIL"), status, "WARN")
        reasons <- c(reasons, sprintf("Singleton rate elevated: %.1f%% > 10%%", singleton_pct))
      }
    }

    # Chimeric
    if (!is.na(chimeric_pct)) {
      if (chimeric_pct > 5) {
        status <- ifelse(status %in% c("SEVERE_FAIL", "FAIL"), status, "FAIL")
        reasons <- c(reasons, sprintf("Chimeric rate high: %.1f%% > 5%%", chimeric_pct))
      } else if (chimeric_pct > 2) {
        status <- ifelse(status %in% c("SEVERE_FAIL", "FAIL"), status, "WARN")
        reasons <- c(reasons, sprintf("Chimeric rate elevated: %.1f%% > 2%%", chimeric_pct))
      }
    }

    # Unexpected fragment length
    if (!is.na(unexpected_fraglen_pct)) {
      if (unexpected_fraglen_pct > 35) {
        status <- ifelse(status %in% c("SEVERE_FAIL", "FAIL"), status, "FAIL")
        reasons <- c(reasons, sprintf("Unexpected fragment length high: %.1f%% > 35%%", unexpected_fraglen_pct))
      } else if (unexpected_fraglen_pct > 20) {
        status <- ifelse(status %in% c("SEVERE_FAIL", "FAIL"), status, "WARN")
        reasons <- c(reasons, sprintf("Unexpected fragment length elevated: %.1f%% > 20%%", unexpected_fraglen_pct))
      }
    }
  } else {
    # Single-end: only check chimeric if available
    if (!is.na(chimeric_pct)) {
      if (chimeric_pct > 5) {
        status <- ifelse(status %in% c("SEVERE_FAIL", "FAIL"), status, "FAIL")
        reasons <- c(reasons, sprintf("Chimeric rate high: %.1f%% > 5%%", chimeric_pct))
      } else if (chimeric_pct > 2) {
        status <- ifelse(status %in% c("SEVERE_FAIL", "FAIL"), status, "WARN")
        reasons <- c(reasons, sprintf("Chimeric rate elevated: %.1f%% > 2%%", chimeric_pct))
      }
    }
  }

  if (length(reasons) == 0 && status == "PASS") {
    reasons <- "All metrics within normal range"
  }

  return(list(
    mapped_pct = mapped_pct,
    properly_paired_pct = properly_paired_pct,
    singleton_pct = singleton_pct,
    chimeric_pct = chimeric_pct,
    unexpected_fraglen_pct = unexpected_fraglen_pct,
    status = status,
    reasons = paste(reasons, collapse = "; ")
  ))
}

# ============================================================================
# Main: Process all log files
# ============================================================================

log_files <- list.files(bam_dir, pattern = "\\.log$", full.names = TRUE)

if (length(log_files) == 0) {
  cat("ERROR: No .log files found in", bam_dir, "\n")
  quit(status = 1)
}

cat("Found", length(log_files), "log files\n")
cat("Processing...\n\n")

results <- list()

for (log_file in log_files) {
  cat("  Processing:", basename(log_file), "\n")

  # Parse log
  metrics <- parse_subread_log(log_file)

  # Compute status
  screening <- compute_status(metrics)

  # Combine (remove status from metrics to avoid overwriting)
  metrics$status <- NULL
  result <- c(metrics, screening)
  results[[length(results) + 1]] <- result
}

# Convert to data frame
results_df <- do.call(rbind, lapply(results, function(x) {
  data.frame(
    sample = x$sample,
    total_fragments = x$total_fragments,
    mapped = x$mapped,
    unmapped = x$unmapped,
    properly_paired = x$properly_paired,
    not_properly_paired = x$not_properly_paired,
    singleton = x$singleton,
    chimeric = x$chimeric,
    unexpected_fraglen = x$unexpected_fraglen,
    is_paired_end = x$is_paired_end,
    mapped_pct = x$mapped_pct,
    properly_paired_pct = x$properly_paired_pct,
    singleton_pct = x$singleton_pct,
    chimeric_pct = x$chimeric_pct,
    unexpected_fraglen_pct = x$unexpected_fraglen_pct,
    status = x$status,
    reasons = x$reasons,
    stringsAsFactors = FALSE
  )
}))

# ============================================================================
# Output: CSV Summary
# ============================================================================

csv_file <- file.path(output_dir, "alignment_qc_summary.csv")
write.csv(results_df, csv_file, row.names = FALSE)
cat("\nSaved CSV summary:", csv_file, "\n")

# ============================================================================
# Output: Text Report
# ============================================================================

report_file <- file.path(output_dir, "alignment_qc_report.txt")
sink(report_file)

cat("================================================================================\n")
cat("ALIGNMENT QC SCREENING REPORT\n")
cat("================================================================================\n\n")

cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Total samples:", nrow(results_df), "\n\n")

cat("IMPORTANT NOTE:\n")
cat("Low mapping rate suggests possible contamination OR wrong reference genome\n")
cat("OR poor library/sequencing quality. This is a screening flag, not definitive\n")
cat("proof of contamination. Review flagged samples carefully.\n\n")

cat("================================================================================\n")
cat("SUMMARY BY STATUS\n")
cat("================================================================================\n\n")

status_counts <- table(results_df$status)
for (s in c("PASS", "WARN", "FAIL", "SEVERE_FAIL", "INSUFFICIENT_DATA")) {
  count <- ifelse(s %in% names(status_counts), status_counts[s], 0)
  cat(sprintf("  %s: %d\n", s, count))
}

cat("\n")

# Report WARN/FAIL samples
flagged <- results_df[results_df$status %in% c("WARN", "FAIL", "SEVERE_FAIL", "INSUFFICIENT_DATA"), ]

if (nrow(flagged) > 0) {
  cat("================================================================================\n")
  cat("FLAGGED SAMPLES (WARN/FAIL/SEVERE_FAIL/INSUFFICIENT_DATA)\n")
  cat("================================================================================\n\n")

  for (i in 1:nrow(flagged)) {
    sample <- flagged[i, ]
    cat(sprintf("Sample: %s\n", sample$sample))
    cat(sprintf("  Status: %s\n", sample$status))
    cat(sprintf("  Total fragments: %s\n", ifelse(is.na(sample$total_fragments), "NA", format(sample$total_fragments, big.mark=","))))
    cat(sprintf("  Mapped: %.1f%%\n", ifelse(is.na(sample$mapped_pct), NA, sample$mapped_pct)))

    if (sample$is_paired_end && !is.na(sample$is_paired_end)) {
      cat(sprintf("  Properly paired: %.1f%%\n", ifelse(is.na(sample$properly_paired_pct), NA, sample$properly_paired_pct)))
      cat(sprintf("  Singleton: %.1f%%\n", ifelse(is.na(sample$singleton_pct), NA, sample$singleton_pct)))
      cat(sprintf("  Chimeric: %.1f%%\n", ifelse(is.na(sample$chimeric_pct), NA, sample$chimeric_pct)))
      cat(sprintf("  Unexpected fragment length: %.1f%%\n", ifelse(is.na(sample$unexpected_fraglen_pct), NA, sample$unexpected_fraglen_pct)))
    } else {
      cat(sprintf("  Chimeric: %.1f%%\n", ifelse(is.na(sample$chimeric_pct), NA, sample$chimeric_pct)))
    }

    cat(sprintf("  Reasons: %s\n", sample$reasons))
    cat("\n")
  }

  cat("RECOMMENDATION:\n")
  if (any(flagged$status == "SEVERE_FAIL")) {
    cat("  ⚠ SEVERE_FAIL detected: Review these samples carefully.\n")
    cat("    Consider removing from analysis or investigating root cause.\n")
  }
  if (any(flagged$status == "FAIL")) {
    cat("  ⚠ FAIL detected: Review QC metrics and consider sample exclusion.\n")
  }
  if (any(flagged$status == "WARN")) {
    cat("  ⚠ WARN detected: Samples are borderline. Proceed with caution.\n")
  }
} else {
  cat("================================================================================\n")
  cat("RESULT: ALL SAMPLES PASS\n")
  cat("================================================================================\n\n")
  cat("✓ All samples have acceptable alignment quality.\n")
}

sink()
cat("Saved text report:", report_file, "\n")

# ============================================================================
# Print summary to console
# ============================================================================

cat("\n")
cat(readLines(report_file), sep = "\n")

# ============================================================================
# Exit with appropriate code
# ============================================================================

cat("\n=== Alignment QC Screening Complete ===\n")

# Remove NA values before checking (defensive)
valid_statuses <- results_df$status[!is.na(results_df$status)]

if (any(valid_statuses %in% c("FAIL", "SEVERE_FAIL"))) {
  cat("Exit code: 2 (FAIL or SEVERE_FAIL detected)\n")
  quit(status = 2)
} else if (any(valid_statuses == "WARN")) {
  cat("Exit code: 1 (WARN detected)\n")
  quit(status = 1)
} else {
  cat("Exit code: 0 (All PASS)\n")
  quit(status = 0)
}
