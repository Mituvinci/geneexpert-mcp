#!/usr/bin/Rscript
# EDITED by Claude Code (2026-01-17):
#   - Modified to match paired version interface (5 args)
#   - Clean column names (basename only, no paths, no .bam extension)
#   - Output as CSV instead of space-delimited text
#   - Fixed output filename to use output_basename
#   - Sort BAMs by control/treatment keywords (controls first, treatments second)
#   - isPairedEnd=FALSE for single-end reads

args = commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
  stop("Usage: featurecounts_unpaired.R annotation output_basename control_keyword treatment_keyword [1.bam 2.bam ...]\n  Example: featurecounts_unpaired.R mm10 DA0036_stroke_vs_control cont ips /path/to/*.bam", call.=FALSE)
}

# Extract parameters
annotation <- args[1]
output_basename <- args[2]
control_keyword <- args[3]
treatment_keyword <- args[4]
# BAM files start from arg 5
bam_files <- args[5:length(args)]

# Sort BAM files: controls first, then treatments
cat("Sorting BAM files by group...\n")
control_bams <- bam_files[grepl(control_keyword, bam_files, ignore.case = TRUE)]
treatment_bams <- bam_files[grepl(treatment_keyword, bam_files, ignore.case = TRUE)]

# Check if sorting worked
if (length(control_bams) == 0) {
  stop("No control BAMs found with keyword: ", control_keyword, "\n")
}
if (length(treatment_bams) == 0) {
  stop("No treatment BAMs found with keyword: ", treatment_keyword, "\n")
}

cat("  Control BAMs (", length(control_bams), "):\n", sep = "")
for (bam in control_bams) {
  cat("    -", basename(bam), "\n")
}
cat("  Treatment BAMs (", length(treatment_bams), "):\n", sep = "")
for (bam in treatment_bams) {
  cat("    -", basename(bam), "\n")
}

# Concatenate: controls first, treatments second
bam_files_sorted <- c(control_bams, treatment_bams)

library("Rsubread")
library("limma")

# Run featureCounts on sorted BAM files (UNPAIRED / single-end)
cat("\nRunning featureCounts (UNPAIRED) on", length(bam_files_sorted), "BAM files...\n")
fc <- featureCounts(bam_files_sorted, annot.inbuilt=annotation, nthreads=40, isPairedEnd=FALSE)

# Clean BAM filenames - remove path and .bam extension
clean_names <- sapply(bam_files_sorted, function(x) {
  # Get basename (removes path)
  bname <- basename(x)
  # Remove .bam extension
  sub("\\.bam$", "", bname)
})

# Create column labels with clean names
lbl <- c("ID", "Length", clean_names)

# Create output data frame
output_df <- data.frame(
  ID = fc$annotation$GeneID,
  Length = fc$annotation$Length,
  fc$counts
)
colnames(output_df) <- lbl

# Use provided output basename
output_file <- paste0(output_basename, ".count.csv")

# Write as CSV instead of space-delimited
write.csv(output_df, file = output_file, row.names = FALSE, quote = FALSE)

cat("\n✓ FeatureCounts (UNPAIRED) complete!\n")
cat("  Output file:", output_file, "\n")
cat("  Total genes:", nrow(output_df), "\n")
cat("  Column order: ID, Length,", length(control_bams), "controls,", length(treatment_bams), "treatments\n")
cat("  Controls first, then treatments (as required for simpleEdger3.R)\n")
