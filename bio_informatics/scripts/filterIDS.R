#!/usr/bin/Rscript
# Modified to read CSV format instead of space-delimited text

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: filterIDS.R counts.csv\n\nMake sure there is a column called 'ID'")
}

# Load bad IDs list (assumes bio_informatics folder structure)
script_dir <- dirname(normalizePath(sys.frame(1)$ofile))
bad_ids_path <- file.path(dirname(script_dir), "reference_data", "shared", "badIDS.txt")
if (!file.exists(bad_ids_path)) {
  stop("ERROR: Cannot find badIDS.txt at: ", bad_ids_path)
}
bad_IDS <- read.table(bad_ids_path, header = T)

# Read CSV input
counts <- read.csv(args[1], header = TRUE)

# Filter out bad IDs
`%notin%` <- Negate(`%in%`)
counts2 <- counts[counts$ID %notin% bad_IDS$GeneID, ]

# Generate output filename
output_file <- sub("\\.count\\.csv$", ".count.filtered.csv", args[1])
if (output_file == args[1]) {
  # If pattern didn't match, just append .filtered.csv
  output_file <- paste0(args[1], ".filtered.csv")
}

# Write CSV output
write.csv(counts2, file = output_file, row.names = FALSE, quote = FALSE)

cat("Filtering complete. Output written to:", output_file, "\n")
cat("Genes before filtering:", nrow(counts), "\n")
cat("Genes after filtering:", nrow(counts2), "\n")
cat("Genes removed:", nrow(counts) - nrow(counts2), "\n")
