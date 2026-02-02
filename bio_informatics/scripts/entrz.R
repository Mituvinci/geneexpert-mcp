#!/usr/bin/Rscript
# Modified to read CSV format instead of space-delimited text

args = commandArgs(trailingOnly = TRUE)

if (length(args)<2) {
  stop("Usage: entrz.R input.rpkm.csv mm10/hg38 [outputfile.csv]", call.=FALSE)
} else if (length(args)==2) {
  # Default output file based on input filename
  output <- sub("\\.rpkm\\.csv$", ".rpkm.entrz.csv", args[1])
  if (output == args[1]) {
    # Fallback if pattern didn't match
    output = "outentrz.csv"
  }
} else {
  output = args[3]
}

# Read CSV input
table <- read.csv(args[1], header = TRUE)

# Read gene symbol / Entrez ID mapping table (relative path from scripts/)
script_dir <- dirname(normalizePath(sys.frame(1)$ofile))
filename <- file.path(script_dir, "..", "reference_data", "shared", paste0(args[2], "_entrzID_GS.txt"))
if (!file.exists(filename)) {
  stop("ERROR: Cannot find reference file: ", filename)
}
ids <- read.table(filename, header = TRUE, sep = "\t", quote = "")

# TODO R builtin? binary search? hash table?
ids$GS = as.character(ids$GS)

# Match Entrez IDs and add gene symbols
index = match(table[,1], ids[,1])

gene_symbol = vector()
for(i in index)
{
   gene_symbol <- append(gene_symbol, ids[i,"GS"])
}

Entrzid = table[,1]
table = cbind(Entrzid, gene_symbol, table[,-1])

cat("Adding gene symbols:\n")
cat("  Input genes:", nrow(table), "\n")
cat("  Matched symbols:", sum(!is.na(gene_symbol)), "\n")
cat("  Unmatched:", sum(is.na(gene_symbol)), "\n")

# EDIT: Write CSV instead of space-delimited, with correct syntax
write.csv(table, file = output, row.names = FALSE, quote = FALSE)

cat("Gene symbol annotation complete. Output written to:", output, "\n")
