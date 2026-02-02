#!/usr/bin/Rscript
# EDITED by Claude Code (2026-01-07):
#   - Read CSV input instead of space-delimited text
#   - Updated filename pattern matching to handle .csv extension

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("\nUsage: edger3xl.R edger3_results.csv\n")
}

library("openxlsx")

wb <- createWorkbook()
addWorksheet(wb, "Sheet 1")

# EDIT: Handle .csv extension instead of .txt
fname <- gsub(".csv", ".xlsx", args[1], fixed = TRUE)
if (fname == args[1]) {
  # Fallback if .csv not found
  fname <- paste0(args[1], ".xlsx")
}

# EDIT: Read CSV instead of space-delimited
edger <- read.csv(args[1], header = TRUE)

writeData(wb, "Sheet 1", edger)

# A   B      C        D     E      F      G
# "X" "id" "logFC" "logCPM" "F" "PValue" "FDR"

# Signficant = |logFC| > cutoff & logCPM > cutoff
s <- seq(1:length(edger[, 1])) + 1
expr <- paste0("=D", s, ">$L$2")
writeData(wb, "Sheet 1", "Expr", startCol = 8, startRow = 1)
writeFormula(wb, "Sheet 1", expr, startCol = 8, startRow = 2)

# Incr or decr: check logFC > 0
chg <- paste0("=IF(AND(C", s, "> $L$3, G", s, " < $L$1, H", s, "),\"incr\",IF(AND(C", s, "<-$L$3, G", s, " < $L$1, H", s, "),\"decr\",\"nchg\"))")
writeData(wb, "Sheet 1", "CG", startCol = 9, startRow = 1)
writeFormula(wb, "Sheet 1", chg, startCol = 9, startRow = 2)

writeData(wb, "Sheet 1", c("FDR", "logCPM", "logFC"), startCol = 11)
writeData(wb, "Sheet 1", c(0.05, 1, 1), startCol = 12)

writeData(wb, "Sheet 1", c("Incr", "Decr"), startCol = 11, startRow = 5)
writeFormula(wb, "Sheet 1", c("=COUNTIFS(I:I,K5)", "=COUNTIFS(I:I,K6)"), startCol = 12, startRow = 5)

saveWorkbook(wb, fname, overwrite = TRUE)

cat("Excel conversion complete. Output written to:", fname, "\n")
