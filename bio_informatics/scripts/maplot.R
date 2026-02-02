#!/usr/bin/Rscript

args <- commandArgs(TRUE)
if (length(args) < 1) {
  stop("\nUsage: maplot.R edger3.txt cpmCutoff fcCutoff fdrCutoff x-range y-range R\nIf no values are supplied the following defaults are use:\nCPM: 1\nFC: 1\nFDR: 0.05\nX: 14\nY: default")
}
suppressPackageStartupMessages({
  library(affy)
})
# t[, 4] is logCPM
# t[, 3] is logFC
fname <- args[1]
cpm <- 1
fc <- 1
fdr <- 0.05
xlimit <- 14
ylimit <- NULL
if (length(args) >= 2) {
  cpm <- as.numeric(args[2])
}
if (length(args) >= 3) {
  fc <- as.numeric(args[3])
}
if (length(args) >= 4) {
  fdr <- as.numeric(args[4])
}
if (length(args) >= 5) {
  xlimit <- as.numeric(args[5])
}

t <- read.csv(args[1], header = TRUE)

ylimit <- range(t[, 3])

if (length(args) >= 6) {
  ylimit <- as.numeric(args[6])
  ylimit <- c(-ylimit, ylimit)
}


if (length(args) >= 7) {
  print("Up and down regulated are reversed.")
  t[, 3] <- -t[, 3]
}
fname <- paste0(args[1], "MaPlot.pdf")
pdf(fname)
smoothScatter(t[, 4], t[, 3], xlab = "logCPM", ylab = "logFC", xlim = c(cpm, xlimit), ylim = ylimit, pch = "")
up <- t[t[, 3] > fc & t[, 4] > cpm & t[, 7] < fdr, ]
dn <- t[t[, 3] < -fc & t[, 4] > cpm & t[, 7] < fdr, ]
points(up[, 4], up[, 3], pch = "+", col = "blue")
points(dn[, 4], dn[, 3], pch = "-", col = "red")

title(args[1])
l <- legend("topright", paste("logCPM >", cpm, "\n|logFC| >", fc, "\nFDR <", fdr, "\nUp-regulated:", length(up[, 1]), "\nDown-regulated:", length(dn[, 1])), bty = "n", cex = 1)
dev.off()
print(paste0("Results are saved as: ", fname))
