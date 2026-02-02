#!/usr/bin/Rscript
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Usage: MDS.R count.txt [labels...]", call. = FALSE)
}

outfile <- paste(args[1], "MDS.pdf", sep = "")

library("Rsubread")
library("edgeR")
library("limma")

fc <- read.table(args[1], header = TRUE, sep = " ")
lbls <- args[2:length(args)]
grp <- factor(lbls)
y <- DGEList(counts = fc[,3:ncol(fc)], genes = fc[, 1], group = grp)
isexpr <- rowSums(cpm(y) > 10) >= 2        #this will exclude genes not expressed
y <- y[isexpr, ]
pdf(outfile)
plotMDS(y, labels = lbls)
dev.off()
