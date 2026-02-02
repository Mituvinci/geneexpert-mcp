#!/usr/bin/Rscript
args = commandArgs(TRUE);
if(length(args) < 1) {
stop("\nUsage: densityplot.R edger3.txt\n")
}
t <- read.table(args[1], quote = "", sep = " ", header = TRUE)
t <- na.omit(t)
pdf(paste0(args[1], "_density.pdf"))
plot(density(t[, 4]))
dev.off()
