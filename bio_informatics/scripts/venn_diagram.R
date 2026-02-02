#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) <= 1) {
  stop("Usage: venn_diagram.R file1 file2 ...\n\nWill make a pdf with the name file1_venn.pdf", call. = FALSE)
}

library(VennDiagram)

lists <- data.frame()

fname <- args[1]
fname <- gsub("\\.[^.]*$", "", fname)
fname <- paste0(fname, "_venn.png")

stuff <- list()
for (i in args)
{
  tmp <- read.table(i, header = F, quote = "", comment.char = "")
  stuff <- append(stuff, tmp)
  print(i)
}
venn.diagram(
  x = stuff,
  filename = fname,
  category.names = args,
  height = 1000,
  width = 1000,
  disable.logging = T,
  output = TRUE
)
