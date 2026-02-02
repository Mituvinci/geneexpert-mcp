#!/usr/bin/env Rscript
library(ggplot2)
library(optparse)

# set up some nice args?
# set title and group names

group_help <- paste0(
  "Example group information format:\n",
  "\t\t\t\"g1,g1,g1,g2,g2,g2,g3,g3,g3\"\n",
  "\t\tYou can replace the \"g\" with whatever you want\n",
  "\t\tOr specify the group as \"auto\" and it will try its best"
)

option_list <- list(
  make_option(c("-i", "--infile"),
    type = "character", default = NULL,
    help = "Input RPKM file", metavar = "character"
  ),
  make_option(c("-g", "--groups"),
    type = "character", default = "auto",
    help = group_help, metavar = "character"
  ),
  make_option(c("-t", "--title"),
    type = "character", default = "PCA",
    help = "Title of the plot"
  ),
  make_option(c("-o", "--outfile"),
    type = "character", default = NULL
  ),
  make_option(c("-e", "--expression_threshold"),
    type = "double", default = 2,
    help = "Min. RPKM for a gene to count as 'expressed'"
  ),
  make_option(c("-n", "--num_expressed"),
    type = "integer", default = 2,
    help = "Min. number of samples a gene needs to be expressed in"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$infile) | is.null(opt$groups)) {
  print_help(OptionParser(option_list = option_list))
  quit()
}

if (is.null(opt$outfile)) {
  opt$outfile <- paste0(opt$infile, ".pdf")
}

fc <- read.table(opt$infile, header = TRUE, sep = " ")

# filter out non-expresssed genes
counts <- fc[, -1]

keep <- rowSums(counts > opt$expression_threshold) >= opt$num_expressed

counts <- counts[keep, ]

# auto generate groups
if (opt$groups == "auto") {
  opt$groups <- substr(colnames(counts), 0, 1)
} else {
  opt$groups <- strsplit(opt$groups, ",")[[1]]
}

print(paste0("Saving PCA plot as: ", opt$outfile))

pca_data <- prcomp(t(counts))
pca_data_perc <- round(100 * pca_data$sdev^2 / sum(pca_data$sdev^2), 1)
df_pca_data <- data.frame(
  PC1 = pca_data$x[, 1],
  PC2 = pca_data$x[, 2], sample = opt$groups
)

# make the legend look a bit neater
Groups <- opt$groups

pdf(opt$outfile)
print(ggplot(
  df_pca_data,
  aes(PC1, PC2, color = Groups)
) +
  geom_point(size = 3) +
  ggtitle(opt$title) +
  labs(
    x = paste0("PC1 (", pca_data_perc[1], "%)"),
    y = paste0("PC2 (", pca_data_perc[2], "%)")
  ))
dev.off()
