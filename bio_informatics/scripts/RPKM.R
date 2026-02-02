#!/usr/bin/Rscript
# EDITED by Claude Code (2026-01-07):
#   - Read CSV input instead of space-delimited text
#   - Write CSV output instead of space-delimited text
#   - Fixed default output filename to be based on input filename

# RPKM = number of reads per kilo bps per million reads
# reads: n
# totals reads in sample S: N
# gene length: L
# output in format log2(RPKM+1)+0.0001

args = commandArgs(trailingOnly = TRUE)

if (length(args)==0) {
  stop("Usage: RPKM.R input.count.filtered.csv [outputfile.csv]", call.=FALSE)
} else if (length(args)==1) {
  # EDIT: Default output file based on input filename
  # Change .count.filtered.csv to .rpkm.csv
  args[2] <- sub("\\.count\\.filtered\\.csv$", ".rpkm.csv", args[1])
  if (args[2] == args[1]) {
    # Fallback if pattern didn't match
    args[2] = "outRPKM.csv"
  }
}

# EDIT: Read CSV instead of space-delimited
table <- read.csv(args[1], header = TRUE)

# calculate total number of reads for each sample
# have to deal with special case of 1 sample...
totalreads <- vector()
if(ncol(table) == 3){
	totalreads <- append(totalreads, sum(table[, 3]))
} else {
	for(i in table[,-1:-2])
	{
		totalreads <- append(totalreads, sum(i))
	}
}

#debugging
cat("Total reads:\n")
print(totalreads)

totalreads <- rep(totalreads, each = length(table[,1]))

# extract gene length
L <- table[,2]	# gene length

# replace read count with log2(n/((L/1000)*(N/1000000))+1)+0.0001 and put in table
table2 <- data.frame(cbind(table[1], log2(table[,-1:-2]/((L/1000.0)*(totalreads/1000000.0))+1)+0.0001))

#disable scientific notation
options(scipen=999)

# EDIT: Write CSV instead of space-delimited
# Keep only ID column (drop Length) from original table
colnames(table2) <- c("ID", colnames(table)[-(1:2)])
write.csv(table2, file = args[2], row.names = FALSE, quote = FALSE)

cat("RPKM normalization complete. Output written to:", args[2], "\n")
