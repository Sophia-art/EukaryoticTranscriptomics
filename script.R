

setwd("/home/vlad/Documents/programming_ls/group-project/")

reads <- read.table("raw_countstdl.tsv", header = TRUE, sep = '\t')
reads_nonzero <- subset(reads, ctl1 != 0 | ctl2 != 0 | ctl3 != 0 | ctl4 != 0 | treat1 != 0 | treat2 != 0 | treat3 != 0 | treat4 != 0)

head(reads_nonzero)
