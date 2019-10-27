
install packages:
```ruby
install.packages("BiocManager")
library(BiocManager)
BiocManager::install("DESeq2")
library("DESeq2") 
```
load data and have a look how it looks like:
```ruby
read.counts=read.table("C:\\Users\\Sophia Schmidt\\Documents\\Uni\\Master\\Programming for Life Science\\raw_countstdl.txt")
head(read.counts, n=5)
names(read.counts)
```
how to do it:
[http://www.genomicdataanalysis.com/genomicdatatype/rna-seq/]
```ruby
library(DESeq2)
library(magrittr)
read_counts <- read.table("results.txt", header = TRUE )# Read in the count data 
row.names(read_counts) <- read_counts$Geneid
read_counts <- read_counts[, -c(1:6)] # Exclude the columns that do not contain read counts
 
# Create a data frame with metadata
# Old command:
sample_info <- data.frame(condition = "control, treatment", row.names = names(read_counts ))

#New command:
sample_info <- data.frame(condition=rep(c("control","treatment"), each=4),row.names=names(read_counts))
sample_info
 
# Generate the DESeqDataSet
DESeq.ds <- DESeqDataSetFromMatrix(countData = readcounts,
                                      colData = sample_info,
                                      design = ~condition)
 
# Check the DESeqDataSet
colData(DESeq.ds) %>% head
assay(DESeq.ds) %>% head
rowRanges(DESeq.ds) %>% head
 
# Remove the genes without any counts
DESeq.ds <- DESeq.ds[rowSums(counts(DESeq.ds)) > 0,]
 
# Normalize the count reads with size factor
DESeq.ds <- estimateSizeFactors(DESeq.ds)
sizeFactors(DESeq.ds)
norm_counts <- counts(DESeq.ds, normalized = TRUE)
 
# Transform normalized read counts to log2 scale using a pseudocount of 1
log.norm.counts <- log2(norm_counts + 1)
 
# Generate a boxplot of log2-transformed read counts
boxplot(log.norm.counts, notch = TRUE ,
             main = "log2-transformed read counts",
              ylab = "log2(read_counts)")
```


