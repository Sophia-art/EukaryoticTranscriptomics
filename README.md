*Make Sure you use the right Browser!!!*
Which, apparently, is not Internet explorer...

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

Question to everyone: It does not recognize DESeqDataSetFromMatrix. But it is a known object. See [https://www.rdocumentation.org/packages/DESeq2/versions/1.12.3/topics/DESeqDataSet-class]
```ruby
DESeq . ds <- DESeqDataSetFromMatrix ( countData = readcounts ,
                                       25 colData = sample _info ,
                                       26 design = ~ condition )
DESeqDataSetFromMatrix
sumOfRow = rowSums(read.counts)
read.counts=cbind(read.counts,sumOfRow)
read.counts.keep=subset(read.counts, read.counts$sumOfRow >0)
head(read.counts.keep, n=5)
```


**Abstract:**<br/>
- read count tables: cDNA raw read counts aligned with STAR for the conditions 4 treatment and 4 control experiments<br/>
- Reference genome: This is one genome consisting of 45706 (expected) genes / exons<br/>
Considering: *STAR (...) use the entire genome as the reference and existing gene annotation as a guide for where to expect large gaps (...) Additionally, lowly expressed isoforms may have very few reads that span their specific splice junctions while, conversely, splice junctions that are supported by few reads are more likely to be false positives. Therefore, novel splice junctions will show a bias towards strongly expressed genes.* [http://chagall.med.cornell.edu/RNASEQcourse/Intro2RNAseq.pdf]
- Output "featureCounts" = raw read counts: Each condition has for each gene an integer number. 0 means "not transcribed" and high integer score means transcribed in a linear dependency:<br/>
*total read count associated with a gene (meta-feature) = the sum of reads associated with each of the exons (feature) that “belong” to that gene*<br/> 
[https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/05_counting_reads.html] 
- task: find genes, that are significant more or less transcribed in the treatment experiments than the control experiments.<br/>
- Strategy: evaluating whether *DESeq2*, *limma/voom*, or *EdgeR* is the most suitable package and applying the best one:
  - [x] analyzing comparisonments of the packages in literature
  [http://chagall.med.cornell.edu/RNASEQcourse/Intro2RNAseq.pdf]
  - [x] analyzing whether the constraints are fulfilled<br/>
  *Instead of using a linear model, DESeq2 and edgeR rely on a negative binomial model to fit the observed read counts to arrive at the estimate for the difference. Originally, read counts had been modeled using the Poisson distribution*<br/>
  - [ ] finding out, which one gives the best results (criteria are p-value, confidential interval, power, (...))<br/>
  *statistical test based on the null hypothesis that the difference is close to zero, which would mean that there is no difference in the gene expression values that could be explained by the conditions.*<br/>
  implement the t-test, ANOVA

**Using DESeq2**<br/>
[http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html]
inputs: 
  - [x] *the number of sequence fragments that have been assigned to each gene.*
  - [x] *in the form of a matrix of integer values. The value in the i-th row and the j-th column of the matrix tells how many reads can be assigned to gene i in sample j.*
  - [x] *The values in the matrix should be un-normalized counts or estimated counts of sequencing reads. The DESeq2 model internally corrects for library size, so transformed or normalized values such as counts scaled by library size should not be used as input.*
- [x] We have already a SummarizedExperiment input. Thus, the next steps would be: 
1. pre-filtering: e.g. minimal pre-filtering to keep only rows that have at least 10 reads total.
```ruby
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
```
2. Note on factor levels: distinguish between control and treatment
```ruby
dds$condition <- factor(dds$condition, levels = c("untreated","treated"))
```
3. Collapsing technical replicates: *DESeq2 provides a function collapseReplicates which can assist in combining the counts from technical replicates into single columns of the count matrix. The term technical replicate implies multiple sequencing runs of the same library.*
4. Differential expression analysis: a single function: DESeq<br/>
*DESeq2’s default method to normalize read counts to account for differences in sequencing depths is imple- mented in estimateSizeFactors() (..) for every gene (= row), determine the geometric mean of its read counts across all samples (yielding the ”pseudo-reference”, i.e. one value per gene), divide every value of the count matrix by the corresponding pseudo-reference value, for every sample (= column), determine the median of these ratios. This is the size factor.* [http://chagall.med.cornell.edu/RNASEQcourse/Intro2RNAseq.pdf]<br/>
Results are the *log2 fold changes, p values and adjusted p values. With no additional arguments to results, the log2 fold change and Wald test p value will be for the last variable in the design formula, and if this is a factor, the comparison will be the last level of this variable over the reference level*<br/>
```ruby
dds <- DESeq(dds)
res <- results(dds)
res
```
5. p value filtering by Independent hypothesis weighting: *A generalization of the idea of p value filtering is to weight hypotheses to optimize power. A Bioconductor package, IHW, is available that implements the method of Independent Hypothesis Weighting*
6. MA-plot to show log2 fold change: * to a given variable over the mean of normalized counts for all the samples in the DESeqDataSet. Points will be colored red if the adjusted p value is less than 0.1.*
(...)
7. Transformation of read counts including variance shrinkage

**Using limma/voom**<br/>
[https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html]
- [x] data input fits
1. loading data
2. preprocessing: Make groups "treatment" and "control"  --> Calculate normalization factors
3. Voom transformation and calculation of variance weights: *Counts are transformed to log2 counts per million reads. A linear model is fitted to the log2 CPM for each gene, and the residuals are calculated. A smoothed curve is fitted to the sqrt(residual standard deviation) by average expression. The smoothed curve is used to obtain weights for each gene and sample that are passed into limma along with the log2 CPMs.*
4. Fitting linear models in limma: Comparisons between groups (log fold-changes) are obtained as contrasts of these fitted linear models: Group "treatment" can be compared with group "control"
Estimate contrast for each gene
```ruby
tmp <- contrasts.fit(fit, contr)
```
Empirical Bayes smoothing of standard errors
```ruby
tmp <- eBayes(tmp)
```
Answer what genes are most differntially expressed with the formula:
```ruby
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
```
(...)

**Using EdgeR for Differential Expression Analysis**<br/>
[https://bioinformatics-core-shared-training.github.io/cruk-bioinf-sschool/Day3/Supplementary-RNAseq-practical.pdf]<br/>
<br/>
<br/>
<br/>
<br/>
**Exploring global read count patterns: Pairwise correlation, hierarchical clustering or Principle Component Analysis**<br/>
[http://chagall.med.cornell.edu/RNASEQcourse/Intro2RNAseq.pdf]<br/>
*technical and biological replicates should show similar expression patterns while the expression patterns of, say, two experimental conditions should be more dissimilar.*<br/>
*The Pearson correlation coefficient, r, is a measure of the strength of the linear relationship between two variables and is often used to assess the similarity of RNA-seq samples in a pair-wise fashion. It is defined as the covariance of two variables divided by the product of their standard deviation.*<br/>
<br/>
*The goal of Principle Component Analysis is to find groups of features (e.g., genes) that have something in common (e.g., certain patterns of expression across different samples), so that the information from thousands of features is captured and represented by a reduced number of groups.The result of PCA are principal components that represent the directions along which the variation in the original multi-dimensional data matrix is maximal.*
```ruby
pc <- prcomp(t(rlog.norm.counts)) 2
plot(pc$x[,1], pc$x[,2],
```

**Create an R environement in Anaconda**
```ruby
conda create -n r_env r-essentials r-base
conda activate r_env
conda install r-BiocManager
```




<br/><br/><br/><br/>
**Code of Dimitri**<br/>
**1. read data**
```ruby
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

%matplotlib inline

df = pd.read_csv('raw_countstdl.sec',sep='\t',index_col=False)
df.rename(columns = {'Unnamed: 0':'Entry'} ,inplace=True)
```

**2. Filtering Data**<br/>
```ruby
data_func=data.apply(lambda x: (x/sum(x))*1000000,axis=0)
```
```ruby
data.apply(lambda x: (x/sum(x))*1000000,axis=0)
```

```ruby
data.apply(lambda x: (x/sum(x))*1000000,axis=0)
df.apply(lambda x: (x*2),axis=0)
df.apply(lambda x: (sum(x))*1000000,axis=0) 
```
```ruby
keep = data_func[data_func[data_func > 100].sum(1) >=4]
```

**3. Normalizing with TMM (trimmed mean of M values as a normalization factor)**<br/>
*Because we assume, that the majority of genes, in both the treated and untreated samples, are not differentially expressed*

**4. Data exploration**<br/>
**5. Estimating the dispersion**<br/>
**6. Differential Expression**<br/>


