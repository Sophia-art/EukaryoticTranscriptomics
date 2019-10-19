**Abstract:**<br/>
- read count tables: cDNA raw read counts aligned with STAR for the conditions 4 treatment and 4 control experiments<br/>
- Reference genome: This is one genome consisting of 45706 genes<br/>
- Single values: Each condition has for each gene an integer number. 0 means "not transcribed" and high integer score means transcribed in a linear dependency:<br/>
*total read count associated with a gene (meta-feature) = the sum of reads associated with each of the exons (feature) that “belong” to that gene* [https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/05_counting_reads.html] 
- task: find genes, that are significant more or less transcribed in the treatment experiments than the control experiments.<br/>
- Strategy: evaluating whether *DESeq2*, *limma/voom*, or *EdgeR* is the most suitable package and applying the best one:
  - [x] analyzing whether the constraints for the usage of the packages are fulfilled in our data 
  - [ ] analyzing comparisonments of the packages in literature
  - [ ] finding out, which one gives the best results (criteria are (?) p-value, variance, graphs, boxplots etc)

**Using DESeq2**<br/>
[http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html]
inputs: 
  - [x] *the number of sequence fragments that have been assigned to each gene.*
  - [x] *in the form of a matrix of integer values. The value in the i-th row and the j-th column of the matrix tells how many reads can be assigned to gene i in sample j.*
  - [x] *The values in the matrix should be un-normalized counts or estimated counts of sequencing reads. The DESeq2 model internally corrects for library size, so transformed or normalized values such as counts scaled by library size should not be used as input.*
- [x] We have already a SummarizedExperiment input. Thus, the next steps would be: 
      1. pre-filtering
      


**Using limma/voom**<br/>

**Using EdgeR for Differential Expression Analysis**<br/>
[https://bioinformatics-core-shared-training.github.io/cruk-bioinf-sschool/Day3/Supplementary-RNAseq-practical.pdf]<br/>

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


