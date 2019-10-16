**What we have:**<br/>
read count tables: cDNA raw read counts aligned with STAR for 4 treatment and 4 control experiments<br/>
Reference genome: considering one gene consisting of different exons<br/>

**Using DESeq2**<br/>

**Using limma/voom **<br/>

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
data.apply(lambda x: (x/sum(x))*1000000,axis=0)```

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


