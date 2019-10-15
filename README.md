This is the **Master branch**. Here, we can write any text using *Markdown* [(https://guides.github.com/features/mastering-markdown/)] A really short interduction to get started with github is the following one: [(https://guides.github.com/activities/hello-world/
)]. 


**This is the code from Dimitri**

**1. read data**
```ruby
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

%matplotlib inline

df = pd.read_csv('raw_countstdl.sec',sep='\t',index_col=False)
df.rename(columns = {'Unnamed: 0':'Entry'} ,inplace=True)
```

**2. Differential Expression Analysis by edgeR**<br/>
[https://bioinformatics-core-shared-training.github.io/cruk-bioinf-sschool/Day3/Supplementary-RNAseq-practical.pdf]<br/>
**2.1. Normalizing data by total count**<br/>
*Why?*
```ruby
data_func=data.apply(lambda x: (x/sum(x))*1000000,axis=0)
```

*I tried to understand and tried this code out in easy examples:*
```data.apply(lambda x: (x/sum(x))*1000000,axis=0)```
, *but I dont understand it: Refers x not to each element in data? Should it not give always the same value within a column? Like in following easy example:*

```ruby
data.apply(lambda x: (x/sum(x))*1000000,axis=0)
df.apply(lambda x: (x*2),axis=0)
df.apply(lambda x: (sum(x))*1000000,axis=0) 
```

*Or is sum(x) not constant but is the sum of all data of the column till the element (iterating from row 0 till the last row)?* 

**2.2. Filtering data**<br/>
Why? *Because we want to kick out data, where there is no count*

```ruby
keep = data_func[data_func[data_func > 100].sum(1) >=4]
```


**2.3. Normalizing with TMM (trimmed mean of M values as a normalization factor)**<br/>
Why? *Because we assume, that the majority of genes, in both the treated and untreated samples, are not differentially expressed*

**2.4. Data exploration**<br/>
**2.5. Estimating the dispersion**

Feature Count
