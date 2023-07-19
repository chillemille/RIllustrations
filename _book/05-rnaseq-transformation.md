
# RNA-Seq Transformation


Branching point: continue to the next chapter if conversion of gene IDs is not necessary. However, if conversion was/is necessary, proceed to Chapter ?? Cluster Profiler.


## Content of this script 

In this file, we will transform the RNA-Seq data to align its characteristics with those of microarray measurements using two approaches: \

 - approach 1: transformation using voom 
 - approach 2: transformation using DESeq2's varianceStabilizingTransformation
 
 
Note that PADOG and GSEA, which are the two tools in this paper that require a (manual) transformation of the RNA-Seq data, require the genes in a different format: 

- PADOG: Entrez gene ID -> requires the gene expression data with converted gene IDs
- web-based tool GSEA: Ensembl gene ID -> required  gene expression data with initial gene IDs 

We here illustrate the procedure for both transformation methods for Entrez ID, while the identical procedure for Ensembl ID is provided in the corresponding .R-file. 

 
### Libraries 

All necessary packages are available on Bioconductor, and should be installed from there if not already available on your machine.


```r
install.packages("BiocManager")
BiocManager::install("tweeDEseqCountData")
BiocManager::install("limma")
BiocManager::install("edgeR")
BiocManager::install("DESeq2")
```


```r
library(tweeDEseqCountData) # for pickrell data set 
library(limma)
library(edgeR)
library(DESeq2)
```

Description of the libraries: 

- tweeDESeqCountData: Both methods for the transformation of the RNA-Seq data require the conditions of the samples of the gene expression data set which, in the case of the Pickrell data set, we obtain by the corresponding library. 

- limma: Since the first method for the transformation of the RNA-Seq data is based on limma's voom, we have to load the corresponding library. 

- edgeR: Normalization is an additional component of the transformation of the RNA-Seq data and the first method employs (by default) a normalization technique provided by the edgeR library. 

- DESeq2: Since the second method for the transformation of the RNA-Seq data is provided by DESeq2, we have to load the corresponding library. 


#### Load data 

We load the pre-filtered gene expression data sets 

1. with the genes identified in the initial Ensembl ID format.
2. with the genes converted to the NCBI (Entrez) ID format.

Note that for the purpose of simplicity, we here assume that the gene expression measurements have been filtered using the simple approach to pre-filtering introduced by DESeq2. 

To obtain the sample conditions of the pickrell data set (which are the respective gender), we also have to load the data. 



```r
# for PADOG:
load("./data/Results_GeneID_Conversion_DuplicatesRemoval/exprdat_filter_conv_DESeq2.Rdata")
# for GSEA 
load("./data/Results_PreFiltering/expression_data_filterDESeq2.Rdata")

# load pickrell data 
data(pickrell)
# assign the conditions of the samples to object "sample_conditions"
sample_conditions <- pickrell.eset$gender
```

For the purpose of readability, we will work with both objects with more neutral names: 


```r
# gene expression data set with gene IDs converted to Entrez ID format: 
expression_data_filt_conv <- exprdat_filter_conv_DESeq2
# gene expression data set with gene IDs in initial Ensembl ID format: 
expression_data_filt <- expression_data_filterDESeq2
```


## Transformation of the gene expression data set with converted (Entrez) IDs 


### Approach 1: Transformation using voom 


#### step 1: Generate a DGEList object from the gene expression data 


```r
counts_voom <-DGEList(counts = expression_data_filt_conv, 
                      group = sample_conditions)
```

Function arguments: 

- counts: corresponds the gene expression measurements 

- group: corresponds to the conditions of the samples 

#### step 2: Normalization

Note that we here use the default normalization method TMM ("Trimmed mean of M-values"). Other choices for normalization can be specified in the argument "method" in function calcNormFactors(). 


```r
counts_voom <- calcNormFactors(counts_voom)
```

#### step 3: Run voom 

Note that we do NOT make use of the precision weights that are an official part of voom. Instead, we only proceed with the cpm-transformed (and normalized) gene expression measurements. 


```r
expression_data_voomtransformed <- voom(counts_voom)
```


#### step 4: Convert resulting data set to data frame

At this point, the transformed gene expression measurements are stored in the object "expression_data_voomtransformed" in the form of an EList object. For the further use in GSA (here specifically: PADOG), we want to convert the transformed gene expression measurements to a data frame. 


```r
expression_data_voomtransformed_Entrez <- as.data.frame(expression_data_voomtransformed)

# inspect (part of) the transformed gene expression data set 
expression_data_voomtransformed_Entrez[1:10, 1:10]
#>         NA18486  NA18498   NA18499  NA18501  NA18502
#> 8813  3.5635974 4.994885 4.0568571 5.047515 4.808758
#> 57147 3.5635974 4.924837 5.4652000 4.994566 4.911852
#> 55732 1.5311760 2.828374 2.1093245 3.462552 2.686139
#> 2268  4.2615689 4.413337 4.0920466 4.319188 4.626761
#> 2519  3.9543874 4.586669 3.7830963 4.230379 3.614586
#> 4800  7.3077585 6.731166 7.1453673 6.716623 7.143454
#> 81887 4.0941122 2.765639 0.8869321 2.340562 2.686139
#> 22875 4.7859899 6.986313 6.1557990 6.356387 5.726296
#> 5893  0.8790993 1.361248 1.4174468 1.953539 1.191375
#> 572   4.6716572 4.622514 4.7504301 3.538501 4.626761
#>        NA18504   NA18505   NA18507  NA18508   NA18510
#> 8813  4.504121 5.3436244 5.4990159 4.767209 4.5153346
#> 57147 4.542089 5.4980707 4.6236036 4.682940 4.5153346
#> 55732 2.362766 2.3995808 2.5490566 2.837111 1.9792817
#> 2268  4.157946 3.3678720 4.6645338 2.271514 4.5153346
#> 2519  3.701568 4.2854098 4.4704468 4.002697 3.7448165
#> 4800  6.688931 6.7738643 6.4559472 6.451025 6.9021139
#> 81887 1.734734 2.8302152 2.5490566 2.913060 3.0265875
#> 22875 5.455212 6.9394139 6.4559472 5.899121 5.4211193
#> 5893  1.445228 0.6833738 0.6745875 1.328097 0.8797461
#> 572   5.120793 4.3517523 4.8175455 4.787529 5.3852741
```

### Approach 2: Transformation using DESeq2's VarianceStabilizingTransformation 


#### step 1: generate a data frame to contain the condition of each sample: 

Just like starting point of differential expression analysis, we here need to generate the input object(s) required when working with DESeq2.  Here, we generate the data frame colnames which contains for each sample (rownames) the corresponding condition (in the form of "treated" vs. "untreated").


```r
# the names of the samples are stored as the row names of the data frame
# -> important: make sure that the order of the conditions in sample_conditions corresponds to the order of the 
#               samples in expression_data
  coldata <- data.frame(sample_conditions, 
                        row.names = colnames(expression_data_filt_conv))
  
# rename the column header to "condition" 
colnames(coldata) <- "condition"
  
# recode the variable condition as a factor 
# rename the sample conditions (in this case from "female" and "male") to "untreated" and "treated"
# note: make sure that the control level in variable condition is coded as the first level (i.e. "untreated")

coldata$condition <- factor(coldata$condition, 
                            labels = c("untreated","treated"))

# inspect resulting data frame: 
head(coldata, n = 10)
#>         condition
#> NA18486   treated
#> NA18498   treated
#> NA18499 untreated
#> NA18501   treated
#> NA18502 untreated
#> NA18504   treated
#> NA18505 untreated
#> NA18507   treated
#> NA18508 untreated
#> NA18510   treated
```

#### step 2: Generate DESeqDataSet

DESeq2 eventually operates on a "DESeqDataSet", which contains the gene expression measurements, the information on the samples (here: the sample conditions) as well as the indication on which variables the count data of each sample depend on. 


```r
dds<-DESeqDataSetFromMatrix(countData = expression_data_filt_conv, 
                            colData = coldata, 
                            design = ~ condition)
```

Relevant arguments in function DESeqDataSetFromMatrix:

- countData: count data from the gene expression data set

- colData: data frame that contains information on the samples (see above) such as the conditions of the samples (required) and possibly further variables to correct for (such as batch effects)

- design: indicates which variables from colData are used for modelling
  + more detailed: the argument design is used to estimate the dispersions and the log2 fold changes of the model 
  + if more than one variable from colData are used in argument design (e.g. a second variable "batch"), the syntax changes to the     following formula: "design ~ batch + condition" 
  
Here, we disregard any possible batch effects and focus on the conditions of the samples. 

#### step 3: Perform DESeq2's varianceStabilizingTransformation 

This transformation results in values with a variance approximately that is approximately constant throughout the range of the mean. Furthermore, the values are normalized for library size.  



```r
# perform variance stabilizing transformation
expression_data_vsttransformed <-vst(dds)
    
# since the data set is now in the format of a DESeqTransform data set, the transformed count data are not directly accessible use function assay() which lets us access the count data:
expression_data_vsttransformed <- assay(expression_data_vsttransformed)
```

#### step 4: Convert to data frame 
  
We convert the transformed gene expression data set to a data frame. 


```r
expression_data_vsttransformed_Entrez <- as.data.frame(expression_data_vsttransformed)

# inspect (part of) the transformed gene expression data set 
expression_data_vsttransformed_Entrez[1:10, 1:10]
#>        NA18486  NA18498  NA18499  NA18501  NA18502  NA18504
#> 8813  5.475282 6.471173 5.786311 6.551076 6.317231 6.107874
#> 57147 5.475282 6.416409 6.857113 6.508884 6.396753 6.135961
#> 55732 4.471560 5.050287 4.706549 5.425412 4.967019 4.813408
#> 2268  5.939797 6.032867 5.810294 5.996914 6.179690 5.859580
#> 2519  5.727813 6.159507 5.605070 5.933449 5.487017 5.554671
#> 4800  8.527528 7.963648 8.355275 7.998985 8.336457 7.943355
#> 81887 5.822793 5.018101 4.256734 4.818117 4.967019 4.540723
#> 22875 6.327592 8.199396 7.451246 7.670505 7.061865 6.857525
#> 5893  4.246067 4.421462 4.433534 4.644343 4.349009 4.429268
#> 572   6.240388 6.186132 6.286555 5.472331 6.179690 6.583569
#>        NA18505  NA18507  NA18508  NA18510
#> 8813  6.751417 6.884505 6.290273 6.117851
#> 57147 6.879356 6.190724 6.226492 6.117851
#> 55732 4.831633 4.910830 5.048319 4.647757
#> 2268  5.343466 6.221367 4.774482 6.117851
#> 2519  5.939749 6.077764 5.742060 5.585061
#> 4800  8.004757 7.718627 7.704235 8.140429
#> 81887 5.044543 4.910830 5.088151 5.156246
#> 22875 8.157830 7.718627 7.215519 6.830071
#> 5893  4.187743 4.199785 4.400327 4.244141
#> 572   5.986951 6.337580 6.305771 6.800351
```


