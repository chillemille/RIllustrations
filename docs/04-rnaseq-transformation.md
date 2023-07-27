
# Transformation of the RNA-Seq data



In this file, we will transform the RNA-Seq data to align its characteristics with those of microarray measurements using two approaches: \

 - **approach 1**: transformation using voom 
 - **approach 2**: transformation using DESeq2's varianceStabilizingTransformation
 
 
Note that PADOG and web-based tool GSEA, which are the two tools in this paper that require a (manual) transformation of the RNA-Seq data, require the genes in a different format: 

- **PADOG**: NCBI (Entrez) ID. Therefore requires the gene expression data with converted gene IDs
- **GSEA**: ENSEMBL ID. Therefore requires gene expression data with initial gene IDs 

We here illustrate the procedure for both transformation methods for **Entrez ID**, while the identical procedure for ENSEMBL ID is provided in the corresponding the R file 'Instructions_RNASeq_Transformation.R'.

 
## Libraries 

All necessary packages are available on Bioconductor, and should be installed from there if not already available on your machine.


```r
install.packages("BiocManager")
BiocManager::install("tweeDEseqCountData")
BiocManager::install("limma")
BiocManager::install("edgeR")
BiocManager::install("DESeq2")
```

**Load libraries**

```r
library(tweeDEseqCountData) 
library(limma)
library(edgeR)
library(DESeq2)
```

Description of the libraries: 

- **tweeDEseqCountData**: Both methods for the transformation of the RNA-Seq data require the conditions of the samples of the gene expression data set which, in the case of the Pickrell data set, we obtain from the library tweeDEseqCountData. 

- **limma**: Since the first method for the transformation of the RNA-Seq data is based on limma's voom, we have to load the corresponding library. 

- **edgeR**: Normalization is an additional component of the transformation of the RNA-Seq data and the first method employs (by default) a normalization technique provided by the edgeR library. 

- **DESeq2**: Since the second method for the transformation of the RNA-Seq data is provided by DESeq2, we have to load the corresponding library. 


#### Load data 

We load the pre-filtered gene expression data sets with the genes identified in the initial NCBI (Entrez) ID format. Note that for the purpose of simplicity, we here assume that the gene expression measurements have been filtered using the simple approach to pre-filtering introduced by DESeq2. 


```r

load("./data/Results_GeneID_Conversion_DuplicatesRemoval/exprdat_filter_conv_DESeq2.Rdata")

# for GSEA, we would load the pre-filtered RNA-Seq data with the genes identified in the ENSEMBL ID format 
#load("./data/Results_PreFiltering/expression_data_filterDESeq2.Rdata")
```

To obtain the sample conditions of the pickrell data set (which are the respective genders), we also have to load the data. 


```r
# load pickrell data 
data(pickrell)

# assign the conditions of the samples to object "sample_conditions"
sample_conditions <- pickrell.eset$gender
```



For the purpose of readability, we will assign the gene expression data set the more neutral name "expression_data_filt_conv".

```r
# gene expression data set with gene IDs converted to Entrez ID format: 
expression_data_filt_conv <- exprdat_filter_conv_DESeq2
```


## Transformation of the gene expression data set with converted (Entrez) IDs 


### Approach 1: Transformation using voom 


**Step 1: Generate a DGEList object from the gene expression data**


```r
counts_voom <-DGEList(counts = expression_data_filt_conv, 
                      group = sample_conditions)
```

Function arguments: 

- **counts**: corresponds the gene expression measurements 

- **group**: corresponds to the conditions of the samples 

**Step 2: Normalization**

Note that we here use the default normalization method TMM ("Trimmed mean of M-values") [@robinson2010scaling]. Other choices for normalization can be specified in the argument "method" in function *calcNormFactors()*. 


```r
counts_voom <- calcNormFactors(counts_voom)
```

**Step 3: Run voom** 

Note that we do NOT make use of the precision weights that are an official part of voom. Instead, we only proceed with the cpm-transformed (and normalized) gene expression measurements. 


```r
expression_data_voomtransformed <- voom(counts_voom)
```


**Step 4: Convert resulting data set to data frame**

At this point, the transformed gene expression measurements are stored in the object "expression_data_voomtransformed" in the form of an EList object. For the further use in GSA (here specifically: PADOG), we want to convert the transformed gene expression measurements to a data frame. 


```r
expression_data_voomtransformed_Entrez <- as.data.frame(expression_data_voomtransformed)
```

Inspect the first few entries of the transformed RNA-Seq data set.

```r
expression_data_voomtransformed_Entrez[1:10, 1:5]
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
```

### Approach 2: Transformation using DESeq2's VarianceStabilizingTransformation [@love2014moderated]

**Step 1: generate a data frame to contain the condition of each sample**

Just like the starting point of differential expression analysis, we here need to generate the input object(s) required when working with DESeq2.  Here, we generate the data frame colnames which contains for each sample (rownames) the corresponding condition (in the form of "treated" vs. "untreated").


```r
# the names of the samples are stored as the row names of the data frame
# -> important: make sure that the order of the conditions in sample_conditions corresponds to the order of the samples in expression_data
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

**Step 2: Generate DESeqDataSet**

DESeq2 eventually operates on a "DESeqDataSet", which contains the gene expression measurements, the information on the samples (here: the sample conditions) as well as the indication on which variables the count data of each sample depend on. 


```r
dds<-DESeqDataSetFromMatrix(countData = expression_data_filt_conv, 
                            colData = coldata, 
                            design = ~ condition)
```

Relevant arguments in function DESeqDataSetFromMatrix:

- **countData**: count data from the gene expression data set

- **colData**: data frame that contains information on the samples (see above) such as the conditions of the samples (required) and possibly further variables to correct for (such as batch effects)

- **design**: indicates which variables from colData are used for modelling
  + more detailed: the argument design is used to estimate the dispersions and the log2 fold changes of the model 
  + if more than one variable from colData are used in argument design (e.g. a second variable "batch"), the syntax changes to the following formula: "design ~ batch + condition" 
  
Here, we disregard any possible batch effects and focus on the conditions of the samples. 

**Step 3: Perform DESeq2's varianceStabilizingTransformation**

This transformation results in values with a variance approximately that is approximately constant throughout the range of the mean. Furthermore, the values are normalized for library size.  



```r
# perform variance stabilizing transformation
expression_data_vsttransformed <-vst(dds)
    
# since the data set is now in the format of a DESeqTransform data set, the transformed count data are not directly accessible use function assay() which lets us access the count data:
expression_data_vsttransformed <- assay(expression_data_vsttransformed)
```

**Step 4: Convert to data frame**
  
We convert the transformed gene expression data set to a data frame. 


```r
expression_data_vsttransformed_Entrez <- as.data.frame(expression_data_vsttransformed)
```

Inspect the first few entries of the transformed gene expression data set 

```r
expression_data_vsttransformed_Entrez[1:10, 1:5]
#>        NA18486  NA18498  NA18499  NA18501  NA18502
#> 8813  5.475282 6.471173 5.786311 6.551076 6.317231
#> 57147 5.475282 6.416409 6.857113 6.508884 6.396753
#> 55732 4.471560 5.050287 4.706549 5.425412 4.967019
#> 2268  5.939797 6.032867 5.810294 5.996914 6.179690
#> 2519  5.727813 6.159507 5.605070 5.933449 5.487017
#> 4800  8.527528 7.963648 8.355275 7.998985 8.336457
#> 81887 5.822793 5.018101 4.256734 4.818117 4.967019
#> 22875 6.327592 8.199396 7.451246 7.670505 7.061865
#> 5893  4.246067 4.421462 4.433534 4.644343 4.349009
#> 572   6.240388 6.186132 6.286555 5.472331 6.179690
```


