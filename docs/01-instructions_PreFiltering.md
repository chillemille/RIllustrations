
# (PART) Common Processing Steps {-}



# Pre-filtering


In this script, we will learn about two options to exclude all genes that do not have a sufficient number of read counts across 
 all samples. We distinguish between these two approaches since the different methods for differential expression analysis (see Chapter 'Differential Expression Analysis') propose different methods for pre-filtering. 

 - **option 1**: pre-filtering using a function provided by edgeR
 - **option 2**: pre-filtering proposed by DESeq2
 
## Libraries 

All necessary packages are available on Bioconductor, and should be installed from there if not already available on your machine.


```r
install.packages("BiocManager")
BiocManager::install("tweeDEseqCountData")
BiocManager::install("edgeR")

```

**Load Libraries** 

```r
library(tweeDEseqCountData)
library(edgeR)
```

**Description of the libraries:** 

- **tweeDESeqCountData**: from this library we obtain the gene expression data set we will use for our illustrations. 

- **edgeR**: offers a function for pre-filtering we will use below. 

## Preparation of RNA-Seq data set used for illustration

For the purpose of simplicity and readability, we store the gene expression measurements and sample conditions from the Pickrell data sets in objects with neutral names, namely in 'expression_data'and 'sample_conditions', respectively.


```r
# load pickrell data set 
data(pickrell)

# access and store gene expression measurements
expression_data <- Biobase::exprs(pickrell.eset)

# access and store sample conditions 
sample_conditions <- pickrell.eset$gender 
```

Take a look at a few entries of the gene expression data set: 

```r
# inspect the read counts of the first 5 genes in the first 5 samples
expression_data[1:5, 1:5]
#>                 NA18486 NA18498 NA18499 NA18501 NA18502
#> ENSG00000000003       0       0       0       0       0
#> ENSG00000000005       0       0       0       0       0
#> ENSG00000000419      22     105      40      55      67
#> ENSG00000000457      22     100     107      53      72
#> ENSG00000000460       5      23      10      18      15
```

Take a look at the sample conditions:

```r
sample_conditions
#>  [1] male   male   female male   female male   female male  
#>  [9] female male   female male   female male   female male  
#> [17] female female male   female male   female female male  
#> [25] female male   female female male   female female male  
#> [33] female male   female female female female male   female
#> [41] male   male   female female male   female female male  
#> [49] female female male   female male   male   female female
#> [57] male   female male   female male   female female male  
#> [65] female female female male   female
#> Levels: female male
```

Inspect the number of genes (rows) and the number of samples (columns) in the gene expression data set:

```r
dim(expression_data)
#> [1] 52580    69
```

## Pre-Filtering

### Option 1: Pre-Filtering using edgeR's builting function

This approach works with the function *filterByExpr()* from the package edgeR and is the proposed method for pre-filtering for the methods for differential expression analysis edgeR and voom/limma. This approach operates on the cpm-transformed count data (cpm: counts-per-million) and excludes all genes that do NOT have a certain number of counts-per-million in a certain number of samples. Note that this approach accounts for differences in the library size between the different samples. 


**step 1: Generate input object required by filterByExpr()**


```r
expression_data_filterByExpr <- DGEList(counts = expression_data, 
                                        group = sample_conditions)
```


Description of function:

- **DGEList()**: object to contain RNA sequencing measurements and additional information

Description of arguments: 

- **counts**: matrix of RNA-Seq data 
- **group**: vector that contains the condition of each sample 

**step 2: Filter out lowly expressed genes**

The function *filterByExpr()* creates an indicator which on which genes do and which do not have a sufficient amount of read counts across all samples. Based on this indicator, the gene expression data set is then filtered.



```r
# (i) for each gene, indicate if it fulfils the requirements to be kept for the subsequent analysis 
indicator_keep <- filterByExpr(expression_data_filterByExpr)

# (ii) filter the gene expression data set such that only those genes are kept which fulfill the requirements
expression_data_filterByExpr <- expression_data_filterByExpr[indicator_keep,, keep.lib.sizes = FALSE]
```
Note that the index **keep.lib.sizes = FALSE** ensures that the library size of each sample is recalculated after pre-filtering. 


**step 3: Obtain final pre-filtered gene expression data set**

At this point, we transform the gene expression measurements back to a data frame. The reason for this is that some subsequent steps (such as the conversion of gene IDs do not work with the DGEList format). 



```r
expression_data_filterByExpr <- as.data.frame(expression_data_filterByExpr$counts)

# inspect the number of genes and the number of samples in the final pre-filtered gene expression data set 
nrow(expression_data_filterByExpr)
#> [1] 6246
```
Observe that the number of genes has been reduced compared to the original (unfiltered) gene expression data set.


### Option 2: Simpler Pre-filtering approach (proposed by DESeq2)  

A simpler approach for pre-filtering has been proposed by the method for differential expression analysis DESeq2. In this approach, only those genes are kept for further analysis that have a pre-specified number of counts $X$ (such as 10) across all samples. A higher value of $X$ thereby leads to more genes being removed. \
Note that DESeq2 proposes a stricter version of pre-filtering in which those genes are kept which have $X$ number of counts in at least $Y$ samples. \
Note that none of the these two "simpler" approaches to pre-filtering take differences in library size into account. 

**step 1: Pre-Filtering**


```r
# indicate which genes have at least 10 read counts across all samples:
indicator_keep <- rowSums(expression_data) >= 10 

# alternative (and more strict) indicator:
# indicator_keep <- rowSums( expression_data >=10) >= 10

# subset gene expression data set accordingly 
expression_data_filterDESEq2 <- expression_data[indicator_keep,]
```

**step 2: Inspect final pre-filtered gene expression data set**


```r
dim(expression_data_filterDESEq2)
#> [1] 10151    69
```
Note that the number of genes in the gene expression data set has decreased compared to the initial (unfiltered) gene expression data set. 




