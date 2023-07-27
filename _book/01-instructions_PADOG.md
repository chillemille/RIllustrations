
# (PART) Common Processing Steps {-}

# Prefiltering


 - In this script, we will learn about two options to exclude all genes that do not have a sufficient number of read counts across 
 all samples. We distinguish between these two approaches since the different methods for differential expression analysis (see file "Instructions_Differential_Expression_Analysis") propose different methods for pre-filtering. 

 - option 1: pre-filtering using a function provided by edgeR
 - option 2: pre-filtering proposed by DESeq2
 
## Libraries 

All necessary packages are available on Bioconductor, and should be installed from there if not already available on your machine.


```r
install.packages("BiocManager")
BiocManager::install("tweeDESeqCountData")
BiocManager::install("edgeR")

```

### Load Libraries 

```r

library(tweeDEseqCountData)
library(edgeR)
```

Description of the packages: 

- tweeDESeqCountData: From this library we obtain the gene expression data set we will use for our illustrations. 

- edgeR: edgeR offers a function for pre-filtering we will use below. 

## Preparation of RNA-Seq data set used for illustration

The RNA-Seq data set we will use for this illustration is provided by Pickrell et al. (2010) and is part of the library tweeDEseqCountData. In this data set, the sample conditions (i.e. phenotypes) of the respective samples correspond to the genders. For the purpose of simplicity and readability, we store the gene expression measurements and sample conditions in objects with neutral names.


```r
# load pickrell data set 
data(pickrell)

# access and store gene expression measurements
expression_data <- Biobase::exprs(pickrell.eset)
# access and store sample conditions 
sample_conditions <- pickrell.eset$gender 


# take a look at the gene expression measurements: 
head(expression_data, n = 5)
#>                 NA18486 NA18498 NA18499 NA18501 NA18502
#> ENSG00000000003       0       0       0       0       0
#> ENSG00000000005       0       0       0       0       0
#> ENSG00000000419      22     105      40      55      67
#> ENSG00000000457      22     100     107      53      72
#> ENSG00000000460       5      23      10      18      15
#>                 NA18504 NA18505 NA18507 NA18508 NA18510
#> ENSG00000000003       0       5       0       0       0
#> ENSG00000000005       0       0       0       0       0
#> ENSG00000000419      37      88     127      70      43
#> ENSG00000000457      38      98      69      66      43
#> ENSG00000000460       8      11      16      18       7
#>                 NA18511 NA18516 NA18517 NA18519 NA18520
#> ENSG00000000003       0       0       0       0       0
#> ENSG00000000005       0       0       0       0       0
#> ENSG00000000419     104      42      66      58      82
#> ENSG00000000457     103      51      91      87      77
#> ENSG00000000460       9       7      19      15      22
#>                 NA18522 NA18523 NA18852 NA18853 NA18855
#> ENSG00000000003       0       0       0       0       0
#> ENSG00000000005       0       0       0       0       0
#> ENSG00000000419      95      40      32      30      34
#> ENSG00000000457      87      39      65      50      43
#> ENSG00000000460      17       7      14       8       5
#>                 NA18856 NA18858 NA18861 NA18862 NA18870
#> ENSG00000000003       0       1       0       0       0
#> ENSG00000000005       0       0       0       0       0
#> ENSG00000000419      44      47     113      79      38
#> ENSG00000000457      48      52     117      73      50
#> ENSG00000000460      13      15      23      15       8
#>                 NA18871 NA18909 NA18912 NA18913 NA18916
#> ENSG00000000003       0       0       0       0       0
#> ENSG00000000005       0       0       0       0       0
#> ENSG00000000419      41      95      45      54      49
#> ENSG00000000457      72      93      44      73      81
#> ENSG00000000460       8      17       2       7       8
#>                 NA19093 NA19098 NA19099 NA19101 NA19102
#> ENSG00000000003       0       0       0       0       0
#> ENSG00000000005       0       0       0       0       0
#> ENSG00000000419      80      86      73      59      36
#> ENSG00000000457     103     158     118      89      15
#> ENSG00000000460      15      15      25      13       2
#>                 NA19108 NA19114 NA19116 NA19119 NA19127
#> ENSG00000000003       0       0       0       0       0
#> ENSG00000000005       0       0       0       0       0
#> ENSG00000000419     124      56      64      25     104
#> ENSG00000000457     108      19     150      16      38
#> ENSG00000000460      16       6      28       2      18
#>                 NA19128 NA19130 NA19131 NA19137 NA19138
#> ENSG00000000003       0       0       0       0       1
#> ENSG00000000005       0       0       0       0       0
#> ENSG00000000419      94     117      62      43      46
#> ENSG00000000457      50     104     113      71     120
#> ENSG00000000460      22      25      19      14      10
#>                 NA19140 NA19143 NA19144 NA19147 NA19152
#> ENSG00000000003       0       0       0       0       0
#> ENSG00000000005       0       0       0       0       0
#> ENSG00000000419      94      82      56      96      66
#> ENSG00000000457      46     106      74     116     102
#> ENSG00000000460      21       6      16      13      24
#>                 NA19153 NA19159 NA19160 NA19171 NA19172
#> ENSG00000000003       0       0       0       0       0
#> ENSG00000000005       0       0       0       0       0
#> ENSG00000000419      78      24      63      34      66
#> ENSG00000000457     119      29      64      55      70
#> ENSG00000000460      21       7       9      15      17
#>                 NA19190 NA19192 NA19193 NA19200 NA19201
#> ENSG00000000003       0       0       0       0       2
#> ENSG00000000005       0       0       1       0       0
#> ENSG00000000419      70      85      43      27      77
#> ENSG00000000457     102      89      31      68      96
#> ENSG00000000460       8      28      10       5      16
#>                 NA19203 NA19204 NA19209 NA19210 NA19222
#> ENSG00000000003       0       0       0       0       0
#> ENSG00000000005       0       0       0       0       0
#> ENSG00000000419      82      57      63      89      60
#> ENSG00000000457      58      86     113      48      71
#> ENSG00000000460      10      19      17      12      12
#>                 NA19225 NA19238 NA19239 NA19257
#> ENSG00000000003       0       0       0       0
#> ENSG00000000005       0       0       0       0
#> ENSG00000000419      76      69      84      76
#> ENSG00000000457      81      73      87      81
#> ENSG00000000460       7      21      35      11
# take a look at the sample conditions:
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

# inspect the number of genes (rows) and the number of samples (columns) in the gene expression data set 
dim(expression_data)
#> [1] 52580    69
```

### Option 1: Pre-Filtering using edgeR's builting function
This approach works with the function filterByExpr() from the package edgeR and is the proposed method for pre-filtering for the methods for differential expression analysis edgeR and voom/limma. This approach operates on the cpm-transformed count data (cpm: counts-per-million), which excludes all genes that do NOT have a certain number of counts-per-million in a certain number of samples. Note that this approach accounts for differences in the library size between the different samples. 


#### step 1: Generate input object required by filterByExpr()


```r
expression_data_filterByExpr <- DGEList(counts = expression_data, 
                                        group = sample_conditions)
```


Function description: 

- DGEList(): object to contain RNA sequencing measurements and additional information

Argument description: 

- counts: matrix of RNA-Seq data 
- group: vector that contains the condition of each sample 

#### step 2: Pre-filtering 

The function filterByExpr() creates an indicator which on which genes do and which do not have a sufficient amount of read counts across all samples. Based on this indicator, the gene expression data set is then filtered.



```r
# (i) for each gene, indicate if it fulfils the requirements to be kept for the subsequent analysis 
indicator_keep <- filterByExpr(expression_data_filterByExpr)

# (ii) filter the gene expression data set such that only those genes are kept which fulfill the requirements
expression_data_filterByExpr <- expression_data_filterByExpr[indicator_keep,, keep.lib.sizes = FALSE]
```
Note that the index "keep.lib.sizes = FALSE" ensures that the library size of each sample is recalculated after pre-filtering. 


#### step 3: Obtain final pre-filtered gene expression data set
At this point, we transform the gene expression measurements back to a data frame. The reason for this is that some subsequent steps (such as the conversion of gene IDs do not work with the DGEList format). 



```r
expression_data_filterByExpr <- as.data.frame(expression_data_filterByExpr$counts)

# inspect the number of genes and the number of samples in the final pre-filtered gene expression data set 
nrow(expression_data_filterByExpr)
#> [1] 6246
```
Observe that the number of genes has been reduced compared to the original (unfiltered) gene expression data set.


### Option 2: Simpler Pre-filtering approach (proposed by DESeq2)  
A simpler approach for pre-filtering has been proposed by the method for differential expression analysis DESeq2 (see file "Instructions_Differential_Expression_Analysis". In this approach, only those genes are kept for further analysis that have a pre-specified number of counts X (such as 10) across all samples. A higher value of X thereby leads to more genes being removed. \
Note that DESeq2 proposes a stricter version of pre-filtering in which those genes are kept which have X number of counts in at least Y samples. \
Note that none of the these two "simpler" approaches to pre-filtering take differences in library size into account. 

#### step 1: Pre-Filtering


```r
# indicate which genes have at least 10 read counts across all samples:
indicator_keep <- rowSums( expression_data ) >= 10 

# alternative (and more strict) indicator:
# indicator_keep <- rowSums( expression_data >=10) >= 10

# subset gene expression data set accordingly 
expression_data_filterDESEq2 <- expression_data[indicator_keep,]
```

#### step 2: Inspect final pre-filtered gene expression data set 


```r
dim(expression_data_filterDESEq2)
#> [1] 10151    69
```
Note that the number of genes in the gene expression data set has decreased compared to the initial (unfiltered) gene expression data set. 




