
# Differential expression analysis



In this script, we will perform differential expression analysis for each of the three parametric methods: 

- voom/limma

- DESeq2 

- edgeR 

Here, we illustrate this process for two gene ID format Entrez IDs, while we additionally present the identical process for the format Ensembl ID in the R script "Instructions_Differential_Expression_Analysis.R".

Note that at at the end of the illustration for each method, we rename some of the columns in the corresponding results tables to unify them across the different methods. While this step is not required for differential expression analysis or for the subsequent use of gene set analysis IN GENERAL, renaming of some columns (such as adjusted p-value) is necessary in this context so that we can use the same code to illustrate further conduct of GSA in the following chapters, independent of the method for differential expression analysis used. 


## Libraries 

All necessary packages are available on Bioconductor, and should be installed from there if not already available on your machine.



```r
install.packages("BiocManager")
BiocManager::install("tweeDESeqCountData")
BiocManager::install("limma")
BiocManager::install("DESeq2")
BiocManager::install("edgeR")
BiocManager::install("dplyr")
```

**Load libraries**


```r
library(tweeDEseqCountData)
library(limma)
library(DESeq2)
library(edgeR)
library(dplyr)
```

Description of the libraries: 

- **tweeDEseqCountData**: In addition the the pre-filtered gene expression data set (with converted gene IDs) that we load in the next step, we need the conditions of the samples which are provided by tweeDEseqCountData. 

- **limma**: First option of methods for differential expression analysis.

- **DESeq2**: Second option of methods for differential expression analysis. 

- **edgeR**: Third option of methods for differential expression analysis. 

- **dplyr**: Used to unify some of the columns in the results table of the three methods for differential expression analysis. 


### Load data 

Note that depending on the method for differential expression analysis, different methods for pre-filtering are proposed: 

- **voom/limma** and **egdeR**: pre-filtering using edgeR's built-in function filterByExpr()
- **DESeq2**: manual/simple pre-filtering

We therefore work with both pre-filtered gene expression data sets generated in Chapter 'Pre-Filtering':


```r
load("./data/Results_GeneID_Conversion_DuplicatesRemoval/exprdat_filter_conv_DESeq2.Rdata")
load("./data/Results_GeneID_Conversion_DuplicatesRemoval/exprdat_filter_conv_filterByExpr.Rdata")
```



**Obtain the sample conditions** 
Obtain the sample conditions (i.e. phenotypes) of pickrell data set and store them in object "sample_conditions":


```r
# load pickrell data set 
data(pickrell)

# the conditions/phenotypes of the samples can be addressed with the following syntax: 
pickrell.eset$gender
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

# store sample conditions 
sample_conditions <- pickrell.eset$gender
```


## Differential Expression Analysis (with Entrez gene ID)

A brief description of the methods for differential expression analysis is provided in the supplement.

### Option 1: Differential expression analysis using limma 

**1. Generate required input object**

At this point, the gene expression data set 'exprdat_filter_conv_filterByExpr' is
a data frame and needs to be converted to a DGEList object. 


```r
# store expression data with corresponding sample conditions in object of class DGEList
y <- DGEList(counts = exprdat_filter_conv_filterByExpr, 
              group = sample_conditions)
```

Required function arguments: 

- **counts**: matrix that contains RNA-Seq data (i.e. count data)

Optional function arguments:

- **group**: indicates condition of each sample/library. This argument must be specified sooner or later (such as in subsequent functions) so we just specify it at this point 

Note: we leave the remaining arguments in their default state since the corresponding info will be added through the following pipeline of functions. 

**2. Normalization**

The following piece of code generates a normalization factor for each sample which accounts for sample-specific effects (such as differences in library sizes and effects of compositionality). If not accounted for, these effects prevent a comparison between the samples. 


```r
y <- calcNormFactors(y)
```

Note that this function does not transform the count data, but rather generates a normalization factor for each sample which is separately incorporated into the subsequent analysis. 


**3. voom transformation** 

Perform the voom [@law2014voom] transformation to make limma, which was initially developed for microarray measurements, available to RNA-Seq data. voom transforms the RNA-Seq data to log counts-per-million and generates a precision weight for each entry in the gene expression data sets to be incorporated into the limma pipeline. 


```r
# (i) generate design matrix (rows correspond to samples, columns indicate which coefficients are to be estimated)
design_matrix <- model.matrix( ~ sample_conditions)

# (ii) voom transformation 
y <- voom(y, design = design_matrix)
```

**4. Test for differential expression** 


```r
# (i) fit a linear model for each gene
y <- lmFit(y)

# (ii) calculate the statistics for the assessment of differential expression 
y <- eBayes(y)

# (iii) Get the result table for each gene whose differential expression was assessed 
DE_results_limma_Entrez <- topTable(fit = y, 
                                    number =  nrow(y))
# note: 'number = nrow(y)' ensures that all genes are displayed in results
```



**5. Rename columns in results of differential expression analysis**

We rename some of the columns in the results table to unify them across the different methods for differential expression analysis. For this, we use the function *rename()* provided by the library dplyr. 



```r
# first, we transform the results table to a data frame so that we see the results table directly when accessing it through the name "res"
DE_results_limma_Entrez <- as.data.frame(DE_results_limma_Entrez)

# rename column that contains adjusted p-values: rename padj to p_adj
DE_results_limma_Entrez <- dplyr::rename(DE_results_limma_Entrez, 
                                         p_adj = "adj.P.Val")
```

**6. Take a look at top 10 genes in the results table generated with limma**

```r
head(DE_results_limma_Entrez, n = 10)
#>             logFC   AveExpr         t      P.Value
#> 6192    9.2105431 1.9863048 47.186499 1.020114e-55
#> 9087    5.0806396 0.5128935 26.668835 6.058859e-39
#> 8228   -0.9240302 5.3334193 -8.779270 5.579187e-13
#> 8226   -0.8544468 2.1626976 -4.415343 3.497072e-05
#> 159013 -0.5508457 3.7163074 -4.251355 6.307821e-05
#> 55787  -0.4042722 5.1160616 -4.045191 1.302619e-04
#> 50486  -1.4703261 5.1922528 -3.855673 2.493894e-04
#> 3725   -0.5880074 8.2731003 -3.920361 2.001855e-04
#> 317     0.4046249 5.9681588  3.822268 2.791390e-04
#> 3552   -1.4855077 3.5435981 -3.666767 4.682471e-04
#>               p_adj           B
#> 6192   6.133944e-52 84.20170055
#> 9087   1.821596e-35 62.68255252
#> 8228   1.118255e-09 19.22260877
#> 8226   5.256973e-02  2.15126317
#> 159013 7.585786e-02  1.55871370
#> 55787  1.305441e-01  0.74145045
#> 50486  1.864959e-01  0.15615850
#> 3725   1.719594e-01  0.11995311
#> 317    1.864959e-01 -0.07537785
#> 3552   2.616375e-01 -0.24387530
```


### option 2: Differential expression analysis using DESeq2 

**1. Generate required input object**

DESeq2 operates on the format DESeqDataSet which contains information on the count data, the conditions of the samples and the design. Note that for gene expression data sets that must be imported to R, additional steps are necessary before the following code can be run. 

**Generate a data frame which contains the condition of each sample**

Here, the information on the sample conditions is the only column in the data frame. However, if further variables (such as batch effects) are to be controlled for, the corresponding variables must additionally be added to coldata. 



```r
# the names of the samples are stored as the row names of the data frame
# -> important: make sure that the order of the conditions in sample_conditions corresponds to the order of the samples in expression_data
coldata <- data.frame(sample_conditions, 
                      row.names = colnames(exprdat_filter_conv_DESeq2))

# rename the column header to "condition" 
colnames(coldata) <- "condition"

# recode the variable condition as a factor 
# rename the sample conditions (in this case from "female" and "male") to "untreated" and "treated"
# note: make sure that the control level in variable condition is coded as the first level (i.e. "untreated")
coldata$condition <- factor(coldata$condition, 
                            labels = c("untreated","treated"))
```

**Generate the DESeqDataSet**


```r
dds_Entrez <- DESeqDataSetFromMatrix(countData = exprdat_filter_conv_DESeq2, 
                                     colData = coldata, 
                                     design = ~ condition)
```

Relevant arguments in function DESeqDataSetFromMatrix:

- **countData**: count data from the gene expression data set

- **colData**: data frame that contains information on the samples (see above) such as the conditions of the samples (required) and possibly further variables to correct for (such as batch effects)

- **design**: indicates which variables from colData are used for modelling

  + more detailed: the argument design is used to estimate the dispersions and the log2 fold changes of the model 
  
  + if more than one variable from colData are used in argument design (e.g. a second variable "batch"), the syntax changes to the following formula: 'design ~ batch + condition'
  
  + make sure that the variable of interest (here: variable that represents conditions of the samples) is placed at the end of the formula


**2. Test for differential expression**

DESeq2 calculates the normalization factor of each sample, estimates the dispersion parameter, fits a negative binomial model. A p-value of differential expression for each gene is then computed using a Wald test. 


```r
# perform default differential expression analysis 
dds_Entrez <- DESeq(dds_Entrez)

# generate results table which provides
DE_results_DESeq2_Entrez <- results(dds_Entrez)
```

Note: function *results()* provides the base mean across all samples, the estimated log2 fold changes, standard errors, test statistics, p-values, and adjusted p-values for each gene. 


**3. Rename columns in results of differential expression analysis**



```r
# first, we transform to results table to a data frame so that we see the results
# table directly when accessing it through the name "res"
DE_results_DESeq2_Entrez <- as.data.frame(DE_results_DESeq2_Entrez)

# rename column that contains adjusted p-values: rename padj to p_adj
DE_results_DESeq2_Entrez <- dplyr::rename(DE_results_DESeq2_Entrez, 
                                          p_adj = padj)
```

**4. Take a look at top 10 genes in the results table generated with DESeq2**

```r
head(DE_results_DESeq2_Entrez, n = 10 )
#>         baseMean log2FoldChange      lfcSE        stat
#> 8813   62.726871    0.039434784 0.11345444  0.34758255
#> 57147  70.453343    0.028533065 0.11131002  0.25633868
#> 55732  13.001060    0.164117365 0.15102073  1.08672080
#> 2268   43.577756    0.113711805 0.19545410  0.58178266
#> 2519   37.963998   -0.006527064 0.10803176 -0.06041801
#> 4800  228.478329    0.051147440 0.09641906  0.53047024
#> 81887  11.534031    0.025191913 0.18501229  0.13616346
#> 22875 175.833526    0.299999462 0.16347952  1.83508896
#> 5893    4.688308    0.280425754 0.19291917  1.45359198
#> 572    66.369975   -0.264896424 0.10233495 -2.58852362
#>            pvalue     p_adj
#> 8813  0.728153713 0.9313348
#> 57147 0.797689329 0.9506228
#> 55732 0.277160217 0.7502690
#> 2268  0.560713083 0.8849136
#> 2519  0.951822712 0.9884098
#> 4800  0.595785935 0.8982377
#> 81887 0.891692065 0.9763145
#> 22875 0.066492509 0.5209390
#> 5893  0.146059463 0.6396912
#> 572   0.009638834 0.3141234
```


### option 3: Differential expression analysis using edgeR 


**1. Generate required input object**

Note that edgeR operates on a DGEList object.



```r
y <- DGEList(counts = exprdat_filter_conv_filterByExpr, 
             group = sample_conditions)
```
Required function arguments: 

- **counts**: matrix that contains RNA-Seq data (i.e. count data)

Optional function arguments: 

- **group**: indicates condition of each sample/library 
  + this argument must be specified sooner or later (such as in subsequent functions) so we just specify it at this point 

Note that we leave the remaining arguments in their default state since the corresponding info will be added through the following pipeline of functions. 

**2. Normalization** 

The following piece of code generates a normalization factor for each sample. This accounts for sample-specific effect (such as differences in library sizes and effects of compositionality). If not accounted for, these effects prevent a comparison between the samples. 


```r
y <- calcNormFactors(y)
```
Note that this function does not transform the count data, but rather generates a normalization factor for each sample which is incorporated into the subsequent analysis separately. 


**3. Estimation of dispersion** 

The following code estimates common and tagwise dispersion, namely the variation of the true abundance of a given gene between different samples. This is required to assess differential expression realistically. 


```r
y <- estimateDisp(y)
```


**4. Test for differential expression**



```r
# test each gene for differential expression: 
DE_results_edgeR_Entrez <- exactTest(y)

# extract pre-specified (n) number of genes
DE_results_edgeR_Entrez <- topTags(DE_results_edgeR_Entrez, 
                                   n = nrow(DE_results_edgeR_Entrez))
```

Note that argument 'n' specifies the number of top differentially expressed genes to be displayed in the results. Setting it to 'n = nrow(DE_results_Entrez)' ensures the results of all genes whose differential expression was assessed are displayed. 


**5. Rename columns in the results of differential expression analysis**



```r
# first, we transform the results table to a data frame so that we see the results
# table directly when accessing it through the name "res"
DE_results_edgeR_Entrez <- as.data.frame(DE_results_edgeR_Entrez)

# rename column that contains adjusted p-values: rename padj to p_adj
DE_results_edgeR_Entrez <- dplyr::rename(DE_results_edgeR_Entrez, 
                                         p_adj = FDR)
```

**4. Take a look at top 10 genes in the results table generated with edgeR**

```r
head(DE_results_edgeR_Entrez, n = 20)
#>                logFC   logCPM        PValue         p_adj
#> 6192      10.1330599 6.154999  0.000000e+00  0.000000e+00
#> 9087       5.5791707 2.519656 6.188599e-121 1.860602e-117
#> 8228      -0.9448255 5.485319  1.821705e-17  3.651303e-14
#> 8226      -0.7768228 2.571002  3.067361e-06  4.611011e-03
#> 6615      -1.6460041 3.288135  4.698096e-06  5.627893e-03
#> 10912     -1.5870378 2.679647  5.615726e-06  5.627893e-03
#> 284486    -1.1257052 2.302952  1.949732e-05  1.674820e-02
#> 1521       1.5019085 3.365855  2.906843e-05  2.184856e-02
#> 55787     -0.4253221 5.207647  3.522694e-05  2.269295e-02
#> 10732     -0.6028374 5.756558  3.773981e-05  2.269295e-02
#> 100191040 -1.1506301 3.121643  4.727503e-05  2.584225e-02
#> 412       -0.7312091 6.861545  1.064460e-04  5.333830e-02
#> 159013    -0.5219775 3.911283  1.225919e-04  5.670349e-02
#> 317        0.4048454 6.062183  2.164067e-04  9.294666e-02
#> 51655     -1.3164084 3.094734  2.570863e-04  1.003474e-01
#> 6662      -0.9774711 8.100796  2.804477e-04  1.003474e-01
#> 9404       0.3361246 9.579044  2.837030e-04  1.003474e-01
#> 171177    -0.6747490 6.350249  3.161790e-04  1.031954e-01
#> 10479      0.3635303 6.142157  3.320976e-04  1.031954e-01
#> 3725      -0.5320863 8.429198  3.593785e-04  1.031954e-01
```



