
# Differential Expression Analysis


## Overview

In this script, we will perform differential expression analysis for each of the three parametric methods: 
- voom/limma
- DESeq2 
- edgeR 

Here, we illustrate this process for two gene ID format Entrez IDs, while we additionally present the identical process for the format Ensembl ID in the R script "Instructions_Differential_Expression_Analysis.R".

Note that at at the end of the illustration for each method, we rename some of the columns in the corresponding results tables to unify them across the different methods. While this step is not required for differential expression analysis or for the subsequent use of gene set analysis IN GENERAL, renaming of some columns (such as adjusted p-value) is necessary in this context so that we can use the same code to illustrate further conduct of GSA in the following R scripts, independent of the method for differential expression analysis used. 


## Libraries 

### Install Libraries

All necessary packages are available on Bioconductor, and should be installed from there if not already available on your machine.



```r
install.packages("BiocManager")
BiocManager::install("tweeDESeqCountData")
BiocManager::install("limma")
BiocManager::install("DESeq2")
BiocManager::install("edgeR")
BiocManager::install("dplyr")
```

### Load libraries


```r
library(tweeDEseqCountData)
library(limma)
library(DESeq2)
library(edgeR)
library(dplyr)
```

Description of the stats libraries: 

<!-- I removed {dplyr} from the explanations because it's already so well-known -->

#### tweeDESeqCountData

Info about the package here...

#### limma

#### DESeq2

#### edgeR



### Load data 

Note that depending on the method for differential expression analysis, different methods for pre-filtering are proposed: 

- voom/limma and egdeR: pre-filtering using edgeR's built-in function filterByExpr()
- DESeq2: manual/simple pre-filtering

 We therefore work with both pre-filtered gene expression data sets generated in "Instructions_PreFiltering": 

#### Previously Generated

```r
# data is not yet in the folder. set to eval=FALSE for now
load("./data/Results_GeneID_Conversion_DuplicatesRemoval/exprdat_filter_conv_DESeq2.Rdata")
load("./data/Results_GeneID_Conversion_DuplicatesRemoval/exprdat_filter_conv_filterByExpr.Rdata")
```


#### {pickrell} Sample Data

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

(Provide brief description of the options here so people know which one to jump to)

### Option 1: Differential expression analysis using limma 

#### step 1: Generate required input object 
As described in "Instructions_Differential_Expression_Analysis.R", the gene expression data set exprdat_filter_conv_filterByExpr is
a data frame and needs to be converted to a DGEList object. 


```r
# store expression data with corresponding sample conditions in object of class DGEList

y <- DGEList(counts = exprdat_filter_conv_filterByExpr, 
              group = sample_conditions)
```

required function arguments: 

- counts: matrix that contains RNA-Seq data (i.e. count data)

optional function arguments:

- group: indicates condition of each sample/library. This argument must be specified sooner or later (such as in subsequent functions) so we just specify it at this point 

Note: we leave the remaining arguments in their default state since the corresponding info will be added through the following pipeline of functions. 

#### step 2: Normalization

The following piece of code generates a normalization factor for each sample which accounts for sample-specific effect (such as differences in library sizes and effects of compositionality). If not accounted for, these effects prevent a comparison between the samples. 


```r
y <- calcNormFactors(y)
```

Note that this function does not transform the count data, but rather generates a normalization factor for each sample which is separately incorporated into the subsequent analysis. 


#### step 3: voom transformation 

Perform the voom transformation to make limma, which was initially developed for microarray measurements, available to RNA-Seq data. voom transforms the RNA-Seq data to log counts-per-million and generates a precision weight for each entry in the gene expression data sets to be incorporated into the limma pipeline. 


```r
# (i) generate design matrix (rows correspond to samples, columns indicate which coefficients are to be estimated)
design_matrix <- model.matrix( ~ sample_conditions)

# (ii) voom transformation 
y <- voom(y, design = design_matrix)
```

#### step 4: Differential expression analysis 


```r
# (i) fit a linear model for each gene
y <- lmFit(y)

# (ii) calculate the statistics for the assessment of differential expression 
y <- eBayes(y)

# (iii) Get the result table for each gene whose differential expression was assessed 
DE_results_limma_Entrez <- topTable(y, number =  nrow(y))
#> Removing intercept from test coefficients
# number = nrow(y) ensures that all genes are displayed in results

```



#### step 5: Rename columns in results of differential expression analysis ######

We rename some of the columns in the results table to unify them across the different methods for differential expression analysis. For this, we use the function rename() provided by the library dplyr. 



```r
# first, we transform the results table to a data frame so that we see the results
# table directly when accessing it through the name "res"
DE_results_limma_Entrez <- as.data.frame(DE_results_limma_Entrez)

# rename column that contains adjusted p-values: rename padj to p_adj
DE_results_limma_Entrez <- dplyr::rename(DE_results_limma_Entrez, p_adj = "adj.P.Val")


# take a look at the results table: 
head(DE_results_limma_Entrez, n = 20)
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
#> 284486 -1.0824842 1.6291231 -3.567807 6.466165e-04
#> 4092   -0.3090285 7.2622871 -3.660091 4.786317e-04
#> 63933  -0.3410894 5.3603234 -3.587703 6.062388e-04
#> 2275   -0.6200465 3.3046963 -3.481030 8.544975e-04
#> 9404    0.3288050 9.5202902  3.575849 6.299981e-04
#> 9697    0.4792856 7.6786140  3.549883 6.851673e-04
#> 10732  -0.5086504 5.5892686 -3.492114 8.247942e-04
#> 4214    0.5165346 4.7408391  3.445856 9.556092e-04
#> 63894   0.1859711 7.0457057  3.483087 8.489115e-04
#> 857    -0.9967185 2.0665254 -3.262960 1.690566e-03
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
#> 284486 2.746607e-01 -0.43082603
#> 4092   2.616375e-01 -0.66791157
#> 63933  2.746607e-01 -0.72462491
#> 2275   2.831770e-01 -0.79305058
#> 9404   2.746607e-01 -0.97053199
#> 9697   2.746607e-01 -1.02933463
#> 10732  2.831770e-01 -1.03194804
#> 4214   2.831770e-01 -1.09133522
#> 63894  2.831770e-01 -1.19771238
#> 857    3.428077e-01 -1.28178906
```

### option 2: Differential expression analysis using DESeq2 

#### step 1: Generate object required by DESeq2 
DESeq2 operates on the format DESeqDataSet which contains information on the count data, the conditions of the samples and the design. Note that for gene expression data sets that must be imported to R, additional steps are necessary before the following code can be run. 

##### Generate a data frame which contains the condition of each sample

Here, the information on the sample conditions is the only column in the data frame. However, if further variables (such as batch effects) are to be controlled for, the corresponding variables must additionally be added to coldata. 



```r
# the names of the samples are stored as the row names of the data frame
# -> important: make sure that the order of the conditions in sample_conditions corresponds to the order of the 
# samples in expression_data
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

##### Generate the DESeqDataSet


```r
dds_Entrez <- DESeqDataSetFromMatrix(countData = exprdat_filter_conv_DESeq2, 
                                     colData = coldata, 
                                     design = ~condition)
```

Relevant arguments in function DESeqDataSetFromMatrix:

- countData: count data from the gene expression data set

- colData: data frame that contains information on the samples (see above) such as the conditions of the samples (required) and possibly further variables to correct for (such as batch effects)

- design: indicates which variables from colData are used for modelling
- -> more detailed: the argument design is used to estimate the dispersions and the log2 fold changes of the model 
- -> if more than one variable from colData are used in argument design (e.g. a second variable "batch"), the syntax changes to the following formula: design ~ batch + condition 
- make sure that the variable of interest (here: variable that represents conditions of the samples) is placed at the end of the formula


#### step 2: Differential expression analysis 

DESeq2 calculates the normalization factor of each sample, estimates the dispersion parameter, fits a negative binomial model. A p-value of differential expression for each gene is then computed using a Wald test. 


```r
# perform default differential expression analysis 
dds_Entrez <- DESeq(dds_Entrez)
#> estimating size factors
#> estimating dispersions
#> gene-wise dispersion estimates
#> mean-dispersion relationship
#> final dispersion estimates
#> fitting model and testing
#> -- replacing outliers and refitting for 117 genes
#> -- DESeq argument 'minReplicatesForReplace' = 7 
#> -- original counts are preserved in counts(dds)
#> estimating dispersions
#> fitting model and testing

# generate results table which provides
DE_results_DESeq2_Entrez <- results(dds_Entrez)
```
Note: function results() provided the base mean across all samples, the estimated log2 fold changes, standard errors, test statistics, p-values, and adjusted p-values for each gene. 


#### step 3: Rename columns in results of differential expression analysis 



```r
# first, we transform to results table to a data frame so that we see the results
# table directly when accessing it through the name "res"
DE_results_DESeq2_Entrez <- as.data.frame(DE_results_DESeq2_Entrez)

# rename column that contains adjusted p-values: rename padj to p_adj
DE_results_DESeq2_Entrez <- dplyr::rename(DE_results_DESeq2_Entrez, p_adj = padj)

# take a look at the results table
head(DE_results_DESeq2_Entrez, n = 20 )
#>           baseMean log2FoldChange      lfcSE         stat
#> 8813    62.7268709   0.0394347840 0.11345444  0.347582551
#> 57147   70.4533427   0.0285330645 0.11131002  0.256338685
#> 55732   13.0010604   0.1641173646 0.15102073  1.086720805
#> 2268    43.5777559   0.1137118045 0.19545410  0.581782664
#> 2519    37.9639983  -0.0065270645 0.10803176 -0.060418014
#> 4800   228.4783293   0.0511474400 0.09641906  0.530470243
#> 81887   11.5340315   0.0251919133 0.18501229  0.136163455
#> 22875  175.8335262   0.2999994619 0.16347952  1.835088963
#> 5893     4.6883080   0.2804257542 0.19291917  1.453591975
#> 572     66.3699746  -0.2648964238 0.10233495 -2.588523618
#> 9957    10.1224447  -0.1994929894 0.31663473 -0.630041397
#> 51384    3.4400718  -0.4871178913 0.50088873 -0.972507196
#> 3927   158.4286342   0.1552048117 0.18239266  0.850937829
#> 29916  843.8205108   0.0009164434 0.11692479  0.007837888
#> 4074  1015.3676532   0.0204941689 0.08518719  0.240578061
#> 90293    0.7530882  -0.1378597175 0.75449224 -0.182718537
#> 79007  924.7932944  -0.4977855669 0.26299623 -1.892747927
#> 6542     0.4013411   0.4437578747 0.73930371  0.600237584
#> 952      5.5336896   0.2786075955 0.32467188  0.858120507
#> 2288   609.7373992  -0.0974512098 0.11711125 -0.832125125
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
#> 9957  0.528667500 0.8707508
#> 51384 0.330798286 0.7871383
#> 3927  0.394803890 0.8165289
#> 29916 0.993746334 0.9984058
#> 4074  0.809882155 0.9537648
#> 90293 0.855018870 0.9662640
#> 79007 0.058391399 0.4994559
#> 6542  0.548347909 0.8785679
#> 952   0.390825925 0.8152442
#> 2288  0.405338321 0.8212059
```

### option 3: Differential expression analysis using edgeR 


#### step 1: generate required input object 

Note that edgeR operates on a DGEList object 



```r
y <- DGEList(counts = exprdat_filter_conv_filterByExpr, 
             group = sample_conditions)
```
Required function arguments: 

- counts: matrix that contains RNA-Seq data (i.e. count data)

Optional function arguments: 

- group: indicates condition of each sample/library 
-  -> this argument must be specified sooner or later (such as in subsequent functions) so we just specify it at this point 

Note thta we leave the remaining arguments in their default state since the corresponding info will be added through the following pipeline of functions. 

#### step 2: Normalization 

The following piece of code generates a normalization factor for each sample. This accounts for sample-specific effect (such as differences in library sizes and effects of compositionality). If not accounted for, these effects prevent a comparison between the samples. 


```r
y <- calcNormFactors(y)
```
Note that this function does not transform the count data, but rather generates a normalization factor for each sample which is incorporated into the subsequent analysis separately. 


#### step 3: Estimation of dispersion 

The following code estimates common and tagwise dispersion, namely the variation of the true abundance of a given gene between different samples. This is required to assess differential expression realistically. 


```r
y <- estimateDisp(y)
#> Using classic mode.
```


#### step 4: Test for differential expression 



```r
# test each gene for differential expression: 
DE_results_edgeR_Entrez <- exactTest(y)

# extract pre-specified (n) number of genes
DE_results_edgeR_Entrez <- topTags(DE_results_edgeR_Entrez, n = nrow(DE_results_edgeR_Entrez))
```

Note that argument n specifies the number of top differentially expressed genes to be displayed in the results. Setting it to n = nrow(DE_results_Entrez) ensures the results of all genes whose differential expression was assessed are displayed. 


#### step 5: Rename columns in the results of differential expression analysis



```r
# first, we transform the results table to a data frame so that we see the results
# table directly when accessing it through the name "res"
DE_results_edgeR_Entrez <- as.data.frame(DE_results_edgeR_Entrez)

# rename column that contains adjusted p-values: rename padj to p_adj
DE_results_edgeR_Entrez <- dplyr::rename(DE_results_edgeR_Entrez, p_adj = FDR)


# take a look at the results table
head(DE_results_edgeR_Entrez, n = 20)
#>                logFC   logCPM        PValue         p_adj
#> 6192      10.1330599 6.154999  0.000000e+00  0.000000e+00
#> 9087       5.5791707 2.519656 6.188598e-121 1.860602e-117
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




