
# Cluster Profiler: GSEA GO


In this script, we will do the following two things: 
- 1. Based on the results of differential expression analysis from voom/limma, DESeq2, and edgeR, we will go through 
    all steps required to run clusterProfiler's GSEA with gene set database GO.
- 2. We will go through all (meaningful) researchers' degrees of freedom.

Note that we provide a separate file for the gene set databases GO and KEGG since clusterProfiler differs in the accepted gene ID formats for them. While for GO, a variety of formats, including Ensembl ID, is accepted, the gene IDs must be converted to NCBI (Entrez) ID when choosing KEGG.  

### Libraries 

All necessary packages are available on Bioconductor, and should be installed from there if not already available on your machine.


```r
install.packages("BiocManager")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("DESeq2")
```

#### Load Libraries 

```r

library(clusterProfiler)
library(org.Hs.eg.db)
library(DESeq2)
```

Description of the packages 

- clusterProfiler: We here use clusterProfiler's implementation of GSEA.

- org.Hs.eg.db: Provides genome-wide annotation for humans. When working with a different organism,  the user has to provide a different package ('Overview_Packages'). Note that this library is required when running clusterProfiler's GSEA with gene set database GO. 

- DESeq2: In contrast detecting the differential expressed genes, which is done solely based on adjusted $p$-vaules, we additionally use the estimated log fold changes when generating the gene ranking. In DESeq2, the estimated log fold change values are shrunken using an additional function, hence loading the library.  



### Load the data 
Note that, as described above, clusterProfiler accepts Ensembl IDs when working with gene set database GO. When generating the ranking of genes with DESeq2, a shrinkage of the log fold change estimates is performed which requires the output of function DESeq(), here named "dds" (see "Instructions_Differential_Expression_Analysis"). 

```r
load("./data/Results_Differential_Expression_Analysis/DE_results_limma_Ensembl.Rdata")
load("./data/Results_Differential_Expression_Analysis/DE_results_DESeq2_Ensembl.Rdata")
load("./data/Results_Differential_Expression_Analysis/DE_results_edgeR_Ensembl.Rdata")

# load dds object 
load("./data/Results_Differential_Expression_Analysis/dds_Ensembl.Rdata")
```

## Run GSEA 


### step 1: Generation of Required Input Object

In the following, we illustrate the creation of the required input object using the three different methods for differential expression analysis, voom/limma, DESeq2, and edgeR, which shall serve as three options. The required input object for clusterProfiler is an "order ranked geneList", i.e. a vector with a ranking value for each gene, named with the respective gene IDs and ordered in a descending manner. 
    

#### 1.1 Option 1: Generate required input using voom/limma


The formula for generating the gene ranking from the results table of differential expression analysis is $(-1) * log_{10}(p\text{-value}) * \text{sign}(\text{log fold change})$. Note that "$p$-value" denotes the non-adjusted p-value. 


```r

# 1. Subset the gene expression data set to those genes that have a p-value (i.e.
# which have been NOT been excluded from differential expression analysis)

# indicate those genes WITH a p-value
ind_nonNA_pvalue_limma_Ensembl <- !is.na(DE_results_limma_Ensembl$P.Value)

# subset gene expression data set to those genes with a p-value
DE_results_noNA_Ensembl <- DE_results_limma_Ensembl[ind_nonNA_pvalue_limma_Ensembl, ]

# 2. create vector that contains the value of the gene-level ranking metric for each gene
rankvec_limma_Ensembl <- sign(DE_results_noNA_Ensembl$logFC)*(-1)*log10(DE_results_noNA_Ensembl$P.Value)

# 3. assign respective gene ID to each value in the vector
names(rankvec_limma_Ensembl) <- rownames(DE_results_noNA_Ensembl)

# 4. sort the vector in descending order
rankvec_limma_Ensembl <- sort(rankvec_limma_Ensembl, decreasing=TRUE)

# inspect the first entries of the ranking vector
head(rankvec_limma_Ensembl , n = 10)
#> ENSG00000129824 ENSG00000099749 ENSG00000154620 
#>       54.973990       50.395302       38.292623 
#> ENSG00000120868 ENSG00000110031 ENSG00000151445 
#>        3.569036        3.262283        3.222725 
#> ENSG00000065308 ENSG00000095015 ENSG00000198689 
#>        3.198940        3.033258        2.949591 
#> ENSG00000102893 
#>        2.915174
```

#### 1.2 Option 2: Generate required input using DESeq2

##### (i) Perform Differential Expression Analysis

This step is performed in the script Instructions_Differential_Expression_Analysis.R. Note that we must exit the DESeq2 workflow from the script a few steps early since there is an additional step required in for the creation of the required input object. We proceed with the object "dds_Ensembl" we've obtained in step 3 of DESeq2's workflow for differential expression analysis.


```r
dds_Ensembl
#> class: DESeqDataSet 
#> dim: 10151 69 
#> metadata(1): version
#> assays(6): counts mu ... replaceCounts replaceCooks
#> rownames(10151): ENSG00000000419 ENSG00000000457 ...
#>   ENSG00000254245 ENSG00000254290
#> rowData names(23): baseMean baseVar ... maxCooks
#>   replace
#> colnames(69): NA18486 NA18498 ... NA19239 NA19257
#> colData names(3): condition sizeFactor replaceable
```



##### (ii) Shrink estimated log fold change values

The intuition behind shrinkage of log fold change values is that, as mentioned in the paper, RNA-Seq data consists (in its raw form) of count data in which is inherently heteroscedastic, i.e. the variance of the count data depends on the mean count of the count data. It is observable that ratios between counts are considerably noisier in low magnitudes of counts compared to higher magnitudes, i.e. the log fold changes between both conditions are higher if the overall magnitude of counts is lower.

DESeq2 addresses this issue by shrinking the estimated log fold changes in the direction of 0 the magnitude of shrinkage is higher if the available information for a gene is lower (which may be because of a low magnitude of counts, a high dispersion
or few degrees of freedom.). A more detailed description is provided in the DESeq2 paper by Love et al. (2014).

Perform shrinkage:

```r
  DE_results_DESeq2_shrink_Ensembl <- lfcShrink(dds_Ensembl,
                                               coef = "condition_treated_vs_untreated",
                                               type="apeglm")
```

Function arguments:

- type: method to perform shrinkage. We opt for the default "apeglm" but you can choose from two alternative options
 as well

- coef: indicate the coefficients to be shrunk. We can obtain the right argument from the following function call:


```r
resultsNames(dds_Ensembl)
#> [1] "Intercept"                     
#> [2] "condition_treated_vs_untreated"
```

This command shows us that we can either shrink the intercept or the "condition_treated_vs_untreated". Since we do not want to shrink the intercept but the estimated log fold changes, we opt for the second option "condition_treated_vs_untreated".


##### (iii) Generate the ranked list of the genes based on the results of differential expression analysis

The formula for generating the gene ranking from the results table of differential expression analysis is $ (-1) * log_{10}(p\text{-value}) * \text{sign}(\text{log fold change})$. Note that "$p$-value" denotes the non-adjusted p-value. 



```r
# 1. Subset the gene expression data set to those genes that have a p-value (i.e.
# which have been NOT been excluded from differential expression analysis)

# indicate those genes WITH a p-value
ind_nonNA_pvalue_Ensembl <- !is.na(DE_results_DESeq2_shrink_Ensembl$pvalue)

# subset gene expression data set to those genes with a p-value
DE_results_noNA_Ensembl <- DE_results_DESeq2_shrink_Ensembl[ind_nonNA_pvalue_Ensembl, ]

# 2. create vector that contains the value of the gene-level ranking metric for each gene
rankvec_DESeq2_Ensembl <- sign(DE_results_noNA_Ensembl$log2FoldChange)*(-1)*log10(DE_results_noNA_Ensembl$pvalue)

# 3. assign respective gene ID to each value in the vector
names(rankvec_DESeq2_Ensembl) <- rownames(DE_results_noNA_Ensembl)

# 4. sort the vector in descending order
rankvec_DESeq2_Ensembl <- sort(rankvec_DESeq2_Ensembl, decreasing=TRUE)


# inspect the first entries of the ranking vector
head(rankvec_DESeq2_Ensembl , n = 10)
#> ENSG00000129824 ENSG00000099749 ENSG00000198692 
#>      204.454993      105.507250       75.109639 
#> ENSG00000154620 ENSG00000157828 ENSG00000183878 
#>       64.872131       39.560339        7.260319 
#> ENSG00000143632 ENSG00000172543 ENSG00000082397 
#>        4.867166        4.469387        3.999177 
#> ENSG00000107731 
#>        3.973429
```


### 1.3 Option 3: Generate required input using edgeR

The formula for generating the gene ranking from the results table of differential expression analysis is $(-1) * log_{10}(p\text{-value}) * \text{sign}(\text{log fold change})$. Note that "$p$-value" denotes the non-adjusted p-value. 



```r
# 1. Subset the gene expression data set to those genes that have a p-value (i.e.
# which have been NOT been excluded from differential expression analysis)

# indicate those genes WITH a p-value
ind_nonNA_pvalue_Ensembl <- !is.na(DE_results_edgeR_Ensembl$PValue)

# subset gene expression data set to those genes with a p-value
DE_results_noNA_Ensembl <- DE_results_edgeR_Ensembl[ind_nonNA_pvalue_Ensembl, ]

# 2. create vector that contains the value of the gene-level ranking metric for each gene
rankvec_edgeR_Ensembl <- sign(DE_results_noNA_Ensembl$logFC)*(-1)*log10(DE_results_noNA_Ensembl$p_adj)

# 3. assign respective gene ID to each value in the vector
names(rankvec_edgeR_Ensembl) <- rownames(DE_results_noNA_Ensembl)

# 4. sort the vector in descending order
rankvec_edgeR_Ensembl <- sort(rankvec_edgeR_Ensembl, decreasing=TRUE)

# inspect ranking vector: 
head(rankvec_edgeR_Ensembl, n = 10)
#> ENSG00000129824 ENSG00000099749 ENSG00000154620 
#>             Inf     192.9498029     118.5930433 
#> ENSG00000172543 ENSG00000120868 ENSG00000110031 
#>       1.7058094       1.0563194       1.0370201 
#> ENSG00000198689 ENSG00000065308 ENSG00000095015 
#>       0.9865108       0.8790009       0.8790009 
#> ENSG00000228253 
#>       0.8654283

# 5. special problem here: gene ENSG00000129824 has a ranking value Inf since its adjusted
# p-value in the results table of differential expression analysis amounts to 0

# here, we deal with this issue by resetting this ranking value to the highest ranking value
# that occurs among the remaining genes
# -> note that there is NO common way of dealing with this issue
rankvec_edgeR_Ensembl[rankvec_edgeR_Ensembl == Inf] <- max(rankvec_edgeR_Ensembl[rankvec_edgeR_Ensembl != Inf])

# inspect ranking vector once more: 
head(rankvec_edgeR_Ensembl, n = 10)
#> ENSG00000129824 ENSG00000099749 ENSG00000154620 
#>     192.9498029     192.9498029     118.5930433 
#> ENSG00000172543 ENSG00000120868 ENSG00000110031 
#>       1.7058094       1.0563194       1.0370201 
#> ENSG00000198689 ENSG00000065308 ENSG00000095015 
#>       0.9865108       0.8790009       0.8790009 
#> ENSG00000228253 
#>       0.8654283
```


### step 2: Run GSEA with gene set database GO

For the purpose of simplicity, we work with the ranking generated using limma. However, you can switch around at your discretion.


```r
rankvec_Ensembl <- rankvec_limma_Ensembl

# alternatively: 
#rankvec_Ensembl <- rankvec_DESeq2_Ensembl
#rankvec_Ensembl <- rankvec_edgeR_Ensembl
```

Run GSEA with gene set database GO. Note that it is important to set a seed for reproducibility.


```r

# important: set seed for reprodubibility
set.seed(1)

GSEA_GO <- gseGO(geneList = rankvec_Ensembl,
                   ont = "BP",
                   OrgDb = org.Hs.eg.db,
                   keyType = "ENSEMBL",
                   seed = TRUE) # set seed for reproducibility


# inspect the corresponding results 
head(GSEA_GO, n = 10 )
#>                    ID
#> GO:0002694 GO:0002694
#> GO:0050865 GO:0050865
#> GO:0045321 GO:0045321
#> GO:0050863 GO:0050863
#> GO:0001775 GO:0001775
#> GO:0051249 GO:0051249
#> GO:0050867 GO:0050867
#> GO:0051251 GO:0051251
#> GO:0002696 GO:0002696
#> GO:0098609 GO:0098609
#>                                             Description
#> GO:0002694           regulation of leukocyte activation
#> GO:0050865                regulation of cell activation
#> GO:0045321                         leukocyte activation
#> GO:0050863              regulation of T cell activation
#> GO:0001775                              cell activation
#> GO:0051249          regulation of lymphocyte activation
#> GO:0050867       positive regulation of cell activation
#> GO:0051251 positive regulation of lymphocyte activation
#> GO:0002696  positive regulation of leukocyte activation
#> GO:0098609                           cell-cell adhesion
#>            setSize enrichmentScore       NES       pvalue
#> GO:0002694     219      -0.4538535 -1.871372 2.478751e-07
#> GO:0050865     232      -0.4464877 -1.843713 1.990043e-07
#> GO:0045321     334      -0.4142950 -1.780532 9.811474e-08
#> GO:0050863     146      -0.5022111 -1.977090 4.610938e-07
#> GO:0001775     377      -0.3951803 -1.717809 4.099930e-07
#> GO:0051249     185      -0.4665593 -1.899849 5.661556e-07
#> GO:0050867     149      -0.4915790 -1.941660 7.600078e-07
#> GO:0051251     130      -0.5062219 -1.963565 8.874654e-07
#> GO:0002696     146      -0.4938384 -1.944128 1.180204e-06
#> GO:0098609     268      -0.4192255 -1.756435 1.120187e-06
#>                p.adjust       qvalue rank
#> GO:0002694 0.0003077782 0.0002683139  841
#> GO:0050865 0.0003077782 0.0002683139  778
#> GO:0045321 0.0003077782 0.0002683139  841
#> GO:0050863 0.0003435149 0.0002994683  841
#> GO:0001775 0.0003435149 0.0002994683  782
#> GO:0051249 0.0003514882 0.0003064193  910
#> GO:0050867 0.0004044327 0.0003525751 1417
#> GO:0051251 0.0004132261 0.0003602409 1375
#> GO:0002696 0.0004396259 0.0003832557 1417
#> GO:0098609 0.0004396259 0.0003832557  778
#>                              leading_edge
#> GO:0002694 tags=26%, list=13%, signal=23%
#> GO:0050865 tags=24%, list=12%, signal=22%
#> GO:0045321 tags=23%, list=13%, signal=21%
#> GO:0050863 tags=29%, list=13%, signal=25%
#> GO:0001775 tags=22%, list=13%, signal=20%
#> GO:0051249 tags=28%, list=15%, signal=24%
#> GO:0050867 tags=40%, list=23%, signal=32%
#> GO:0051251 tags=42%, list=22%, signal=33%
#> GO:0002696 tags=41%, list=23%, signal=33%
#> GO:0098609 tags=24%, list=12%, signal=22%
#>                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            core_enrichment
#> GO:0002694                                                                                                                                                                                                                                                                                                                                                                                                                                 ENSG00000121594/ENSG00000125735/ENSG00000095059/ENSG00000116670/ENSG00000126353/ENSG00000183655/ENSG00000204287/ENSG00000171791/ENSG00000237541/ENSG00000172977/ENSG00000091972/ENSG00000179388/ENSG00000232629/ENSG00000198053/ENSG00000164430/ENSG00000127666/ENSG00000112149/ENSG00000104951/ENSG00000183624/ENSG00000126368/ENSG00000131981/ENSG00000119508/ENSG00000090776/ENSG00000125726/ENSG00000108861/ENSG00000117090/ENSG00000197329/ENSG00000204592/ENSG00000107338/ENSG00000213689/ENSG00000196839/ENSG00000159958/ENSG00000163874/ENSG00000198435/ENSG00000167775/ENSG00000160223/ENSG00000113302/ENSG00000176170/ENSG00000232810/ENSG00000223496/ENSG00000167604/ENSG00000130164/ENSG00000101017/ENSG00000110944/ENSG00000105246/ENSG00000002330/ENSG00000143507/ENSG00000156127/ENSG00000125657/ENSG00000185338/ENSG00000171223/ENSG00000130522/ENSG00000172216/ENSG00000105974/ENSG00000101665/ENSG00000115008
#> GO:0050865                                                                                                                                                                                                                                                                                                                                                                                                                                 ENSG00000095059/ENSG00000116670/ENSG00000126353/ENSG00000183655/ENSG00000154127/ENSG00000204287/ENSG00000171791/ENSG00000237541/ENSG00000172977/ENSG00000091972/ENSG00000179388/ENSG00000232629/ENSG00000198053/ENSG00000164430/ENSG00000127666/ENSG00000112149/ENSG00000104951/ENSG00000183624/ENSG00000126368/ENSG00000131981/ENSG00000119508/ENSG00000090776/ENSG00000125726/ENSG00000108861/ENSG00000117090/ENSG00000197329/ENSG00000204592/ENSG00000107338/ENSG00000213689/ENSG00000010278/ENSG00000196839/ENSG00000159958/ENSG00000163874/ENSG00000198435/ENSG00000167775/ENSG00000160223/ENSG00000113302/ENSG00000176170/ENSG00000232810/ENSG00000223496/ENSG00000167604/ENSG00000130164/ENSG00000101017/ENSG00000110944/ENSG00000105246/ENSG00000002330/ENSG00000143507/ENSG00000156127/ENSG00000125657/ENSG00000185338/ENSG00000171223/ENSG00000130522/ENSG00000172216/ENSG00000105974/ENSG00000101665/ENSG00000115008
#> GO:0045321                                                                                 ENSG00000121594/ENSG00000125735/ENSG00000177697/ENSG00000120738/ENSG00000095059/ENSG00000116670/ENSG00000126353/ENSG00000183655/ENSG00000204287/ENSG00000197405/ENSG00000171791/ENSG00000237541/ENSG00000172977/ENSG00000091972/ENSG00000215021/ENSG00000179388/ENSG00000232629/ENSG00000155760/ENSG00000198053/ENSG00000164430/ENSG00000127666/ENSG00000112149/ENSG00000104951/ENSG00000183624/ENSG00000168264/ENSG00000012061/ENSG00000205220/ENSG00000126368/ENSG00000131981/ENSG00000119508/ENSG00000090776/ENSG00000090339/ENSG00000125726/ENSG00000153879/ENSG00000108861/ENSG00000117090/ENSG00000197329/ENSG00000105438/ENSG00000204592/ENSG00000107338/ENSG00000213689/ENSG00000168447/ENSG00000196839/ENSG00000071564/ENSG00000184371/ENSG00000159958/ENSG00000131188/ENSG00000163874/ENSG00000198435/ENSG00000105374/ENSG00000167775/ENSG00000160223/ENSG00000113302/ENSG00000176170/ENSG00000232810/ENSG00000103522/ENSG00000223496/ENSG00000167604/ENSG00000130164/ENSG00000078902/ENSG00000140968/ENSG00000101017/ENSG00000110944/ENSG00000105246/ENSG00000002330/ENSG00000143507/ENSG00000156127/ENSG00000125657/ENSG00000160683/ENSG00000185338/ENSG00000171223/ENSG00000130522/ENSG00000172216/ENSG00000105974/ENSG00000101665/ENSG00000115008/ENSG00000177606
#> GO:0050863                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 ENSG00000121594/ENSG00000125735/ENSG00000095059/ENSG00000126353/ENSG00000183655/ENSG00000204287/ENSG00000237541/ENSG00000172977/ENSG00000179388/ENSG00000232629/ENSG00000198053/ENSG00000164430/ENSG00000112149/ENSG00000104951/ENSG00000131981/ENSG00000090776/ENSG00000125726/ENSG00000108861/ENSG00000117090/ENSG00000197329/ENSG00000204592/ENSG00000107338/ENSG00000213689/ENSG00000196839/ENSG00000159958/ENSG00000163874/ENSG00000198435/ENSG00000160223/ENSG00000113302/ENSG00000167604/ENSG00000110944/ENSG00000105246/ENSG00000002330/ENSG00000143507/ENSG00000156127/ENSG00000125657/ENSG00000185338/ENSG00000171223/ENSG00000172216/ENSG00000105974/ENSG00000101665/ENSG00000115008
#> GO:0001775 ENSG00000121594/ENSG00000125735/ENSG00000177697/ENSG00000120738/ENSG00000095059/ENSG00000116670/ENSG00000126353/ENSG00000183655/ENSG00000154127/ENSG00000204287/ENSG00000197405/ENSG00000171791/ENSG00000237541/ENSG00000172977/ENSG00000091972/ENSG00000215021/ENSG00000179388/ENSG00000232629/ENSG00000155760/ENSG00000198053/ENSG00000164430/ENSG00000127666/ENSG00000112149/ENSG00000135124/ENSG00000104951/ENSG00000183624/ENSG00000168264/ENSG00000012061/ENSG00000205220/ENSG00000126368/ENSG00000131981/ENSG00000119508/ENSG00000090776/ENSG00000090339/ENSG00000125726/ENSG00000153879/ENSG00000108861/ENSG00000117090/ENSG00000197329/ENSG00000105438/ENSG00000204592/ENSG00000107338/ENSG00000213689/ENSG00000010278/ENSG00000168447/ENSG00000196839/ENSG00000071564/ENSG00000184371/ENSG00000159958/ENSG00000131188/ENSG00000163874/ENSG00000198435/ENSG00000105374/ENSG00000167775/ENSG00000160223/ENSG00000113302/ENSG00000176170/ENSG00000232810/ENSG00000103522/ENSG00000223496/ENSG00000167604/ENSG00000130164/ENSG00000078902/ENSG00000140968/ENSG00000101017/ENSG00000110944/ENSG00000105246/ENSG00000002330/ENSG00000143507/ENSG00000156127/ENSG00000125657/ENSG00000186222/ENSG00000160683/ENSG00000106211/ENSG00000185338/ENSG00000171223/ENSG00000130522/ENSG00000172216/ENSG00000105974/ENSG00000101665/ENSG00000115008/ENSG00000177606
#> GO:0051249                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 ENSG00000196126/ENSG00000189114/ENSG00000121594/ENSG00000125735/ENSG00000095059/ENSG00000116670/ENSG00000126353/ENSG00000183655/ENSG00000204287/ENSG00000171791/ENSG00000237541/ENSG00000172977/ENSG00000179388/ENSG00000232629/ENSG00000198053/ENSG00000164430/ENSG00000127666/ENSG00000112149/ENSG00000104951/ENSG00000183624/ENSG00000131981/ENSG00000090776/ENSG00000125726/ENSG00000108861/ENSG00000117090/ENSG00000197329/ENSG00000204592/ENSG00000107338/ENSG00000213689/ENSG00000196839/ENSG00000159958/ENSG00000163874/ENSG00000198435/ENSG00000167775/ENSG00000160223/ENSG00000113302/ENSG00000223496/ENSG00000167604/ENSG00000101017/ENSG00000110944/ENSG00000105246/ENSG00000002330/ENSG00000143507/ENSG00000156127/ENSG00000125657/ENSG00000185338/ENSG00000171223/ENSG00000172216/ENSG00000105974/ENSG00000101665/ENSG00000115008
#> GO:0050867                                                                                                                                                                                                                                                                                                                                                                 ENSG00000113811/ENSG00000174130/ENSG00000186827/ENSG00000130475/ENSG00000103653/ENSG00000185950/ENSG00000019582/ENSG00000154096/ENSG00000107968/ENSG00000149273/ENSG00000135426/ENSG00000033327/ENSG00000198502/ENSG00000105639/ENSG00000135334/ENSG00000168040/ENSG00000196126/ENSG00000189114/ENSG00000121594/ENSG00000125735/ENSG00000095059/ENSG00000116670/ENSG00000126353/ENSG00000183655/ENSG00000204287/ENSG00000171791/ENSG00000237541/ENSG00000172977/ENSG00000179388/ENSG00000232629/ENSG00000198053/ENSG00000127666/ENSG00000112149/ENSG00000104951/ENSG00000183624/ENSG00000119508/ENSG00000090776/ENSG00000125726/ENSG00000117090/ENSG00000197329/ENSG00000204592/ENSG00000107338/ENSG00000196839/ENSG00000159958/ENSG00000167775/ENSG00000160223/ENSG00000113302/ENSG00000232810/ENSG00000223496/ENSG00000167604/ENSG00000101017/ENSG00000110944/ENSG00000105246/ENSG00000002330/ENSG00000143507/ENSG00000125657/ENSG00000185338/ENSG00000130522/ENSG00000105974/ENSG00000115008
#> GO:0051251                                                                                                                                                                                                                                                                                                                                                                                                                                                                 ENSG00000186827/ENSG00000130475/ENSG00000103653/ENSG00000185950/ENSG00000019582/ENSG00000154096/ENSG00000107968/ENSG00000149273/ENSG00000135426/ENSG00000198502/ENSG00000105639/ENSG00000135334/ENSG00000168040/ENSG00000196126/ENSG00000189114/ENSG00000121594/ENSG00000125735/ENSG00000095059/ENSG00000116670/ENSG00000126353/ENSG00000183655/ENSG00000204287/ENSG00000171791/ENSG00000237541/ENSG00000172977/ENSG00000179388/ENSG00000232629/ENSG00000198053/ENSG00000127666/ENSG00000112149/ENSG00000104951/ENSG00000183624/ENSG00000090776/ENSG00000125726/ENSG00000117090/ENSG00000197329/ENSG00000204592/ENSG00000107338/ENSG00000196839/ENSG00000159958/ENSG00000167775/ENSG00000160223/ENSG00000113302/ENSG00000223496/ENSG00000167604/ENSG00000101017/ENSG00000110944/ENSG00000105246/ENSG00000002330/ENSG00000143507/ENSG00000125657/ENSG00000185338/ENSG00000105974/ENSG00000115008
#> GO:0002696                                                                                                                                                                                                                                                                                                                                                                 ENSG00000113811/ENSG00000174130/ENSG00000186827/ENSG00000130475/ENSG00000103653/ENSG00000185950/ENSG00000019582/ENSG00000154096/ENSG00000107968/ENSG00000149273/ENSG00000135426/ENSG00000033327/ENSG00000198502/ENSG00000105639/ENSG00000135334/ENSG00000168040/ENSG00000196126/ENSG00000189114/ENSG00000121594/ENSG00000125735/ENSG00000095059/ENSG00000116670/ENSG00000126353/ENSG00000183655/ENSG00000204287/ENSG00000171791/ENSG00000237541/ENSG00000172977/ENSG00000179388/ENSG00000232629/ENSG00000198053/ENSG00000127666/ENSG00000112149/ENSG00000104951/ENSG00000183624/ENSG00000119508/ENSG00000090776/ENSG00000125726/ENSG00000117090/ENSG00000197329/ENSG00000204592/ENSG00000107338/ENSG00000196839/ENSG00000159958/ENSG00000167775/ENSG00000160223/ENSG00000113302/ENSG00000232810/ENSG00000223496/ENSG00000167604/ENSG00000101017/ENSG00000110944/ENSG00000105246/ENSG00000002330/ENSG00000143507/ENSG00000125657/ENSG00000185338/ENSG00000130522/ENSG00000105974/ENSG00000115008
#> GO:0098609                                                                                                                                                                                                                                                                                                 ENSG00000095059/ENSG00000116670/ENSG00000126353/ENSG00000148948/ENSG00000183655/ENSG00000154127/ENSG00000204287/ENSG00000171791/ENSG00000237541/ENSG00000140678/ENSG00000172977/ENSG00000172270/ENSG00000171236/ENSG00000091972/ENSG00000179388/ENSG00000232629/ENSG00000198053/ENSG00000154134/ENSG00000108622/ENSG00000112149/ENSG00000104951/ENSG00000129474/ENSG00000105376/ENSG00000161638/ENSG00000131981/ENSG00000119508/ENSG00000090776/ENSG00000090339/ENSG00000125726/ENSG00000108861/ENSG00000164442/ENSG00000117090/ENSG00000197329/ENSG00000130669/ENSG00000204592/ENSG00000179820/ENSG00000107338/ENSG00000010278/ENSG00000196839/ENSG00000159958/ENSG00000163874/ENSG00000198435/ENSG00000253958/ENSG00000160223/ENSG00000113302/ENSG00000167642/ENSG00000232810/ENSG00000167604/ENSG00000169992/ENSG00000076826/ENSG00000110944/ENSG00000105246/ENSG00000002330/ENSG00000125398/ENSG00000143507/ENSG00000125657/ENSG00000186222/ENSG00000171840/ENSG00000106211/ENSG00000185338/ENSG00000172216/ENSG00000105974/ENSG00000101665/ENSG00000115008
```

Arguments of function gseGO():

- geneList: vector of gene ranking (generated in manner of step 1)

- ont: subontology of GO (must be one of "BP", "CC", "MF")

- OrgDb: indicates organism from which the gene expression measurements are taken

- keyTypes: for our purposes here, argument keyType indicates the gene ID format of the gene names in vector rankvec
- available keytypes can be found in the following: 



```r
keytypes(org.Hs.eg.db)
#>  [1] "ACCNUM"       "ALIAS"        "ENSEMBL"     
#>  [4] "ENSEMBLPROT"  "ENSEMBLTRANS" "ENTREZID"    
#>  [7] "ENZYME"       "EVIDENCE"     "EVIDENCEALL" 
#> [10] "GENENAME"     "GENETYPE"     "GO"          
#> [13] "GOALL"        "IPI"          "MAP"         
#> [16] "OMIM"         "ONTOLOGY"     "ONTOLOGYALL" 
#> [19] "PATH"         "PFAM"         "PMID"        
#> [22] "PROSITE"      "REFSEQ"       "SYMBOL"      
#> [25] "UCSCKG"       "UNIPROT"
```

- seed = TRUE: set the seed you have indicated above (here: seed = 1 )


- additional arguments (these do not appear in the above code lines since we keep the default options)

- exponent 1: in the calculation of the enrichment score, each gene is weighted by its absolute value of the ranking metric

- eps = 1e-10: each p-value that is smaller than 1e-10 is indicated as 1e-10

- pvalueCutoff = 0.05: only those gene sets with a p-value < 0.05 are indicated in the results. Note that if you want to inspect ALL gene sets in the results, set pvalueCutoff = 1

- pAdjustMethod = "BH": adjustment for multiple testing using the method by Benjamini and Hochberg

### step 3: Interpretation of the results 

Columns in the results table:

- ID: ID of gene set
- Description: Description of Gene Set

- setSize: Size of the gene set

- enrichmentScore: Enrichment Score
--> note: this enrichment score has not been normalized for gene set size
--> means: larger gene sets automatically have a bigger (absolute) enrichment score,
-->        independent of actual differential enrichment
--> the raw enrichment score is therefore not comparable between different gene sets

- NES: Normalized version of column enrichmentScore
-> NES can be compared between different gene sets

- pvalue: p-value of enrichment of given gene set
-> note: this raw p-value has not been adjusted for multiple testing
-> therefore: CANNOT be used to assess differential enrichment of a given gene set

- p.adjust: ADJUSTED p-value of a given gene set. Note that this p-value now has been adjusted for multiple testing and can therefore be used to assess differential enrichment. For instance, detect all genes with p.adjust < 0.05 as differentially enriched.

- qvalue: ADJUSTED p-value of a given gene set. Note that a qvalue is the analog to p.adjust, but adjusted for multiple testing
         using a different method.

->> you can either use the column p.adjust or qvalue to assess whether a gene set is differentially enriched or not.

- rank: position in the ranked list of the genes in which the maximum difference between the two sums occurs, i.e. the rank at which the enrichment score is extracted.

- leading_edge:
  + tags: The percentage of genes  in the ranked list that are members of gene set before (for positive enrichment score)
         or after (for negative enrichment score) the position from which the enrichment score is extracted.
  + list:  The percentage of genes in the ranked gene list before (for positive enrichment
           score) or after (for negative enrichment score) the position from which the enrichment score is extracted.
  + signal: enrichment signal strength
           combines the statistics "tags" and "list"

  + core_enrichment: genes that contribute most to the enrichment results

Note that a more detailed information on leading_edge and core_enrichment can be found in the user manual provided for GSEA:
 https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideFrame.html : 
 
- search for Detailed Enrichment Results (section "Leading Edge")
 
- search for CORE ENRICHMENT




## Researchers' degrees of freedom 

In this part, we will go through all parameters that can be adapted in the GOSeq workflow. It is important to note that the intention behind going through the researchers' degrees of freedom is to give you an understanding of what you can do to adapt the given (parameter) setting to the research question. It is even more important to keep in mind that the intention behind going through these flexible parameters is NOT to change them in order to help you obtain the most preferable results by systematically changing these parameters since such behaviour would correspond to "cherry-picking". Any changes in the parameter choice must be documented transparently. 


### change 1: change exponent ###

Note that, in contrast to the web-based tool GSEA, clusterProfiler does not make any indication as to which alternative exponent values are allowed or sensible. oOne could stick with the web-based tool GSEA which suggests, in addition to the default exponent value of 1, the alternative values 0, 1.5, and 2.

Here, we change the exponent value to 2


```r
# important: set seed for reproducibility
set.seed(1)

GSEA_GO_exponent <- gseGO(geneList = rankvec_Ensembl,
                            ont = "BP",
                            OrgDb = org.Hs.eg.db,
                            keyType = "ENSEMBL",
                            exp = 2,
                            seed = TRUE)

# inspect the results
head(GSEA_GO_exponent, n = 10) 
#> [1] ID              Description     setSize        
#> [4] enrichmentScore NES             pvalue         
#> [7] p.adjust        qvalue         
#> <0 rows> (or 0-length row.names)
```






