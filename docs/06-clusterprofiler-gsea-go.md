
# clusterProfiler's GSEA (with gene set database GO)



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
BiocManager::install("apeglm")
```

**Load Libraries **

```r
library(clusterProfiler)
library(org.Hs.eg.db)
library(DESeq2)
```

Description of the packages 

- **clusterProfiler**: We here use clusterProfiler's implementation of GSEA.

- **org.Hs.eg.db**: Provides genome-wide annotation for humans. When working with a different organism,  the user has to provide a different package (see chapter 'About'). Note that this library is required when running clusterProfiler's GSEA with gene set database GO. 

- **DESeq2**: In contrast detecting the differential expressed genes, which is done solely based on adjusted $p$-vaules, we additionally use the estimated log fold changes when generating the gene ranking. In DESeq2, the estimated log fold change values are shrunken using an additional function, hence loading the library.  



## Load data 
Note that, as described above, clusterProfiler accepts Ensembl IDs when working with gene set database GO. When generating the ranking of genes with DESeq2, a shrinkage of the log fold change estimates is performed which requires the output of function *DESeq()*, here named **'dds'** (see Chapter 'Differential expression analysis'). 


```r
load("./data/Results_Differential_Expression_Analysis/DE_results_limma_Ensembl.Rdata")
load("./data/Results_Differential_Expression_Analysis/DE_results_DESeq2_Ensembl.Rdata")
load("./data/Results_Differential_Expression_Analysis/DE_results_edgeR_Ensembl.Rdata")

# load dds object 
load("./data/Results_Differential_Expression_Analysis/dds_Ensembl.Rdata")
```

## Run GSEA 


### Step 1: Generation of required input object

In the following, we illustrate the creation of the required input object using the three different methods for differential expression analysis, voom/limma, DESeq2, and edgeR, which shall serve as three options. The required input object for clusterProfiler is an "order ranked geneList", i.e. a vector with a ranking value for each gene, named with the respective gene IDs and ordered in a descending manner. \
Here we choose the following gene-level statistic based on which we can transform the results of differential expression analysis to the gene ranking: 
$$(-1) * \log_{10}(p\text{-value}) * \text{sign}(\text{log fold change})$$. \
Note that "$p$-value" denotes the non-adjusted p-value. 
    

#### Option 1: Generate required input using voom/limma

Before starting, we want to have a short look at the results table of differential expression analysis generated with limma so that we can identify the relevant columns. 


```r
head(DE_results_limma_Ensembl, n = 10)
#>                      logFC   AveExpr         t      P.Value
#> ENSG00000129824  9.2133770 1.9601987 47.132083 1.061719e-55
#> ENSG00000099749  6.1709577 0.4153409 40.465602 4.024370e-51
#> ENSG00000154620  5.0848991 0.4867874 26.729729 5.097737e-39
#> ENSG00000006757 -0.9214380 5.3073131 -8.805751 4.963544e-13
#> ENSG00000130021 -0.8516173 2.1365915 -4.388951 3.846411e-05
#> ENSG00000185753 -0.5479423 3.6902013 -4.217130 7.121081e-05
#> ENSG00000086712 -0.4017993 5.0899555 -4.023482 1.403927e-04
#> ENSG00000123689 -1.4679601 5.1661467 -3.850433 2.537689e-04
#> ENSG00000177606 -0.5844304 8.2469942 -3.910485 2.069770e-04
#> ENSG00000120868  0.4076940 5.9420527  3.832337 2.697515e-04
#>                        p_adj           B
#> ENSG00000129824 6.631497e-52 94.00193786
#> ENSG00000099749 1.256811e-47 86.16860821
#> ENSG00000154620 1.061349e-35 68.53453635
#> ENSG00000006757 7.750575e-10 19.36978341
#> ENSG00000130021 4.804937e-02  2.04582318
#> ENSG00000185753 7.413045e-02  1.30583651
#> ENSG00000086712 1.252704e-01  0.48313770
#> ENSG00000123689 1.684868e-01 -0.04984613
#> ENSG00000177606 1.615973e-01 -0.13323415
#> ENSG00000120868 1.684868e-01 -0.25078873
```
Generate gene ranking: 

```r
# 1. Subset the gene expression data set to those genes that have a p-value (i.e. which have been NOT been excluded from differential expression analysis)

# indicate those genes with an existing p-value
ind_nonNA_pvalue_limma_Ensembl <- !is.na(DE_results_limma_Ensembl$P.Value)

# subset gene expression data set to those genes with an existing p-value
DE_results_noNA_Ensembl <- DE_results_limma_Ensembl[ind_nonNA_pvalue_limma_Ensembl, ]

# 2. create vector that contains the value of the gene-level ranking metric for each gene
rankvec_limma_Ensembl <- sign(DE_results_noNA_Ensembl$logFC)*(-1)*log10(DE_results_noNA_Ensembl$P.Value)

# 3. assign respective gene ID to each value in the vector
names(rankvec_limma_Ensembl) <- rownames(DE_results_noNA_Ensembl)

# 4. sort the vector in descending order
rankvec_limma_Ensembl <- sort(rankvec_limma_Ensembl, decreasing=TRUE)
```

Inspect the first entries of the ranking vector

```r
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


#### Option 2: Generate required input using DESeq2

**Step 1: Perform differential expression analysis**

This step is performed in Chapter 'Differential expression analysis'. Note that we must exit the DESeq2 workflow from the script a few steps early since there is an additional step required in for the creation of the required input object. We proceed with the object "dds_Ensembl" we've obtained in step 2 of DESeq2's workflow for differential expression analysis.


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



**Step 2: Shrink estimated log fold change values**

The intuition behind shrinkage of log fold change values is that, as mentioned in the paper, RNA-Seq data consists (in its raw form) of count data which is inherently heteroscedastic, i.e. the variance of the count data depends on the mean count of the count data. It is observable that ratios between counts are considerably noisier in low magnitudes of counts compared to higher magnitudes, i.e. the log fold changes between both conditions are higher if the overall magnitude of counts is lower. \
DESeq2 addresses this issue by shrinking the estimated log fold changes in the direction of 0 the magnitude of shrinkage is higher if the available information for a gene is lower (which may be because of a low magnitude of counts, a high dispersion or few degrees of freedom.). A more detailed description is provided in the DESeq2 paper by Love et al. (2014) [@love2014moderated].

Perform shrinkage:

```r
DE_results_DESeq2_shrink_Ensembl <- lfcShrink(dds_Ensembl,
                                               coef = "condition_treated_vs_untreated",
                                               type="apeglm")
```

Function arguments:

- **type**: method to perform shrinkage. We opt for the default 'type = apeglm', however, there are two alternative options the user can choose from. 

- **coef**: indicate the coefficients to be shrunk. We can obtain the right argument from the following function call:


```r
resultsNames(dds_Ensembl)
#> [1] "Intercept"                     
#> [2] "condition_treated_vs_untreated"
```

This command shows us that we can either shrink the intercept or the "condition_treated_vs_untreated". Since we do not want to shrink the intercept but the estimated log fold changes, we opt for the second option "condition_treated_vs_untreated".

Before generating the ranking, we want to inspect the results table generated with DESeq2 (incl. shrinkage) so that we can identify the relevant columns


```r
head(DE_results_DESeq2_shrink_Ensembl, n = 10)
#> log2 fold change (MAP): condition treated vs untreated 
#> Wald test p-value: condition treated vs untreated 
#> DataFrame with 10 rows and 5 columns
#>                  baseMean log2FoldChange      lfcSE
#>                 <numeric>      <numeric>  <numeric>
#> ENSG00000000419  62.71639   -1.75964e-05 0.00144263
#> ENSG00000000457  70.46817    2.47225e-06 0.00144258
#> ENSG00000000460  13.00098    7.54395e-06 0.00144264
#> ENSG00000000938  43.56413    7.35510e-06 0.00144267
#> ENSG00000001036  37.95191    1.68099e-06 0.00144257
#> ENSG00000001167 228.48025    6.71863e-07 0.00144253
#> ENSG00000001497  11.52425    7.90345e-07 0.00144265
#> ENSG00000001561 175.93013    1.18172e-05 0.00144266
#> ENSG00000002016   4.68668    7.85799e-06 0.00144267
#> ENSG00000002330  66.34023   -2.62299e-05 0.00144267
#>                     pvalue      padj
#>                  <numeric> <numeric>
#> ENSG00000000419 0.72425162  0.971027
#> ENSG00000000457 0.79118089  0.981324
#> ENSG00000000460 0.27497228  0.820697
#> ENSG00000000938 0.55928481  0.945493
#> ENSG00000001036 0.95622491  0.997305
#> ENSG00000001167 0.59046443  0.953105
#> ENSG00000001497 0.88864354  0.992301
#> ENSG00000001561 0.06610429  0.575449
#> ENSG00000002016 0.14450244  0.703418
#> ENSG00000002330 0.00949353  0.350398
```

Now we proceed to generating the ranking. 

**Step 3: Generate the gene ranking based on the results of differential expression analysis**



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
```

Inspect the first entries of the ranking vector

```r
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


#### Option 3: Generate required input using edgeR

Before generating the ranking, we want inspect the results table of differential expression analysis generated with edgeR to identify the relevant columns.


```r
head(DE_results_edgeR_Ensembl, n = 10)
#>                      logFC   logCPM        PValue
#> ENSG00000129824 10.1353289 6.130860  0.000000e+00
#> ENSG00000099749  7.8142855 2.949875 3.594389e-197
#> ENSG00000154620  5.5800954 2.493788 1.225959e-122
#> ENSG00000006757 -0.9412925 5.457916  1.516865e-17
#> ENSG00000130021 -0.7752345 2.545089  3.626637e-06
#> ENSG00000124216 -1.6423626 3.261210  4.708108e-06
#> ENSG00000130222 -1.5833419 2.652436  5.638958e-06
#> ENSG00000196407 -1.1211411 2.276703  1.959274e-05
#> ENSG00000172543  1.5008703 3.336720  2.836816e-05
#> ENSG00000086712 -0.4218274 5.181013  3.922485e-05
#>                         p_adj
#> ENSG00000129824  0.000000e+00
#> ENSG00000099749 1.122528e-193
#> ENSG00000154620 2.552447e-119
#> ENSG00000006757  2.368584e-14
#> ENSG00000130021  4.530395e-03
#> ENSG00000124216  4.901141e-03
#> ENSG00000130222  5.031562e-03
#> ENSG00000196407  1.529703e-02
#> ENSG00000172543  1.968750e-02
#> ENSG00000086712  2.449984e-02
```

Now we proceed to generating the ranking.

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
```

Inspect the first entries of the ranking vector: 

```r
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


### Step 2: Run GSEA with gene set database GO

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
```

Arguments of function *gseGO()*:

- **geneList**: vector of gene ranking (generated in manner of step 1)

- **ont**: subontology of GO (must be one of "BP", "CC", "MF")

- **OrgDb**: indicates organism from which the gene expression measurements are taken

- **keyTypes**: for our purposes here, argument keyType indicates the gene ID format of the gene names in vector rankvec
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

- **seed = TRUE**: set the seed you have indicated above (here: seed = 1 )


- additional arguments (these do not appear in the above code lines since we keep the default options)

- **exponent = 1**: in the calculation of the enrichment score, each gene is weighted by its absolute value of the ranking metric

- **eps = 1e-10**: each p-value that is smaller than 1e-10 is indicated as 1e-10

- **pvalueCutoff = 0.05**: only those gene sets with a p-value < 0.05 are indicated in the results. Note that if you want to inspect ALL gene sets in the results, set pvalueCutoff = 1

- **pAdjustMethod = "BH"**: adjustment for multiple testing using the method by Benjamini and Hochberg

### Step 3: Interpretation of the results 

**Inspect the results table**

```r
head(GSEA_GO, n = 10 )
#>                    ID
#> GO:0002694 GO:0002694
#> GO:0045321 GO:0045321
#> GO:0001775 GO:0001775
#> GO:0050865 GO:0050865
#> GO:0050863 GO:0050863
#> GO:0051249 GO:0051249
#> GO:0002696 GO:0002696
#> GO:0050867 GO:0050867
#> GO:0046649 GO:0046649
#> GO:0051240 GO:0051240
#>                                                        Description
#> GO:0002694                      regulation of leukocyte activation
#> GO:0045321                                    leukocyte activation
#> GO:0001775                                         cell activation
#> GO:0050865                           regulation of cell activation
#> GO:0050863                         regulation of T cell activation
#> GO:0051249                     regulation of lymphocyte activation
#> GO:0002696             positive regulation of leukocyte activation
#> GO:0050867                  positive regulation of cell activation
#> GO:0046649                                   lymphocyte activation
#> GO:0051240 positive regulation of multicellular organismal process
#>            setSize enrichmentScore       NES       pvalue
#> GO:0002694     218      -0.4550306 -1.885682 1.443823e-07
#> GO:0045321     335      -0.4102866 -1.776603 1.721317e-07
#> GO:0001775     379      -0.3913563 -1.715310 1.627551e-07
#> GO:0050865     231      -0.4475860 -1.862707 3.379430e-07
#> GO:0050863     145      -0.5024266 -1.973006 6.812416e-07
#> GO:0051249     185      -0.4654907 -1.896787 6.286296e-07
#> GO:0002696     145      -0.4964230 -1.949430 1.090979e-06
#> GO:0050867     148      -0.4941648 -1.945621 1.030633e-06
#> GO:0046649     282      -0.4097538 -1.745495 1.526809e-06
#> GO:0051240     490      -0.3558447 -1.593477 1.476326e-06
#>                p.adjust       qvalue rank
#> GO:0002694 0.0002129843 0.0001866874  841
#> GO:0045321 0.0002129843 0.0001866874  841
#> GO:0001775 0.0002129843 0.0001866874  782
#> GO:0050865 0.0003136111 0.0002748900  778
#> GO:0050863 0.0004214615 0.0003694242  841
#> GO:0051249 0.0004214615 0.0003694242  910
#> GO:0002696 0.0005062142 0.0004437126 1417
#> GO:0050867 0.0005062142 0.0004437126 1417
#> GO:0046649 0.0005667516 0.0004967755  910
#> GO:0051240 0.0005667516 0.0004967755  841
#>                              leading_edge
#> GO:0002694 tags=26%, list=13%, signal=23%
#> GO:0045321 tags=23%, list=13%, signal=21%
#> GO:0001775 tags=22%, list=13%, signal=20%
#> GO:0050865 tags=24%, list=12%, signal=22%
#> GO:0050863 tags=29%, list=13%, signal=26%
#> GO:0051249 tags=28%, list=15%, signal=24%
#> GO:0002696 tags=42%, list=23%, signal=33%
#> GO:0050867 tags=41%, list=23%, signal=33%
#> GO:0046649 tags=24%, list=15%, signal=22%
#> GO:0051240 tags=24%, list=13%, signal=23%
#>                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            core_enrichment
#> GO:0002694                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 ENSG00000121594/ENSG00000125735/ENSG00000095059/ENSG00000116670/ENSG00000126353/ENSG00000183655/ENSG00000204287/ENSG00000171791/ENSG00000237541/ENSG00000172977/ENSG00000091972/ENSG00000179388/ENSG00000232629/ENSG00000198053/ENSG00000164430/ENSG00000127666/ENSG00000112149/ENSG00000104951/ENSG00000183624/ENSG00000126368/ENSG00000131981/ENSG00000119508/ENSG00000090776/ENSG00000125726/ENSG00000108861/ENSG00000117090/ENSG00000197329/ENSG00000204592/ENSG00000107338/ENSG00000213689/ENSG00000196839/ENSG00000159958/ENSG00000163874/ENSG00000198435/ENSG00000167775/ENSG00000160223/ENSG00000113302/ENSG00000176170/ENSG00000232810/ENSG00000223496/ENSG00000167604/ENSG00000130164/ENSG00000101017/ENSG00000110944/ENSG00000105246/ENSG00000002330/ENSG00000143507/ENSG00000156127/ENSG00000125657/ENSG00000185338/ENSG00000171223/ENSG00000130522/ENSG00000172216/ENSG00000105974/ENSG00000101665/ENSG00000115008
#> GO:0045321                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 ENSG00000121594/ENSG00000125735/ENSG00000177697/ENSG00000120738/ENSG00000095059/ENSG00000116670/ENSG00000126353/ENSG00000183655/ENSG00000204287/ENSG00000197405/ENSG00000171791/ENSG00000237541/ENSG00000172977/ENSG00000091972/ENSG00000215021/ENSG00000179388/ENSG00000232629/ENSG00000155760/ENSG00000198053/ENSG00000164430/ENSG00000127666/ENSG00000112149/ENSG00000104951/ENSG00000183624/ENSG00000168264/ENSG00000012061/ENSG00000205220/ENSG00000126368/ENSG00000131981/ENSG00000119508/ENSG00000090776/ENSG00000090339/ENSG00000125726/ENSG00000153879/ENSG00000108861/ENSG00000117090/ENSG00000197329/ENSG00000105438/ENSG00000204592/ENSG00000107338/ENSG00000213689/ENSG00000168447/ENSG00000196839/ENSG00000071564/ENSG00000184371/ENSG00000159958/ENSG00000131188/ENSG00000163874/ENSG00000198435/ENSG00000105374/ENSG00000167775/ENSG00000160223/ENSG00000113302/ENSG00000176170/ENSG00000232810/ENSG00000103522/ENSG00000223496/ENSG00000167604/ENSG00000130164/ENSG00000078902/ENSG00000140968/ENSG00000101017/ENSG00000110944/ENSG00000105246/ENSG00000002330/ENSG00000143507/ENSG00000156127/ENSG00000125657/ENSG00000160683/ENSG00000185338/ENSG00000171223/ENSG00000130522/ENSG00000172216/ENSG00000105974/ENSG00000101665/ENSG00000115008/ENSG00000177606
#> GO:0001775                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 ENSG00000121594/ENSG00000125735/ENSG00000177697/ENSG00000120738/ENSG00000095059/ENSG00000116670/ENSG00000126353/ENSG00000183655/ENSG00000154127/ENSG00000204287/ENSG00000197405/ENSG00000171791/ENSG00000237541/ENSG00000172977/ENSG00000091972/ENSG00000215021/ENSG00000179388/ENSG00000232629/ENSG00000155760/ENSG00000198053/ENSG00000164430/ENSG00000127666/ENSG00000112149/ENSG00000135124/ENSG00000104951/ENSG00000183624/ENSG00000168264/ENSG00000012061/ENSG00000205220/ENSG00000126368/ENSG00000131981/ENSG00000119508/ENSG00000090776/ENSG00000090339/ENSG00000125726/ENSG00000153879/ENSG00000108861/ENSG00000117090/ENSG00000197329/ENSG00000105438/ENSG00000204592/ENSG00000107338/ENSG00000213689/ENSG00000010278/ENSG00000168447/ENSG00000196839/ENSG00000071564/ENSG00000184371/ENSG00000159958/ENSG00000131188/ENSG00000163874/ENSG00000198435/ENSG00000105374/ENSG00000167775/ENSG00000160223/ENSG00000113302/ENSG00000176170/ENSG00000232810/ENSG00000103522/ENSG00000223496/ENSG00000167604/ENSG00000130164/ENSG00000078902/ENSG00000140968/ENSG00000101017/ENSG00000110944/ENSG00000105246/ENSG00000002330/ENSG00000143507/ENSG00000156127/ENSG00000125657/ENSG00000186222/ENSG00000160683/ENSG00000106211/ENSG00000185338/ENSG00000171223/ENSG00000130522/ENSG00000172216/ENSG00000105974/ENSG00000101665/ENSG00000115008/ENSG00000177606
#> GO:0050865                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 ENSG00000095059/ENSG00000116670/ENSG00000126353/ENSG00000183655/ENSG00000154127/ENSG00000204287/ENSG00000171791/ENSG00000237541/ENSG00000172977/ENSG00000091972/ENSG00000179388/ENSG00000232629/ENSG00000198053/ENSG00000164430/ENSG00000127666/ENSG00000112149/ENSG00000104951/ENSG00000183624/ENSG00000126368/ENSG00000131981/ENSG00000119508/ENSG00000090776/ENSG00000125726/ENSG00000108861/ENSG00000117090/ENSG00000197329/ENSG00000204592/ENSG00000107338/ENSG00000213689/ENSG00000010278/ENSG00000196839/ENSG00000159958/ENSG00000163874/ENSG00000198435/ENSG00000167775/ENSG00000160223/ENSG00000113302/ENSG00000176170/ENSG00000232810/ENSG00000223496/ENSG00000167604/ENSG00000130164/ENSG00000101017/ENSG00000110944/ENSG00000105246/ENSG00000002330/ENSG00000143507/ENSG00000156127/ENSG00000125657/ENSG00000185338/ENSG00000171223/ENSG00000130522/ENSG00000172216/ENSG00000105974/ENSG00000101665/ENSG00000115008
#> GO:0050863                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 ENSG00000121594/ENSG00000125735/ENSG00000095059/ENSG00000126353/ENSG00000183655/ENSG00000204287/ENSG00000237541/ENSG00000172977/ENSG00000179388/ENSG00000232629/ENSG00000198053/ENSG00000164430/ENSG00000112149/ENSG00000104951/ENSG00000131981/ENSG00000090776/ENSG00000125726/ENSG00000108861/ENSG00000117090/ENSG00000197329/ENSG00000204592/ENSG00000107338/ENSG00000213689/ENSG00000196839/ENSG00000159958/ENSG00000163874/ENSG00000198435/ENSG00000160223/ENSG00000113302/ENSG00000167604/ENSG00000110944/ENSG00000105246/ENSG00000002330/ENSG00000143507/ENSG00000156127/ENSG00000125657/ENSG00000185338/ENSG00000171223/ENSG00000172216/ENSG00000105974/ENSG00000101665/ENSG00000115008
#> GO:0051249                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 ENSG00000196126/ENSG00000189114/ENSG00000121594/ENSG00000125735/ENSG00000095059/ENSG00000116670/ENSG00000126353/ENSG00000183655/ENSG00000204287/ENSG00000171791/ENSG00000237541/ENSG00000172977/ENSG00000179388/ENSG00000232629/ENSG00000198053/ENSG00000164430/ENSG00000127666/ENSG00000112149/ENSG00000104951/ENSG00000183624/ENSG00000131981/ENSG00000090776/ENSG00000125726/ENSG00000108861/ENSG00000117090/ENSG00000197329/ENSG00000204592/ENSG00000107338/ENSG00000213689/ENSG00000196839/ENSG00000159958/ENSG00000163874/ENSG00000198435/ENSG00000167775/ENSG00000160223/ENSG00000113302/ENSG00000223496/ENSG00000167604/ENSG00000101017/ENSG00000110944/ENSG00000105246/ENSG00000002330/ENSG00000143507/ENSG00000156127/ENSG00000125657/ENSG00000185338/ENSG00000171223/ENSG00000172216/ENSG00000105974/ENSG00000101665/ENSG00000115008
#> GO:0002696                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 ENSG00000113811/ENSG00000174130/ENSG00000186827/ENSG00000130475/ENSG00000103653/ENSG00000185950/ENSG00000019582/ENSG00000154096/ENSG00000066336/ENSG00000107968/ENSG00000149273/ENSG00000135426/ENSG00000033327/ENSG00000198502/ENSG00000105639/ENSG00000135334/ENSG00000168040/ENSG00000196126/ENSG00000189114/ENSG00000121594/ENSG00000125735/ENSG00000095059/ENSG00000116670/ENSG00000126353/ENSG00000183655/ENSG00000204287/ENSG00000171791/ENSG00000237541/ENSG00000172977/ENSG00000179388/ENSG00000232629/ENSG00000198053/ENSG00000127666/ENSG00000112149/ENSG00000104951/ENSG00000183624/ENSG00000119508/ENSG00000090776/ENSG00000125726/ENSG00000117090/ENSG00000197329/ENSG00000204592/ENSG00000107338/ENSG00000196839/ENSG00000159958/ENSG00000167775/ENSG00000160223/ENSG00000113302/ENSG00000232810/ENSG00000223496/ENSG00000167604/ENSG00000101017/ENSG00000110944/ENSG00000105246/ENSG00000002330/ENSG00000143507/ENSG00000125657/ENSG00000185338/ENSG00000130522/ENSG00000105974/ENSG00000115008
#> GO:0050867                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 ENSG00000113811/ENSG00000174130/ENSG00000186827/ENSG00000130475/ENSG00000103653/ENSG00000185950/ENSG00000019582/ENSG00000154096/ENSG00000066336/ENSG00000107968/ENSG00000149273/ENSG00000135426/ENSG00000033327/ENSG00000198502/ENSG00000105639/ENSG00000135334/ENSG00000168040/ENSG00000196126/ENSG00000189114/ENSG00000121594/ENSG00000125735/ENSG00000095059/ENSG00000116670/ENSG00000126353/ENSG00000183655/ENSG00000204287/ENSG00000171791/ENSG00000237541/ENSG00000172977/ENSG00000179388/ENSG00000232629/ENSG00000198053/ENSG00000127666/ENSG00000112149/ENSG00000104951/ENSG00000183624/ENSG00000119508/ENSG00000090776/ENSG00000125726/ENSG00000117090/ENSG00000197329/ENSG00000204592/ENSG00000107338/ENSG00000196839/ENSG00000159958/ENSG00000167775/ENSG00000160223/ENSG00000113302/ENSG00000232810/ENSG00000223496/ENSG00000167604/ENSG00000101017/ENSG00000110944/ENSG00000105246/ENSG00000002330/ENSG00000143507/ENSG00000125657/ENSG00000185338/ENSG00000130522/ENSG00000105974/ENSG00000115008
#> GO:0046649                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 ENSG00000196126/ENSG00000189114/ENSG00000114383/ENSG00000121594/ENSG00000125735/ENSG00000177697/ENSG00000120738/ENSG00000095059/ENSG00000116670/ENSG00000126353/ENSG00000183655/ENSG00000204287/ENSG00000171791/ENSG00000237541/ENSG00000172977/ENSG00000215021/ENSG00000179388/ENSG00000232629/ENSG00000155760/ENSG00000198053/ENSG00000164430/ENSG00000127666/ENSG00000112149/ENSG00000104951/ENSG00000183624/ENSG00000168264/ENSG00000012061/ENSG00000205220/ENSG00000131981/ENSG00000090776/ENSG00000090339/ENSG00000125726/ENSG00000153879/ENSG00000108861/ENSG00000117090/ENSG00000197329/ENSG00000105438/ENSG00000204592/ENSG00000107338/ENSG00000213689/ENSG00000196839/ENSG00000071564/ENSG00000159958/ENSG00000131188/ENSG00000163874/ENSG00000198435/ENSG00000105374/ENSG00000167775/ENSG00000160223/ENSG00000113302/ENSG00000103522/ENSG00000223496/ENSG00000167604/ENSG00000140968/ENSG00000101017/ENSG00000110944/ENSG00000105246/ENSG00000002330/ENSG00000143507/ENSG00000156127/ENSG00000125657/ENSG00000160683/ENSG00000185338/ENSG00000171223/ENSG00000172216/ENSG00000105974/ENSG00000101665/ENSG00000115008
#> GO:0051240 ENSG00000128604/ENSG00000111012/ENSG00000064393/ENSG00000089094/ENSG00000092820/ENSG00000196126/ENSG00000189114/ENSG00000114383/ENSG00000136997/ENSG00000121594/ENSG00000125735/ENSG00000126262/ENSG00000095739/ENSG00000186350/ENSG00000121966/ENSG00000120738/ENSG00000073150/ENSG00000095059/ENSG00000131979/ENSG00000116670/ENSG00000126353/ENSG00000148834/ENSG00000183655/ENSG00000102103/ENSG00000103811/ENSG00000136717/ENSG00000204287/ENSG00000148926/ENSG00000197405/ENSG00000171791/ENSG00000237541/ENSG00000140678/ENSG00000133321/ENSG00000172977/ENSG00000213015/ENSG00000111725/ENSG00000171236/ENSG00000091972/ENSG00000185507/ENSG00000179388/ENSG00000187109/ENSG00000139146/ENSG00000232629/ENSG00000213626/ENSG00000124920/ENSG00000198053/ENSG00000110092/ENSG00000164430/ENSG00000127666/ENSG00000112149/ENSG00000104951/ENSG00000183624/ENSG00000176788/ENSG00000060138/ENSG00000161638/ENSG00000111424/ENSG00000119508/ENSG00000140995/ENSG00000090776/ENSG00000158470/ENSG00000125726/ENSG00000153879/ENSG00000164442/ENSG00000117090/ENSG00000105372/ENSG00000197329/ENSG00000108561/ENSG00000130669/ENSG00000153487/ENSG00000204592/ENSG00000107338/ENSG00000173511/ENSG00000196839/ENSG00000184371/ENSG00000159958/ENSG00000111087/ENSG00000170345/ENSG00000163794/ENSG00000163874/ENSG00000073712/ENSG00000103257/ENSG00000167775/ENSG00000160223/ENSG00000113302/ENSG00000168056/ENSG00000176170/ENSG00000232810/ENSG00000168779/ENSG00000186652/ENSG00000223496/ENSG00000101670/ENSG00000215301/ENSG00000131408/ENSG00000167604/ENSG00000124216/ENSG00000169992/ENSG00000140564/ENSG00000140968/ENSG00000105245/ENSG00000101017/ENSG00000110944/ENSG00000105246/ENSG00000002330/ENSG00000125398/ENSG00000143507/ENSG00000156127/ENSG00000196843/ENSG00000125657/ENSG00000106211/ENSG00000130222/ENSG00000185338/ENSG00000143878/ENSG00000130522/ENSG00000172216/ENSG00000105974/ENSG00000101665/ENSG00000115008/ENSG00000123689
```


Columns in the results table:

- **ID**: ID of gene set

- **Description**: Description of Gene Set

- **setSize**: Size of the gene set

- **enrichmentScore**: Enrichment Score
  - note: this enrichment score has not been normalized for gene set size
  - means: larger gene sets automatically have a bigger (absolute) enrichment score, independent of actual differential enrichment
  - the raw enrichment score is therefore not comparable between different gene sets

- **NES**: Normalized version of column enrichmentScore
  - NES can be compared between different gene sets

- **pvalue**: p-value of enrichment of given gene set
  - note: this raw p-value has not been adjusted for multiple testing
  - therefore: CANNOT be used to assess differential enrichment of a given gene set

- **p.adjust**: **adjusted** p-value of a given gene set. 
  - Note that this p-value now has been adjusted for multiple testing and can therefore be used to assess differential enrichment. For instance, detect all genes with p.adjust < 0.05 as differentially enriched.

- **qvalue**: **adjusted** p-value of a given gene set. Note that a qvalue is the analog to p.adjust, but adjusted for multiple testing using a different method.

Note that you can either use the column p.adjust or qvalue to assess whether a gene set is differentially enriched or not.

- **rank**: position in the ranked list of the genes in which the maximum difference between the two sums occurs, i.e. the rank at which the enrichment score is extracted.

- **leading_edge**:
  + **tags**: The percentage of genes  in the ranked list that are members of gene set before (for positive enrichment score) or after (for negative enrichment score) the position from which the enrichment score is extracted.
  + **list**:  The percentage of genes in the ranked gene list before (for positive enrichment
           score) or after (for negative enrichment score) the position from which the enrichment score is extracted.
  + **signal**: enrichment signal strength
           combines the statistics "tags" and "list"

  + **core_enrichment**: genes that contribute most to the enrichment results

Note that a more detailed information on leading_edge and core_enrichment can be found in the user manual provided for GSEA:
 (https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideFrame.html): 
 
- search for **Detailed Enrichment Results** (section "Leading Edge")
 
- search for **CORE ENRICHMENT**




## Researchers' degrees of freedom 

In this part, we will go through all parameters that can be adapted in the GOSeq workflow. It is important to note that the intention behind going through the researchers' degrees of freedom is to give you an understanding of what you can do to adapt the given (parameter) setting to the research question. It is even more important to keep in mind that the intention behind going through these flexible parameters is NOT to change them in order to help you obtain the most preferable results by systematically changing these parameters since such behaviour would correspond to "cherry-picking". Any changes in the parameter choice must be documented transparently. 


### Change 1: change exponent 

Note that, in contrast to the web-based tool GSEA, clusterProfiler does not make any indication as to which alternative exponent values are allowed or sensible. One can stick with the web-based tool GSEA which suggests, in addition to the default exponent value of 1, the alternative values 0, 1.5, and 2.

Here, we change the exponent value to 1.5.


```r
# important: set seed for reproducibility
set.seed(1)

GSEA_GO_exponent <- gseGO(geneList = rankvec_Ensembl,
                            ont = "BP",
                            OrgDb = org.Hs.eg.db,
                            keyType = "ENSEMBL",
                            exp = 1.5,
                            seed = TRUE)
```

**Inspect the results**

```r
head(GSEA_GO_exponent, n = 10) 
#>                    ID                         Description
#> GO:0042572 GO:0042572           retinol metabolic process
#> GO:0045321 GO:0045321                leukocyte activation
#> GO:0001775 GO:0001775                     cell activation
#> GO:0006721 GO:0006721         terpenoid metabolic process
#> GO:0001523 GO:0001523          retinoid metabolic process
#> GO:0002694 GO:0002694  regulation of leukocyte activation
#> GO:0050865 GO:0050865       regulation of cell activation
#> GO:0055088 GO:0055088                   lipid homeostasis
#> GO:0051249 GO:0051249 regulation of lymphocyte activation
#> GO:0006629 GO:0006629             lipid metabolic process
#>            setSize enrichmentScore       NES       pvalue
#> GO:0042572      16      -0.9375514 -1.999582 1.516074e-06
#> GO:0045321     335      -0.5323438 -1.717534 8.735337e-06
#> GO:0001775     379      -0.5155165 -1.681090 8.893550e-06
#> GO:0006721      24      -0.8807910 -2.037460 1.655703e-05
#> GO:0001523      20      -0.9263665 -2.057425 3.114786e-05
#> GO:0002694     218      -0.5788244 -1.797112 2.686939e-05
#> GO:0050865     231      -0.5719874 -1.781633 3.022575e-05
#> GO:0055088      52      -0.7733441 -2.026849 5.967564e-05
#> GO:0051249     185      -0.5863351 -1.798671 7.530031e-05
#> GO:0006629     441      -0.4774551 -1.573053 6.777393e-05
#>               p.adjust      qvalue rank
#> GO:0042572 0.005627665 0.005122733   22
#> GO:0045321 0.011004286 0.010016946  782
#> GO:0001775 0.011004286 0.010016946  782
#> GO:0006721 0.015364923 0.013986333   22
#> GO:0001523 0.016517266 0.015035284   22
#> GO:0002694 0.016517266 0.015035284  778
#> GO:0050865 0.016517266 0.015035284  778
#> GO:0055088 0.027689495 0.025205104  285
#> GO:0051249 0.027951474 0.025443577  778
#> GO:0006629 0.027951474 0.025443577  836
#>                              leading_edge
#> GO:0042572  tags=12%, list=0%, signal=12%
#> GO:0045321 tags=22%, list=13%, signal=20%
#> GO:0001775 tags=21%, list=13%, signal=19%
#> GO:0006721  tags=12%, list=0%, signal=13%
#> GO:0001523  tags=10%, list=0%, signal=10%
#> GO:0002694 tags=25%, list=12%, signal=22%
#> GO:0050865 tags=24%, list=12%, signal=22%
#> GO:0055088  tags=29%, list=5%, signal=28%
#> GO:0051249 tags=25%, list=12%, signal=23%
#> GO:0006629 tags=19%, list=13%, signal=18%
#>                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            core_enrichment
#> GO:0042572                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 ENSG00000062282/ENSG00000006757
#> GO:0045321                                                                                                                                                 ENSG00000120738/ENSG00000095059/ENSG00000116670/ENSG00000126353/ENSG00000183655/ENSG00000204287/ENSG00000197405/ENSG00000171791/ENSG00000237541/ENSG00000172977/ENSG00000091972/ENSG00000215021/ENSG00000179388/ENSG00000232629/ENSG00000155760/ENSG00000198053/ENSG00000164430/ENSG00000127666/ENSG00000112149/ENSG00000104951/ENSG00000183624/ENSG00000168264/ENSG00000012061/ENSG00000205220/ENSG00000126368/ENSG00000131981/ENSG00000119508/ENSG00000090776/ENSG00000090339/ENSG00000125726/ENSG00000153879/ENSG00000108861/ENSG00000117090/ENSG00000197329/ENSG00000105438/ENSG00000204592/ENSG00000107338/ENSG00000213689/ENSG00000168447/ENSG00000196839/ENSG00000071564/ENSG00000184371/ENSG00000159958/ENSG00000131188/ENSG00000163874/ENSG00000198435/ENSG00000105374/ENSG00000167775/ENSG00000160223/ENSG00000113302/ENSG00000176170/ENSG00000232810/ENSG00000103522/ENSG00000223496/ENSG00000167604/ENSG00000130164/ENSG00000078902/ENSG00000140968/ENSG00000101017/ENSG00000110944/ENSG00000105246/ENSG00000002330/ENSG00000143507/ENSG00000156127/ENSG00000125657/ENSG00000160683/ENSG00000185338/ENSG00000171223/ENSG00000130522/ENSG00000172216/ENSG00000105974/ENSG00000101665/ENSG00000115008/ENSG00000177606
#> GO:0001775                                                                 ENSG00000120738/ENSG00000095059/ENSG00000116670/ENSG00000126353/ENSG00000183655/ENSG00000154127/ENSG00000204287/ENSG00000197405/ENSG00000171791/ENSG00000237541/ENSG00000172977/ENSG00000091972/ENSG00000215021/ENSG00000179388/ENSG00000232629/ENSG00000155760/ENSG00000198053/ENSG00000164430/ENSG00000127666/ENSG00000112149/ENSG00000135124/ENSG00000104951/ENSG00000183624/ENSG00000168264/ENSG00000012061/ENSG00000205220/ENSG00000126368/ENSG00000131981/ENSG00000119508/ENSG00000090776/ENSG00000090339/ENSG00000125726/ENSG00000153879/ENSG00000108861/ENSG00000117090/ENSG00000197329/ENSG00000105438/ENSG00000204592/ENSG00000107338/ENSG00000213689/ENSG00000010278/ENSG00000168447/ENSG00000196839/ENSG00000071564/ENSG00000184371/ENSG00000159958/ENSG00000131188/ENSG00000163874/ENSG00000198435/ENSG00000105374/ENSG00000167775/ENSG00000160223/ENSG00000113302/ENSG00000176170/ENSG00000232810/ENSG00000103522/ENSG00000223496/ENSG00000167604/ENSG00000130164/ENSG00000078902/ENSG00000140968/ENSG00000101017/ENSG00000110944/ENSG00000105246/ENSG00000002330/ENSG00000143507/ENSG00000156127/ENSG00000125657/ENSG00000186222/ENSG00000160683/ENSG00000106211/ENSG00000185338/ENSG00000171223/ENSG00000130522/ENSG00000172216/ENSG00000105974/ENSG00000101665/ENSG00000115008/ENSG00000177606
#> GO:0006721                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 ENSG00000145545/ENSG00000062282/ENSG00000006757
#> GO:0001523                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 ENSG00000062282/ENSG00000006757
#> GO:0002694                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 ENSG00000095059/ENSG00000116670/ENSG00000126353/ENSG00000183655/ENSG00000204287/ENSG00000171791/ENSG00000237541/ENSG00000172977/ENSG00000091972/ENSG00000179388/ENSG00000232629/ENSG00000198053/ENSG00000164430/ENSG00000127666/ENSG00000112149/ENSG00000104951/ENSG00000183624/ENSG00000126368/ENSG00000131981/ENSG00000119508/ENSG00000090776/ENSG00000125726/ENSG00000108861/ENSG00000117090/ENSG00000197329/ENSG00000204592/ENSG00000107338/ENSG00000213689/ENSG00000196839/ENSG00000159958/ENSG00000163874/ENSG00000198435/ENSG00000167775/ENSG00000160223/ENSG00000113302/ENSG00000176170/ENSG00000232810/ENSG00000223496/ENSG00000167604/ENSG00000130164/ENSG00000101017/ENSG00000110944/ENSG00000105246/ENSG00000002330/ENSG00000143507/ENSG00000156127/ENSG00000125657/ENSG00000185338/ENSG00000171223/ENSG00000130522/ENSG00000172216/ENSG00000105974/ENSG00000101665/ENSG00000115008
#> GO:0050865                                                                                                                                                                                                                                                                                                                                                                                                                                                 ENSG00000095059/ENSG00000116670/ENSG00000126353/ENSG00000183655/ENSG00000154127/ENSG00000204287/ENSG00000171791/ENSG00000237541/ENSG00000172977/ENSG00000091972/ENSG00000179388/ENSG00000232629/ENSG00000198053/ENSG00000164430/ENSG00000127666/ENSG00000112149/ENSG00000104951/ENSG00000183624/ENSG00000126368/ENSG00000131981/ENSG00000119508/ENSG00000090776/ENSG00000125726/ENSG00000108861/ENSG00000117090/ENSG00000197329/ENSG00000204592/ENSG00000107338/ENSG00000213689/ENSG00000010278/ENSG00000196839/ENSG00000159958/ENSG00000163874/ENSG00000198435/ENSG00000167775/ENSG00000160223/ENSG00000113302/ENSG00000176170/ENSG00000232810/ENSG00000223496/ENSG00000167604/ENSG00000130164/ENSG00000101017/ENSG00000110944/ENSG00000105246/ENSG00000002330/ENSG00000143507/ENSG00000156127/ENSG00000125657/ENSG00000185338/ENSG00000171223/ENSG00000130522/ENSG00000172216/ENSG00000105974/ENSG00000101665/ENSG00000115008
#> GO:0055088                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 ENSG00000126368/ENSG00000109084/ENSG00000073060/ENSG00000167695/ENSG00000119655/ENSG00000087237/ENSG00000105698/ENSG00000204444/ENSG00000101670/ENSG00000215301/ENSG00000131408/ENSG00000130164/ENSG00000062282/ENSG00000105974/ENSG00000006757
#> GO:0051249                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 ENSG00000095059/ENSG00000116670/ENSG00000126353/ENSG00000183655/ENSG00000204287/ENSG00000171791/ENSG00000237541/ENSG00000172977/ENSG00000179388/ENSG00000232629/ENSG00000198053/ENSG00000164430/ENSG00000127666/ENSG00000112149/ENSG00000104951/ENSG00000183624/ENSG00000131981/ENSG00000090776/ENSG00000125726/ENSG00000108861/ENSG00000117090/ENSG00000197329/ENSG00000204592/ENSG00000107338/ENSG00000213689/ENSG00000196839/ENSG00000159958/ENSG00000163874/ENSG00000198435/ENSG00000167775/ENSG00000160223/ENSG00000113302/ENSG00000223496/ENSG00000167604/ENSG00000101017/ENSG00000110944/ENSG00000105246/ENSG00000002330/ENSG00000143507/ENSG00000156127/ENSG00000125657/ENSG00000185338/ENSG00000171223/ENSG00000172216/ENSG00000105974/ENSG00000101665/ENSG00000115008
#> GO:0006629 ENSG00000101255/ENSG00000100889/ENSG00000141858/ENSG00000184207/ENSG00000099377/ENSG00000221968/ENSG00000079459/ENSG00000166908/ENSG00000120738/ENSG00000074660/ENSG00000126353/ENSG00000183260/ENSG00000183655/ENSG00000149084/ENSG00000130766/ENSG00000148926/ENSG00000133321/ENSG00000172977/ENSG00000163964/ENSG00000242372/ENSG00000110536/ENSG00000166394/ENSG00000111725/ENSG00000165175/ENSG00000174915/ENSG00000215021/ENSG00000103502/ENSG00000181513/ENSG00000169710/ENSG00000112531/ENSG00000121769/ENSG00000129474/ENSG00000100979/ENSG00000204386/ENSG00000133328/ENSG00000103642/ENSG00000126368/ENSG00000119508/ENSG00000125741/ENSG00000158470/ENSG00000145545/ENSG00000123179/ENSG00000213689/ENSG00000168447/ENSG00000073060/ENSG00000158715/ENSG00000119655/ENSG00000167468/ENSG00000165782/ENSG00000087237/ENSG00000083807/ENSG00000100294/ENSG00000176170/ENSG00000232810/ENSG00000196743/ENSG00000179598/ENSG00000167969/ENSG00000101670/ENSG00000167130/ENSG00000131408/ENSG00000173599/ENSG00000124216/ENSG00000120833/ENSG00000130164/ENSG00000086544/ENSG00000167508/ENSG00000089163/ENSG00000125398/ENSG00000087076/ENSG00000111666/ENSG00000110921/ENSG00000174326/ENSG00000100600/ENSG00000147383/ENSG00000185338/ENSG00000125652/ENSG00000062282/ENSG00000105974/ENSG00000172893/ENSG00000101846/ENSG00000196407/ENSG00000115008/ENSG00000006757
```

When switching the exponent to 2, the number of differentially enriched gene sets has decreased to two. 



