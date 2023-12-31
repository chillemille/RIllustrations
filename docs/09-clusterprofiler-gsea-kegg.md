
# (PART) With Gene ID Conversion {-}

# clusterProfiler's GSEA (with gene set database KEGG)



In this script, we will do the following two things:

- 1. Based on the results of differential expression analysis from voom/limma, DESeq2, and edgeR, we will go through all steps required to run clusterProfiler's GSEA with gene set database GO.

- 2. We will go through all (meaningful) researchers' degrees of freedom.

Note that we provide a separate file for the gene set databases GO and KEGG since clusterProfiler differs in the accepted gene ID formats for them. While for GO, a variety of formats, including Entrez ID, is accepted, the gene IDs must be converted to NCBI (Entrez) ID when choosing KEGG.  

## Libraries 

All necessary packages are available on Bioconductor, and should be installed from there if not already available on your machine.


```r
install.packages("BiocManager")
BiocManager::install("clusterProfiler")
BiocManager::install("DESeq2")
```

**Load Libraries** 

```r
library(clusterProfiler)
library(DESeq2)
```

Description of the packages 

- **clusterProfiler**: We here use clusterProfiler's implementation of GSEA. 

- **DESeq2**: In contrast detecting the differential expressed genes, which is done solely based on adjusted $p$-vaules, we additionally use the estimated log fold changes when generating the gene ranking. In DESeq2, the estimated log fold change values are shrunken using an additional function, hence loading the library.  


## Load data 
Note that, as described above, clusterProfiler accepts Entrez IDs when working with gene set database GO. When generating the ranking of the genes with DESeq2, a shrinkage of the log fold change estimates is performed which requires the output of function *DESeq()*, here named **dds** (see Chapter 'Differential Expression Analysis).


```r
load("./data/Results_Differential_Expression_Analysis/DE_results_limma_Entrez.Rdata")
load("./data/Results_Differential_Expression_Analysis/DE_results_DESeq2_Entrez.Rdata")
load("./data/Results_Differential_Expression_Analysis/DE_results_edgeR_Entrez.Rdata")

# load dds object 
load("./data/Results_Differential_Expression_Analysis/dds_Entrez.Rdata")
```

## Run GSEA 


### Step 1: Generation of Required Input Object

In the following, we illustrate the creation of the required input object based on the results table of the three different methods for differential expression analysis, voom/limma, DESeq2, and edgeR, which shall serve as three options. The required input object for clusterProfiler is an "order ranked geneList", i.e. a vector with a ranking value for each gene, named with the respective gene IDs and ordered in a descending manner. \
    
The formula for generating the gene ranking from the results table of differential expression analysis is 
$$(-1) * \log_{10}(p\text{-value}) * \text{sign}(\text{log fold change})$$. \
Note that "$p$-value" denotes the non-adjusted p-value.     

#### Option 1: Generate required input using voom/limma



```r
# 1. Subset the gene expression data set to those genes that have a p-value (i.e.
# which have been NOT been excluded from differential expression analysis)

# indicate those genes WITH a p-value
ind_nonNA_pvalue_limma_Entrez <- !is.na(DE_results_limma_Entrez$P.Value)

# subset gene expression data set to those genes with a p-value
DE_results_noNA_Entrez <- DE_results_limma_Entrez[ind_nonNA_pvalue_limma_Entrez, ]

# 2. create vector that contains the value of the gene-level ranking metric for each gene
rankvec_limma_Entrez <- sign(DE_results_noNA_Entrez$logFC)*(-1)*log10(DE_results_noNA_Entrez$P.Value)

# 3. assign respective gene ID to each value in the vector
names(rankvec_limma_Entrez) <- rownames(DE_results_noNA_Entrez)

# 4. sort the vector in descending order
rankvec_limma_Entrez <- sort(rankvec_limma_Entrez, decreasing=TRUE)

# inspect the first entries of the ranking vector
head(rankvec_limma_Entrez , n = 10)
#>      6192      9087       317      9404      9697     63894 
#> 54.991351 38.217609  3.554179  3.200661  3.164203  3.071138 
#>      4214     10479      5257     29915 
#>  3.019720  2.924894  2.875024  2.814519
```

#### Option 2: Generate required input using DESeq2

**(i) Perform Differential Expression Analysis**

This step is performed in the script Instructions_Differential_Expression_Analysis.R. Note that we must exit the DESeq2 workflow from the script a few steps early since there is an additional step required in for the creation of the required input object. We proceed with the object "dds_Entrez" we've obtained in step 3 of DESeq2's workflow for differential expression analysis.


```r
dds_Entrez
#> class: DESeqDataSet 
#> dim: 9492 69 
#> metadata(1): version
#> assays(6): counts mu ... replaceCounts replaceCooks
#> rownames(9492): 8813 57147 ... 56104 56112
#> rowData names(23): baseMean baseVar ... maxCooks
#>   replace
#> colnames(69): NA18486 NA18498 ... NA19239 NA19257
#> colData names(3): condition sizeFactor replaceable
```



**(ii) Shrink estimated log fold change values**

The intuition behind shrinkage of log fold change values is that, as mentioned in the paper, RNA-Seq data consists (in its raw form) of count data in which is inherently heteroscedastic, i.e. the variance of the count data depends on the mean count of the count data. It is observable that ratios between counts are considerably noisier in low magnitudes of counts compared to higher magnitudes, i.e. the log fold changes between both conditions are higher if the overall magnitude of counts is lower. \
DESeq2 addresses this issue by shrinking the estimated log fold changes in the direction of 0 the magnitude of shrinkage is higher if the available information for a gene is lower (which may be because of a low magnitude of counts, a high dispersion or few degrees of freedom.). A more detailed description is provided in the DESeq2 paper by Love et al. (2014) [@love2014moderated].

Perform shrinkage:

```r
DE_results_DESeq2_shrink_Entrez <- lfcShrink(dds_Entrez,
                                              coef = "condition_treated_vs_untreated",
                                              type="apeglm")
```

Function arguments:

- **type**: method to perform shrinkage. We opt for the default "apeglm" but you can choose from two alternative options as well

- **coef**: indicate the coefficients to be shrunk. We can obtain the right argument from the following function call:


```r
resultsNames(dds_Entrez)
#> [1] "Intercept"                     
#> [2] "condition_treated_vs_untreated"
```

This command shows us that we can either shrink the intercept or the "condition_treated_vs_untreated". Since we do not want to shrink the intercept but the estimated log fold changes, we opt for the second option "condition_treated_vs_untreated".


**(iii) Generate the ranked list of the genes**



```r
# 1. Subset the DESeq2 results to those genes that have a p-value (i.e.
# which have been NOT been excluded from differential expression analysis)

# indicate those genes with an existing p-value
ind_nonNA_pvalue_Entrez <- !is.na(DE_results_DESeq2_shrink_Entrez$pvalue)

# subset the DESeq2 results to those genes with an existing p-value
DE_results_noNA_Entrez <- DE_results_DESeq2_shrink_Entrez[ind_nonNA_pvalue_Entrez, ]

# 2. create vector that contains the value of the gene-level ranking metric for each gene
rankvec_DESeq2_Entrez <- sign(DE_results_noNA_Entrez$log2FoldChange)*(-1)*log10(DE_results_noNA_Entrez$pvalue)

# 3. assign respective gene ID to each value in the vector
names(rankvec_DESeq2_Entrez) <- rownames(DE_results_noNA_Entrez)

# 4. sort the vector in descending order
rankvec_DESeq2_Entrez <- sort(rankvec_DESeq2_Entrez, decreasing=TRUE)


# inspect the first entries of the ranking vector
head(rankvec_DESeq2_Entrez , n = 10)
#>       6192       9086       9087       7404         58 
#> 204.374390  75.096224  64.815184   7.205688   4.854231 
#>       1521      23136     219699     126129        317 
#>   4.467907   3.996186   3.966326   3.918029   3.744345
```


#### Option 3: Generate required input using edgeR


```r
# 1. Subset the gene expression data set to those genes that have a p-value (i.e.
# which have been NOT been excluded from differential expression analysis)

# indicate those genes WITH a p-value
ind_nonNA_pvalue_Entrez <- !is.na(DE_results_edgeR_Entrez$PValue)

# subset gene expression data set to those genes with a p-value
DE_results_noNA_Entrez <- DE_results_edgeR_Entrez[ind_nonNA_pvalue_Entrez, ]

# 2. create vector that contains the value of the gene-level ranking metric for each gene
rankvec_edgeR_Entrez <- sign(DE_results_noNA_Entrez$logFC)*(-1)*log10(DE_results_noNA_Entrez$p_adj)

# 3. assign respective gene ID to each value in the vector
names(rankvec_edgeR_Entrez) <- rownames(DE_results_noNA_Entrez)

# 4. sort the vector in descending order
rankvec_edgeR_Entrez <- sort(rankvec_edgeR_Entrez, decreasing=TRUE)

# inspect ranking vector
head(rankvec_edgeR_Entrez, n = 10)
#>        6192        9087        1521         317        9404 
#>         Inf 116.7303466   1.6605772   1.0317662   0.9984938 
#>       10479        4214        9697        4509        9645 
#>   0.9863396   0.8622382   0.8467278   0.8434267   0.8251273

# 5. special problem here: gene 6192 has a ranking value Inf since its adjusted
# p-value in the results table of differential expression analysis amounts to 0

# here, we deal with this issue by resetting this ranking value to the highest ranking value
# that occurs among the remaining genes
# -> note that there is NO common way of dealing with this issue
rankvec_edgeR_Entrez[rankvec_edgeR_Entrez == Inf] <- max(rankvec_edgeR_Entrez[rankvec_edgeR_Entrez != Inf])
 
# inspect ranking vector once more 
head(rankvec_edgeR_Entrez, n = 10)
#>        6192        9087        1521         317        9404 
#> 116.7303466 116.7303466   1.6605772   1.0317662   0.9984938 
#>       10479        4214        9697        4509        9645 
#>   0.9863396   0.8622382   0.8467278   0.8434267   0.8251273
```


### Step 2: Run GSEA with gene set database KEGG

For the purpose of simplicity, we work with the ranking generated using limma. However, you can switch around at your discretion.



```r
rankvec_Entrez <- rankvec_limma_Entrez

# alternatively: 
# rankvec_Entrez <- rankvec_DESeq2_Entrez
# rankvec_Entrez <- rankvec_edgeR_Entrez
```

Run clusterProfiler with gene set database KEGG 


```r
# important: set seed for reproducibility 
set.seed(1)
  
GSEA_KEGG <- gseKEGG(geneList = rankvec_Entrez, 
                     seed = TRUE)  
```

 additional arguments (these do not appear in the above code lines since we keep the default options)

- **'exponent = 1'**: in the calculation of the enrichment score, each gene is weighted by its absolute value of the ranking metric 

- **'organism = "hsa"'**: human organism. Note that the list of accepted organisms is available in the following link: \
https://www.genome.jp/kegg/catalog/org_list.html

- **'seed = TRUE'**: set the seed you have indicated above (here: seed = 1 )

- **'eps = 1e-10'**: each (adjusted) p-value that is smaller than 1e-10 is indicated as 1e-10

- **'pvalueCutoff'** = 0.05: only those gene sets with a p-value < 0.05 are indicated in the results 
 
- **note**: if you want to inspect ALL gene sets in the results, set pvalueCutoff = 1 

- **pAdjustMethod** = "BH": adjustment for multiple testing using Benjamini and Hochberg 
 
 
### step 3: Interpretation of the results 

We inspect the results for the first 10 gene sets 

```r
#head(GSEA_KEGG, n = 10 )
```


Columns in Result Tables: 

- **ID**: ID of the gene set 

- **Description**: Description of Gene Set 

- **setSize**: Size of the gene set

- **enrichmentScore**: Enrichment Score
  + note: this enrichment score has not been normalized for gene set size
  + means: larger gene sets automatically have a bigger (absolute) enrichment score independent of actual differential enrichment
  + the raw enrichment score is therefore not comparable between different gene sets

- **NES**: Normalized version of column enrichmentScore
  + NES can be compared between different gene sets 

- **pvalue**: p-value of enrichment of given gene set
  + note: this raw p-value has **not** been adjusted for multiple testing
  + therefore: **cannot** be used to assess differential enrichment of a given gene set 
 
- **p.adjust**: **adjusted** p-value of a given gene set 
  + this p-value now has been adjusted for multiple testing and can therefore 
  + can therefore be used to assess differential enrichment 
  + example: detect all genes with p.adjust < 0.05 as differentially enriched 
 
- **qvalue**: ADJUSTED p-value of a given gene set 
  + note: a qvalue is the analog to p.adjust, but adjusted for multiple testing 
         using a different method 

You can either use the column p.adjust or qvalue to assess whether a gene set 
     is differentially enriched or not 

- **rank**: position in the ranked list of the genes in which the maximum difference between the two sums occurs, i.e. the rank at which the enrichment score is extracted 

- **leading_edge**: 
  + **tags**: The percentage of genes  in the ranked list that are members of gene set before (for positive enrichment score) or after (for negative enrichment score) the position from which the enrichment score is extracted. 
  + **list**: The percentage of genes in the ranked gene list before (for positive enrichment 
           score) or after (for negative enrichment score) the position from which the enrichment score is extracted. 
  + **signal**: enrichment signal strength combines the statistics "tags" and "list"

- **core_enrichment**: genes that contribute most to the enrichment results 

- Note that a more detailed information on leading_edge and core_enrichment can be found in the user 
 manual provided for GSEA: 
 https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideFrame.html
  + search for **Detailed Enrichment Results** (section "Leading Edge")
  + search for **CORE ENRICHMENT** 




