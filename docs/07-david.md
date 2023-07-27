
# DAVID



In this script, we will

1. preprocess and export the required input object AND

2. generate the alternative universe as part of the researchers' degrees of freedom

for the web-based tool DAVID [@huang2009systematic; @huang2009bioinformatics].

The required input object for DAVID consists of a list of differentially expressed genes.\ 
As an alternative universe, we choose all genes whose differential expression was assessed in the experiment. Note that it is not required to provide an alternative universe. We illustrate its use as part of the researchers' degrees of freedom.




## Libraries 

For this part of the instructions, there are no further libraries required. 


## Load data 

Both objects, the list of differentially expressed genes as well as the alternative universe, can be extracted from the results table of differential expression analysis. Note that DAVID works on a variety of gene ID formats, including Ensembl ID. 

Here, we will illustrate the process for the limma results only since the process is identical for all methods for differential exression analysis (due to the unification of the column names that refer to the adjusted p-values). 


```r
load("./data/Results_Differential_Expression_Analysis/DE_results_limma_Ensembl.Rdata")
```

Alternatively, you can proceed with the results generated with DESeq2 or edgeR:

```r
# load("./Results_Differential_Expression_Analysis/DE_results_DESeq2_Ensembl.Rdata")
# load("./Results_Differential_Expression_Analysis/DE_results_edgeR_Ensembl.Rdata")
```

For the purpose of redability, we store limma's differential expression analysis results in an object with a more neutral name:


```r
DE_results <- DE_results_limma_Ensembl
```

Take a look at the results table:


```r
head(DE_results, n  = 10)
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

## Prepare and export required input oject 

**Step 1: Generate required input object**

From the results table of differential expression analysis, we generate a list of differentially expressed genes. 


```r
# (i) indicate which genes have an (existing!) adjusted p-value < 0.05.
ind_DE <- (!is.na(DE_results$p_adj)) & (DE_results$p_adj < 0.05)

# get overview of number of differentially expressed genes
table(ind_DE)
#> ind_DE
#> FALSE  TRUE 
#>  6241     5

# (ii) use the indicator to obtain the list of differentially expressed genes
DEG_vec <- rownames(DE_results)[ind_DE]
```

We take a look at the first genes in the vector of differentially expressed genes and also want to look at the number of differentially expressed genes.


```r
# inspect the first genes from the vector
head(DEG_vec, n = 10)
#> [1] "ENSG00000129824" "ENSG00000099749" "ENSG00000154620"
#> [4] "ENSG00000006757" "ENSG00000130021"

# get the number of differentially expressed genes
length(DEG_vec)
#> [1] 5
```

Note that in this specific example, there are very few differentially expressed genes. 


**Step 2: Export the vector of differentially expressed genes**

We store the vector of differentially expressed genes in the file "DEG_vec.txt" in the folder "Input_Objects_DAVID" (which is a subfolder of the folder "data"). 


```r
# the path indicates the location (folder(s) and file) of the vector to be stored in
path <- "./data/Input_Objects_DAVID/DEG_vec.txt"

# export 
 write.table(DEG_vec,
            file = path,
            quote=FALSE,
            row.names=FALSE,
            col.names = FALSE)
```

The file "DEG_vec.txt" should now have appeared in the folder "Input_Objects_DAVID". 


**Step 3: Upload to DAVID**

The resulting .txt file can be directly uploaded to the DAVID website: https://david.ncifcrf.gov/



## Researchers' degrees of freedom 


In this part, we will illustrate how to create an alternative universe to DAVID, which is a researchers' degree of freedom. It is important to note that the intention behind going through the researchers' degrees of freedom is to give you an understanding of what you can do to adapt the given (parameter) setting to the research question. It is even more important to keep in mind that the intention behind going through these flexible parameters is NOT to change them in order to help you obtain the most preferable results by systematically changing these parameters since such behaviour would correspond to "cherry-picking". Any changes in the parameter choice must be documented transparently.


### Change 1: Alternative universe

**Step 1: Generate alterative universe**

For the alternative universe, we want to consider all genes whose differential expression was measured in the experiment For some methods for differential expression analysis, such as DESeq2, the adjusted p-values of some genes are set to NA which means that these genes cannot be detected as neither differentially expressed nor not differentially expressed. We therefore want to remove these from the universe.



```r
# indicate which genes have an adjusted p-value in the results of differential expression analysis
ind_nona_p <- !is.na(DE_results$p_adj)

# filter the list of genes to those with an existing adjusted p-value
alternative_universe <- rownames(DE_results)[ind_nona_p]
```

Inspect the first few genes from the universe


```r
head(alternative_universe, n = 10)
#>  [1] "ENSG00000129824" "ENSG00000099749" "ENSG00000154620"
#>  [4] "ENSG00000006757" "ENSG00000130021" "ENSG00000185753"
#>  [7] "ENSG00000086712" "ENSG00000123689" "ENSG00000177606"
#> [10] "ENSG00000120868"
```

**Step 2: Export alternative universe**

We store the vector of differentially expressed genes in the file "alternative_universe.txt" in the folder "Input_Objects_DAVID" (which is a subfolder of the folder "data"). 



```r
path_alt_universe <- "./data/Input_Objects_DAVID/alternative_universe.txt"

# export 
write.table(alternative_universe,
            file = path_alt_universe,
            quote=FALSE,
            row.names=FALSE,
            col.names = FALSE)
```

The file "alternative_universe.txt" should now have appeared in the folder "Input_Objects_DAVID". 


