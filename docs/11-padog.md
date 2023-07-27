
# PADOG



In this script, we will go through the process of running PADOG. Note that PADOG requires the transformed RNA-Seq data and the genes in the Entrez ID format.

## Libraries

All necessary packages are available on Bioconductor, and should be installed from there if not already available on your machine.


```r
install.packages("BiocManager")
BiocManager::install("PADOG")
BiocManager::install("tweeDEseqCountData")
BiocManager::install("KEGGREST")
```

**Load libraries**


```r
library(PADOG)
library(tweeDEseqCountData)
library(KEGGREST)
```

Summary of the packages: 

- **PADOG**: Provides the implementation of the method PADOG 

- **tweeDEseqCountData**: In addition to the preprocessed Pickrell data set, which we load in the following step, PADOG requires the conditions of the samples 

- **KEGGREST**: To obtain the character value for a variety of organisms the user is required to provide in the function *padog()*. While the human organism is indicated by default, it might become necessary to consult the KEGGREST package if working with a different organism

## Load Data


```r
# we load the voom-transformed Pickrell data set 
load("data/Results_RNASeq_Transformation/expression_data_voomtransformed_Entrez.Rdata")

# alternatively: load the gene expression measurements that have been transformed using 
# load("data/expression_data_vsttransformed_Entrez.Rdata")

# additionally, we load the pickrell data set so that we can access the sample conditions
data(pickrell)
```

The sample conditions (i.e. phenotype labels) of the pickrell data set can be accessed using


```r
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
```

We proceed with the voom-transformed pickrell data set and the corresponding phenotype labels.


```r
# gene expression measurements (transformed)
# note: you can also proceed with the vst-transformed gene expression measurements 
expression_data_transformed <- expression_data_voomtransformed_Entrez

# sample conditions
sample_conditions <- pickrell.eset$gender
```
## PADOG

### Step 1: Prepare Sample Conditions

First, we inspect the form of the initial (raw) sample conditions.


```r
# look at the class: 
class(sample_conditions)
#> [1] "factor"
# -> the sample labels are already coded as factors

# the current levels are:
levels(sample_conditions)
#> [1] "female" "male"
```

PADOG requires a character vector with the class labels of the samples. It can only contain "c" for control samples or "d" for disease samples.


```r
# prepare sample conditions
# we want to convert 
# (i) "female" to "c"
# (ii) "male" to "d"
sample_conditions_prep <- factor(sample_conditions, 
                                levels=c("female","male"), 
                                labels=c("c","d"))
```

Inspect the prepared sample conditions: 

```r
sample_conditions_prep
#>  [1] d d c d c d c d c d c d c d c d c c d c d c c d c d c c
#> [29] d c c d c d c c c c d c d d c c d c c d c c d c d d c c
#> [57] d c d c d c c d c c c d c
#> Levels: c d
```


### Step 2: Run PADOG

**It is recommended to set a seed to ensure exact reproducibility of the results if the code is run at multiple time points**

You can specify any integer number as the seed. It is VERY IMPORTANT to choose the seed arbitrarily and WITHOUT INSPECTING the results. The seed should NEVER be specified based on which value yields the most preferable results.


```r
# run PADOG: 
PADOG_results <- padog(esetm = as.matrix(expression_data_transformed), 
                        group = sample_conditions_prep, 
                        dseed = 1)
```

Function arguments:

-   **esetm**: matrix that contains the expression measurements
    - note: since the expression data is initially stored in a data frame, we transform it to a matrix when running PADOG
    
- **group**: sample conditions (has values "c" and "d")

- **dseed** : seed for random number generation (used in the process of phenotype permutation)

additional arguments (we do not use here):

- **paired**: indicates whether the samples in both groups are paired

- **block**: if the samples are paired (i.e. argument paired = TRUE), then the paired samples must have the same block value

- **gslist**: gives instructions on how to cluster the genes into gene sets

  - **'gslist = "KEGGRESTpathway"'** (default): gene sets correspond to KEGG pathways
    
  - alternative: provide a user-defined list of gene sets
    
- **organism**: organism from which the gene expression measurements are taken

  - for human, set organism = "hsa"
    
  - the required character value for other organisms can be extracted from the            KEGGREST package:
    
  - obtain required organisms from column organism 
    
- **annotation**: required if gslist is set to "KEGGRESTpathway" and the rownames of esetm are probe IDs 
    
  - can be set to NULL of gslist is set to "KEGGRESTpathway" and the rownames of esetm are in the Entrez gene ID format 
    
  - if rownames are other gene IDs, then sett annotation = NULL and make sure that the rownames are elements of gslist (and unique!)
    
- **gs.names**: contains names of gene sets -> character vector
    
  - must have the same length as gslist 
    
- **NI**: number of phenotype permutations employed in the assessment of the significance of a given gene set 


We want to take a first look at the results table. Note, however, that that we still have to add the adjusted p-values (see next step), based on which differential enrichment is eventually assessed.


```r
head(PADOG_results , n = 10)
#>       Name    ID Size meanAbsT0   padog0 PmeanAbsT Ppadog
#> 04380 <NA> 04380   66  2.622104 3.832812     0.005  0.002
#> 04810 <NA> 04810  111  2.326027 5.476170     0.004  0.004
#> 03010 <NA> 03010   89  5.378671 9.659220     0.014  0.010
#> 05213 <NA> 05213   23  1.528150 1.461810     0.027  0.010
#> 05171 <NA> 05171  118  4.788746 8.020661     0.013  0.013
#> 05417 <NA> 05417  112  2.094534 3.096663     0.015  0.017
#> 05162 <NA> 05162   73  2.089532 3.042040     0.010  0.020
#> 00760 <NA> 00760   17  1.871113 3.583774     0.026  0.023
#> 04064 <NA> 04064   59  1.931306 3.552921     0.019  0.025
#> 04010 <NA> 04010  169  1.987648 4.770146     0.018  0.028
```


### Step 3: Adjust for Multiple Testing

PADOG does **not** perform multiple testing adjustment internally so that must be performed by the user. We here work with a function from the R package base which uses the method by Benjamini and Hochberg by default. 



```r
# add adjusted p-value in column Ppadog_adjusted
PADOG_results$Ppadog_adjusted <- p.adjust(PADOG_results$Ppadog)
```


### Step 4: Interpretation of the results 
We want to inspect (a part of) the results and interpret the columns provided in the results table. 


```r
head(PADOG_results, n = 10)
#>       Name    ID Size meanAbsT0   padog0 PmeanAbsT Ppadog
#> 04380 <NA> 04380   66  2.622104 3.832812     0.005  0.002
#> 04810 <NA> 04810  111  2.326027 5.476170     0.004  0.004
#> 03010 <NA> 03010   89  5.378671 9.659220     0.014  0.010
#> 05213 <NA> 05213   23  1.528150 1.461810     0.027  0.010
#> 05171 <NA> 05171  118  4.788746 8.020661     0.013  0.013
#> 05417 <NA> 05417  112  2.094534 3.096663     0.015  0.017
#> 05162 <NA> 05162   73  2.089532 3.042040     0.010  0.020
#> 00760 <NA> 00760   17  1.871113 3.583774     0.026  0.023
#> 04064 <NA> 04064   59  1.931306 3.552921     0.019  0.025
#> 04010 <NA> 04010  169  1.987648 4.770146     0.018  0.028
#>       Ppadog_adjusted
#> 04380           0.696
#> 04810           1.000
#> 03010           1.000
#> 05213           1.000
#> 05171           1.000
#> 05417           1.000
#> 05162           1.000
#> 00760           1.000
#> 04064           1.000
#> 04010           1.000
```

Differential enrichment of a given gene set can now be assessed based on the adjusted p-value in column **Ppadog_adjusted**. For instance: detect all gene sets with Ppadog_adjusted < 0.05 as differentially enriched.

Additional columns: 

- **Name**: Name of gene set 

- **ID:** Identifier of gene set 

- **Size**: number of genes in gene set 

- **meanAbsT0**: Mean of absolute (moderated) t-statistic of all genes that are a member of the given gene set 

- **padog0**: Mean of abolute weighted moderate t-statistic of all genes that are a member of the given gene set 

- **PmeanAbsT**: significance of of meanAbsT0

- **Ppadog**: significance of padog0
