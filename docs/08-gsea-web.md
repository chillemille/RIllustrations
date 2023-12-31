
# GSEA (web-based tool)



In this script, we will 

1. export the transformed RNA-Seq measurements 

2. prepare and export sample conditions 

Note that some further preprocessing steps are required in Excel. For these further preprocessing steps, click on the following link: \
http://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats. 


## Libraries

All necessary packages are available on Bioconductor, and should be installed from there if not already available on your machine.


```r
install.packages("BiocManager")
BiocManager::install("tweeDEseqCountData")
```

**Load libraries**: 


```r
library(tweeDEseqCountData) 
```

Description of the library: 

- **tweeDEseqCountData**: Here, we need the conditions of the samples of the pickrell data set because these are (in a processed form) required for the web-based tool GSEA. 


## Load data

For the web-based tool, we work with the pre-filtered and transformed gene expression measurements. Since GSEA accepts a wide variety of gene ID formats, we can input the gene IDs in the Ensembl ID format. \
Note that for the purpose of simplicity, we here work with the voom-transformed data. However, you can easily switch to the gene expression data set transformed using DESeq2's varianceStabilizingTransformation. 


```r
load("./data/Results_RNASeq_Transformation/expression_data_voomtransformed_Ensembl.Rdata")

# or alternatively:
# load("./data/Results_RNASeq_Transformation/expression_data_vsttransformed_Ensembl.Rdata")
```

For an easier readability, we store the gene expression data set in a new object with a neutral name. 


```r
expression_data_transformed <- expression_data_voomtransformed_Ensembl
```

As mentioned above, GSEA requires (a preprocessed version of) the sample conditions. In the case of the pickrell data set used for these illustrations, we obtain these using the following commands: 


```r
# load pickrell data 
data(pickrell)

# store sample conditions
sample_conditions <- pickrell.eset$gender
```


### step 1: Export the (transformed) gene expression measurements

We export the transformed gene expression measurements to a .txt file and into the pre-specified folder "Input_Objects_GSEA_web". 


```r
# 1. Generate the path that indicates that we will store the transformed gene expression measurements in the object "expression_data_transformed.txt, which is located in the file "Input_Objects_GSEA_web"
path_measurements  <- "./data/Input_Objects_GSEA_web/expression_data_transformed.txt"

# export the gene expression measurements 
write.table(expression_data_transformed,
            file = path_measurements,
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)
```

Arguments of function *write.table()*:

 - **'quote = FALSE'**: ensures that none of the characters (in this case the gene and sample
 identifiers, i.e. row and column names) are surrounded by double quotes
 
 - **'row.names = TRUE'**: ensures that the gene IDs are included in the export 
 
 - **'col.names = TRUE'**: ensures that the samples IDs are included in the export 
 
 
### step 2: Prepare and export the sample conditions

**Prepare the sample conditions**
 
GSEA accepts the sample conditions in the (binary) format with values 0 and 1. For the pickrell data set, we inspect the current values in which the sample conditions are stored:


```r
# inspect the raw sample conditions:
head(sample_conditions, n = 10)
#>  [1] male   male   female male   female male   female male  
#>  [9] female male  
#> Levels: female male

# inspect the levels of the sample conditions:
levels(sample_conditions)
#> [1] "female" "male"
```

Currently, the sample conditions are coded as "female" and "male". We now want to convert both levels to 0 and 1:



```r
# (i) create vector to contain the sample conditions in the right format
sample_conditions_prepared <- c()

# (ii) assign all "females", which is the first level of the factor, the value 0
sample_conditions_prepared[ pickrell.eset$gender == levels(pickrell.eset$gender)[[1]]] <- 0

# (iii) assign all "males", which is the second level of the factor, the value 1
sample_conditions_prepared[ pickrell.eset$gender == levels(pickrell.eset$gender)[[2]]] <- 1
```

Inspect the (first few) sample conditions: 

```r
head(sample_conditions_prepared, n = 10)
#>  [1] 1 1 0 1 0 1 0 1 0 1
```


**Export the sample conditions**: 


```r
# the following path indicates that we store the prepared sample conditions in the object "sample_conditions_prepared.txt" in the file "Input_Objects_GSEA_web". 
path_conditions <- "./data/Input_Objects_GSEA_web/sample_conditions_prepared.txt"

# export 
write.table(x = sample_conditions_prepared,
            file = path_conditions,
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)
```

Arguments of function *write.table()*:

 - **'quote = FALSE'**: ensures that none of the characters (in this case the gene and sample
 identifiers, i.e. row and column names) are surrounded by double quotes
 
 - **'row.names = TRUE'** ensures that the gene IDs are included in the export 
 
 - **'col.names = TRUE'** ensures that the samples IDs are included in the export 
 
 
 
### step 3: Further preparation of input object in Excel

For instructions of the further preparation of the text file in Excel, open the following link: \
http://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats \

Follow the instructions in section 

- **'GCT: Gene Cluster Text file format (.gct)'**: for further preparation of expression_data_transformed

- **'CLS: Categorical (e.g tumor vs. normal) class file format (.cls)'**: for further preparation of sample_conditions_prepared



