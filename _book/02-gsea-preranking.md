
# GSEA Preranking

Skip to Chapter 4 Differentiation Expression Analysis if conversion of gene IDs is not needed.


In this script, we will go through the process of preparing the required input for GSEAPreranked (which is part of the web-based tool GSEA) \
Note that in contrast to the remaining tools considered in this work, it is recom-
mended for GSEAPreranked to provide the input with the genes ID HGNC (HUGO) gene
symbols. The reason for this is described in the supplement of this work.

Therefore, our starting point is the pre-filtered gene expression data set with the
gene IDs in the initial (Ensembl) ID format.  

We proceed in the following order:

1. Conversion of the gene IDs to HGNC (HUGO) symbols and removal of resulting duplicated gene IDs

2. Differential expression analysis using voom/limma, DESeq2, and edgeR

3. Generation of the gene ranking from each results table generated in (ii)

4. Export of gene rankings to text file

5. Link to information on how to further process the gene rankings in Excel


## Libraries 

All necessary packages are available on Bioconductor, and should be installed from there if not already available on your machine.


```r
install.packages("BiocManager")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("tweeDEseqCountData")
BiocManager::install("limma")
BiocManager::install("edgeR")
BiocManager::install("DESeq2")
BiocManager::install("dplyr")
```


Load libraries: 


```r
library(clusterProfiler)
library(org.Hs.eg.db) 
library(tweeDEseqCountData)
library(limma)
library(edgeR)
library(DESeq2)
library(dplyr)
```

Description of the libraries: 

- clusterProfiler: To convert Ensembl ID format to HGNC (HUGO) gene symbols

- org.Hs.eg.db: To indicate that we work with human gene expression data (must be adapted when working with a different organism)

- tweeDEseqCountData: To obtain the conditions of the samples in the gene expression data set

- limma: For differential expression analysis (based on which the gene ranking is generated)

- edgeR: For differential expression analysis (based on which the gene ranking is generated)

- DESeq2: For differential expression analysis (based on which the gene ranking is generated)

- dplyr: To unify the relevant columns in the results of differential expression analysis across the three methods


### Load data 

We will proceed with the gene expression data set which was pre-filtered using edgeR's filterByExpr. 


```r
load("./data/Results_PreFiltering/expression_data_filterByExpr.Rdata")
```

Alternatively, you can also load the gene expression data set which was manually pre-filtered (as proposed by DESeq2):



```r
# load("./data/Results_PreFiltering/expression_data_filterDESeq2.Rdata")
```

For a simplified readability, we will proceed with the pre-filtered gene expression data set under a different name:


```r
expression_data_filtered <- expression_data_filterByExpr
```

## Step 1: Convert Ensembl IDs to HGNC symbols

#### (i) Obtain the mapping of initial (Ensembl) gene IDs to required HGNC gene symbols:



```r
bitr_EnsToSymb <- bitr(geneID = rownames(expression_data_filtered),
                       fromType = "ENSEMBL",
                       toType = "SYMBOL",
                       OrgDb = org.Hs.eg.db)
#> 'select()' returned 1:many mapping between keys and
#> columns
#> Warning in bitr(geneID =
#> rownames(expression_data_filtered), fromType = "ENSEMBL", :
#> 3.73% of input gene IDs are fail to map...
```
Function arguments: 

- gene IDs: gene IDs to be converted (stored in the rownames of the gene expression data set)

- fromType: initial gene ID format that is to be converted

- toType: gene ID format to be converted to 

- OrgDb: annotation data base of organisms from which the gene expression measurements originate. The corresponding argument must be loaded as a library (see above)

Note that the arguments "fromType" and "toType" must be set as one of the following, depending on the given and the required gene ID format: \
ACCNUM, ALIAS, Ensembl, EnsemblPROT, EnsemblTRANS, ENTREZID, ENZYME, EVIDENCE, EVIDENCEALL, GENENAME, GENETYPE, GO, GOALL, IPI, MAP, OMIM, ONTOLOGY, ONTOLOGYALL, PATH, PFAM, PMID, PROSITE, REFSEQ, SYMBOL, UCSCKG, UNIPROT.


#### (ii) Concatenate the initial gene expression data set with the mapping from initial (Ensembl) to required HGNC symbols

Merge by row names of expression data set and ENSEMBL ID of gene ID mapping. 



```r
expression_data_merged <- merge(x = expression_data_filtered,
                                             y = bitr_EnsToSymb,
                                             by.x = 0, 
                                             by.y = "ENSEMBL",
                                             sort=TRUE)
```

Description of function arguments: 

- by.x and by.y specify by which columns expression_data_filterByExpr and bitr_EnsToEntr are concatenated: 

- by.x = 0: use row names of expression_data_filterByExpr (the row names contain the gene IDs)

- by.y = "ENSEMBL": use the column that contains the Ensembl gene IDs 


As we can see from the document "Instructions_GeneID_Conversion_DuplicatesRemoval.R", we work with three separate ways to remove the duplicate gene IDs that result from
gene ID conversion. We will apply all three of these approaches to removal. 

## Step 2: Remove the duplicated gene IDs

### (i) Take a closer look at the duplicates

#### Case 1: single ENSEMBL IDs are mapped to multiple HGNC symbols

In this case, we denote the corresponding ENSEMBL IDs as the duplicated IDs


```r
# obtain number of cases in which an ENSEMBL gene ID was converted to several HGNC symbols, i.e. the number of
# times an Ensembl ID appears more than once in the mapping
sum(duplicated(bitr_EnsToSymb$ENSEMBL))
#> [1] 28
# determine all duplicated ENSEMBL gene IDS (i.e. all ENSEMBL gene IDs that were mapped to multiple distinct HGNC gene symbols):
dupl_ensembl<-unique(bitr_EnsToSymb$ENSEMBL[duplicated(bitr_EnsToSymb$ENSEMBL)])
# number of ENSEMBL IDs that have at least one duplicate
length(dupl_ensembl)
#> [1] 27
```

In the following, we can inspect the conversion pattern for each Ensembl ID that was mapped to multiple distinct HGNC symbols: 


```r
# get subset conversion pattern of all duplicated Ensembl IDs 
duplicated_conversion_ens<-bitr_EnsToSymb[bitr_EnsToSymb$ENSEMBL %in% dupl_ensembl,]
# take a look at conversion pattern for some of the duplicated Ensembl gene IDs:
head(duplicated_conversion_ens, n = 6)
#>              ENSEMBL          SYMBOL
#> 303  ENSG00000063587          ZNF275
#> 304  ENSG00000063587    LOC105373378
#> 615  ENSG00000088298           EDEM2
#> 616  ENSG00000088298 MMP24-AS1-EDEM2
#> 1521 ENSG00000111215            PRH1
#> 1522 ENSG00000111215            PRR4
```
We can see perfectly that these Ensembl ID are each mapped to multiple individual HGNC symbols. 

#### CASE 2: multiple ENSEMBL IDs are mapped to the same HGNC symbol

In this case, we denote the corresponding HGNC symbols as the duplicated IDs. 

Obtain number of cases in which multiple distinct Ensembl IDs are converted to the same HGNC gene symbol

```r
# in these cases, the corresponding HGNC symbols appear repeatedly in the mapping:
sum(duplicated(bitr_EnsToSymb$SYMBOL))
#> [1] 2
# determine all HGNC gene symbols affected by this duplication
dupl_symbol<-unique(bitr_EnsToSymb$SYMBOL[duplicated(bitr_EnsToSymb$SYMBOL)])
# number of HGNC symbols that have at least one duplicate
length(dupl_symbol)
#> [1] 2
```


In the following, we want to take a look at the conversion pattern of some of the duplicated HGNC symbols. 

```r
# subset the conversion pattern to that of the duplicated HGNC symbols
duplicated_conversion<-bitr_EnsToSymb[bitr_EnsToSymb$SYMBOL %in% dupl_symbol,]

# for visibility: order the conversion pattern by the HGNC symbols
duplicated_conversion <- duplicated_conversion[order(duplicated_conversion$SYMBOL),]

# take a look at conversion pattern of the first few duplicated HGNC symbols:
head(duplicated_conversion, n = 6 )
#>              ENSEMBL SYMBOL
#> 2940 ENSG00000140543   DET1
#> 4664 ENSG00000173867   DET1
#> 250  ENSG00000054277   OPN3
#> 5869 ENSG00000203668   OPN3
```

### (ii) Remove the duplicated gene IDs 

For the subsequent analysis (whether differential expression analysis of gene set analysis), we need to deal with the duplicated gene IDs and remove them in a suitable way such that only one unique gene expression measurement among the duplicates remains. There is no recommended way to proceed, i.e. no common approach presented in an official scientific puplication, but instead several approaches suggested by users in corresponding user platforms, generally, the manner of duplicate gene ID removal is performed at the discretion of the user. We therefore present three approaches to duplicate gene ID removal. 


#### Option 1: Keep the duplicate with the lowest subscript 

This is the simplest approach among the three. 


```r
#1. remove duplicated HGNC gene symbols
exprdat_symbol_dupl1<-expression_data_merged[!duplicated(expression_data_merged$SYMBOL),]


#2. remove duplicated ENSEMBL gene IDs
exprdat_symbol_dupl1<-exprdat_symbol_dupl1[!duplicated(exprdat_symbol_dupl1$Row.names),]
dim(exprdat_symbol_dupl1)
#> [1] 6012   71

#3. set HGNC gene symbols as rownames
rownames(exprdat_symbol_dupl1)<-exprdat_symbol_dupl1$SYMBOL
#Remove columns containing ENSEMBL and HGNC symbols
exprdat_symbol_dupl1<-subset(exprdat_symbol_dupl1, select=-c(Row.names,SYMBOL))
dim(exprdat_symbol_dupl1)
#> [1] 6012   69
# we particularly see that the number of columns is now at its initial number
```

We now want to inspect a small part of the gene expression data set: 

```r
exprdat_symbol_dupl1[1:10, 1:10]
#>          NA18486 NA18498 NA18499 NA18501 NA18502 NA18504
#> DPM1          22     105      40      55      67      37
#> SCYL3         22     100     107      53      72      38
#> C1orf112       5      23      10      18      15       8
#> FGR           36      70      41      33      59      29
#> FUCA2         29      79      33      31      29      21
#> NFYA         301     351     344     176     340     170
#> LAS1L         32      22       4       8      15       5
#> ENPP4         52     419     173     137     127      72
#> BAD           48      81      65      19      59      57
#> HS3ST1         8      15       8      18       4      27
#>          NA18505 NA18507 NA18508 NA18510
#> DPM1          88     127      70      43
#> SCYL3         98      69      66      43
#> C1orf112      11      16      18       7
#> FGR           22      71      12      43
#> FUCA2         42      62      41      25
#> NFYA         238     247     226     227
#> LAS1L         15      16      19      15
#> ENPP4        267     247     154      81
#> BAD           44      79      71      79
#> HS3ST1         3      11       5       3
```


#### Option 2: Keep the (rounded) mean expression value of all duplicated gene IDs

Here, we remove duplicated HGNC gene symbols (case 2) before removing the duplicated Ensembl gene IDs (case 1). The reason for this is elaborated below. 

##### (i) Remove the duplicated HGNC gene symbols (case 2)


```r
#generate matrix to contain (rounded) mean expression values of all rows that have same HGNC gene symbol
#ncol=ncol(expression_data_filterByExpr)-2 since data set contains 2 columns with IDs at this point
mean_symbol <- matrix(nrow=0, ncol=ncol(expression_data_merged)-2)


# for each duplicated HGNC gene symbol separately, we gather all rows with the corresponding gene expression data and then extract the (rounded) mean expression value of all rows
for(i in 1:length(dupl_symbol)){

# go through each HGNC symbols which occurs multiple times
# determine all rows whose HGNC symbols correspond to currently considered HGNC symbol
counts_dupl <- expression_data_merged[expression_data_merged$SYMBOL %in% unique(dupl_symbol)[i],]

#compute the mean expression value of all rows that contain to the given HGNC gene symbol
dupl_id <- round(colMeans(counts_dupl[,c(2:(ncol(expression_data_merged)-1))]))

# store rounded mean expression value in matrix
# this matrix is extended by a single row of gene expression data which corresponds to the (rounded) mean expression data that corresponds to the given HGNC gene symbol
mean_symbol <- rbind(mean_symbol,dupl_id)


}

# after completing the for-loop, mean_symbol contains the mean expression measurements of each HGNC gene symbol which contains duplicates resulting from gene ID conversion
# We want to take a look at a part of the data set 
mean_symbol[, 1:9]
#>         NA18486 NA18498 NA18499 NA18501 NA18502 NA18504
#> dupl_id      37     129      64      54      70      54
#> dupl_id      50     182      82     118      60      52
#>         NA18505 NA18507 NA18508
#> dupl_id      62     123      70
#> dupl_id     102      80      78

#set corresponding HGNC gene symbols as rownames
rownames(mean_symbol) <- unique(dupl_symbol)


# test whether the number of rows in mean_symbol corresponds to the number HGNC symbols
# that occur more than once
# result should be TRUE
nrow(mean_symbol) == length(dupl_symbol)
#> [1] TRUE

# remove all rows from the expression data whose HGNC symbol has at least one duplicate
# intuition: we have just dealt with the corresponding rows and do not want them to be considered in the second step (which deals with case 2

exprdat_symbol_dupl2 <- expression_data_merged[!expression_data_merged$SYMBOL %in% dupl_symbol,]

# test whether number of rows in resulting data set equals nrow of inital data set minus number of genes with at least one duplicate
nrow(exprdat_symbol_dupl2) == nrow(expression_data_merged)-nrow(duplicated_conversion)
#> [1] TRUE
dim(exprdat_symbol_dupl2)
#> [1] 6037   71
```

##### (ii) Remove the duplicated ENSEMBL IDs

Caution: Single ENSEMBL IDs that are mapped to multiple HGNC symbol naturally generate identical count data for all corresponding HGNC symbols. It is therefore pointless to compute mean expression values. This is verifiable by looking at data set only containing those ENSEMBL IDs that are mapped by multiple HGNC symbols:


```r
test_dupl_ensembl <- exprdat_symbol_dupl2[exprdat_symbol_dupl2$Row.names %in% dupl_ensembl,]

# take a look at a few entries of the data set:
test_dupl_ensembl[1:6, 1:6]
#>            Row.names NA18486 NA18498 NA18499 NA18501
#> 299  ENSG00000063587     143     116     244      76
#> 300  ENSG00000063587     143     116     244      76
#> 610  ENSG00000088298     202     384     368     207
#> 611  ENSG00000088298     202     384     368     207
#> 1502 ENSG00000111215      17      12       9       2
#> 1503 ENSG00000111215      17      12       9       2
#>      NA18502
#> 299      157
#> 300      157
#> 610      397
#> 611      397
#> 1502       8
#> 1503       8
```
We see that those rows that originate from the same Ensembl ID contain identical count data. \

We therefore proceed as in option 1 and use the HGNC symbol that occurs first while we remove the rest. 

```r
# Keep HGNC symbol that occurs first
exprdat_symbol_dupl2<-exprdat_symbol_dupl2[!duplicated(exprdat_symbol_dupl2$Row.names),]

# inspect dimension
dim(exprdat_symbol_dupl2)
#> [1] 6010   71

#set HGNC symbol as rownames
rownames(exprdat_symbol_dupl2)<-exprdat_symbol_dupl2$SYMBOL

#remove any columns containing IDs
exprdat_symbol_dupl2<-subset(exprdat_symbol_dupl2,select= -c(Row.names,SYMBOL))

#add those rows to data set that contain mean expression values of duplicate HGNC symbols
exprdat_symbol_dupl2 <- rbind(exprdat_symbol_dupl2,mean_symbol)

#inspect the dimension of the remaining gene expression data set
dim(exprdat_symbol_dupl2)
#> [1] 6012   69
```

#### Option 3: Among the duplicates, keep row with highest overall expression values (i.e highest counts across all samples)

The intuituin behind this approach is that the row with highest counts values has  highest power of being detected as differentially expressed. As in option 2, this applies only to the duplicates that result from multiple ENSEMBL IDs being mapped to the same HGNC symbol. 

##### (i) Remove duplicated HGNC symbols


```r
# generate matrix to later contain row with highest count values among ID duplicates this data set is to be filled gradually and with each iteration of the follwing for-loop
highest_count_symbol<-matrix(, nrow=0, ncol=ncol(expression_data_filterByExpr))


# for each duplicated HGNC gene symbol separately, we gather all rows with the corresponding gene expression data and then extract the row with the highest overall magnitude of counts

for(i in 1:length(dupl_symbol)){

  # go through each HGNC symbols which occurs multiple times
  # determine all rows whose HGNC symbols correspond to currently considered HGNC symbol
  counts_dupl<-expression_data_merged[expression_data_merged$SYMBOL %in% unique(dupl_symbol)[i],]

  # order rows in decreasing manner by their number of read counts across all samples
  order_rowsums<-order(rowSums(counts_dupl[,2:(ncol(counts_dupl)-1)]),decreasing=TRUE)
  
  #detect row with highest number of read counts across all samples (i.e. row with rank 1)
  dupl_id<-counts_dupl[order_rowsums==1,]
  
  #store corresponding expression
  highest_count_symbol<-rbind(highest_count_symbol,dupl_id)
  
}


#Remove all initial values with HGNC  duplicates from the dataset initial gene expression data set
exprdat_symbol_dupl3<-expression_data_merged[! expression_data_merged$SYMBOL %in% unique(dupl_symbol),]
```


##### Case 2: Remove duplicated ENSEMBL IDs 

As in option 2, it is pointless to detect the row with highest count values as all rows
that correspond to the same ENSEMBL ID as these rows naturally contain identical count data. We therefore remove the duplicate ENSEMBL ID that occurs first (such as in option 1). 


```r
# Keep the corresponding ENSEMBL ID that occurs first
exprdat_symbol_dupl3<-exprdat_symbol_dupl3[!duplicated(exprdat_symbol_dupl3$Row.names),]

# Add the gene expression rows of all HGNC gene symbols that were initially duplicated
exprdat_symbol_dupl3<-rbind(exprdat_symbol_dupl3,highest_count_symbol )

#Set HGNC symbols as rownames
rownames(exprdat_symbol_dupl3)<-exprdat_symbol_dupl3$SYMBOL
#Remove any column that contains information on gene IDs
exprdat_symbol_dupl3<-subset(exprdat_symbol_dupl3, select=-c(Row.names,SYMBOL))

# Inspect the dimension of the resulting gene expression data set 
dim(exprdat_symbol_dupl3)
#> [1] 6012   69
# we see that the sample size (number of columns) is now back at its initial number
```

## Intermediate step: Choose which converted gene expression data sets to proceed with 

In this illustration, we will proceed with only one of the three gene expression data sets that result from the conversion of the gene IDs and the removal of the duplications. Here, we will proceed with the first approach to
the removal of duplicated gene IDs. However, you can easily switch to another approach at your discretion. 

We will proceed with the gene expression data set stored as object "expression_data_filt_symbol". 

```r
expression_data_filt_symbol <- exprdat_symbol_dupl1
```

## Step 3: Differential expression analysis 

From the corresponding results table of differential expression analysis, we will later transform certain metrics (i.e. columns) into a gene-level statistic for each gene based on which the gene ranking is generated.

We will proceed analogously to the script "Instructions_Differential_Expression_Analysis", in which we perform differential expression analysis using three established methods, namely

1. voom/limma
2. DESeq2
3. edgeR

While we already have the count data at hand that we will be working with, we additionally need the conditions of the samples. Since we work with the Pickrell data set in these illustrations, we have to load the corresponding data.


```r
# Load pickrell data 
data(pickrell)

# Store the sample conditions in the object sample_conditions 
sample_conditions <- pickrell.eset$gender
```


### Option 1: Differential expression analysis using limma 

Note that the code illustrations are based on the following user manual: \ https://bioconductor.org/packages/devel/bioc/vignettes/limma/inst/doc/usersguide.pdf


#### (i) Generate the required input object for limma 

As described in "Instructions_Differential_Expression_Analysis", the gene expression data set expression_data_filt_symbol is currently a data frame and needs to be converted to a DGEList object. 


```r
# store expression data with corresponding sample conditions in object of class DGEList
y <- DGEList(counts = expression_data_filt_symbol,
             group = sample_conditions)
```

Required function arguments: 

- counts: matrix that contains RNA-Seq data (i.e. count data)

Optional function arguments: 

- group: indicates the condition of each sample/library
  + this argument must be specified sooner or later (such as in subsequent functions) so we just specify it at this point.

Note that we leave the remaining arguments in their default state since the corresponding info will be added through the following pipeline of functions. 


#### (ii) Normalization 

The following piece of code generates a normalization factor for each sample. This accounts for sample-specific effects (such as differences in library sizes and effects of compositionality). If not accounted for, these effects prevent a comparison between the samples. 


```r
y <- calcNormFactors(y)
```

Note that this function does not transform the count data, but rather generates a normalization factor for each sample which is incorporated into the subsequent analysis separately. 


#### (iii) voom transformation 

Note that the voom transformation facilitates the use of the subsequent functions that were initially developed for microarray data. 



```r
# design matrix (rows correspond to samples, columns indicate which coefficients are to be estimated)
design_matrix <- model.matrix(~sample_conditions)

# voom transformation (transformation of count data to log-cpm values, computation of observation weight for each
# gene based on mean-variance relationship)
y <- voom(y, design = design_matrix)
```

#### (iv) Differential expression analysis 


```r
# Fit a linear model for each gene
y <- lmFit(y)

# Calculate the statistics for the assessment of differential expression
y <- eBayes(y)

# Get the results table for each gene whose differential expression was assessed
DE_results_limma <- topTable(y, number =  nrow(y))
#> Removing intercept from test coefficients
# Note that number = nrow(y) ensures that all genes are displayed in results, not just a subset
```

#### (v) Rename columns in results table

This step is not required for differential expression analysis or for the subsequent use of gene set analysis IN GENERAL. 

A renaming of some columns (such as adjusted p-value) is necessary in the context of these illustrations since the different methods for differential expression analysis typically differ in the column names of the results tables. A unification of the column names is required so that we can use the same code to illustrate further conduct of GSA in the subsequent code. 


```r
# First, we transform to results table to a data frame so that we see the results table directly when accessing it through the name "res"
DE_results_limma <- as.data.frame(DE_results_limma)

# Rename the column that contains the adjusted p-values: rename padj to p_adj
DE_results_limma <- dplyr::rename(DE_results_limma, p_adj = `adj.P.Val`)
```

Inspect the first 10 genes in the results table: 

```r
head(DE_results_limma, n = 10)
#>              logFC   AveExpr         t      P.Value
#> RPS4Y1   9.2108222 1.9863775 47.179852 1.030921e-55
#> TMSB4Y   5.0808942 0.5129662 26.670268 6.040369e-39
#> PNPLA4  -0.9237796 5.3334920 -8.775861 5.661613e-13
#> PUDP    -0.8541771 2.1627703 -4.413691 3.518147e-05
#> CXorf38 -0.5505387 3.7163801 -4.249159 6.357396e-05
#> TXLNG   -0.4040253 5.1161343 -4.042236 1.316055e-04
#> G0S2    -1.4699085 5.1923255 -3.854451 2.504228e-04
#> JUN     -0.5877351 8.2731730 -3.918457 2.014918e-04
#> APAF1    0.4049046 5.9682315  3.825437 2.761789e-04
#> IL1A    -1.4852480 3.5436708 -3.665932 4.695371e-04
#>                p_adj           B
#> RPS4Y1  6.197899e-52 84.19116319
#> TMSB4Y  1.815735e-35 62.68059267
#> PNPLA4  1.134587e-09 19.20851972
#> PUDP    5.287775e-02  2.14585449
#> CXorf38 7.644132e-02  1.55146299
#> TXLNG   1.318687e-01  0.73174212
#> G0S2    1.844875e-01  0.15225728
#> JUN     1.730527e-01  0.11384668
#> APAF1   1.844875e-01 -0.06540026
#> IL1A    2.645036e-01 -0.24633830
```

### Option 2: Differential expression analysis using DESeq2 

The official DESeq2 vignette can be found through the following link: \
http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html


#### (i) Generate the required input object for DESeq2 

DESeq2 operates on the format DESeqDataSet which contains information on the count data, the conditions of the samples and the design (for further information see below). Note that for gene expression data sets that must be imported to R, additional steps are necessary before the following code can be run. 


##### Generate a data frame which contains the condition of each sample

Here, the information on the sample conditions is the only column in the data frame. However, if further variables (such as batch effects) are to be controlled for, the corresponding variables must additionally be added to coldata. 



```r
# Here, the information on the sample conditions is the only column in the data frame
# However, if further variables (such as batch effects) are to be controlled for, the corresponding variables must additionally be added to coldata

# The names of the samples are stored as the row names of the data frame
# Important: make sure that the order of the conditions in sample_conditions corresponds to the order of the
# samples in expression_data
coldata <- data.frame(sample_conditions,
                      row.names = colnames(expression_data_filt_symbol))

# Rename the column header to "condition"
colnames(coldata) <- "condition"

# Recode the variable "condition" as a factor
# Rename the sample conditions (in this case from "female" and "male") to "untreated" and "treated"
# Note: make sure that the control level in variable condition is coded as the first level (i.e. "untreated")
coldata$condition <- factor(coldata$condition,
                            labels = c("untreated","treated"))
```

##### Generate the DESeqDataSet 



```r
dds <- DESeqDataSetFromMatrix(countData = expression_data_filt_symbol,
                            colData = coldata,
                            design = ~ condition)
```

Relevant function arguments: 

- countData: count data from gene expression data set

- colData: data frame that contains information on the samples (see above)
  + conditions of the samples (required) and possibly further variables to correct for (such as batch effects)

- design: indicates which variables from colData are used for modelling
  + more detailed: the argument design is used to estimate the dispersions and the log2 fold changes of the model
  + if more than one variable from colData are used in argument design (e.g. a second variable "batch"), the syntax changes to the following formula: design ~ batch + condition
  + make sure that the variable of interest (here: variable that represents conditions of the samples) is placed at the end of the formula


#### (ii) Differential expression analysis 


```r
# 1. Perform default differential expression analysis
dds <- DESeq(dds)
# 2. Generate results table which provides
DE_results_DESeq2 <- results(dds)
```
Description of the functions: 

1. DESeq(): Estimation normalization factors, estimation of dispersions, fitting of generalized linear model, Wald statistics

2. results(): Provides base mean across all samples, estimated log2 fold changes, standard errors, test statistics, p-values, adjusted p-values

#### (iii) Rename columns in results table

This step is not required for differential expression analysis or for the subsequent use of gene set analysis IN GENERAL. 


A renaming of some columns (such as adjusted p-value) is necessary in the context of these illustrations since the different methods for differential expression analysis typically differ in the column names of the results tables. A unification of the column names is required so that we can use the same code to illustrate further conduct of GSA in the subsequent code. 



```r
# First, we transform to results table to a data frame so that we see the results table directly when accessing it through the name "res"
DE_results_DESeq2 <- as.data.frame(DE_results_DESeq2)

# Rename the column that contains the adjusted p-values: rename padj to p_adj
DE_results_DESeq2 <- dplyr::rename(DE_results_DESeq2, p_adj = "padj")
```

Inspect the first 10 genes in the results table: 

```r
head(DE_results_DESeq2, n = 10)
#>           baseMean log2FoldChange      lfcSE        stat
#> DPM1      62.75469    0.039750916 0.11391821  0.34894259
#> SCYL3     70.48070    0.029018887 0.11169067  0.25981478
#> C1orf112  13.00584    0.165462332 0.15395569  1.07473997
#> FGR       43.59294    0.114160979 0.19557799  0.58371076
#> FUCA2     37.97983   -0.006418336 0.10904211 -0.05886108
#> NFYA     228.55761    0.051517676 0.09634003  0.53474838
#> LAS1L     11.53808    0.025279778 0.18769732  0.13468375
#> ENPP4    175.90798    0.300378628 0.16292558  1.84365538
#> BAD       66.39161   -0.264365320 0.10275897 -2.57267378
#> HS3ST1    10.12732   -0.199017932 0.31888446 -0.62410671
#>              pvalue     p_adj
#> DPM1     0.72713241 0.9221794
#> SCYL3    0.79500665 0.9425284
#> C1orf112 0.28249116 0.7202446
#> FGR      0.55941492 0.8668048
#> FUCA2    0.95306276 0.9863983
#> NFYA     0.59282386 0.8789008
#> LAS1L    0.89286192 0.9709726
#> ENPP4    0.06523338 0.4823900
#> BAD      0.01009163 0.3001405
#> HS3ST1   0.53255751 0.8522562
```


### Option 3: Differential expression analysis using edgeR 

The instructions for edgeR are based on the following vignette: \ https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf

#### (i) Generate the required input object for edgeR 

Store the expression data with corresponding sample conditions in object of class DGEList


```r
y <- DGEList(counts = expression_data_filt_symbol,
             group = sample_conditions)
```
Required function arguments: 

- counts: matrix that contains RNA-Seq data (i.e. count data)

Optional function arguments: 

- group: indicates condition of each sample/library
  + This argument must be specified sooner or later (such as in subsequent functions) so we just specify it at this point

Note that we leave the remaining arguments in their default state since the corresponding info will be added through the subsequent pipeline of functions.


#### (ii) Normalization 

The following piece of code generates a normalization factor for each sample. This accounts for sample-specific effects (such as differences in library sizes and effects of compositionality). If not accounted for, these effects prevent a comparison between the samples



```r
y <- calcNormFactors(y)
```

Note: this function does not transform the count data, but rather generates a normalization factor for each sample which is incorporated into the subsequent analysis separately. 

#### (iii) Estimation of the dispersion 

Estimate the common and tagwise dispersion which quantify the variation of the true abundance of a given gene between different samples and is required to assess differential expression realistically. 


```r
y <- estimateDisp(y)
#> Using classic mode.
```

#### (iv) Differential expression analysis 


```r
# Test each gene for differential expression:
DE_results_edgeR <- exactTest(y)

# Extract pre-specified number (n) of genes
DE_results_edgeR <- topTags(DE_results_edgeR, n = nrow(DE_results_edgeR))

# Note: argument n specifies the number of top differentially expressed genes to be displayed in the results
# n = nrow(DE_results) ensures the results of all genes whose differential expression was assessed are displayed
```

#### (v) Rename columns in results table 


This step is not required for differential expression analysis or for the subsequent use of gene set analysis IN GENERAL. 

A renaming of some columns (such as adjusted p-value) is necessary in the context of these illustrations since the different methods for differential expression analysis typically differ in the column names of the results tables. A unification of the column names is required so that we can use the same code to illustrate further conduct of GSA in the subsequent code. 


```r
# first, we transform to results table to a data frame so that we see the results
# table directly when accessing it through the name "res"
DE_results_edgeR <- as.data.frame(DE_results_edgeR)

# rename column that contains adjusted p-values: rename padj to p_adj
DE_results_edgeR <- dplyr::rename(DE_results_edgeR, p_adj = "FDR")
```

Inspect the first 10 genes in the results table:

```r
head(DE_results_edgeR, n = 10)
#>              logFC   logCPM        PValue         p_adj
#> RPS4Y1  10.1333123 6.155249  0.000000e+00  0.000000e+00
#> TMSB4Y   5.5794301 2.519825 6.410059e-121 1.926864e-117
#> PNPLA4  -0.9445304 5.485374  1.886604e-17  3.780754e-14
#> PUDP    -0.7765765 2.571056  3.068012e-06  4.611222e-03
#> SNAI1   -1.6458976 3.288273  4.706338e-06  5.476681e-03
#> GADD45G -1.5869489 2.679748  5.465750e-06  5.476681e-03
#> THEM5   -1.1254688 2.303001  1.950080e-05  1.674840e-02
#> CTSW     1.5021901 3.365996  2.908369e-05  2.185640e-02
#> TXLNG   -0.4250226 5.207714  3.574779e-05  2.284273e-02
#> TCFL5   -0.6026382 5.756625  3.799523e-05  2.284273e-02
```
Note that one gene has an adjusted p-value of 0. We will deal with this issue when generating the gene ranking from edgeR. 


## Step 4: Generate gene ranking 

We will now apply a gene-level statistic to the results table of differential expression analysis to generate a gene ranking which will serve as input to GSEAPreranked. 

As we have performed differential expression analysis with
- voom/limma
- DESeq2
- edgeR,
we will also generate three gene rankings, one for each of these three methods. 


The formula of the gene-level ranking metric is \ 
$ -1 * log10(p-value) * sign(log fold change)$. \
Note that by p-value, we refer to the non-adjusted p-value. 


### Option 1: Generate gene ranking from limma/voom results 


```r
# 1. Subset the gene expression data set to those genes that have a p-value (i.e.
# which have been NOT been excluded from differential expression analysis)

# indicate those genes WITH a p-value
ind_nonNA_pvalue_limma <- !is.na(DE_results_limma$P.Value)

# subset gene expression data set to those genes with a p-value
DE_results_noNA <- DE_results_limma[ind_nonNA_pvalue_limma, ]

# 2. create vector that contains the value of the gene-level ranking metric for each gene
rankvec_limma <- sign(DE_results_noNA$logFC)*(-1)*log10(DE_results_noNA$P.Value)

# 3. assign respective gene ID to each value in the vector
names(rankvec_limma) <- rownames(DE_results_noNA)

# 4. sort the vector in descending order
rankvec_limma <- sort(rankvec_limma, decreasing=TRUE)
```


Inspect the first 10 genes from the gene ranking 

```r
head(rankvec_limma, n = 10)
#>    RPS4Y1    TMSB4Y     APAF1      LPXN     TRAM2   VIPAS39 
#> 54.986774 38.218937  3.558810  3.205103  3.168007  3.079378 
#>    MAP3K1    SLC9A6      PHKB     HCFC2 
#>  3.021735  2.928611  2.878186  2.820192
```

### Option 2: Generate gene ranking from DESeq2 results 

##### (i) Shrink estimated log fold change values

Before generating the gene ranking from the DESeq2 results, we need to perform an additional shrinkage that is specific to DESeq2. The intuition behind the shrinkage of log fold change values is as follows: As mentioned in the paper, RNA-Seq data consists (in its rare form) of count data in which is inherently heteroscedastic, i.e. the variance of the count data depends on the mean count of the count data. It is observable that ratios between counts are considerably noisier in low magnitudes of counts compared to higher magnitudes, i.e.the log fold changes between both conditions are higher if the overall magnitude of counts is lower, independent of the acual extent to which the gene is differentially expressed between both samples. DESeq2 addresses this issue by shrinking the estimated log fold changes in the direction of 0. The magnitude of shrinkage is higher if the available information for a gene is lower (which may be because of a low magnitude of counts, a high dispersion or few degrees of freedom.) A more detailed description is provided in the DESeq2 paper by Love et al. (2014). \

The shrinkage is performed on object dds, which is a result of the function DESeq(): 

```r
if (!any(installed.packages() %in% "apeglm")) {
  BiocManager::install("apeglm", update = FALSE, ask = FALSE)
}
DE_results_DESeq2_shrink <- lfcShrink(dds,
                                      coef = "condition_treated_vs_untreated",
                                      type="apeglm")
```

Function arguments: 

- type: method to perform shrinkage
  + we opt for the default "apeglm" but you can choose from two alternative options as well

- coef: indicate the coefficients to be shrunk
  + we can obtain the right argument from the following function call:


```r
resultsNames(dds)
#> [1] "Intercept"                     
#> [2] "condition_treated_vs_untreated"
```

This shows us that we can either shrink the intercept or the "condition_treated_vs_untreated". Since we do not want to shrink the intercept but the log fold changes, we opt for the second option "condition_treated_vs_untreated". 

Finish up the results table: 

```r
# transform results table to data frame
DE_results_DESeq2_shrink <- as.data.frame(DE_results_DESeq2_shrink)

# Inspect the first 10 genes from the results table: 
head(DE_results_DESeq2_shrink, n = 10)
#>           baseMean log2FoldChange      lfcSE     pvalue
#> DPM1      62.75469    0.007239544 0.04886163 0.72713241
#> SCYL3     70.48070    0.005442790 0.04843766 0.79500665
#> C1orf112  13.00584    0.018952417 0.05474525 0.28249116
#> FGR       43.59294    0.008111253 0.05228672 0.55941492
#> FUCA2     37.97983   -0.001170671 0.04792965 0.95306276
#> NFYA     228.55761    0.004565422 0.04684065 0.59282386
#> LAS1L     11.53808    0.001892317 0.05136660 0.89286192
#> ENPP4    175.90798    0.035048816 0.06668445 0.06523338
#> BAD       66.39161   -0.153096248 0.13831883 0.01009163
#> HS3ST1    10.12732   -0.005332940 0.05301329 0.53255751
#>               padj
#> DPM1     0.9221794
#> SCYL3    0.9425284
#> C1orf112 0.7202446
#> FGR      0.8668048
#> FUCA2    0.9863983
#> NFYA     0.8789008
#> LAS1L    0.9709726
#> ENPP4    0.4823900
#> BAD      0.3001405
#> HS3ST1   0.8522562
```
Note that the estimated log fold change values (column log2FoldChange) are smaller (in absolute terms!) compared to the regular DESeq2 results, while the (adjusted and non-adjusted) p-value remain unchanged. \
Also we have not performed shrinkage in Instructions_Differential_Expression_Analysis as there, the goal was to get the list of differentially expressed genes which are detected solely based on the adjusted p-value. 

##### (ii) Generate the gene ranking


```r
# 1. Subset the gene expression data set to those genes that have a p-value (i.e.
# which have been NOT been excluded from differential expression analysis)

# indicate those genes WITH a p-value
ind_nonNA_pvalue <- !is.na(DE_results_DESeq2_shrink$pvalue)

# subset gene expression data set to those genes with a p-value
DE_results_noNA <- DE_results_DESeq2_shrink[ind_nonNA_pvalue, ]

# 2. create vector that contains the value of the gene-level ranking metric for each gene
rankvec_DESeq2 <- sign(DE_results_noNA$log2FoldChange)*(-1)*log10(DE_results_noNA$pvalue)

# 3. assign respective gene ID to each value in the vector
names(rankvec_DESeq2) <- rownames(DE_results_noNA)

# 4. sort the vector in descending order
rankvec_DESeq2 <- sort(rankvec_DESeq2, decreasing=TRUE)
```

Inspect the first 10 genes from the gene ranking 

```r
head(rankvec_DESeq2, n = 10)
#>     RPS4Y1     TMSB4Y       CTSW      APAF1     SLC9A6 
#> 204.389769  63.911362   4.463494   3.750030   3.568346 
#>       LPXN      HCFC2     MAP3K1    VIPAS39      TRAM2 
#>   3.559146   3.342299   3.281639   3.273736   3.174903
```



### Option 3: Generate gene ranking from edgeR results 



```r
# 1. Subset the gene expression data set to those genes that have a p-value (i.e. 
# which have been NOT been excluded from differential expression analysis)

# indicate those genes WITH a p-value
ind_nonNA_pvalue <- !is.na(DE_results_edgeR$PValue)

# subset gene expression data set to those genes with a p-value 
DE_results_noNA <- DE_results_edgeR[ind_nonNA_pvalue, ]

# 2. create vector that contains the value of the gene-level ranking metric for each gene 
rankvec_edgeR <- sign(DE_results_noNA$logFC)*(-1)*log10(DE_results_noNA$p_adj)

# 3. assign respective gene ID to each value in the vector 
names(rankvec_edgeR) <- rownames(DE_results_noNA)

# 4. sort the vector in descending order 
rankvec_edgeR <- sort(rankvec_edgeR, decreasing = TRUE)
```

Inspect the first 10 genes from the gene ranking 

```r
head(rankvec_edgeR, n = 10)
#>      RPS4Y1      TMSB4Y        CTSW       APAF1        LPXN 
#>         Inf 116.7151490   1.6604215   1.0364270   1.0022766 
#>      SLC9A6      MAP3K1       TRAM2        ATP8      MICAL2 
#>   0.9838143   0.8656713   0.8517879   0.8455509   0.8230068
```

Note that gene RPS4Y1 has a ranking value of Inf (infinity). The reason behind this is that the gene has an adjusted p-value of 0 in the edgeR results of differential expression analysis (as already observed above). We have to solve this issue, however, there is no common default approach. For the illustration purposes here, we replace the value Inf by the second biggest ranking value in the ranking. 


```r
# replace value Inf by second biggest value in the ranking 
rankvec_edgeR[rankvec_edgeR == Inf] <- max(rankvec_edgeR[rankvec_edgeR != Inf])

# inspect ranking again
head(rankvec_edgeR, n = 10)
#>      RPS4Y1      TMSB4Y        CTSW       APAF1        LPXN 
#> 116.7151490 116.7151490   1.6604215   1.0364270   1.0022766 
#>      SLC9A6      MAP3K1       TRAM2        ATP8      MICAL2 
#>   0.9838143   0.8656713   0.8517879   0.8455509   0.8230068
```


## Step 5: Export gene ranking to text file

For the purpose of simplicity, we only export one of the gene rankings to a text file. Here, we export the ranking generated with limma. 


```r
# this path indicates that we store the gene ranking as the text file 
# "gene_ranking.txt" in folder "Input_Objects_GSEAPreranked"
path_conditions <- "./data/Input_Objects_GSEAPreranked/gene_ranking.txt"

write.table(x = rankvec_limma,
            file = path_conditions,
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)
```

Function arguments: 

- x: object to be exported to a text file 

- file: location (i.e. folder and file) in which the object is to be stored 

 - quote = FALSE ensures that none of the characters (in this case gene and sample identifiers) are surrounded by double quotes
 
- row.names = TRUE ensures that the gene IDs are included in the export 

- col.names = TRUE ensures that no gene information is in the first row since the web-based tool ignores whatever is in the first row


## Step 6: Further preparation in Excel 

For this step, follow the instructions on the very bottom of the following link: \
http://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats

Then, scroll to section "RNK: Ranked list file format (*.rnk)". 
