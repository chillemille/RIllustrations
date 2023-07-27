--- 
title: "A Practical Guide through Running Gene Set Analysis in R"
author: "Milena Wünsch and Pat Callahan"
date: "2023-07-27"
site: bookdown::bookdown_site
documentclass: book
bibliography: [book.bib, packages.bib]
# url: your book url like https://bookdown.org/yihui/bookdown
# cover-image: path to the social sharing image like images/cover.jpg
description: |
  This is a book about...
biblio-style: apalike
csl: chicago-fullnote-bibliography.csl
---


```r
knitr::opts_chunk$set(cache = TRUE)
```

# About {-}

This book provides a practical illustration of all steps required to run the GSA tools from the selection in this paper. Even though the files build upon each other, the illustrations are structured in a way that each step can be viewed independently. \
The RNA-Seq data set we will use for this illustration is provided by Pickrell et al. (2010) and is part of the library tweeDEseqCountData. The data set contains the RNA-Samples extracted from the lymphoplastoid cell lines of 69 independent Nigerian individuals. The sample conditions (i.e. phenotypes) of the respective samples correspond to the gender of the individuals. Note that in this gene expression data set, the gene IDs are initially in the form of Ensembl IDs [@cunningham2022ensembl].


To run the entire pipeline in this illustration, from pre-filtering to the actual GSA tool, several R packages are required. We here provide an overview of all packages that must be installed (if not already available) and loaded. In the sections that correspond to the individual steps, we give an overview of the locally required libraries. 

## Required R packages

### Available on CRAN: 

- **BiocManager** [@huber2015orchestrating]: offers functionalities to install and manage the packages that are a part of the Bioconductor project. Here, we particularly use it to install the required Bioconductor packages.

- **dplyr**: provides various functionalities for data manipulation. \
We use it here to unify the results tables of the different methods for differential expression analysis. While this is not necessary in general, it simplifies our illustrations as we can use the same GSA pipelines following differential expression analysis, independent of the method for differential expression analysis used.


### Available on Bioconductor 

To install the following libraries, the library **BiocManager** must be installed and loaded first.  

- **tweeDEseqCountData** [@tweeDEseq_package]: provides count data from three RNA-Seq experiments together with information such as the conditions of the samples. We work work the gene expression data set provided by Pickrell et al. (2010) [@pickrell2010understanding] in our illustrations. 

- **Biobase** [@huber2015orchestrating]: contains standardized data structures to represent genomic data. We need this package to access the count data of the pickrell data set from the library 'tweeDEseqCountData' using function *exprs()*.

 - **clusterProfiler** [@wu2021clusterprofiler]: offers functionalities to explore the functional characteristics of high-throughput genomic data. Here, we use its functions to perform ORA and GSEA and for the conversion of the gene IDs. 

- **org.Hs.eg.db** [@orgHs_package]: provides the genome-wide annotation for the human organism. Note that if your gene expression measurements are extracted from a different organism, you must adapt this library accordingly. For instance, the organism mouse requires the library org.Mm.eg.db. \
This library is often required when working with the library clusterProfiler. 

- **limma** [@ritchie2015limma]: provides functionalities for analysing data from gene expression experiments. limma allows the analysis of microarray as well as RNA-Seq data for which the analysis pipelines are similar. We here use limma for differential expression analysis and make use of a transformation of the RNA-Seq data to align its characteristics for microarray data. 

- **edgeR** [@robinson2010edger]: provides functionalities to analyse data from RNA sequencing technologies and focuses on differential expression analysis with biological (as opposed to technical) replicates. 

- **DESeq2**[@love2014moderated]: offers functionalities for visualization and differential expression analysis of RNA-Seq data.  

- **goseq** [@young2010gene]: implements the GSA method GOSeq.

- **PADOG** [@tarca2013comparison]: implements the GSA method PADOG.



