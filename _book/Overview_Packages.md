---
title: "Overview and description of all required libraries"
author: "Milena Wünsch"
date: "2023-07-17"
output: html_document
---



To run the entire pipeline in this illustration, from pre-filtering to the actual GSA tool, several R packages are required. We here provide an overview of all packages that must be installed (if not already available) and loaded. In the sections that correspond to the individual steps, we give an overview of the locally required libraries. 

### Overview of all libraries required throughout the illustrations 

- tweeDEseqCountData: provides count data from three RNA-Seq experiments together with information such as the conditions of the samples. We work work the gene expression data set provided by Pickrell et al. (2010) in our illustrations. 

- BiocManager: offers functionalities to install and manage the packages that are a part of the Bioconductor project. Here, we particularly use it to install the required Bioconductor packages. 

- clusterProfiler: offers functionalities to explore the functional characteristics of high-throughput genomic data. Here, we use its functions to perform ORA and GSEA and for the conversion of the gene IDs. 

- org.Hs.eg.db: provides the genome-wide annotation for the human organism. Note that if your gene expression measurements are extracted from a different organism, you must adapt this library accordingly. For instance, the organism mouse requires the library org.Mm.eg.db. \
This library is often required when working with the library clusterProfiler. 

- limma: provides functionalities for analysing data from gene expression experiments. limma allows the analysis of microarray as well as RNA-Seq data for which the analysis pipelines are similar. We here use limma for differential expression analysis and make use of a transformation of the RNA-Seq data to align its characteristics for microarray data. 

- edgeR: provides functionalities to analyse data from RNA sequencing technologies and focuses on differential expression analysis with biological (as opposed to technical) replicates. 

- DESeq2: offers functionalities for visualization and differential expression analysis of RNA-Seq data.  

- dplyr: provides various functionalities for data manipulation. \
We use it here to unify the results tables of the different methods for differential expression analysis. While this is not necessary in general, it simplifies our illustrations as we can use the same GSA pipelines following differential expression analysis, independent of the method for differential expression analysis used.

- goseq: implements the GSA method GOSeq.

- PADOG: implements the GSA method PADOG.
