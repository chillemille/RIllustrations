################################################################################
### Instructions: clusterProfiler's GSEA  for gene set database GO #############
################################################################################

# empty environment
rm(list =ls())


################################################################################
### Content of this script #####################################################
################################################################################

# In this script, we will do the following two things:
# 1. Based on the results of differential expression analysis, we will go through
#    all steps required to run clusterProfiler's GSEA
# 2. We will go through all (meaningful) researchers' degrees of freedom

################################################################################



# set working directory: File in which all resulting data sets are stored
setwd("/nfsmb/koll/milena.wuensch/Dokumente/GSA_Review/Data")

# load required libraries
library(clusterProfiler)


# In this script, we want to perform clusterProfiler's GSEA
# We want generate the required input object based on the results of differential
# expression analysis

# here, for differential expression analysis we want to utilize the results of the
# parametric methods voom/limma, DESeq2 and edgeR:

# IMPORTANT: the gene ID format(s) accepted by clusterProfiler depends on the chosen gene set database:
# (i) for KEGG, clusterProfiler accepts the gene IDs in the NCBI (Entrez) gene ID format
# (ii) for GO, clusterProfiler accepts the gene IDs in a variety of formats , such as
#      Ensembl gene ID format, Entrez gene ID format and HGNC gene symbols
#     -> for the list of available gene IDs, check
keytypes(org.Hs.eg.db)
# for other organisms, adjust argument accordingly

# since here, we choose gene set database GO: proceed with Ensembl ID

load("./Results_Differential_Expression_Analysis/DE_results_limma_Ensembl.Rdata")
load("./Results_Differential_Expression_Analysis/DE_results_DESeq2_Ensembl.Rdata")
load("./Results_Differential_Expression_Analysis/DE_results_edgeR_Ensembl.Rdata")

# for DESEq2, we additionally need to load the object dds_Ensembl from
# "Instructions_Differential_Expression_Analysis.R":
load("./Results_Differential_Expression_Analysis/dds_Ensembl.Rdata")



# for the purpose of simplicity: shorten term for each result
DE_results_limma <- DE_results_limma_Ensembl
DE_results_DESeq2 <- DE_results_DESeq2_Ensembl
DE_results_edgeR <- DE_results_edgeR_Ensembl



  ################################################################################
  ### step 1: Generation of Required Input Object ################################
  ################################################################################

  # In the following, we illustrate the creation of the required input object
  # using the three different methods for differential expression analysis which
  # shall serve as three options

  # required input object for DESeq2: "order ranked geneList",
  # i.e. a vector with a ranking value for each gene, named with
  # the respective gene IDs and ordered in a descending manner


  ################################################################################
  ### 1.1 Option 1: Generate required input using limma/voom #####################
  ################################################################################

  # load library
  library(limma)


  #########
  # step 1: Generate the ranked list of the genes based on the results of differential
  ######### expression analysis

  # formula for gene-level ranking metric: -1 * log10(p-value) * sign(log fold change)

  # note: by p-value, we mean the non-adjusted p-value

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



  ################################################################################
  ### 1.2 Option 2: Generate required input using DESeq2 #########################
  ################################################################################

  # load library
  library(DESeq2)

  #########
  # step 1: Perform Differential Expression Analysis
  #########

  # this step is performed in the script Instructions_Differential_Expression_Analysis.R
  # note: we must exit the DESeq2 workflow from the script a few steps early since there
  # is an additional step required in for the creation of the required input object

  # we proceed with the object dds_Ensembl we've obtained in step 3 of DESeq2's workflow
  # for differential expression analysis

  dds_Ensembl

  #########
  # step 2: Shrink estimated log fold change values
  #########

  # intuition behind shrinkage of log fold change values:

  # As mentioned in the paper, RNA-Seq data consists (in its rare form) of count
  # data in which is inherently heteroscedastic, i.e. the variance of the count
  # data depends on the mean count of the count data

  # It is observable that ratios between counts are considerably noisier in
  # low magnitudes of counts compared to higher magnitudes, i.e.
  # the log fold changes between both conditions are higher if the overall
  # magnitude of counts is lower.

  # DESeq2 addresses this issue by shrinking the estimated log fold changes in
  # the direction of 0
  # the magnitude of shrinkage is higher if the available information for a gene
  # is lower (which may be because of a low magnitude of counts, a high dispersion
  # or few degrees of freedom.)
  # a more detailed description is provided in the DESeq2 paper by Love et al. (2014)

  # shrinkage:
  DE_results_DESeq2_shrink_Ensembl <- lfcShrink(dds_Ensembl,
                                               coef = "condition_treated_vs_untreated",
                                               type="apeglm")


  # function arguments:

  # type: method to perform shrinkage
  # we opt for the default "apeglm" but you can choose from two alternative options
  # as well

  # coef: indicate the coefficients to be shrunk
  # we can obtain the right argument from the following function call:

  resultsNames(dds_Ensembl)

  # this shows us that we can either shrink the intercept or the "condition_treated_vs_untreated"
  # since we do not want to shrink the intercept but the log fold changes, we opt for
  # the second option "condition_treated_vs_untreated"

  # transform results table to data frame
  DE_results_DESeq2_shrink_Ensembl <- as.data.frame(DE_results_DESeq2_shrink_Ensembl)


  # inspect results table:
  DE_results_DESeq2_shrink_Ensembl

  # note :
  # adjusted and non-adjusted p-values remain the same between DE_results_DESeq2_shrink
  # and DE_results_DESeq2
  # we can see that compared to DE_results_DESeq2, the values in column log2FoldChange in
  # DE_results_DESeq2_shrink have become smaller (on the absolute range)

  # note: we have not performed shrinkage in the script Instructions_Differential_Expression_Analysis.R
  # as there, the goal was to get the list of differentially expressed genes which is independent
  # of the log fold change values


  #########
  # step 2: Generate the ranked list of the genes based on the results of differential
  ######### expression analysis

  # formula for gene-level ranking metric: -1 * log10(p-value) * sign(log fold change)

  # note: by p-value, we mean the non-adjusted p-value

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



  ################################################################################
  ### 1.3 Option 3: Generate required input using edgeR ##########################
  ################################################################################

  # load library
  library(edgeR)


  #########
  # step 1: Generate the ranked list of the genes based on the results of differential
  ######### expression analysis

  # formula for gene-level ranking metric: -1 * log10(p-value) * sign(log fold change)

  # note: by p-value, we mean the non-adjusted p-value

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

  # 5. special problem here: gene ENSG00000129824 has a ranking value Inf since its adjusted
  # p-value in the results table of differential expression analysis amounts to 0

  # here, we deal with this issue by resetting this ranking value to the highest ranking value
  # that occurs among the remaining genes
  # -> note that there is NO common way of dealing with this issue
  rankvec_edgeR_Ensembl[rankvec_edgeR_Ensembl == Inf] <- max(rankvec_edgeR_Ensembl[rankvec_edgeR_Ensembl != Inf])



  ################################################################################
  ### step 2: Run GSEA with gene set database GO  ##############################
  ################################################################################

  # for the purpose of simplicity, we work with the ranking generated using limma here
  # however: you can switch around at your discretion

  rankvec_Ensembl <- rankvec_limma_Ensembl


  # load Genome wide annotation for Human (this is a required argument for function gseGO())
  library(org.Hs.eg.db)
  # note: for a genome that is not human, a corresponding annotation package must be installed and loaded
  # e.g. for mouse:
  # library(org.Mm.eg.db)
  # -> the argument in gseGO() must then be set to OrgDb = org.Mm.eg.db

  # important: set seed for reprodubibility
  set.seed(1)

  GSEA_GO <- gseGO(geneList = rankvec_Ensembl,
                   ont = "BP",
                   OrgDb = org.Hs.eg.db,
                   keyType = "ENSEMBL",
                   seed = TRUE) # set seed for reproducibility

  # arguments:

  # geneList: vector of gene ranking (generated in manner of step 1)

  # ont: subontology of GO (must be one of "BP", "CC", "MF")

  # OrgDb: indicates organism from which the gene expression measurements are taken

  # keyTypes: for our purposes here, argument keyType indicates the gene ID format of
  # the gene names in vector rankvec
  # available keytypes can be found in the following:
  keytypes(org.Hs.eg.db)

  # seed = TRUE: set the seed you have indicated above (here: seed = 1 )


  # additional arguments (these do not appear in the above code lines since we keep the
  # default options)

  # - exponent 1: in the calculation of the enrichment score, each gene is weighted
  # by its absolute value of the ranking metric

  # eps = 1e-10: each p-value that is smaller than 1e-10 is indicated as 1e-10

  # pvalueCutoff = 0.05: only those gene sets with a p-value < 0.05 are indicated in the
  # results
  # - note: if you want to inspect ALL gene sets in the results, set pvalueCutoff = 1

  # pAdjustMethod = "BH": adjustment for multiple testing using Benjamini and Hochberg
  # method



  ################################################################################
  ### Interpretation of Results ##################################################
  ################################################################################


  # inspect results via
  # View(as.data.frame(GSEA_GO))

  # Columns in Result Tables:

  # ID: ID of gene set
  # - Description: Description of Gene Set
  #
  # - setSize: Size of the gene set
  #
  # - enrichmentScore: Enrichment Score
  # -> note: this enrichment score has not been normalized for gene set size
  # -> means: larger gene sets automatically have a bigger (absolute) enrichment score,
  # ->        independent of actual differential enrichment
  # -> the raw enrichment score is therefore not comparable between different gene sets
  #
  # - NES: Normalized version of column enrichmentScore
  # -> NES can be compared between different gene sets
  #
  # - pvalue: p-value of enrichment of given gene set
  # -> note: this raw p-value has not been adjusted for multiple testing
  # -> therefore: CANNOT be used to assess differential enrichment of a given gene set
  #
  # - p.adjust: ADJUSTED p-value of a given gene set
  # -> this p-value now has been adjusted for multiple testing and can therefore
  # -> can therefore be used to assess differential enrichment
  # -> example: detect all genes with p.adjust < 0.05 as differentially enriched
  #
  # - qvalue: ADJUSTED p-value of a given gene set
  # -> note: a qvalue is the analog to p.adjust, but adjusted for multiple testing
  #         using a different method
  #
  # ->> you can either use the column p.adjust or qvalue to assess whether a gene set
  #     is differentially enriched or not
  #
  # - rank: position in the ranked list of the genes in which the maximum difference
  #       between the two sums occurs, i.e. the rank at which the enrichment score
  #       is extracted
  #
  # - leading_edge:
  # * tags: The percentage of genes  in the ranked list that are members of gene set before (for positive enrichment score)
  #         or after (for negative enrichment score) the position from which the enrichment score is extracted.
  # * list:  The percentage of genes in the ranked gene list before (for positive enrichment
  #           score) or after (for negative enrichment score) the position from which the enrichment score is extracted.
  # * signal: enrichment signal strength
  #           combines the statistics "tags" and "list"
  #
  # - core_enrichment: genes that contribute most to the enrichment results

  # note: more detailed information on leading_edge and core_enrichment can be found in the user
  # manual provided for GSEA:
  # https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideFrame.html
  # -> search for Detailed Enrichment Results (section "Leading Edge")
  # -> search for CORE ENRICHMENT








################################################################################
# Researchers's degrees of freedom #############################################
################################################################################

# IMPORTANT NOTE: the intention behind going through the researchers' degrees of freedom
# is to give you an understanding of what you can do to adapt the given (parameter)
# setting to the research question

# ->> MORE IMPORTANTLY: the intention behind going through these flexible parameters is
# NOT to change them in order to help you obtain the most preferable results by systematically
# changing these parameters
# -> such behaviour would correspond to "cherry-picking"

# Any changes in the parameter choice should be documented transparently


# Since we've already covered how to choose between the gene set databases KEGG
# and GO, we will now cover how to change the value of the exponent


###############################
# change 1: change exponent ###
###############################

# note that, in contrast to the web-based tool GSEA, clusterProfiler does not make
# any indication as to which alternative exponent values are allowed or sensible

# one could stick with the web-based tool GSEA which suggests, in addition to the
# default exponent value of 1, the alternative values 0, 1.5, and 2

# here, we change the exponent value to 2


  # important: set seed for reproducibility
  set.seed(1)

  GSEA_GO_exponent <- gseGO(geneList = rankvec_Ensembl,
                            ont = "BP",
                            OrgDb = org.Hs.eg.db,
                            keyType = "ENSEMBL",
                            exp = 2,
                            seed = TRUE)











