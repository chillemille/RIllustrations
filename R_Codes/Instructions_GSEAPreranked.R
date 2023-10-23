################################################################################
### Instructions: GSEAPreranked  ###############################################
################################################################################

# empty environment
rm(list =ls())

# set working directory
# you will have to adapt the path within setwd() to your computer
setwd("/nfsmb/koll/milena.wuensch/Dokumente/GSA_Review")

# create folder to eventually contain the objects generated in this script
dir.create("./data/Input_Objects_GSEAPreranked")

################################################################################
### Content of this script #####################################################
################################################################################

# In this script, we will go through the process of preparing the required input
# for GSEAPreranked (which is part of the web-based tool GSEA)

# Note that in contrast to the remaining tools considered in this work, it is recom-
# mended for GSEAPreranked to provide the input with the genes ID HGNC (HUGO) gene
# symbols. The reason for this is described in the supplement of this work.

# Therefore, we start with the pre-filtered gene expression data set with the
# gene IDs in the initial (Ensembl) ID format

# we proceed as follows:
# (i) conversion of the gene IDs to HGNC (HUGO) symbols and removal of resulting
#     duplicated gene IDs
# (ii) differential expression analysis using voom/limma, DESeq2, and edgeR
# (iii) generation of the gene ranking from each results table generated in (ii)
# (iv) export of gene rankings to text file
# (v) link to information on how to further process the gene rankings in Excel

################################################################################

# we will proceed with the gene expression data set which was pre-filtered using
# edgeR's filterByExpr
load("./data/Results_PreFiltering/expression_data_filterByExpr.Rdata")
# alternatively, you can also load the gene expression data set
# which was manually pre-filtered (as proposed by DESeq2):
# load("./Results_PreFiltering/expression_data_filterDESeq2.Rdata")

# for the purpose of simplicity, we will proceed with the pre-filtered gene
# expression data set under a different name:
expression_data_filtered <- expression_data_filterByExpr



################################################################################
### step 1: Convert Ensembl IDs to HGNC symbols ################################
################################################################################

# Here, we work with the function provided by clusterProfiler.

# load libraries
library(clusterProfiler)
library(org.Hs.eg.db) # for human

# obtain mapping of initial (Ensembl) gene IDs to required HGNC gene symbols:
bitr_EnsToSymb <- bitr(rownames(expression_data_filtered), # the rownames of the gene expression data set correspond to the gene IDs
                       fromType = "ENSEMBL",
                       toType = "SYMBOL",
                       OrgDb = org.Hs.eg.db)



# concatenate initial gene expression data set with the mapping from initial (Ensembl) to required HGNC symbols
# merge by row names of expression data set and ENSEMBL ID of gene ID mapping
expression_data_merged <- merge(x = expression_data_filtered, # initial (pre-filtered) gene expression data set
                                             y = bitr_EnsToSymb, # mapping from HGNC gene symbols to HGNC symbols
                                             by.x = 0, # merge by the rownames (i.e. Ensembl IDs) in gene expression data set
                                             by.y = "ENSEMBL", # merge by Ensembl IDs in mapping
                                             sort=TRUE)



################################################################################
### step 2: Remove duplicated gene IDs  ########################################
################################################################################

# As we can see from the document "Instructions_GeneID_Conversion_DuplicatesRemoval.R",
# we work with three separate ways to remove the duplicate gene IDs that result from
# gene ID conversion.

# We will apply all three of these approaches to removal

###########################################
# (1) Take a closer look at the duplicates
###########################################

########
# CASE 1: single ENSEMBL IDs are mapped to multiple HGNC symbols
########

# obtain number of cases in which an ENSEMBL gene ID was converted to several HGNC symbols, i.e. the number of
# times an Ensembl ID appears more than once in the mapping
sum(duplicated(bitr_EnsToSymb$ENSEMBL))
# determine all duplicated ENSEMBL gene IDS (i.e. all ENSEMBL gene IDs that were mapped to multiple distinct HGNC gene symbols):
dupl_ensembl<-unique(bitr_EnsToSymb$ENSEMBL[duplicated(bitr_EnsToSymb$ENSEMBL)])
# number of ENSEMBL IDs that have at least one duplicate
length(dupl_ensembl)



# in the following, we can inspect the conversion pattern for each Ensembl ID that was mapped to
# multiple distinct HGNC symbols
duplicated_conversion_ens<-bitr_EnsToSymb[bitr_EnsToSymb$ENSEMBL %in% dupl_ensembl,]
# take a look at conversion pattern for duplicated Ensembl gene IDs:
duplicated_conversion_ens
# -> we can see perfectly that these Ensembl ID are each mapped to multiple individual HGNC symbols


########
# CASE 2: multiple ENSEMBL IDs are mapped to the same HGNC symbol
########

# obtain number of cases in which multiple distinct Ensembl IDs are converted to the same HGNC gene symbol
# in these cases, the corresponding HGNC symbols appear repeatedly in the mapping:
sum(duplicated(bitr_EnsToSymb$SYMBOL))
# determine all HGNC gene symbols affected by this duplication
dupl_symbol<-unique(bitr_EnsToSymb$SYMBOL[duplicated(bitr_EnsToSymb$SYMBOL)])
# number of HGNC symbols that have at least one duplicate
length(dupl_symbol)
# display of conversion pattern of duplicated HGNC symbols
duplicated_conversion<-bitr_EnsToSymb[bitr_EnsToSymb$SYMBOL %in% dupl_symbol,]
# take a look at conversion pattern:
duplicated_conversion
# if ordered by HGNC symbols, we can see perfectly that for each of these HGNC gene symbols, there are multiple
# corresponding Ensembl gene IDs
dim(duplicated_conversion)


##################################
# (ii) remove duplicated gene IDs
##################################


# for subsequent analysis (whether differential expression analysis of gene set analysis), we need to deal
# with the duplicated gene IDs and remove them in a suitable way such that only one unique gene expression
# measurement among the duplicates remains

# there is no recommended way to proceed, i.e. no common approach presented in an official scientific puplication,
# but instead several approaches suggested by users in corresponding user platforms,
# generally, the manner of duplicate gene ID removal is performed at the discretion of the user

# we therefore present three approaches to duplicate gene ID removal



############
### option 1: keep first subscript among duplicates
############

# this is the simplest approach among the three


#1. remove duplicated HGNC gene symbols
exprdat_symbol_dupl1<-expression_data_merged[!duplicated(expression_data_merged$SYMBOL),]


#2. remove duplicated ENSEMBL gene IDs
exprdat_symbol_dupl1<-exprdat_symbol_dupl1[!duplicated(exprdat_symbol_dupl1$Row.names),]
dim(exprdat_symbol_dupl1)

#3. set HGNC gene symbols as rownames
rownames(exprdat_symbol_dupl1)<-exprdat_symbol_dupl1$SYMBOL
#Remove columns containing ENSEMBL and HGNC symbols
exprdat_symbol_dupl1<-subset(exprdat_symbol_dupl1, select=-c(Row.names,SYMBOL))
dim(exprdat_symbol_dupl1)
# we particularly see that the number of columns is now at its initial number


#############################################
### final converted gene expression data set:
#############################################

exprdat_symbol_dupl1


#############################################
#############################################




############
### option 2: keep (rounded) mean expression value of all duplicated gene IDs
############


# here, we switch order and remove duplicated HGNC gene symbols (case 2) before removing
# duplicated Ensembl gene IDs (case 1)
# the reason for this is elaborated below

#1.remove duplicated HGNC gene symbols (case 2)
#i.e. multiple different ENSEMBL IDs that are mapped to the same single HGNC symbol

#generate matrix to contain (rounded) mean expression values of all rows that
#have same HGNC gene symbol
#ncol=ncol(expression_data_filterByExpr)-2 since data set contains 2 columns with IDs at this point
mean_symbol<-matrix(, nrow=0, ncol=ncol(expression_data_merged)-2)


# for each duplicated HGNC gene symbol separately, we gather all rows with the corresponding gene expression data
# and then extract the (rounded) mean expression value of all rows
for(i in 1:length(dupl_symbol)){

  #go through each HGNC symbols which occurs multiple times
  #determine all rows whose HGNC symbols correspond to currently considered HGNC symbol
  counts_dupl<-expression_data_merged[expression_data_merged$SYMBOL %in% unique(dupl_symbol)[i],]

  #compute the mean expression value of all rows that contain to the given HGNC gene symbol
  dupl_id<-round(colMeans(counts_dupl[,c(2:(ncol(expression_data_merged)-1))]))
  #store rounded mean expression value in matrix
  # this matrix is extended by a single row of gene expression data which corresponds to the
  # (rounded) mean expression data that corresponds to the given HGNC gene symbol
  mean_symbol<-rbind(mean_symbol,dupl_id)


}

# after completing the for-loop, mean_symbol contains the mean expression measurements of each HGNC gene symbol
# which contains duplicates resulting from gene ID conversion
mean_symbol

#set corresponding HGNC gene symbols as rownames
rownames(mean_symbol)<-unique(dupl_symbol)


# test whether the number of rows in mean_symbol corresponds to the number HGNC symbols
# that occur more than once
# result should be TRUE
nrow(mean_symbol)==length(dupl_symbol)

# remove all rows from the expression data whose HGNC symbol has at least one duplicate
# intuition: we have just dealt with the corresponding rows and do not want them to be considered
# in the second step (which deals with case 2

exprdat_symbol_dupl2<-expression_data_merged[!expression_data_merged$SYMBOL %in% dupl_symbol,]

# test whether number of rows in resulting data set equals nrow of inital data set
# minus number of genes with at least one duplicate
nrow(exprdat_symbol_dupl2)==nrow(expression_data_merged)-nrow(duplicated_conversion)
dim(exprdat_symbol_dupl2)



#2. remove duplicated ENSEMBL IDs
#caution: single ENSEMBL IDs that are mapped to multiple HGNC symbol naturally generate
#identical count data for all corresponding HGNC symbols
#->pointless to compute mean expression values
#verifiable by looking at data set only containing those ENSEMBL IDs that are
#mapped by multiple HGNC symbols:
#test_dupl_ensembl<-expression_data_filterByExpr[expression_data_filterByExpr$Row.names %in% dupl_ensembl,]
#View(test_dupl_ensembl)

#therefore: proceed as in option 1 and use HGNC symbol that occurs first, remove the rest
exprdat_symbol_dupl2<-exprdat_symbol_dupl2[!duplicated(exprdat_symbol_dupl2$Row.names),]
dim(exprdat_symbol_dupl2)
#set HGNC symbol as rownames
rownames(exprdat_symbol_dupl2)<-exprdat_symbol_dupl2$SYMBOL
#remove any columns containing IDs
exprdat_symbol_dupl2<-subset(exprdat_symbol_dupl2,select= -c(Row.names,SYMBOL))
#add rows to data set that contain mean expression values of duplicate HGNC symbols
exprdat_symbol_dupl2<-rbind(exprdat_symbol_dupl2,mean_symbol)
#dimension of remaining expression data set:
#dim(exprdat_symbol_dupl2)


#############################################
### final converted gene expression data set:
#############################################

exprdat_symbol_dupl2



#############################################
############################################




############
### option 3: among duplicates, keep row with highest overall expression values (i.e highest counts across all samples)
############


#intuition: row with highest counts values has  highest power of detecting
#differential expression
#as in option 2, this applies only to duplicates that result from multiple ENSEMBL IDs
#that are mapped to the same HGNC symbol


#case 2: (case 1 below) multiple ENSEMBL IDs that are converted to the same single HGNC symbol

# generate matrix to later contain row with highest count values among ID duplicates
# this data set is to be filled gradually and with each iteration of the follwing for-loop
highest_count_symbol<-matrix(, nrow=0, ncol=ncol(expression_data_filterByExpr))




# for each duplicated HGNC gene symbol separately, we gather all rows with the corresponding gene expression data
# and then extract the row with the highest overall magnitude of counts

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
  #View(highest_count_symbol)
  #remove rows in counts_dupl from count data set successively
}


#Remove all initial values with entrez duplicates from the dataset initial gene expression data set
exprdat_symbol_dupl3<-expression_data_merged[! expression_data_merged$SYMBOL %in% unique(dupl_symbol),]


#as in option 2, pointless to detect row with highest count values as all rows
#corresponding to the same ENSEMBL ID naturally contain identical count data
#therefore: remove duplicate ENSEMBL ID that occurs first
exprdat_symbol_dupl3<-exprdat_symbol_dupl3[!duplicated(exprdat_symbol_dupl3$Row.names),]

# add gene expression rows of all HGNC gene symbols that were initially duplicated
exprdat_symbol_dupl3<-rbind(exprdat_symbol_dupl3,highest_count_symbol )

#Set HGNC symbols as rownames
rownames(exprdat_symbol_dupl3)<-exprdat_symbol_dupl3$SYMBOL
#Remove any column that contains information on gene IDs
exprdat_symbol_dupl3<-subset(exprdat_symbol_dupl3, select=-c(Row.names,SYMBOL))

dim(exprdat_symbol_dupl3)
# we see that the sample size (number of columns) is now back at its initial number


#############################################
### final converted gene expression data set:
#############################################

exprdat_symbol_dupl3




################################################################################
# ! intermediate step !: Choose which converted gene expression data set ###########
# to proceed with ##############################################################

# for the purpose of simplicity, we will proceed with the first approach to
# the removal of duplicated gene IDs
# however, you can easily switch to another approach at your discretion


expression_data_filt_symbol <- exprdat_symbol_dupl1




################################################################################
### step 3: Differential expression analysis with HGNC gene symbols ############
################################################################################





# We proceed with differential expression analysis
# from the corresponding results table, we later transform certain metrics
# into a gene-level statistic for each gene based on which the gene ranking
# is generated

# We will proceed analogously to the script "Instructions_Differential_Expression_Analysis.R",
# in which we perform differential expression analysis using three established methods,
# namely
# (i) voom/limma
# (ii) DESeq2
# (iii) edgeR

# note: we already have the count data at hand that we will be working with
# we additionally need the conditions of the samples
# reminder: in this illustration, we work with the pickrell data set

# load library from whith we then load the pickrell data set
library(tweeDEseqCountData)
data(pickrell)

sample_conditions <- pickrell.eset$gender

################################################################################
### Instructions: limma ########################################################
################################################################################

library(limma)
library(edgeR) # required since limma utilizes a variety of functionalities from the edgeR

# code illustrations are based on the following user manual:
# https://bioconductor.org/packages/devel/bioc/vignettes/limma/inst/doc/usersguide.pdf

###################################
### 1. Generate required object ###
###################################

# As described in "Instructions_Differential_Expression_Analysis.R", expression_data_filt_symbol is
# a data frame and needs to be converted to a DGEList object

# store expression data with corresponding sample conditions in object of class DGEList

y <- DGEList(counts = expression_data_filt_symbol,
             group = sample_conditions)



# required arguments:

# counts: matrix that contains RNA-Seq data (i.e. count data)

# optional arguments:

# group: indicates condition of each sample/library
# -> this argument must be specified sooner or later (such as in subsequent functions) so we just specify it
#    at this point

# note: we leave the remaining arguments in their default state since the corresponding info will be added
# through the following pipeline of functions



#########################
### 2. Normalization  ###
#########################

# the following piece of code generates a normalization factor for each sample
# this accounts for sample-specific effects (such as differences in library sizes and effects of compositionality)
# -> if not accounted for, these effects prevent a comparison between the samples

y <- calcNormFactors(y)

# note: this function does not transform the count data, but rather generates a normalization factor for each sample
# which is incorporated into the subsequent analysis separately


##############################
### 3. voom transformation ###
##############################

# (i) design matrix (rows correspond to samples, columns indicate which coefficients are to be estimated)
design_matrix <- model.matrix(~sample_conditions)

# (ii) voom transformation (transformation of count data to log-cpm values, computation of observation weight for each
# gene based on mean-variance relationship)
y <- voom(y, design = design_matrix)

# note: voom-transformation facilitates the use of the following functions that were initially developed for microarray
# data



###########################################
### 4. differential expression analysis ###
###########################################


# (i) fit linear model for each gene
y <- lmFit(y)

# (ii) calculate statistics for the assessment of differential expression
y <- eBayes(y)

# (iii) Get Results table for each gene whose differential expression was assessed
DE_results_limma <- topTable(y, number =  nrow(y))
# number = nrow(y) ensures that all genes are displayed in results



###########################################################################
### 5. Rename Columns in Results of Differential Expression Analysis ######
###########################################################################

library(dplyr) # contains function rename()

# this step is not required for differential expression analysis or for the subsequent use
# of gene set analysis IN GENERAL

# a renaming of some columns (such as adjusted p-value) is necessary in this context
# since the different methods for differential expression analysis typically differ
# in the column names of the results tables. A unification of the column names is
# required so that we can use the same code to illustrate further conduct of GSA in the
# following R scripts

# first, we transform to results table to a data frame so that we see the results
# table directly when accessing it through the name "res"
DE_results_limma <- as.data.frame(DE_results_limma)

# rename column that contains adjusted p-values: rename padj to p_adj
DE_results_limma <- rename(DE_results_limma, p_adj = adj.P.Val)









################################################################################
### Instructions: DESeq2 #######################################################
################################################################################

# the official DESeq2 vignette can be found through the following link:
# http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

# load library
library(DESeq2)

###########################################################################
### 1. Format required for DESeq2 #########################################
###########################################################################

# DESeq2 operates on the format DESeqDataSet which contains information on the count data, the conditions
# of the samples and the design (for further information see below)

# note: for gene expression data sets that must be imported to R, additional steps are necessary before
# the following code can be run

#########################################################################
# (I) generate a data frame which contains the condition of each sample:
#########################################################################


# here: the information on the sample conditions is the only column in the data frame
# however: if further variables (such as batch effects) are to be controlled for, the corresponding
# variables must additionally be added to coldata

# the names of the samples are stored as the row names of the data frame
# -> important: make sure that the order of the conditions in sample_conditions corresponds to the order of the
# samples in expression_data
coldata <- data.frame(sample_conditions,
                      row.names = colnames(expression_data_filt_symbol))

# rename the column header to "condition"
colnames(coldata) <- "condition"

# recode the variable condition as a factor
# rename the sample conditions (in this case from "female" and "male") to "untreated" and "treated"
# note: make sure that the control level in variable condition is coded as the first level (i.e. "untreated")
coldata$condition <- factor(coldata$condition,
                            labels = c("untreated","treated"))

#############################
# (II) generate DESeqDataSet
#############################

dds <- DESeqDataSetFromMatrix(countData = expression_data_filt_symbol,
                            colData = coldata,
                            design = ~ condition)


########################################################
# relevant arguments in function DESeqDataSetFromMatrix:
########################################################

# (i) countData: count data from gene expression data set

# (ii) colData: data frame that contains information on the samples (see above)
# -> conditions of the samples (required) and possibly further variables to correct for (such as batch effects)

# (iii) design: indicates which variables from colData are used for modelling
# -> more detailed: the argument design is used to estimate the dispersions and the log2 fold changes of the model
# -> if more than one variable from colData are used in argument design (e.g. a second variable "batch"), the syntax
# changes to the following formula: design ~ batch + condition
# ->> make sure that the variable of interest (here: variable that represents conditions of the samples)
# is placed at the end of the formula


###########################################################################
### 2. Differential Expression Analysis ###################################
###########################################################################



# (I) perform default differential expression analysis
dds <- DESeq(dds)
# (II) generate results table which provides
DE_results_DESeq2 <- results(dds)


# more detailed explanations:
# (I) DESeq(): Estimation of normalization factors, estimation of dispersions, fitting of generalized
# linear model, Wald statistics

# (II) results(): base mean across all samples, estimated log2 fold changes, standard errors, test statistics,
# p-values, adjusted p-values



###########################################################################
### 3. Rename Columns in Results of Differential Expression Analysis ######
###########################################################################

library(dplyr) # contains function rename()

# this step is not required for differential expression analysis or for the subsequent use
# of gene set analysis

# a renaming of some columns (such as adjusted p-value) is necessary in this context
# since the different methods for differential expression analysis typically differ
# in the column names of the results tables. A unification of the column names is
# required so that we can use the same code to illustrate further conduct of GSA in the
# following R scripts

# first, we transform to results table to a data frame so that we see the results
# table directly when accessing it through the name "res"
DE_results_DESeq2 <- as.data.frame(DE_results_DESeq2)

# rename column that contains adjusted p-values: rename padj to p_adj
DE_results_DESeq2 <- rename(DE_results_DESeq2, p_adj = padj)





################################################################################
### Instructions: edgeR ########################################################
################################################################################

# the instructions for edgeR are based on the following vignette:
# https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf

# load library
library(edgeR)


###################################
### 1. Generate required object ###
###################################

# store expression data with corresponding sample conditions in object of class DGEList

y <- DGEList(counts = expression_data_filt_symbol,
             group = sample_conditions)



# required arguments:

# counts: matrix that contains RNA-Seq data (i.e. count data)

# optional arguments:

# group: indicates condition of each sample/library
# -> this argument must be specified sooner or later (such as in subsequent functions) so we just specify it
#    at this point

# note: we leave the remaining arguments in their default state since the corresponding info will be added
# through the following pipeline of functions




#########################
### 2. Normalization  ###
#########################

# the following piece of code generates a normalization factor for each sample
# this accounts for sample-specific effect (such as differences in library sizes and effects of compositionality)
# -> if not accounted for, these effects prevent a comparison between the samples

y <- calcNormFactors(y)

# note: this function does not transform the count data, but rather generates a normalization factor for each sample
# which is incorporated into the subsequent analysis separately




###################################
### 3. Estimation of Dispersion ###
###################################

# estimates common and tagwise dispersion: variationof the true abundance of a given gene between
# different samples
# -> required to assess differential expression realistically

y <- estimateDisp(y)


###########################################
### 4. Test for Differential Expression ###
###########################################

# test each gene for differential expression:
DE_results_edgeR <- exactTest(y)

# extract pre-specified (n) number of genes
DE_results_edgeR <- topTags(DE_results_edgeR, n = nrow(DE_results_edgeR))

# argument n specifies the number of top differentially expressed genes to be displayed in the results
# n = nrow(DE_results) ensures the results of all genes whose differential expression was assessed
# are displayed






########################################################################
### 5. Rename Columns in Results of Differential Expression Analysis ###
########################################################################

library(dplyr) # contains function rename()

# this step is not required for differential expression analysis or for the subsequent use
# of gene set analysis

# a renaming of some columns (such as adjusted p-value) is necessary in this context
# since the different methods for differential expression analysis typically differ
# in the column names of the results tables. A unification of the column names is
# required so that we can use the same code to illustrate further conduct of GSA in the
# following R scripts

# first, we transform to results table to a data frame so that we see the results
# table directly when accessing it through the name "res"
DE_results_edgeR <- as.data.frame(DE_results_edgeR)

# rename column that contains adjusted p-values: rename padj to p_adj
DE_results_edgeR <- rename(DE_results_edgeR, p_adj = FDR)











################################################################################
### step 4: Generate gene ranking ##############################################
################################################################################

# we will now apply a gene-level statistic to the results table of differential
# expression analysis to generate a gene ranking which will serve as input to
# GSEAPreranked

# as we have performed differential expression analysis with
# (i) voom/limma
# (ii) DESeq2
# (iii) edgeR,
# we will also generate three gene rankings, one for each of these three methods


# formula for gene-level ranking metric: -1 * log10(p-value) * sign(log fold change)

# note: by p-value, we mean the non-adjusted p-value



################################################################################
### 4.1 Option 1: Generate required input using limma/voom #####################
################################################################################


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



################################################################################
### 4.2 Option 2: Generate required input using DESeq2 #########################
################################################################################

# load library
library(DESeq2)


#########
# step 1: Shrink estimated log fold change values
#########

# Before generating the gene ranking from the DESeq2 results, we need to perform an additional
# shrinkage that is specific to DESeq2.

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

# The shrinkage is performed on object dds, which is a result of the function DESeq()

# shrinkage:
DE_results_DESeq2_shrink <- lfcShrink(dds,
                                      coef = "condition_treated_vs_untreated",
                                      type="apeglm")


# function arguments:

# type: method to perform shrinkage
# we opt for the default "apeglm" but you can choose from two alternative options
# as well

# coef: indicate the coefficients to be shrunk
# we can obtain the right argument from the following function call:

resultsNames(dds)

# this shows us that we can either shrink the intercept or the "condition_treated_vs_untreated"
# since we do not want to shrink the intercept but the log fold changes, we opt for
# the second option "condition_treated_vs_untreated"

# transform results table to data frame
DE_results_DESeq2_shrink <- as.data.frame(DE_results_DESeq2_shrink)


# inspect results table:
DE_results_DESeq2_shrink

# note :
# adjusted and non-adjusted p-values remain the same between DE_results_DESeq2_shrink
# and DE_results_DESeq2
# we can see that compared to DE_results_DESeq2, the values in column log2FoldChange in
# DE_results_DESeq2_shrink have become smaller (on the absolute range)

# note: we have not performed shrinkage in the script Instructions_Differential_Expression_Analysis.R
# as there, the goal was to get the list of differentially expressed genes which is independent
# of the log fold change values


#########
# step 2: Generate the gene ranking

# formula for gene-level ranking metric: -1 * log10(p-value) * sign(log fold change)

# note: by p-value, we mean the non-adjusted p-value

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



################################################################################
### 4.3 Option 3: Generate required input using edgeR ##########################
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
ind_nonNA_pvalue <- !is.na(DE_results_edgeR$PValue)

# subset gene expression data set to those genes with a p-value
DE_results_noNA <- DE_results_edgeR[ind_nonNA_pvalue, ]

# 2. create vector that contains the value of the gene-level ranking metric for each gene
rankvec_edgeR <- sign(DE_results_noNA$logFC)*(-1)*log10(DE_results_noNA$p_adj)

# 3. assign respective gene ID to each value in the vector
names(rankvec_edgeR) <- rownames(DE_results_noNA)

# 4. sort the vector in descending order
rankvec_edgeR <- sort(rankvec_edgeR, decreasing=TRUE)

# 5. special problem here: gene ENSG00000129824 has a ranking value Inf since its adjusted
# p-value in the results table of differential expression analysis amounts to 0

# here, we deal with this issue by resetting this ranking value to the highest ranking value
# that occurs among the remaining genes
# -> note that there is NO common way of dealing with this issue
rankvec_edgeR[rankvec_edgeR == Inf] <- max(rankvec_edgeR[rankvec_edgeR != Inf])



################################################################################
### step 5: export gene ranking to text file  ##################################
################################################################################

# for the purpose of simplicity, we only export one of the gene rankings to
# a text file

# here, we export the ranking generated with limma


# this path indicates that we store the gene ranking as the text file
# "gene_ranking.txt" in folder "Input_Objects_GSEAPreranked"
path_conditions <- "./data/Input_Objects_GSEAPreranked/gene_ranking.txt"

write.table(x = rankvec_limma,
            file = path_conditions,
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)


# note:
# - quote = FALSE ensures that none of the characters (in this case gene and sample
# identifiers) are surrounded by double quotes
# - row.names = TRUE ensures that the gene IDs are included in the export
# - col.names = TRUE ensures that no gene information is in the first row since
#               the web-based tool ignores whatever is in the first row


################################################################################
### step 6: Further preparation of input object in Excel #######################
################################################################################

# For this step, follow the instructions on the very bottom of the following
# link:
# http://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats

# -> section "RNK: Ranked list file format (*.rnk)"





################################################################################
### step 7: Clean-up ###########################################################
################################################################################

# remove all objects from the environment apart from the relevant final objects:

rm(list=ls()[!ls() %in% c("rankvec_limma", "rankvec_DESeq2", "rankvec_edgeR")])


















