################################################################################
### Instructions: Differential Expression Analysis #############################
################################################################################

# set working directory: File in which all resulting data sets are stored 
setwd("/nfsmb/koll/milena.wuensch/Dokumente/GSA_Review/Data")

# in the case that the GSA method of choice requires a different (e.g. Entrez gene ID format):
# obtain the pre-filtered gene expression data sets that result from gene ID conversion and 
# duplicate gene ID removal
load("./Results_GeneID_Conversion_DuplicatesRemoval/exprdat_filter_conv_DESeq2.Rdata")
load("./Results_GeneID_Conversion_DuplicatesRemoval/exprdat_filter_conv_filterByExpr.Rdata")

# in the case that the GSA method of choice accepts the gene IDs in the Ensembl gene ID format, 
# obtain the pre-filtered gene expression data sets with initial (Ensembl) gene IDs 
load("./Results_PreFiltering/expression_data_filterByExpr.Rdata")
load("./Results_PreFiltering/expression_data_filterDESeq2.Rdata")


# For illustrative purposes, we use a dataset from the following library:
library(tweeDEseqCountData) 
# We load the following data set
data(pickrell)

# the gene expression data set can be addressed with the following syntax:
# Biobase::exprs(pickrell.eset)
# the conditions/phenotypes of the samples can be addressed with the following syntax: 
pickrell.eset$gender


# for the purpose of simplicity (of syntax), we store the count data
# and the conditions of the samples of the gene expression data set in new objects: 

# conditions of the samples
sample_conditions <- pickrell.eset$gender


# Depending on the method for differential expression analyis, different methods of pre-filtering are proposed 

# voom/limma and egdeR: pre-filtering using edgeR's builtin function filterByExpr()
# DESeq2: manual/simple pre-filtering
# NOISeq: performs pre-filtering internally by default

# for this reason, we have generated separate gene expression data sets for 
# (i) limma/voom and edgeR: exprdat_filter_conv_filterByExpr
# for DESeq2: exprdat_filter_conv_filterByExpr
# for NOISeq: 





################################################################################
### (I) Differential Expression Analysis with Entrez gene IDs ##################
################################################################################




################################################################################
### Instructions: limma ########################################################
################################################################################

library(limma)
library(edgeR) # required since limma utilizes a variety of functionalities from the edgeR

# code illustrations are based on the following user manual: 
# https://bioconductor.org/packages/devel/bioc/vignettes/limma/inst/doc/usersguide.pdf

########################################
### step 1: Generate required object ###
########################################

# As described in "Instructions_Differential_Expression_Analysis.R", exprdat_filter_conv_filterByExpr is
# a data frame and needs to be converted to a DGEList object 

# store expression data with corresponding sample conditions in object of class DGEList

 y <- DGEList(counts = exprdat_filter_conv_filterByExpr, 
              group = sample_conditions)



# required arguments: 

# counts: matrix that contains RNA-Seq data (i.e. count data)

# optional arguments:

# group: indicates condition of each sample/library 
# -> this argument must be specified sooner or later (such as in subsequent functions) so we just specify it 
#    at this point 

# note: we leave the remaining arguments in their default state since the corresponding info will be added 
# through the following pipeline of functions



##############################
### step 2: Normalization  ###
##############################

# the following piece of code generates a normalization factor for each sample
# this accounts for sample-specific effect (such as differences in library sizes and effects of compositionality)
# -> if not accounted for, these effects prevent a comparison between the samples 

y <- calcNormFactors(y)

# note: this function does not transform the count data, but rather generates a normalization factor for each sample
# which is incorporated into the subsequent analysis separately 


###################################
### step 3: voom transformation ###
###################################

# (i) design matrix (rows correspond to samples, columns indicate which coefficients are to be estimated)
design_matrix <- model.matrix(~sample_conditions)

# (ii) voom transformation (transformation of count data to log-cpm values, computation of observation weight for each 
# gene based on mean-variance relationship)
y <- voom(y, design = design_matrix)

# note: voom-transformation facilitates the use of the following functions that were initially developed for microarray
# data 



################################################
### step 4: differential expression analysis ###
################################################


# (i) fit linear model for each gene
y <- lmFit(y)

# (ii) calculate statistics for the assessment of differential expression 
y <- eBayes(y)

# (iii) Get Results table for each gene whose differential expression was assessed 
DE_results_limma_Entrez <- topTable(y, number =  nrow(y))
# number = nrow(y) ensures that all genes are displayed in results



################################################################################
### step 5: Rename Columns in Results of Differential Expression Analysis ######
################################################################################

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
DE_results_limma_Entrez <- as.data.frame(DE_results_limma_Entrez)

# rename column that contains adjusted p-values: rename padj to p_adj
DE_results_limma_Entrez <- rename(DE_results_limma_Entrez, p_adj = adj.P.Val)









################################################################################
### Instructions: DESeq2 #######################################################
################################################################################

# the official DESeq2 vignette can be found through the following link: 
# http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

# load library
library(DESeq2)

################################################################################
### step 1: Format required for DESeq2 #########################################
################################################################################

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
                      row.names = colnames(exprdat_filter_conv_DESeq2))

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

dds_Entrez<-DESeqDataSetFromMatrix(countData = exprdat_filter_conv_DESeq2, 
                            colData = coldata, 
                            design = ~condition)


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


################################################################################
### step 2: Differential Expression Analysis ###################################
################################################################################



# (I) perform default differential expression analysis 
dds_Entrez <- DESeq(dds_Entrez)
# (II) generate results table which provides
DE_results_DESeq2_Entrez <- results(dds_Entrez)


# more detailed explanations: 
# (I) DESeq(): Estimation of normalization factors, estimation of dispersions, fitting of generalized
# linear model, Wald statistics

# (II) results(): base mean across all samples, estimated log2 fold changes, standard errors, test statistics,
# p-values, adjusted p-values 



################################################################################
### step 3: Rename Columns in Results of Differential Expression Analysis ######
################################################################################

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
DE_results_DESeq2_Entrez <- as.data.frame(DE_results_DESeq2_Entrez)

# rename column that contains adjusted p-values: rename padj to p_adj
DE_results_DESeq2_Entrez <- rename(DE_results_DESeq2_Entrez, p_adj = padj)





################################################################################
### Instructions: edgeR ########################################################
################################################################################

# the instructions for edgeR are based on the following vignette:
# https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf

# load library
library(edgeR)


########################################
### step 1: Generate required object ###
########################################

# store expression data with corresponding sample conditions in object of class DGEList

y <- DGEList(counts = exprdat_filter_conv_filterByExpr, 
             group = sample_conditions)



# required arguments: 

# counts: matrix that contains RNA-Seq data (i.e. count data)

# optional arguments:

# group: indicates condition of each sample/library 
# -> this argument must be specified sooner or later (such as in subsequent functions) so we just specify it 
#    at this point 

# note: we leave the remaining arguments in their default state since the corresponding info will be added 
# through the following pipeline of functions




##############################
### step 2: Normalization  ###
##############################

# the following piece of code generates a normalization factor for each sample
# this accounts for sample-specific effect (such as differences in library sizes and effects of compositionality)
# -> if not accounted for, these effects prevent a comparison between the samples 

y <- calcNormFactors(y)

# note: this function does not transform the count data, but rather generates a normalization factor for each sample
# which is incorporated into the subsequent analysis separately 




########################################
### step 3: Estimation of Dispersion ###
########################################

# estimates common and tagwise dispersion: variationof the true abundance of a given gene between 
# different samples 
# -> required to assess differential expression realistically 

y <- estimateDisp(y)


################################################
### step 4: Test for Differential Expression ###
################################################

# test each gene for differential expression: 
DE_results_edgeR_Entrez <- exactTest(y)

# extract pre-specified (n) number of genes 
DE_results_edgeR_Entrez <- topTags(DE_results_edgeR_Entrez, n = nrow(DE_results_edgeR_Entrez))

# argument n specifies the number of top differentially expressed genes to be displayed in the results
# n = nrow(DE_results_Entrez) ensures the results of all genes whose differential expression was assessed 
# are displayed






################################################################################
### step 5: Rename Columns in Results of Differential Expression Analysis ######
################################################################################

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
DE_results_edgeR_Entrez <- as.data.frame(DE_results_edgeR_Entrez)

# rename column that contains adjusted p-values: rename padj to p_adj
DE_results_edgeR_Entrez <- rename(DE_results_edgeR_Entrez, p_adj = FDR)









################################################################################
### Instructions: NOISeq #######################################################
################################################################################

# # load NOISeq library
# library(NOISeq)
# 
# # NOTE: - NOISeq is considerably less popular compared to voom/limma, DESeq2 and 
# #         edgeR
# #       - In contrast to these three methods, NOISeq is non-parametric which means
# #         that it does not make any distributional assumptions on the underlying 
# #         data structure
# 
# # We want to illustrate its application for the purpose of completeness: 
# 
# 
# ########################################
# ### step 1: Prepare Sample Conditions ##
# ########################################
# 
# # NOISeq requires the sample conditions to be stored in a data frame: 
# 
# sample_conditions_df <- data.frame(condition = sample_conditions) 
# 
# 
# ############################################
# ### step 2: Prepare Required Input Object ##
# ############################################
# 
# 
# expression_data_prep <- readData(data= expression_data, 
#                                  factors = sample_conditions_df)
# 
# 
# # arguments: 
# # - data: matrix/data frame that contains the counts 
# # - factors: data frame that contains the sample conditions (see step 1)
# 
# 
# ######################################
# ### step 3: Run NOISeq (NOISeqbio) ###
# ######################################
# 
# # Important note: NOISeqbio was specifically developed for the case when there 
# #                 are biological replicates (as opposed to technical replicates)
# 
# DE_results_Entrez_NOISeq <- noiseqbio(input = mydata, 
#                                factor = "condition")
# 
# 
# # arguments: 
# # - input: object resulting from step 2 (i.e. from function readData())
# # - factor: indicator of column name of data frame sample_conditions_df which 
# #           contains the conditions to be compared 
# # -> note: in step 1, we set the column name of sample_conditions_prep to "condition"
# # -> Correspondingly, in step 3, we assign argument factor the character value "condition"
# 
# 
# # additional arguments: 
# # - k: value by which zero counts are replaced (k = 0.5 means that all values in the 
# #      count data which are equal to 0 are replaced by 0.5)
# # 
# # - norm: method for normalization (default "rpkm", i.e. "reads per kilobase per million")
# # -> alternatives are: * norm = "tmm", i.e. Trimmed Mean of M-Values
# #                      * norm = "uqua", i.e. upper quartile normalization 
# #
# # - conditions: only needed IF two conditions are to be compared by in the process of 
# #               differential expression and IF factor contains more than two conditions
# # 
# # - random.seed: ensures that the results of random sampling are identical with each run,
# #                leading to identical DE results with each run
# # -> by default, the argument is set to random.seed = 12345
# # 
# # - filter: pre-filtering method to filter out lowly expressed genes prior to differential
# #           expression analysis
# # -> options: more information see ?filtered.data
# #            * filter = 0: No pre-filtering
# #            * filter = 1: CPM method for pre-filtering
# #            * filter = 2: Wilcoxon Test for pre-filtering (not recommended if there are 
# #                          less than 5 samples per condition)
# #            * filter = 3: Proportion Test 
# #
# 
# 
# #######################################################
# ### step 3: Retrieve Differentially Expressed Genes ###
# #######################################################
# 
# DEGs_NOISeq <- degenes(object = DE_results_Entrez_NOISeq, 
#                        q = 0.95, 
#                        M = NULL)
# 
# 
# # arguments: 
# # - object: Results of Differential Expression Analysis generated in step 3
# # - q: threshold of probability of differential expression 
# # -> all genes with a probability of differential expression > q are displayed
# #    in the results 
# # - M: indicates whether only up-regulated (M = "up"), down-regulated (M = "down")  
# #      or all genes (M = NULL) with a probability of differential expression > q 
# #      are to be displayed in the results
# 
# 
# #############################################
# ### step 4: Interpretation of the Results ###
# #############################################
# 
# # Differentially expressed genes are selected based on column prob
# # -> based on argument q from step 3, those genes with prob > 0.95
# #    are selected as differentially expressed 
# 
# # Important note: the probability of differential expression prob does NOT correspond
# # to 1 - pvalue, but rather to 1-FDR, where FDR can be considered an adjusted p-value
# # -> this justifies the choice of q = 0.05 in step 3, since in limma/voom, DESeq2, 
# #    edgeR, ... typically those genes with an adjusted p-value < 0.05 are detected 
# #    as differentially expressed 
# 
# 
# # additional columns: 
# # - female_mean: 
# # - male_mean
# 
# # -> note: "female" in female_mean and "male" in male_mean are the levels of the factor in 
# #          sample_conditions_df
# # -> for a different gene expression data sets, the column names will be different
# 
# 
# 
# 
# rowMeans(expression_data[rownames(expression_data) == "9087", pickrell.eset$gender == "female"])
# 











################################################################################
### Save Results ###############################################################
################################################################################

# Save results of differential expression analysis for subsequent analyses 

save(DE_results_limma_Entrez, file = "./Results_Differential_Expression_Analysis/DE_results_limma_Entrez.Rdata")
save(DE_results_DESeq2_Entrez, file = "./Results_Differential_Expression_Analysis/DE_results_DESeq2_Entrez.Rdata")
save(DE_results_edgeR_Entrez, file = "./Results_Differential_Expression_Analysis/DE_results_edgeR_Entrez.Rdata")


# note: we need the object dds (from the DESeq2 workflow) for clusterProfiler's GSEA 
# -> there, we will perform an additional shrinkage of the estimated log fola changes,
# which required the object dds 

save(dds_Entrez, file = "./Results_Differential_Expression_Analysis/dds_Entrez.Rdata")











################################################################################
### (II) Differential Expression Analysis with Ensembl gene IDs ################
################################################################################




################################################################################
### Instructions: limma ########################################################
################################################################################

library(limma)
library(edgeR) # required since limma utilizes a variety of functionalities from the edgeR

# code illustrations are based on the following user manual: 
# https://bioconductor.org/packages/devel/bioc/vignettes/limma/inst/doc/usersguide.pdf

########################################
### step 1: Generate required object ###
########################################

# As described in "Instructions_Differential_Expression_Analysis.R", exprdat_filter_conv_filterByExpr is
# a data frame and needs to be converted to a DGEList object 

# store expression data with corresponding sample conditions in object of class DGEList

y <- DGEList(counts = expression_data_filterByExpr, 
             group = sample_conditions)



# required arguments: 

# counts: matrix that contains RNA-Seq data (i.e. count data)

# optional arguments:

# group: indicates condition of each sample/library 
# -> this argument must be specified sooner or later (such as in subsequent functions) so we just specify it 
#    at this point 

# note: we leave the remaining arguments in their default state since the corresponding info will be added 
# through the following pipeline of functions



##############################
### step 2: Normalization  ###
##############################

# the following piece of code generates a normalization factor for each sample
# this accounts for sample-specific effect (such as differences in library sizes and effects of compositionality)
# -> if not accounted for, these effects prevent a comparison between the samples 

y <- calcNormFactors(y)

# note: this function does not transform the count data, but rather generates a normalization factor for each sample
# which is incorporated into the subsequent analysis separately 


###################################
### step 3: voom transformation ###
###################################

# (i) design matrix (rows correspond to samples, columns indicate which coefficients are to be estimated)
design_matrix <- model.matrix(~sample_conditions)

# (ii) voom transformation (transformation of count data to log-cpm values, computation of observation weight for each 
# gene based on mean-variance relationship)
y <- voom(y, design = design_matrix)

# note: voom-transformation facilitates the use of the following functions that were initially developed for microarray
# data 



################################################
### step 4: differential expression analysis ###
################################################


# (i) fit linear model for each gene
y <- lmFit(y)

# (ii) calculate statistics for the assessment of differential expression 
y <- eBayes(y)

# (iii) Get Results table for each gene whose differential expression was assessed 
DE_results_limma_Ensembl <- topTable(y, number =  nrow(y))
# number = nrow(y) ensures that all genes are displayed in results



################################################################################
### step 5: Rename Columns in Results of Differential Expression Analysis ######
################################################################################

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
DE_results_limma_Ensembl <- as.data.frame(DE_results_limma_Ensembl)

# rename column that contains adjusted p-values: rename padj to p_adj
DE_results_limma_Ensembl <- dplyr::rename(DE_results_limma_Ensembl, p_adj = adj.P.Val)









################################################################################
### Instructions: DESeq2 #######################################################
################################################################################

# the official DESeq2 vignette can be found through the following link: 
# http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

# load library
library(DESeq2)

################################################################################
### step 1: Format required for DESeq2 #########################################
################################################################################

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
                      row.names = colnames(expression_data_filterDESeq2))

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

dds_Ensembl<-DESeqDataSetFromMatrix(countData = expression_data_filterDESeq2, 
                            colData = coldata, 
                            design = ~condition)


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


################################################################################
### step 2: Differential Expression Analysis ###################################
################################################################################



# (I) perform default differential expression analysis 
dds_Ensembl <- DESeq(dds_Ensembl)
# (II) generate results table which provides
DE_results_DESeq2_Ensembl <- results(dds_Ensembl)


# more detailed explanations: 
# (I) DESeq(): Estimation of normalization factors, estimation of dispersions, fitting of generalized
# linear model, Wald statistics

# (II) results(): base mean across all samples, estimated log2 fold changes, standard errors, test statistics,
# p-values, adjusted p-values 



################################################################################
### step 3: Rename Columns in Results of Differential Expression Analysis ######
################################################################################

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
DE_results_DESeq2_Ensembl <- as.data.frame(DE_results_DESeq2_Ensembl)

# rename column that contains adjusted p-values: rename padj to p_adj
DE_results_DESeq2_Ensembl <- dplyr::rename(DE_results_DESeq2_Ensembl, p_adj = padj)





################################################################################
### Instructions: edgeR ########################################################
################################################################################

# the instructions for edgeR are based on the following vignette:
# https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf

# load library
library(edgeR)


########################################
### step 1: Generate required object ###
########################################

# store expression data with corresponding sample conditions in object of class DGEList

y <- DGEList(counts = expression_data_filterByExpr, 
             group = sample_conditions)



# required arguments: 

# counts: matrix that contains RNA-Seq data (i.e. count data)

# optional arguments:

# group: indicates condition of each sample/library 
# -> this argument must be specified sooner or later (such as in subsequent functions) so we just specify it 
#    at this point 

# note: we leave the remaining arguments in their default state since the corresponding info will be added 
# through the following pipeline of functions




##############################
### step 2: Normalization  ###
##############################

# the following piece of code generates a normalization factor for each sample
# this accounts for sample-specific effect (such as differences in library sizes and effects of compositionality)
# -> if not accounted for, these effects prevent a comparison between the samples 

y <- calcNormFactors(y)

# note: this function does not transform the count data, but rather generates a normalization factor for each sample
# which is incorporated into the subsequent analysis separately 




########################################
### step 3: Estimation of Dispersion ###
########################################

# estimates common and tagwise dispersion: variationof the true abundance of a given gene between 
# different samples 
# -> required to assess differential expression realistically 

y <- estimateDisp(y)


################################################
### step 4: Test for Differential Expression ###
################################################

# test each gene for differential expression: 
DE_results_edgeR_Ensembl <- exactTest(y)

# extract pre-specified (n) number of genes 
DE_results_edgeR_Ensembl <- topTags(DE_results_edgeR_Ensembl, n = nrow(DE_results_edgeR_Ensembl))

# argument n specifies the number of top differentially expressed genes to be displayed in the results
# n = nrow(DE_results_Entrez) ensures the results of all genes whose differential expression was assessed 
# are displayed






################################################################################
### step 5: Rename Columns in Results of Differential Expression Analysis ######
################################################################################

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
DE_results_edgeR_Ensembl <- as.data.frame(DE_results_edgeR_Ensembl)

# rename column that contains adjusted p-values: rename padj to p_adj
DE_results_edgeR_Ensembl <- dplyr::rename(DE_results_edgeR_Ensembl, p_adj = FDR)









################################################################################
### Instructions: NOISeq #######################################################
################################################################################

# # load NOISeq library
# library(NOISeq)
# 
# # NOTE: - NOISeq is considerably less popular compared to voom/limma, DESeq2 and 
# #         edgeR
# #       - In contrast to these three methods, NOISeq is non-parametric which means
# #         that it does not make any distributional assumptions on the underlying 
# #         data structure
# 
# # We want to illustrate its application for the purpose of completeness: 
# 
# 
# ########################################
# ### step 1: Prepare Sample Conditions ##
# ########################################
# 
# # NOISeq requires the sample conditions to be stored in a data frame: 
# 
# sample_conditions_df <- data.frame(condition = sample_conditions) 
# 
# 
# ############################################
# ### step 2: Prepare Required Input Object ##
# ############################################
# 
# 
# expression_data_prep <- readData(data= expression_data, 
#                                  factors = sample_conditions_df)
# 
# 
# # arguments: 
# # - data: matrix/data frame that contains the counts 
# # - factors: data frame that contains the sample conditions (see step 1)
# 
# 
# ######################################
# ### step 3: Run NOISeq (NOISeqbio) ###
# ######################################
# 
# # Important note: NOISeqbio was specifically developed for the case when there 
# #                 are biological replicates (as opposed to technical replicates)
# 
# DE_results_Entrez_NOISeq <- noiseqbio(input = mydata, 
#                                       factor = "condition")
# 
# 
# # arguments: 
# # - input: object resulting from step 2 (i.e. from function readData())
# # - factor: indicator of column name of data frame sample_conditions_df which 
# #           contains the conditions to be compared 
# # -> note: in step 1, we set the column name of sample_conditions_prep to "condition"
# # -> Correspondingly, in step 3, we assign argument factor the character value "condition"
# 
# 
# # additional arguments: 
# # - k: value by which zero counts are replaced (k = 0.5 means that all values in the 
# #      count data which are equal to 0 are replaced by 0.5)
# # 
# # - norm: method for normalization (default "rpkm", i.e. "reads per kilobase per million")
# # -> alternatives are: * norm = "tmm", i.e. Trimmed Mean of M-Values
# #                      * norm = "uqua", i.e. upper quartile normalization 
# #
# # - conditions: only needed IF two conditions are to be compared by in the process of 
# #               differential expression and IF factor contains more than two conditions
# # 
# # - random.seed: ensures that the results of random sampling are identical with each run,
# #                leading to identical DE results with each run
# # -> by default, the argument is set to random.seed = 12345
# # 
# # - filter: pre-filtering method to filter out lowly expressed genes prior to differential
# #           expression analysis
# # -> options: more information see ?filtered.data
# #            * filter = 0: No pre-filtering
# #            * filter = 1: CPM method for pre-filtering
# #            * filter = 2: Wilcoxon Test for pre-filtering (not recommended if there are 
# #                          less than 5 samples per condition)
# #            * filter = 3: Proportion Test 
# #
# 
# 
# #######################################################
# ### step 3: Retrieve Differentially Expressed Genes ###
# #######################################################
# 
# DEGs_NOISeq <- degenes(object = DE_results_Entrez_NOISeq, 
#                        q = 0.95, 
#                        M = NULL)
# 
# 
# # arguments: 
# # - object: Results of Differential Expression Analysis generated in step 3
# # - q: threshold of probability of differential expression 
# # -> all genes with a probability of differential expression > q are displayed
# #    in the results 
# # - M: indicates whether only up-regulated (M = "up"), down-regulated (M = "down")  
# #      or all genes (M = NULL) with a probability of differential expression > q 
# #      are to be displayed in the results


#############################################
### step 4: Interpretation of the Results ###
#############################################

# Differentially expressed genes are selected based on column prob
# -> based on argument q from step 3, those genes with prob > 0.95
#    are selected as differentially expressed 

# Important note: the probability of differential expression prob does NOT correspond
# to 1 - pvalue, but rather to 1-FDR, where FDR can be considered an adjusted p-value
# -> this justifies the choice of q = 0.05 in step 3, since in limma/voom, DESeq2, 
#    edgeR, ... typically those genes with an adjusted p-value < 0.05 are detected 
#    as differentially expressed 


# additional columns: 
# - female_mean: 
# - male_mean

# -> note: "female" in female_mean and "male" in male_mean are the levels of the factor in 
#          sample_conditions_df
# -> for a different gene expression data sets, the column names will be different




#rowMeans(expression_data[rownames(expression_data) == "9087", pickrell.eset$gender == "female"])












################################################################################
### Save Results ###############################################################
################################################################################

# Save results of differential expression analysis for subsequent analyses 


save(DE_results_limma_Ensembl, file = "./Results_Differential_Expression_Analysis/DE_results_limma_Ensembl.Rdata")
save(DE_results_DESeq2_Ensembl, file = "./Results_Differential_Expression_Analysis/DE_results_DESeq2_Ensembl.Rdata")
save(DE_results_edgeR_Ensembl, file = "./Results_Differential_Expression_Analysis/DE_results_edgeR_Ensembl.Rdata")


# note: we need the object dds (from the DESeq2 workflow) for clusterProfiler's GSEA 
# -> there, we will perform an additional shrinkage of the estimated log fola changes,
# which required the object dds 

save(dds_Ensembl, file = "./Results_Differential_Expression_Analysis/dds_Ensembl.Rdata")




