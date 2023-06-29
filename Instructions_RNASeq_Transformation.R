################################################################################
### Instructions: Transformation (and Normalization) of RNA-Seq data ###########
################################################################################

# empty environment
rm(list =ls())


################################################################################
### Content of this script #####################################################
################################################################################

# Transformation of the RNA-Seq data to align its characteristics with those of #
# microarray measurements using two approaches: 
# approach 1: transformation using voom 
# approach 2: transformation using DESeq2's varianceStabilizingTransformation

# we repeat this procedure for two gene ID formats:
# (i) Entrez gene IDs 
# (ii) Ensembl gene IDs 

################################################################################



# set working directory: File in which all resulting data sets are stored 
setwd("/nfsmb/koll/milena.wuensch/Dokumente/GSA_Review/Data")

# note: we assume that at this point, the gene expression data has already been 

# 1. pre-filtered (if applicable/necessary)
# 2. converted in terms of gene IDs (if necessary)

# for the purpose of simplicity, we here work with the pickrell example data set provided by
library(tweeDEseqCountData)


# note: we generate a new objects
# (i) expression_data in which we store the pickrell data set (WITH CONVERTED GENE IDs) AND 
# (ii) phenotype labels in which we store the phenotypes/conditions of the samples
# for the purpose of illustration



# BUT: if you have (an already pre-filtered and/or converted; if applicable) gene 
# expression data set (along with the sample conditions), you can easily substitute 
# both objects with your at hand

# as the expression data set we choose the data set resulting from the first manner of 
# duplicate gene ID removal and manual pre-filtering proposed by the authors of DESeq2 

data(pickrell)
sample_conditions <- pickrell.eset$gender



# Note: PADOG and GSEA, which are the two tools in this paper that require 
# a (manual) transformation of the RNA-Seq data, require the genes in a different
# format: 
# (i) PADOG: Entrez gene ID -> need gene expression data with converted gene IDs
# (ii) GSEA: Ensembl gene ID -> need gene expression data with initial gene IDs 

# -> that's why: we perform a transformation for BOTH gene expression data sets: 

# for (i): load PRE-FILTERED gene expression data set with CONVERTED gene IDs 
# here, we work with the approach to pre-filtering proposed by DESeq2, 
# but you can also swith to pre-filtering using filterByExpr()
load("./Results_GeneID_Conversion_DuplicatesRemoval/exprdat_filter_conv_DESeq2.Rdata")

# for (ii): Load PRE-FILTERED gene expression data set with INITIAL gene IDs 
# here, we work with the approach to pre-filtering proposed by DESeq2, 
# but you can also swith to pre-filtering using filterByExpr()
load("./Results_PreFiltering/expression_data_filterDESeq2.Rdata")


# In the following, we present two approaches to RNA-Seq transformation
# We apply both approaches to the gene expression data set with 
# (i) converted (Entrez) gene IDs (for PADOG)
expression_data_filt_conv <- exprdat_filter_conv_DESeq2
# (ii) initial (Ensembl) gene IDs (for GSEA)
expression_data_filt <- expression_data_filterDESeq2



################################################################################
################################################################################
### (i) Transformation of gene expression data set with converted (Entrez) IDs #
################################################################################
################################################################################




################################################################################
### Approach 1: Transformation of RNA-Seq data using voom (limma) ##############
################################################################################

# load required libraries: 
library(limma)
# voom makes use of certain functionalities provided by edgeR: 
library(edgeR)





  ########
  # step 1: generate DGEList object from the expression data 
  ########

  counts_voom <-DGEList(expression_data_filt_conv, group=sample_conditions)

  #######
  # step 2: perform normalization (calculate normalization/scaling factor)
  #######

  # note: we here use the default normalization method TMM
  # other choices for normalization can be made using argument method
  # in function calcNormFactors()

  counts_voom <- calcNormFactors(counts_voom)

  ########
  # step 3: perform voom function (do NOT make use of precision weights)
  ########

  expression_data_voomtransformed <- voom(counts_voom)

  ########
  # step 4: convert data set (from EList class) to data frame 
  ########

  expression_data_voomtransformed_Entrez <- as.data.frame(expression_data_voomtransformed)



###########################################################################################################
### Approach 2: Transformation of RNA-Seq data using DESeq2's Variance Stabilizing Transformation #########
###########################################################################################################

# load DESeq2 library
library(DESeq2)

  ########
  # step 1: generate a data frame which contains the condition of each sample: 
  ########

  # note: this corresponds to the format required by DESeq2


  # the names of the samples are stored as the row names of the data frame
  # -> important: make sure that the order of the conditions in sample_conditions corresponds to the order of the 
  # samples in expression_data
  coldata <- data.frame(sample_conditions, 
                        row.names = colnames(expression_data_filt_conv))
  
  # rename the column header to "condition" 
  colnames(coldata) <- "condition"
  
  # recode the variable condition as a factor 
  # rename the sample conditions (in this case from "female" and "male") to "untreated" and "treated"
  # note: make sure that the control level in variable condition is coded as the first level (i.e. "untreated")
  coldata$condition <- factor(coldata$condition, 
                              labels = c("untreated","treated"))
  
  
  ########
  # step 2: generate DESeq2 data set
  ########
  
  
  dds<-DESeqDataSetFromMatrix(countData = expression_data_filt_conv, 
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
  
  
  
  ########
  # step 3: perform DESeq2's varianceStabilizingInformation
  ########
  
    
        # transform expression data set  
        expression_data_vsttransformed <-vst(dds)
    
        # since the data set is now in the format of a DESeqTransform
        # data set, the transformed count data are not directly accessible 
        # use function assay() which lets us access the count data:
  
        expression_data_vsttransformed <- assay(expression_data_vsttransformed)
  
    
  
  ########
  # step 4: convert data set to data frame 
  ########
  
  expression_data_vsttransformed_Entrez <- as.data.frame(expression_data_vsttransformed)
        
        
        
        
        
################################################################################
################################################################################
### (ii) Transformation of gene expression data set with initial (Ensembl) IDs #
################################################################################
################################################################################
        
        
  ################################################################################
  ### Approach 1: Transformation of RNA-Seq data using voom (limma) ##############
  ################################################################################
        
  # load required libraries: 
  library(limma)
  # voom makes use of certain functionalities provided by edgeR: 
  library(edgeR)
        
        
        
        
        
  ########
  # step 1: generate DGEList object from the expression data 
  ########
        
  counts_voom <-DGEList(expression_data_filt, group=sample_conditions)
        
  #######
  # step 2: perform normalization (calculate normalization/scaling factor)
  #######
        
  # note: we here use the default normalization method TMM
  # other choices for normalization can be made using argument method
  # in function calcNormFactors()
        
  counts_voom <- calcNormFactors(counts_voom)
        
  ########
  # step 3: perform voom function (do NOT make use of precision weights)
  ########
        
  expression_data_voomtransformed <- voom(counts_voom)
        
  ########
  # step 4: convert data set (from EList class) to data frame 
  ########
        
  expression_data_voomtransformed_Ensembl <- as.data.frame(expression_data_voomtransformed)
        
        
        
  ###########################################################################################################
  ### Approach 2: Transformation of RNA-Seq data using DESeq2's Variance Stabilizing Transformation #########
  ###########################################################################################################
        
  # load DESeq2 library
  library(DESeq2)
        
  ########
  # step 1: generate a data frame which contains the condition of each sample: 
  ########
        
  # note: this corresponds to the format required by DESeq2
        
        
  # the names of the samples are stored as the row names of the data frame
  # -> important: make sure that the order of the conditions in sample_conditions corresponds to the order of the 
  # samples in expression_data
  coldata <- data.frame(sample_conditions, 
            row.names = colnames(expression_data_filt_conv))
        
  # rename the column header to "condition" 
  colnames(coldata) <- "condition"
        
  # recode the variable condition as a factor 
  # rename the sample conditions (in this case from "female" and "male") to "untreated" and "treated"
  # note: make sure that the control level in variable condition is coded as the first level (i.e. "untreated")
  coldata$condition <- factor(coldata$condition, 
                              labels = c("untreated","treated"))
        
  ########
  # step 2: generate DESeq2 data set
  ########
        
        
  dds<-DESeqDataSetFromMatrix(countData = expression_data_filt_conv, 
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
        
        
        
  ########
  # step 3: perform DESeq2's varianceStabilizingInformation
  ########
        
  # transform expression data set  
  expression_data_vsttransformed <-vst(dds)
  
  # since the data set is now in the format of a DESeqTransform
  # data set, the transformed count data are not directly accessible 
  # use function assay() which lets us access the count data:
        
  expression_data_vsttransformed <- assay(expression_data_vsttransformed)
        
        
        
  ########
  # step 4: convert data set to data frame 
  ########
        
  expression_data_vsttransformed_Ensembl <- as.data.frame(expression_data_vsttransformed)

        
        
        
################################################################################
### Save Results ###############################################################
################################################################################
        
# save voom-transformed RNA-Seq data 
 
save(expression_data_voomtransformed_Entrez, file = "./Results_RNASeq_Transformation/expression_data_voomtransformed_Entrez.Rdata")     
save(expression_data_voomtransformed_Ensembl, file = "./Results_RNASeq_Transformation/expression_data_voomtransformed_Ensembl.Rdata")        

# save variance transformed RNA-Seq data 
save(expression_data_vsttransformed_Entrez, file = "./Results_RNASeq_Transformation/expression_data_vsttransformed_Entrez.Rdata")     
save(expression_data_vsttransformed_Ensembl, file = "./Results_RNASeq_Transformation/expression_data_vsttransformed_Ensembl.Rdata")        


    






  









