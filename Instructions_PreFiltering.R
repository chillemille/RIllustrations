################################################################################
### Instructions: Pre-Filtering  ###############################################
################################################################################

# set working directory: File in which all resulting data sets are stored 
setwd("/nfsmb/koll/milena.wuensch/Dokumente/GSA_Review/Data")


# intuition: exclude all genes from subsequent analysis that do not have a sufficient number of read counts across 
# all samples

# load gene expression data set for illustration purposes 
# the following library contains the gene expression data: 
library(tweeDEseqCountData)
# load gene expression data: 
data(pickrell)

# for illustration purposes, we create the objects expression_data and sample_conditions which contain the gene 
# expression measurements and conditions of the samples

# gene expression measurements 
expression_data <- Biobase::exprs(pickrell.eset)
# sample conditions: 
sample_conditions <- pickrell.eset$gender 




################################################################################
### option 1: edgeR's builtin function #########################################
################################################################################

library(edgeR)


# note: this approach to pre-filtering is used for voom/limma in addition to edgeR


########################################
### step 1: Generate required object ###
########################################

# store expression data with corresponding sample conditions in object of class DGEList

expression_data_filterByExpr <- DGEList(counts = expression_data, 
             group = sample_conditions)



# required arguments: 

# counts: matrix that contains RNA-Seq data (i.e. count data)

# optional arguments:

# group: indicates condition of each sample/library 
# -> this argument must be specified sooner or later (such as in subsequent functions) so we just specify it 
#    at this point 

# note: we leave the remaining arguments in their default state since the corresponding info will be added 
# through the following pipeline of functions

########################################
### step 2: Pre-Filtering ##############
########################################

# egdeR's builtin function filterByExpr() operates on the cpm-transformed count data (cpm: counts-per-million), 
# which excludes all genes that do NOT have a certain number of counts-per-million in a certain number of samples

# note: filtering based on cpm-transformed values accounts for differences in library sizes between the samples 

# (i) for each gene, indicate for each gene if it contains sufficient cpm in enough samples
keep <- filterByExpr(expression_data_filterByExpr)

# (ii) filter the gene expression data set such that only those genes are kept which fulfill the requirements
expression_data_filterByExpr <- expression_data_filterByExpr[keep,, keep.lib.sizes = FALSE]


###################################################
### Final pre-filtered gene expression data set ###
###################################################

# note that at this point, expression_data_filterByExpr is still a DGEList object 
# however: at this point, we want to proceed with the count data as a simple data frame

# the reason for this is that we do not only use the pre-filtered gene expression data set(s)
# for differential expression analysis, in which we cannot work with a DGEList object 
# -> in the course of differential expression analysis, we will then make it a DGEList object again 

expression_data_filterByExpr <- as.data.frame(expression_data_filterByExpr$counts)


# note: this object is of the class DGEList:
class(expression_data_filterByExpr)




################################################################################
### option 2: Pre-Filtering proposed by DESeq2 #################################
################################################################################

library(DESeq2)


# DESeq2 operates on the format DESeqDataSet which contains information on the count data, the conditions
# of the samples and the design (for further information see below)

# note: for gene expression data sets that must be imported to R, additional steps are necessary before
# the following code can be run 


################################################################################
### step 1: Pre-Filtering ######################################################
################################################################################

# DESeq2 employs the pre-filtering of lowly expressed genes to ... 
# (i) reduce the memory size of dds
# (ii) increase the speed of the DESeq2 pipeline


######################################
# (I) generate pre-filtering indicator 
######################################

# to assess which row (i.e. gene) from the gene expression data set fulfills the minimal 
#read count requirements and which genes are removed for subsequent analysis

# approach 1: keep those rows (i.e. genes) with at least X read counts across all samples 
# the user is free in the choice of the pre-filtering threshold, but keep in mind: 
# a higher X leads to the exclusion of more genes from further analysis, while:
# a higher X causes less genes to be removed 

# approach 1: quite minimal pre-filtering -> set X = 10
indicator_keep1 <- rowSums( expression_data ) >= 10 


# approach 2: keep those rows (i.e. genes) that have at least X read counts in at least Y samples 
# here: keep those genes with at least 10 read counts in at least 10 samples 
indicator_keep2 <- rowSums( expression_data >=10) >= 10

# comment the code of the pre-filtering indicator you NO NOT choose


# note: 
# (i) with function counts(), we can access the count data of the dds object, i.e. the data initially stored
# in object expression_data
# (ii) this "simple" manner of pre-filtering operated on the raw counts data and therefore does NOT take into 
# account the different library sizes of the samples 

############################################
# (II) subset gene expression data set (dds) 
############################################
# to those genes that fulfill the minimum requirement 

# option 1: 
expression_data_filterDESEq2_1 <- expression_data[indicator_keep1,]


# option 2: 
expression_data_filterDESEq2_2 <- expression_data[indicator_keep2,]





######################################################
### Final pre-filtered gene expression data set(s) ###
######################################################

# inspect dimensions of both data sets
dim(expression_data_filterDESEq2_1)
dim(expression_data_filterDESEq2_2)

# -> we see that pre-filtering using option 2 is more stringent than option 1





################################################################################
### Save Results  ##############################################################
################################################################################

### save results for subsequent GSA workflow 

# note: for purpose of simplicity, we only proceed with option 1 of DESeq2's pre-filtering
# however, you can alternatively proceed with the second approach to pre-filtering

expression_data_filterDESeq2 <- as.data.frame(expression_data_filterDESEq2_1)



# save gene expression data set filtered with DESeq2's approach
save(expression_data_filterDESeq2, file = "./Results_PreFiltering/expression_data_filterDESeq2.Rdata")
# save gene expression data set filtered with edgeR's approach
save(expression_data_filterByExpr, file = "./Results_PreFiltering/expression_data_filterByExpr.Rdata")





