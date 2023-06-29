################################################################################
### Instructions: GOSeq ########################################################
################################################################################

# empty environment
rm(list =ls())


################################################################################
### Content of this script #####################################################
################################################################################

# In this script, we will do the following two things: 
# 1. Based on the results of differential expression analysis, we will go through 
#    all steps required to run GOSeq 
# 2. We will go through all (meaningful) researchers' degrees of freedom 

################################################################################







# set working directory: File in which all resulting data sets are stored 
setwd("/nfsmb/koll/milena.wuensch/Dokumente/GSA_Review/Data")

library(goseq)
library(tweeDEseqCountData) # for pickrell data set 




# important note: it is required to run Instructions_Differential_Expression_Analysis
# with a gene expression data set before running this script:
# note that GOSeq accepts the gene IDs in the Ensembl gene ID format which is the initial format of the 
# pickrell data set 

# -> starting point: DE_results_limma, DE_results_edgeR, DE_results_DESeq2 (in Ensembl gene ID format)
load("./Results_Differential_Expression_Analysis/DE_results_limma_Ensembl.Rdata")
load("./Results_Differential_Expression_Analysis/DE_results_DESeq2_Ensembl.Rdata")
load("./Results_Differential_Expression_Analysis/DE_results_edgeR_Ensembl.Rdata")

# we arbitrarily set the object DE_results we will work with in the following code
DE_results <- DE_results_DESeq2_Ensembl

# however: you can easily change to the results generated with limma or edgeR:

# DE_results <- DE_results_limma_Ensembl
# DE_results <- DE_results_DESeq2_Ensembl


# generally required: Differential expression results table which contains an adjusted p-value 
# for each gene 
# we also assume that the column that contains the adjusted p-values is named "p_adj"



# ADDITIONALLY: in addition to the results of differential expression analysis, 
#               we also need the raw count data for the case that we want to 
#               adjust the bias to account for in the analysis (see step 4 of 
#               'Researcher's degrees of freedom at the end of this script)

# load pickrell data set 
data(pickrell)
# for the purpose of simplicity, we will work with this data set under the name "expression_data"
expression_data <- Biobase::exprs(pickrell.eset)


#################
# IMPORTANT NOTE: GOSeq accepts genes in the format ENSEMBL gene ID as well as Entrez gene ID 
#################

# -> the type of gene ID must be indicated in both functions run in the course of GOSeq 


################################################################################
### step 1: Preparation of required input object ###############################
################################################################################

# important note: GOSeq's required input object differs from the list of differentially
# expressed genes as it is typically required by ORA

# INSTEAD: 

# GOSeq requires as input a "named binary vector" with ...
# (i) values 0 (for not differentially expressed genes) and 1 (for differentially
# expressed genes)
# (ii) with names that correspond to the given gene IDs 

# -> this vector contains a value (0 or 1) for each gene from the experiment
# with each value named by the corresponding gene ID 

# note: some methods for differential expression analysis, such as DESeq2, 
# set the adjusted p-values of some genes to NA
# this means that these genes are actually not tested for differential expression
# consequence: these genes cannot be assigned neither 1 nor 0 

# approach here: remove genes with an adjusted p-value set to NA from the
# creation of the named binary vector


### Generation of named binary vector: ###


# 1. for each gene from the experiment, indicate those that have an EXISTING (i.e.
# NON-MISSING) p-value
indicator_nonNA_p <- !is.na(DE_results$p_adj) 


# 2. using the indicator from 1., subset the results of differential expression 
# analysis to only those genes with an EXISTING (i.e. NON-MISSING) p-value
DE_results_noNA <- DE_results[indicator_nonNA_p,]

# 3. create vector with binary value (0 or 1) for each gene in DE_results_noNA
# for this, we focus on the column p_adj which contains the adjusted p-values
bin_vec_DEG <- ifelse(DE_results_noNA$p_adj < 0.05 , 1,0)

# note: function ifelse() goes through each adjusted p-value in the results
# of differential expression analysis and in a new vector, assigns value 1 
# if the adjusted p-value is < 0.05 and else value 0 

# 4. assign each value of the vector the corresponding gene ID 
names(bin_vec_DEG) <- rownames(DE_results_noNA)

# 5. results: named binary vector bin_vec_DEG
bin_vec_DEG


# for illustration purposes: compare the number of differentially expressed genes 
# from the results table of differential expression analysis to the number of 
# differentially expressed genes from the named binary vector 

# initial results table of differential expression analysis 
table(DE_results$p_adj < 0.05)

# results table of differential expression analysis subsetted to those genes 
# with a non-missing adjusted p-value
table(DE_results_noNA$p_adj < 0.05)

# named binary vector
table(bin_vec_DEG)

# note: bin_vec_DEG stands for BINary VECtor of Differentially Expressed Genes 

################################################################################
### step 2: Calculate Value of Probability Weighting Function PWF ##############
################################################################################

# GOSeq works with function nullp() which provides one row of data for each 
# gene with the following information: 
# - DEgenes: indicates status of differential expression (0/1) 
# -> this information can be directly extracted from the input object 
# - bias.data: numeric value of bias concerning the detection of differential expression
# -> usually corresponds to the length of the gene 
# pwf: value of the probability weighting function of the gene 

# note: the values of bias.data depend on the function argument bias.data 
# in function nullp() (further information below)


ProbabilityWeighting<- nullp(DEgenes = bin_vec_DEG, 
             genome= "hg19", 
             id= "ensGene")

# function arguments: 
# DEgenes: named binary vector generated in step 1
# genome: indicates organism from which gene expression is measured 
# -> specify organism = "hg19" for human 
# -> for other organisms, the argument must be adapted accordingly (such as "mm9" for mouse)
#
# id: indicates gene identifier
# -> specify id = "knownGene" for Entrez gene ID
# -> specify id = "ensGene" for Ensembl gene ID (which is the case here )
# -> specify id = "geneSymbol" for HGNC gene symbols

# for an overview of available gene identifiers and supported organisms see: 
# supportedGeneIDs()
# for better overview: 
# View(supportedGeneIDs())

# note 1: argument bias.data: specifies data on which the detection of differential expression 
# might be dependent on 
# default value bias.data = NULL: the length of each gene is retrieved from the UCSC
# genome browser
# alternative: the user can provide a vector which contains an entry for each gene 
# contained by bin_vec_DEG

# -> an alternative suggested by the authors of GOSeq would be to account for the 
# read count bias of each gene, i.e. for the total number of read counts across 
# all samples that are assigned to each individual gene 


# note 2: by default, the probability weight function is plotted when nullp() is 
# run
# -> this can be deactivated using adding the function argument plot.fit = FALSE


################################################################################
### step 3: Test for Differential Enrichment ###################################
################################################################################


# test for differential enrichment 
GOSeq_results <- goseq(pwf = ProbabilityWeighting, 
                       genome = "hg19", 
                       id = "ensGene", 
                       test.cats = "GO:BP")

# function arguments: 
# pwf: typically the output of the nullp() function (step 2), i.e. probability weighting function
# genome: indicates organism from which gene expression is measured (see also step 2)
# -> specify organism = "hg19" for human 
# -> for other organisms, the argument must be adapted accordingly (such as "mm9" for mouse)
# id: indicates gene identifier (see also step 2)
# -> specify id = "knownGene" for Entrez gene ID
# -> specify id = "ensGene" for Ensembl gene ID 
# -> specify id = "geneSymbol" for HGNC gene symbols


# test.cats: gene set database
# -> "GO:BP", "GO:MF", "GO:CC" : geneset database GO with respective subontology
# -> "KEGG": geneset database KEGG



################################################################################
### step 4: Multiple Test Adjustment  ##########################################
################################################################################

# important note: GOSeq does not perform adjustment for multiple testing internally
# so that is MUST be performed by the user manually 

# approach: add new column that contains adjusted p-values 
# note: relevant for over-representation test are the p-values from column "over_represented_pvalue"


GOSeq_results$p_adj_overrepresented <- p.adjust(GOSeq_results$over_represented_pvalue) 

# note: by default, p.adjust performs multiple test adjustment based on Benjamini
# and Hochberg 


################################################################################
### step 5: Interpretation of Results Table  ###################################
################################################################################

# final results table: 
GOSeq_results

# Category: provides the ID of the gene set (based on the choice of the gene set
# database in step 4)

# over_represented_pvalue: p-value of over-representation
# IMPORTANT NOTE: DO NOT use these p-values to assess differential enrichment 
# of a given gene set, since these p-value have not been adjusted for multiple
# testing 

# under_represented_p-value: p-value of under-representation
# -> not relevant if you want to test for OVER-REPRESENTATION
# IMPORTANT NOTE: if you are interested in testing for UNDER-REPRESENTATION, 
# you must adjust the p-values for multiple testing: 
# GOSeq_results$p_adj_underrepresented <- p.adjust(GOSeq_results$under_represented_pvalue) 
 

# numDEInCat: gives the number of differentially expressed genes from the input
# that are members of the given gene set 

# numInCat: number of genes that are members of the given gene set 

# term: description of the gene set 
# -> NOTE: this column is only provided for the geneset database GO 

# ontology: subontology of geneset database GO (based on choice in step 3)
# -> this column is only provided for the geneset database GO 



# p_adj_overrepresented: p-value of over-representation that has been ADJUSTED 
# for multiple testing
# -> based on these adjusted p-values, differential enrichment (i.e. significant
# over-representation) can be adjust
# typically: detect those gene sets as differentially enriched with a value of 
# p_adj_overrepresented < 0.05 








################################################################################
### Researchers' Degrees of Freedom in GOSeq ###################################
################################################################################

# In this part, we want to go through all parameters that can be adapted in the 
# GOSeq workflow 

# IMPORTANT NOTE: the intention behind going through the researchers' degrees of freedom
# is to give you an understanding of what you can do to adapt the given (parameter)
# setting to the research question 

# ->> MORE IMPORTANTLY: the intention behind going through these flexible parameters is 
# NOT to change them in order to help you obtain the most preferable results by systematically 
# changing these parameters 
# -> such behaviour would correspond to "cherry-picking" 

# Any changes in the parameter choice should be documented transparently 


#######################################
# change 1: change geneset database ###
#######################################

# the gene set database can be adapted in function goseq() in argument "test.cats"

# example: change gene set database to KEGG 
GOSeq_results_database <- goseq(pwf = ProbabilityWeighting, 
                       genome = "hg19", 
                       id = "ensGene", 
                       test.cats = "KEGG")

# -> note: when the gene set database KEGG is specified, then no column "term"
#          is provided to give a description of the corresponding gene set 

# other gene set databases can be specified by setting the argument 
# test.cats = "GO:CC" -> GO with subontology Cellular Component
# test.cats = "GO:MF" -> GO with subontology Molecular Function 



###############################################################################
# change 2: Include genes that are no nember of any gene set in calculation ###
########### of p-value ########################################################
###############################################################################

# the inclusion of these additional genes can be specified in function goseq()
# in argument "use_genes_without_cat"

GOSeq_results_method <- goseq(pwf = ProbabilityWeighting, 
                                    genome = "hg19", 
                                    id = "ensGene", 
                                    test.cats = "GO:BP", 
                                    use_genes_without_cat = TRUE)



################################################################################
# change 3: change method for calculation of p-value ###########################
################################################################################

# the method for the calculation of the p-value can be adapted in function goseq()
# in argument "method"


# example: change the method for the calculation of the p-value to the computationally 
# expensive resampling 
GOSeq_results_method <- goseq(pwf = ProbabilityWeighting, 
                              genome = "hg19", 
                              id = "ensGene", 
                              test.cats = "GO:BP", 
                              method = "Sampling") 

# IMPORTANT NOTE: by default, 2000 samples are calculated in the course of resampling
#                 -> see argument "repcnt"
#                 -> this argument should NEVER be played around with to generate 
#                    preferable results
#                 -> instead: it is advisable to keep this parameter in its default 
#                    configuration


# note: theoretically, GOSeq additionally offers the use of the standard hypergeometric 
# distribution as the method to calcualte a p-value of over-representation
# -> however: users are explicitly advised against using this option since the 
#             standard hypergeometric distribution does not adjust for ANY biases
#             that might be present in the data 



################################################################################
# change 4: change the bias to account for in the analysis #####################
################################################################################

# example: as an example taken from the user manual of GOSeq, we now want to adjust 
# for total number of counts of each gene. 

# the idea behind this is that in the context of RNA-Seq data, the magnitude of 
# counts of a given gene sets reflects its overall expression level. In RNA-Seq data,
# which takes the form of count data, a higher magnitude of counts leads to an increased
# statistical power (to detect a gene as differentially expressed).

# IMPORTANT NOTE: accounting for the total number of counts might also remove the 
# bias that results from actual differential expression between both conditions 

# Account for total read counts: 

#########
# step 1: create a vector that contains the total number of read counts for each gene 
########  from the experiment that is a part of the input object 

# 1. Compute the sum of each gene across all samples from the gene expression in the gene 
# expression data set 
countbias <- rowSums(expression_data)

# 2. Subset the vector countbias to all genes that are part of the input object 
#    bin_vec_DEG 

# (i) indicate which genes are a part of the input object: 
ind_input <- names(countbias) %in% names(bin_vec_DEG)

# (ii) subset vector countbias 
countbias <- countbias[ind_input]

#########
# step 2: Provide countbias to function nullp() using argument "bias.data"
#########


ProbabilityWeighting_countbias <- nullp(DEgenes = bin_vec_DEG, 
                                        genome= "hg19", 
                                        id= "ensGene", 
                                        bias.data = countbias) # countbias 

#########
# step 3: Run function GOSeq with adjusted probability weightings 
########


GOSeq_results_countbias <- goseq(pwf = ProbabilityWeighting_countbias, # adjusted probability weightings
                                 genome = "hg19", 
                                 id = "ensGene", 
                                 test.cats = "GO:BP")



#########
# step 4: Multiple test adjustment 
########


GOSeq_results_countbias$p_adj_overrepresented <- p.adjust(GOSeq_results_countbias$over_represented_pvalue) 































