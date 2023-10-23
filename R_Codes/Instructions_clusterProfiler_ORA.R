################################################################################
### Instructions: clusterProfiler's ORA ########################################
################################################################################

# empty environment
rm(list =ls())


# set working directory
setwd("/nfsmb/koll/milena.wuensch/Dokumente/GSA_Review")

################################################################################
### Content of this script #####################################################
################################################################################

# In this script, we will 
# (i) run clusterProfiler's ORA tool for the gene set databases 
# - KEGG 
# - GO (with default subontology "MF")
# (ii) go through all meaningful researchers' degrees of freedom 


################################################################################


library(clusterProfiler)
# the following library is required for ORA in combination with gene set database GO
# AND for human gene expression measurements (package contains genome-wide annotation for humans)
library(org.Hs.eg.db) 
# note: if other organisms are investigated, other annotation packages must be loaded, e.g. for mice: 
# library(org.Mm.eg.db)


# important note: it is required to run Instructions_Differential_Expression_Analysis.R

load("./data/Results_Differential_Expression_Analysis/DE_results_limma_Entrez.Rdata")
load("./data/Results_Differential_Expression_Analysis/DE_results_DESeq2_Entrez.Rdata")
load("./data/Results_Differential_Expression_Analysis/DE_results_edgeR_Entrez.Rdata")

# with a gene expression data set with converted gene IDs before running this script

# -> starting point: DE_results_limma, DE_results_edgeR, DE_results_DESeq2 

# we arbitrarily set the object DE_results we will work with in the following code
DE_results <- DE_results_DESeq2_Entrez



# generally required: Differential expression results table which contains an adjusted p-value 
# for each gene 
# we also assume that the column that contains the adjusted p-values is named "p_adj"


################################################################################
### step 1: Preparation of required input object ###############################
################################################################################


# clusterProfiler's ORA requires as input a list of differentially expressed genes 
# we extract such list from the results table of differential expression analysis 

# for each gene from the results of differential expression analysis, we indicate
# whether it is differentially expressed (TRUE) or not differentially expressed (FALSE)
# based on the following two criteria: 
# (i)  it was tested for differential expression, i.e. has a non-missing adjusted p-value AND
# (ii) it has an adjusted p-value < 0.05

ind_differentially_expressed <- ((!is.na(DE_results$p_adj)) & (DE_results$p_adj<0.05))

# using this indicator, we extract the vector of differentially expressed genes from the results
# of differential expression analysis 

DEG_vec <- rownames(DE_results[ind_differentially_expressed,])




################################################################################
### step 2: Run clusterProfiler's ORA ##########################################
################################################################################

# ORA can be run with the common geneset databases KEGG and GO as well as user-
# defined geneset databases 


###########
# option 1: geneset database GO (with subontology Molecular Function by default)
###########


ORA_results_GO <- enrichGO(gene = DEG_vec, 
                        OrgDb = org.Hs.eg.db, # for humans
                        ont = "MF") # "MF" is default, alternatives are "BP" and "CC"

# arguments: 
# gene: vector of differentially expressed genes 
# OrgDb: annotation package for organism at hand (here: human)
# ont: subontology ("MF" by default, alternatives: "BP" and "CC")





###########
# option 2: geneset database KEGG 
###########

ORA_results_KEGG <- enrichKEGG(gene = DEG_vec, 
                          organism = "hsa")

# arguments: 
# gene: vector of differentially expressed genes 
# organism: organism from which gene expression measurements are obtained
# -> default: "hsa"
# must be adapted for other organisms (such as organism = "mmu" for mice)



# description of columns in results table: 

# GeneRatio: number of genes from the input list that are members of the given gene  
#            set divided by the number of genes from the input list that are NOT members
#            of the given gene set 
#
# BgRatio: number of genes from the universe that are members of the gene set divided by 
#          the total number of genes in the universe 
#
# pvalue: p-value of enrichment of the given gene set 
#
# p.adjust: p-value of enrichment ADJUSTED for multiple testing 
# qvalue: p-value of enrichment ADJUSTED for multiple testing 
# -> note: p.adjust and qvalue are adjusted using slightly different approaches 
#
# geneID: list of genes from the input list that are members of the given gene set 
# count: number of genes from the input list that are members of the given gene set 





################################################################################
### Researchers' Degrees of Freedom in clusterProfiler's ORA ###################
################################################################################

# In this part, we want to go through all parameters that can be adapted in  
# clusterProfiler's ORA workflow 

# IMPORTANT NOTE: the intention behind going through the researchers' degrees of freedom
# is to give you an understanding of what you can do to adapt the given (parameter)
# setting to the research question 

# ->> MORE IMPORTANTLY: the intention behind going through these flexible parameters is 
# NOT to change them in order to help you obtain the most preferable results by systematically 
# changing these parameters 
# -> such behaviour would correspond to "cherry-picking" 

# Any changes in the parameter choice should be documented transparently 


###############################
# change 1: change universe ###
###############################

# here, we change the universe to all genes measured in the experiment 

# note 1: we here do not change the universe to the interception between all 
# genes from the experiment and the list of genes annotated to the given gene 
# set database since we found no way to obtain the latter list of genes 

# note 2: we want to restrict ourselves to all genes in the experiment that 
# HAVE an adjusted p-value (i.e. whose expression was indeed measured)
# -> intuition: e.g. for DESeq2, some genes are filtered out internally and 
# therefore do not have an adjusted p-value. These genes therefore neither be  
# detected as differentially expressed or not differentially expressed so it would
# not be feasible to include them in the universe. 

########
# step 1: setup alternative universe 
########

# (i) indicate which genes have an adjusted p-value
ind_adjp <- !is.na(DE_results$p_adj)

# (ii) filter the genes from the experiment to those genes that do have an adjusted 
# p-value 
alternative_universe <- rownames(DE_results)[ind_adjp]


########
# step 2: add alternative universe as a parameter to ORA 
########

# (a) gene set database GO: specify parameter "universe"

ORA_results_GO_universe  <- enrichGO(gene = DEG_vec, 
                           OrgDb = org.Hs.eg.db, 
                           universe = alternative_universe)


# (b) gene set database KEGG: specify parameter "universe"

ORA_results_KEGG_universe <- enrichKEGG(gene = DEG_vec, 
                               organism = "hsa",
                               universe = alternative_universe)









