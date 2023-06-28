################################################################################
### Instructions: clusterProfiler's ORA ########################################
################################################################################


# set working directory: File in which all resulting data sets are stored 
setwd("/nfsmb/koll/milena.wuensch/Dokumente/GSA_Review/Data")



library(clusterProfiler)
# the following library is required for ORA in combination with geneset database GO
# AND for human gene expression measurements (package contains genome-wide annotation for humans)
library(org.Hs.eg.db) 
# note: if other organisms are investigated, other annotation packages must be loaded, e.g. for mice: 
# library(org.Mm.eg.db)


# important note: it is required to run Instructions_Differential_Expression_Analysis.R

load("./Results_Differential_Expression_Analysis/DE_results_limma_Entrez.Rdata")
load("./Results_Differential_Expression_Analysis/DE_results_DESeq2_Entrez.Rdata")
load("./Results_Differential_Expression_Analysis/DE_results_edgeR_Entrez.Rdata")

# with a gene expression data set with converted gene IDs before running this script

# -> starting point: DE_results_limma, DE_results_edgeR, DE_results_DESeq2 

# we arbitrarily set the object DE_results we will work with in the following code
DE_results <- DE_results_edgeR_Entrez



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


ORA_results <- enrichGO(gene = DEG_vec, 
                        OrgDb = org.Hs.eg.db, # for humans
                        ont = "MF") # "MF" is default, alternatives are "BP" and "CC"

# arguments: 
# gene: vector of differentially expressed genes 
# OrgDb: annotation package for organism at hand (here: human)
# ont: subontology ("MF" by default, alternatives: "BP" and "CC")





###########
# option 2: geneset database KEGG 
###########

ORA_results <- enrichKEGG(gene = DEG_vec, 
                          organism = "hsa")

# arguments: 
# gene: vector of differentially expressed genes 
# organism: organism from which gene expression measurements are obtained
# -> default: "hsa"
# must be adapted for other organisms (such as organism = "mmu" for mice)


