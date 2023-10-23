################################################################################
### Instructions: DAVID ########################################################
################################################################################

# empty environment
rm(list =ls())

# set working directory
setwd("/nfsmb/koll/milena.wuensch/Dokumente/GSA_Review")

# create the folder "Input_Objects_DAVID" to eventually contain the objects 
# generated in this script
dir.create("./data/Input_Objects_DAVID")

################################################################################
### Content of this script #####################################################
################################################################################

# In this script, we will

# (i) preprocess and export the required input object AND
# (ii) generate the alternative universe
# for the web-based tool DAVID



################################################################################


# Required input object: a list of differentially expressed genes
# Alternative universe: all genes whose differential expression was assessed in the experiment
# -> Note that it is not required to provide an alternative universe. We illustrate its use as
# part of the researchers' degrees of freedom

# Note that DAVID works on Ensembl IDs.


# Load data

# Both objects, the list of differentially expressed genes as well as the alternative universe
# can be extracted from the results table of differential expression analysis.

# Here, we will illustrate the process for the limma results only since the process is identical
# for all methods for differential exression analysis (due to the unification of the column names
# that refer to the adjusted p-values).

load("./data/Results_Differential_Expression_Analysis/DE_results_limma_Ensembl.Rdata")

# Alternatively, you can proceed with the results generated with DESeq2 or edgeR:
# load("./Results_Differential_Expression_Analysis/DE_results_DESeq2_Ensembl.Rdata")
# load("./Results_Differential_Expression_Analysis/DE_results_edgeR_Ensembl.Rdata")

# For the purpose of redability, we store limma's differential expression analysis results
# in an object with a more neutral name:

DE_results <- DE_results_limma_Ensembl

# Take a look at the results table:
head(DE_results, n  = 10)


################################################################################
### step 1: Generate required input object #####################################
################################################################################

# Create the vector of differentially expressed genes

# (i) indicate which genes have an (existing!) adjusted p-value < 0.05.
ind_DE <- (!is.na(DE_results$p_adj)) & (DE_results$p_adj < 0.05)

# get overview of number of differentially expressed genes
table(ind_DE)

# (ii) use the indicator to obtain the list of differentially expressed genes
DEG_vec <- rownames(DE_results)[ind_DE]

# take a look at the vector
head(DEG_vec, n = 10 )


################################################################################
### step 2: Export vector of differentially expressed genes ####################
################################################################################


# create path
path <- "./data/Input_Objects_DAVID/DEG_vec.txt"
# this path indicates that the vector is stored in the text document "DEG_vec" in folder "Input_Objects_DAVID".

# export
write.table(DEG_vec,
            file = path,
            quote=FALSE,
            row.names=FALSE,
            col.names = FALSE)


################################################################################
### step 3: Upload to DAVID ####################################################
################################################################################

# The resulting .txt file can be directly uploaded to the DAVID website:

#  https://david.ncifcrf.gov/





################################################################################
### Researchers' Degrees of freedom: Alternative universe ######################
################################################################################

# In this part, we will illustrate how to create an alternative universe to DAVID, which is
# a researchers' degree of freedom. It is important to note that the intention behind going
# through the researchers' degrees of freedom is to give you an understanding of what you can
# do to adapt the given (parameter) setting to the research question. It is even more important
# to keep in mind that the intention behind going through these flexible parameters is NOT to change
# them in order to help you obtain the most preferable results by systematically changing these parameters
# since such behaviour would correspond to "cherry-picking". Any changes in the parameter choice must be
# documented transparently.


#########
# step 1: Generate alternative universe
#########

# For the alternative universe, we want to consider all genes whose differential expression was measured in the experiment
# For some methods for differential expression analysis, such as DESeq2, the adjusted p-values of some genes are set to NA
# which means that these genes cannot be detected as neither differentially expressed nor not differentially expressed. We
# therefore want to remove these from the universe



# indicate which genes have an adjusted p-value in the results of differential expression analysis
ind_nona_p <- !is.na(DE_results$p_adj)

# filter the list of genes to those with an existing adjusted p-value
alternative_universe <- rownames(DE_results)[ind_nona_p]


# inspect the first few genes from the universe
head(alternative_universe, n = 10)


#########
# step 2: Export alternative universe
########


# create path
path_alt_universe <- "./data/Input_Objects_DAVID/alternative_universe.txt"
# this path indicates that the alternative universe is stored in the text file "alternative universe" in the folder
# "Input_Objects_DAVID"

write.table(alternative_universe,
            file = path_alt_universe,
            quote=FALSE,
            row.names=FALSE,
            col.names = FALSE)

# For a description of the arguments see above.





