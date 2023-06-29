################################################################################
### Instructions: PADOG ########################################################
################################################################################

# empty environment
rm(list =ls())

# set working directory 
setwd("/nfsmb/koll/milena.wuensch/Dokumente/GSA_Review/Data")




################################################################################
### Content of this script #####################################################
################################################################################

# in this script, we will go through all steps required to run PADOG

################################################################################



# load required libraries 
library(PADOG)


# note: 
# PADOG was initially developed for Microarray Data (which differ strongly 
#   in their characteristics compared to RNA-Seq data)

# ->> we therefore need to perform a suitable transformation of the RNA-Seq data
# to (approximately) align its distributional characteristics with Microarray Data 

# we load the voom-transformed Pickrell data set 
load("./Results_RNASeq_Transformation/expression_data_voomtransformed_Entrez.Rdata")
# alternatively: load the gene expression measurements that have been transformed using 
# DESeq2's varianceStabilizingTransformation 
load("./Results_RNASeq_Transformation/expression_data_vsttransformed_Entrez.Rdata")



# additionally, we load the pickrell data set so that we can access the sample conditions
library(tweeDEseqCountData)
data(pickrell)

# the sample conditions (i.e. phenotype labels) of the pickrell data set can be accessed using
pickrell.eset$gender 


# we proceed with the voom-transformed pickrell data set and the corresponding phenotype labels

# gene expression measurements (transformed)
# note: you can also proceed with the vst-transformed gene expression measurements 
expression_data_transformed <- expression_data_voomtransformed_Entrez
# sample conditions
sample_conditions <- pickrell.eset$gender 



################################################################################
### step 1: Prepare Sample Conditions ##########################################
################################################################################


# first, we inspect the form of the initial (raw) sample conditions

# look at the class: 
class(sample_conditions)
# -> the sample labels are already coded as factor

# the current levels are:
levels(sample_conditions)


# PADOG requires character vector with class labels of the samples
#-> can only contain "c" for control samples or "d" for disease samples


  # prepare sample conditions
  # we want to convert 
  # (i) "female" to "c"
  # (ii) "male" to "d"
  sample_conditions_prep <- factor(sample_conditions, 
                                  levels=c("female","male"), 
                                  labels=c("c","d"))
  

################################################################################
### step 2: Run PADOG ##########################################################
################################################################################
  
  
  # NOTE: It is recommended to set a seed to ensure exact reproducibility 
  #       of the results if the code is run at multiple time points 
  

  
  # note: you can specify any integer number as the seed. It is VERY IMPORTANT
  # to choose the seed arbitrarily and WITHOUT INSPECTING the results 
  # -> the seed should NEVER be specified based on which value yields the 
  # most preferable results 
  
  # run PADOG: 
   PADOG_results <- padog(esetm = as.matrix(expression_data_transformed), 
                          group = sample_conditions_prep, 
                          dseed = 1)
   
   # arguments: 
   #
   # - esetm: matrix that contains the expression measurements 
   # -> note: since the expression data is initially stored in a data frame, 
   #          we transform it to a matrix when running PADOG 
   #
   # - group: sample conditions (has values "c" and "d")
   # - dseed:seed for random number generation (used in the process of phenotype
   #         permutation)
  
   # additional arguments: 
   #
   # - paired: indicates whether the samples in both groups are paired
   #
   # - block: if the samples are paired (i.e. argument paired = TRUE), then the
   #          paired samples must have the same block value 
   #
   # - gslist: gives instructions on how to cluster the genes into gene sets 
   # -> gslist = "KEGGRESTpathway": gene sets correspond to KEGG pathways
   # -> alternative: provide a user-defined list of gene sets 
   # 
   # - organism: organism from which the gene expression measurements are taken
   # -> for human, set organism = "hsa"
   # -> the required character value for other organisms can be extracted 
   #    from the KEGGREST package: 
  
   # if not already installed: 
   # BiocManager::install("KEGGREST")
    library(KEGGREST)
    organisms <- keggList("organism")
   # -> obtain required organisms from column organism 
    
   # annotation: required if gslist is set to "KEGGRESTpathway" and
   #             the rownames of esetm are probe IDs 
   # -> can be set to NULL of gslist is set to "KEGGRESTpathway" and 
   #    the rownames of esetm are in the Entrez gene ID format 
   # -> if rownames are other gene IDs, then sett annotation = NULL and 
   #    make sure that the rownames are elements of gslist (and unique!)
   #
   # gs.names: contains names of gene sets -> character vector 
   # -> must have the same length as gslist 
   #
   # NI: number of phenotype permutations employed in the assessment of 
   #     the significance of a given gene set 
    
    
    
################################################################################
### step 3: Adjust for Multiple Testing ########################################
################################################################################
  
# IMPORTANT: PADOG does not perform multiple test adjustment internally so that 
#            must be performed by the user 
    
# In the PADOG results, the raw (UN-ADJUSTED) p-values are stored in column Ppadog
    
# add adjusted p-value in column Ppadog_adjusted
PADOG_results$Ppadog_adjusted <- p.adjust(PADOG_results$Ppadog)

# by default, p.adjust performs multiple test adjustment using Benjamini and 
# Hochberg method 
  
    

################################################################################
### step 4: Interpretation of Results ##########################################
################################################################################
    
# Differential enrichment of a given gene set can now be assessed based on the
# adjusted p-value in column Ppadog_adjusted 
    
# for instance: detect all gene sets with Ppadog_adjusted < 0.05 as differentially
# enriched 
    
# additional columns: 
    
  # - Name: Name of gene set 
  # - ID: Identifier of gene set 
  # - Size: number of genes in gene set 
  # - meanAbsT0: Mean of absolute (moderated) t-statistic of all genes that are 
  #              a member of the given gene set 
  # - padog0: Mean of abolute weighted moderate t-statistic of all genes that are
  #           a member of the given gene set 
  # - PmeanAbsT: significance of of meanAbsT0
  # - Ppadog: significance padog0
    

  
  
  
  
  
  
  

  




