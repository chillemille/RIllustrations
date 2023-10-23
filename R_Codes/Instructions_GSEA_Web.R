################################################################################
### Instructions: GSEA (web application) #######################################
################################################################################

# empty environment
rm(list =ls())


# set working directory 
# you will have to adapt the path within setwd() to your computer 
setwd("/nfsmb/koll/milena.wuensch/Dokumente/GSA_Review")

# create folder to eventually store the objects generated in this script
dir.create("./data/Input_Objects_GSEA_web")


################################################################################
### Content of this script #####################################################
################################################################################

# In this script, we will 
# (i) export the transformed RNA-Seq measurements 
# (ii) prepare and export sample conditions 

# note: some further preprocessing steps are required in Excel (link with
# instructions are provided below


# http://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats


#############################
### required input: #########
#############################


### 1. required input object: (transformed) gene expression measurements ###

# transformed (and normalized) gene expression data set 
# as gene ID format, GSEA accepts a quite wide variety, thereunder 
# (i) Ensembl gene ID ("Human_Ensembl_Gene_ID_MSigDB.v...")
# (ii) Entrez gene ID ("Human_NCBI_Gene_ID_MSigDB.v...")
# ...

# note: if the tool accepts the genes in the format they are initially identified in 
# in the gene expression data set, it does make sense to work with with exactly this 
# format 
# reason: gene ID conversion leads to an immediate reduction of the gene expressiond data
# set (and therefore also of information)

# We therefore proceed with the pre-filtered and transformed gene expression data set 
# with the genes idenfitied in the ENSEMBL gene ID format 
load("./data/Results_RNASeq_Transformation/expression_data_voomtransformed_Ensembl.Rdata")

# we create a new object which contains the transformed gene expression data set for 
# the purpose of clarity and so that it is easier for you to insert your own (transformed)
# RNA-Seq data set if desired 
expression_data_transformed <- expression_data_voomtransformed_Ensembl


### 2. required input object: assignments of samples to the conditions ###

# load pickrell data set (which is the RNA-Seq data set on which our entire analysis)
# here is based
library(tweeDEseqCountData)
data(pickrell)
# for the pickrell data set used in these illustrations, we obtain the assignments 
# of the conditions of the samples and assign it to a new object, again, for the 
# purpose of clarity: 


sample_conditions <- pickrell.eset$gender





########################################################################
### step 1: Export of the (transformed) gene expression measurements ###
########################################################################

# we export the transformed gene expression measurements to a .txt file and into
# the pre-specified folder "Input_Objects_GSEA_web"

# generate path 
# this path indicates that we will store the transformed gene expression measurements 
# in the object "expression_data_transformed.txt, which is located in the file 
# "Input_Objects_GSEA_web"
path_measurements  <- "./data/Input_Objects_GSEA_web/expression_data_transformed.txt"

### export 
write.table(expression_data_transformed,
            file = path_measurements,
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)

# note: 
# - quote = FALSE ensures that none of the characters (in this case gene and sample
# identifiers) are surrounded by double quotes
# - row.names = TRUE ensures that the gene IDs are included in the export 
# - col.names = TRUE ensures that the samples IDs are included in the export 


########################################################
### step 2: Prepare and export the sample conditions ###
########################################################

# GSEA accepts the sample conditions in the (binary) format with values 
# 0 and 1 

# for the pickrell data set, we inspect the current values in which the 
# sample conditions are stored:

levels(sample_conditions)
# currently, the sample conditions are coded as "male" and "female"


# convert levels ("female", "male") to 0 and 1:
# (i) create vector to contain the sample conditions in the right format
sample_conditions_prepared <- c()
# (ii) assign all "female"s, which is the first level of the factor, the value 0
sample_conditions_prepared[ pickrell.eset$gender == levels(pickrell.eset$gender)[[1]]] <- 0
# (iii) assign all "male"s, which is the second level of the factor, the value 1
sample_conditions_prepared[ pickrell.eset$gender == levels(pickrell.eset$gender)[[2]]] <- 1


# this path indicates that we store the prepared sample conditions in the object 
# sample_conditions_prepared.txt
path_conditions <- "./data/Input_Objects_GSEA_web/sample_conditions_prepared.txt"

write.table(x = sample_conditions_prepared,
            file = path_conditions,
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)


# note: 
# - quote = FALSE ensures that none of the characters (in this case gene and sample
# identifiers) are surrounded by double quotes
# - row.names = TRUE ensures that the gene IDs are included in the export 
# - col.names = TRUE ensures that the samples IDs are included in the export 


################################################################################
### step 3: Further preparation of input object in Excel #######################
################################################################################

# For instructions of the further preparation of the text file  in Excel open 
# the following link: 
# http://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats

# and follow the instructions in section 

# (i) "GCT: Gene Cluster Text file format (*.gct)"
# -> for further preparation of expression_data_transformed

# (ii) "CLS: Categorical (e.g tumor vs. normal) class file format (*.cls) "
# -> for further preparation of sample_conditions_prepared"





