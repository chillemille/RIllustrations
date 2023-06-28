################################################################################
### Convert gene IDs from Ensembl gene ID to Entrez gene ID ####################
################################################################################

library(tweeDEseqCountData)
data(pickrell)

expression_data <- Biobase::exprs(pickrell.eset)

# inspect dimension of expression_data: 
dim(expression_data)

# -> in the case of the pickrell gene expression data set, the gene expression of 
# 52580 genes is measured 


################################################################################
### Some preprocessing: (simple) pre-filtering #################################
################################################################################

# In alignment with the pre-processing workflow proposed in the main part of the paper, 
# we perform pre-filtering before we convert the gene IDs 
# for the purpose of simplicity, we use the very simple pre-filtering approach proposed 
# alongside with DESeq2: 

# keep all genes that have at least 10 read counts across all samples: 


# filtering indicator which indicates whether the condition is true for each sample: 
indicator_keep <- rowSums(expression_data) >= 10 

# filter based on pre-filtering indicator
expression_data <- expression_data[indicator_keep,]

# inspect how many genes are left after pre-filtering: 
dim(expression_data)


################################################################################
### option 1: ##################################################################
################################################################################


library(org.Hs.eg.db)


# mapIds() gets the mapped IDs (column) for a set of keys that are of a particular
# keytype
conversion1 <- mapIds(org.Hs.eg.db, 
       keys = rownames(expression_data), 
       column = 'ENTREZID', 
       keytype = 'ENSEMBL')


# -> this function returns a vector that contains as values the corresponding Entrez gene ID
# for each Ensembl gene ID (Ensembl gene IDs are the names of the vector (see names(conversion1)))

# Inspection of Conversion Pattern: 

# create data frame with two columns that contain Ensembl IDs and corresponding Entrez IDs 
conversion_df <- data.frame(Entrez = conversion1, 
                                Ensembl = names(conversion1)) 

# inspect dimenstion:
dim(conversion_df)

# at this point, the data frame has the same number of observations as the data frame 
# expression data. This makes it as if each Ensembl ID could be converted to a corresponding
# Entrez gene ID. However, if we inspect the number of NAs in the column Entrez, 
# we see the following: 

sum(is.na(conversion_df$Entrez))

# this means that there is a relatively high number of Ensembl IDs that could not be 
# converted to a corresponding Entrez gene ID 

# We now remove these NAs from the data set: 
conversion_df <- conversion_df[!is.na(conversion_df$Entrez),]

# resulting dimension: 
dim(conversion_df)

# We now want to inspect the number of duplicates 
sum(duplicated(conversion_df$Entrez))

# conversion duplicates
duplicates <- conversion_df[duplicated(conversion_df$Entrez),]
# -> by ordering this data frame by Entrez IDs, we can see that there are 
# several cases in which multiple distinct Ensembl gene IDs were mapped to 
# the same Entrez gene ID 

################################################################################
### option 2: clusterProfiler ##################################################
################################################################################

# load required libraries:
library(clusterProfiler)
# load annotation database to indicate the organism from which the gene expression
# measurements are taken 
# for human: 
library(org.Hs.eg.db)
# e.g. for mouse, the following library must be loaded: 
# library(org.Mm.eg.db)

# get mapping pattern from ENSEMBL gene IDs to Entrez gene IDs 
bitr_mapping <- bitr(geneID = rownames(expression_data), 
                     fromType = "ENSEMBL", 
                     toType = "ENTREZID", 
                     OrgDb = org.Hs.eg.db)

# note: bitr_mapping is a data frame 

# arguments: 
# - geneID: vector of gene IDs in the initial format to be converted to the accepted format
# - fromType: initial format of gene IDs 
# - toType: accepted/required format of gene IDs 
# - OrgDb: annotation database to indiacte organism (corresponding library must be loaded
#         prior to running bitr() (see above))


# 'fromType' and 'toType' should be one of ACCNUM, ALIAS, ENSEMBL, ENSEMBLPROT, ENSEMBLTRANS, ENTREZID, 
# ENZYME, EVIDENCE, EVIDENCEALL, GENENAME, GENETYPE, GO, GOALL, IPI, MAP, OMIM, ONTOLOGY, 
# ONTOLOGYALL, PATH, PFAM, PMID, PROSITE, REFSEQ, SYMBOL, UCSCKG, UNIPROT.


### inspection of data frame bitr_mapping 

dim(bitr_mapping)
# -> we see a direct reduction of the genes (from 52580) to 31336 genes 
# this means that for 52580 - 31336 = 21244 genes in the initial gene ID format,
# there exists no corresponding gene ID in the accepted/required format 

### inspection of duplicates: 
 

# (i) case 1: a single initial (Ensembl) gene IDs is mapped to multiple 
# accepted (Entrez) gene ID
sum(duplicated(bitr_mapping$ENSEMBL))

# (ii) case 2: several distinct initial (Ensembl) gene IDs are mapped to a single 
# accepted (Entrez) gene ID 
sum(duplicated(bitr_mapping$ENTREZID))






################################################################################
### option 3: biomaRt ##########################################################
################################################################################

library(biomaRt)


mart <- useDataset(dataset = "hsapiens_gene_ensembl", 
                   mart = useMart("ensembl"))

# update the Mart object  
# - dataset: select mart database 
#           -> corresponds to the type of data you are interested in 
#           -> get an overview over all possible data sets: listDatasets(mart = mart)
# - mart: select mart database 
#         -> corresponds to the species you are interested in and need to collect 
#            data from 



genes <-  rownames(expression_data)


gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","entrezgene_id"),
                  values = genes, mart= mart)

