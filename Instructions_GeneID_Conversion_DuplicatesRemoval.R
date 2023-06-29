################################################################################
### Guide: Conversion of gene IDs and removal of resulting duplicates ##########
################################################################################

# empty environment
rm(list =ls())


# set working directory: File in which all resulting data sets are stored 
setwd("/nfsmb/koll/milena.wuensch/Dokumente/GSA_Review/Data")



################################################################################
### Content of this script #####################################################
################################################################################

# In this script, we want to go through the process of 
# (1.) converting the given gene IDs in the initial format to a different format 
# (2.) removing the resulting duplicated gene IDs 

# we repeat this process for the two pre-filtered gene expression data sets created
# in script "Instructions_PreFiltering.R":
# (i) expression_data_filterByExpr: gene expression data set that has been pre-filtered
#                                   using edgeR's builtin function filterByExpr()
# (ii) expression_data_filterDESeq2: gene expression data set that has been pre-filtered
#                                    using the approach proposed by DESeq2 

################################################################################





# Most GSA methods/tools accept the genes from the gene expression data set in 
# (a number of) specific formats which might not always coincide with the gene 
# ID format of the given gene expression data set 

# example: you dataset at hand identifies genes in the Ensembl gene ID format, 
# whereas your method of choice only accepts Entrez gene ID formats 

# a conversion of the given to the required/accepted gene ID format is therefore
# necessary 



# In this example, we assume that 
# 1. gene expression measurements originate from humans
# 2. the genes in the given gene expression data set are provided in the Ensembl gene ID format, 
# 3. whereas the gene set analysis tool requires the genes to be identified in the Entrez gene ID format

# 4. the gene expression data set has already been pre-filtered 
# based on script "Instructions_RNA_Seq_Transformation.R"
load("./Results_PreFiltering/expression_data_filterDESeq2.Rdata")
load("./Results_PreFiltering/expression_data_filterByExpr.Rdata")

# IMPORTANT NOTE: We proceed as follows
# (i) We perform the workflow of gene ID conversion and duplicate gene ID removal for the 
#     pre-filtered gene expression data set expression_data_filterByExpr
# (ii) We perform the EXACT SAME workflow for the second pre-filtered gene expression data set
#      expression_data_DESeq2 afterwards

# the reason behind this approach is that for the different methods of differential expression analysis, 
# which we will illustrate in "Instructions_Differential_Expression_Analysis.R"
# different methods of pre-filtering are proposed. 

# In the following, we work with the functionalities provided by clusterProfiler
library(clusterProfiler)



################################################################################
# (i) expression_data_filterByExpr #############################################
################################################################################

# inspect data set: 
# View(expression_data_filterByExpr)
dim(expression_data_filterByExpr)

# ->> this expression data set contains gene expression measurements of 6246 genes and 69 samples 

# ->> goal: convert Ensembl gene IDs to Entrez gene IDs 

# here, we work with functionalities provided by clusterProfiler 



#############################
# (I) Conversion of gene IDs
#############################




#######################################################################################
### step 1: indicate which organism the gene expression measurements originate from ###
#######################################################################################



# assumption 1: gene expression measurements originate from humans 
# -> in function clusterProfiler's function bitr(), this corresponds to argument OrgDb = org.Hs.eg.db

# IMPORTANT: library "org.Hs.eg.db", which contains the corresponding genome wide annotation,  
# must be loaded if the given organism is human: 
library(org.Hs.eg.db)

# note: if, for instance, the organism is mouse, then the library "org.Mm.eg.db" must be loaded via 
# library(org.Mm.eg.db)
# and the argument OrgDb must be set to OrgDb = org.Mm.eg.db

# assumption 2: genes originate from Ensembl gene ID format 
# -> this corresponds to the argument fromType = "ENSEMBL"

# assumption 3: gene set analysis tool requires genes to be in Entrez gene ID format 
# -> this corresponds to argument toType = "ENTREZID" 


# NOTE: arguments "fromType" and "toType" must be set as one of the following, depending on the given and the 
# required gene ID format
  # ACCNUM, ALIAS, ENSEMBL, ENSEMBLPROT, ENSEMBLTRANS, ENTREZID, ENZYME, EVIDENCE, 
  # EVIDENCEALL, GENENAME, GENETYPE, GO, GOALL, IPI, MAP, OMIM, ONTOLOGY, ONTOLOGYALL, PATH, PFAM, PMID, PROSITE, 
  # REFSEQ, SYMBOL, UCSCKG, UNIPROT.
  
  
  
  
##############################################################
### step 2: Gene ID conversion via clusterProfiler::bitr() ###
##############################################################
  
  
  
  # obtain mapping of initial (Ensembl) gene IDs to required (Entrez) gene IDs: 
  bitr_EnsToEntr <- bitr(rownames(expression_data_filterByExpr), # the rownames of the gene expression data set correspond to the gene IDs 
                         fromType = "ENSEMBL", 
                         toType = "ENTREZID", 
                         OrgDb = org.Hs.eg.db)
  #note: not all ENSEMBL IDs could be converted to a corresponding ENTREZ Gene ID: 
  dim(bitr_EnsToEntr)
  
  #results: 
  #- not all Ensembl gene IDs can be mapped to a corresponding Entrez gene ID
  #- some individual Ensembl gene IDs were mapped to multiple distinct Entrez gene IDs
  #- some distinct Ensembl gene IDs were mapped to an identical Entrez gene ID 
  
  ### merge
  
  # concatenate initial gene expression data set with the mapping from initial (Ensembl) to required (Entrez) gene ID format
  # merge by row names of expression data set and ENSEMBL ID of gene ID mapping
  counts_expression_data_filterByExpr <- merge(x = expression_data_filterByExpr, 
                                               y = bitr_EnsToEntr, 
                                               by.x=0, 
                                               by.y="ENSEMBL", 
                                               all.y=TRUE, sort=TRUE)
  
  
  # note on function arguments: 
  # by.x and by.y specify by which columns expression_data_filterByExpr and bitr_EnsToEntr are concatenated: 
  # by.x = 0: use row names of expression_data_filterByExpr (the row names contain the gene IDs)
  # by.y = "ENSEMBL": use the column that contains the Ensembl gene IDs 
  
  
  
  # take a look at dimension of resulting data set
  dim(counts_expression_data_filterByExpr)
  # - number of genes in counts_expression_data_filterByExpr has decreased to 6040
  # -> this is directly caused by the circumstance that typically, not all initial (Ensembl) gene IDs can be 
  # converted to the required (Entrez) gene ID 
  # - the expression data set has been extended by two columns:
  # (1.) the very last column (here called "ENTREZID" containing the converted (Entrez) gene ID
  # (2.) the very first column (here called Row.names) with the initial (Ensembl) gene IDs that were originally stored in the row names
  
  # note: in the resulting data set, the row names now correspond to the row number 
  
  
##############################################
### step 3: take closer look at duplicates ###
##############################################
  
  ########
  # CASE 1: single ENSEMBL IDs are mapped to multiple ENTREZ IDs
  ########
  
  # obtain number of cases in which an ENSEMBL gene ID was converted to several ENTREZ IDs, i.e. the number of 
  # times an Ensembl ID appears more than once in the mapping
  sum(duplicated(bitr_EnsToEntr$ENSEMBL)) 
  # determine all duplicated ENSEMBL gene IDS (i.e. all ENSEMBL gene IDs that were mapped to multiple distinct Entrez gene IDs):
  dupl_ensembl<-unique(bitr_EnsToEntr$ENSEMBL[duplicated(bitr_EnsToEntr$ENSEMBL)])
  # number of ENSEMBL IDs that have at least one duplicate
  length(dupl_ensembl)
  
  
  
  # in the following, we can inspect the conversion pattern for each Ensembl ID that was mapped to 
  # multiple distinct Entrez IDs 
  duplicated_conversion_ens<-bitr_EnsToEntr[bitr_EnsToEntr$ENSEMBL %in% dupl_ensembl,]
  # take a look at conversion pattern for duplicated Ensembl gene IDs:
  duplicated_conversion_ens
  # -> we can see perfectly that these Ensembl ID are each mapped to multiple individual Entrez IDs 

  
  ########
  # CASE 2: multiple ENSEMBL IDs are mapped to single Entrez ID
  ########
  
  # obtain number of cases in which multiple distinct Ensembl IDs are converted to the same Entrez gene ID 
  # in these cases, the corresponding Entrez IDs appear repeatedly in the mapping:
  sum(duplicated(bitr_EnsToEntr$ENTREZ)) 
  # determine all Entrez gene IDs affected by this duplication
  dupl_entrez<-unique(bitr_EnsToEntr$ENTREZID[duplicated(bitr_EnsToEntr$ENTREZID)])
  # number of ENTREZ IDs that have at least one duplicate
  length(dupl_entrez)
  # display of conversion pattern of duplicated ENTREZ IDs
  duplicated_conversion_entrez<-bitr_EnsToEntr[bitr_EnsToEntr$ENTREZID %in% dupl_entrez,]
  # take a look at conversion pattern: 
  duplicated_conversion_entrez
  # if ordered by Entrez IDs, we can see perfectly that for each of these Entrez gene IDs, there are multiple
  # corresponding Ensembl gene IDs 
  dim(duplicated_conversion_entrez)
  
  
  
  ##############################################
  ### step 4: remove duplicated gene IDs  ######
  ##############################################
  
  
  # for subsequent analysis (whether differential expression analysis of gene set analysis), we need to deal 
  # with the duplicated gene IDs and remove them in a suitable way such that only one unique gene expression 
  # measurement among the duplicates remains
  
  # there is no recommended way to proceed, i.e. no common approach presented in an official scientific puplication, 
  # but instead several approaches suggested by users in corresponding user platforms, 
  # generally, the manner of duplicate gene ID removal is performed at the discretion of the user 
  
  # we therefore present three approaches to duplicate gene ID removal 
  
  
  
  ############
  ### option 1: keep first subscript among duplicates 
  ############
  
  # this is the simplest approach among the three 
  

    #1. remove duplicated ENTREZ gene IDs
    exprdat_filterByExpr_dupl1<-counts_expression_data_filterByExpr[!duplicated(counts_expression_data_filterByExpr$ENTREZID),]

    
    #2. remove duplicated ENSEMBL gene IDs
    exprdat_filterByExpr_dupl1<-exprdat_filterByExpr_dupl1[!duplicated(exprdat_filterByExpr_dupl1$Row.names),]
    dim(exprdat_filterByExpr_dupl1)
    
    #3. set Entrez gene IDs as rownames 
    rownames(exprdat_filterByExpr_dupl1)<-exprdat_filterByExpr_dupl1$ENTREZID
    #Remove columns containing ENSEMBL and ENTREZ IDs
    exprdat_filterByExpr_dupl1<-subset(exprdat_filterByExpr_dupl1, select=-c(Row.names,ENTREZID))
    dim(exprdat_filterByExpr_dupl1)
    # we particularly see that the number of columns is now at its initial number 
    
    
    #############################################
    ### final converted gene expression data set: 
    #############################################
    
    exprdat_filterByExpr_dupl1
    
    # note: exprdat_filterByExpr_dupl1 is NOT in the form that is required by edgeR and limma: 
    class(exprdat_filterByExpr_dupl1)
    # wherease edgeR and limma require a DGEList object 
    
    # We will therefore generate a DGEList object in the R script Instructions_Differential_Expression_Analysis.R
    
    #############################################
    #############################################
    
    
    
    
    ############
    ### option 2: keep (rounded) mean expression value of all duplicated gene IDs
    ############
    
    
    # here, we switch order and remove duplicated Entrez gene IDs (case 2) before removing 
    # duplicated Ensembl gene IDs (case 1)
    # the reason for this is elaborated below 
    
    #1.remove duplicated ENTREZ gene IDs (case 2)
    #i.e. multiple different ENSEMBL IDs that are mapped to the same single ENTREZ ID
    
    #generate matrix to contain (rounded) mean expression values of all rows that 
    #have same ENTREZ gene ID
    #ncol=ncol(expression_data_filterByExpr)-2 since data set contains 2 columns with IDs at this point
    mean_entrez<-matrix(, nrow=0, ncol=ncol(counts_expression_data_filterByExpr)-2)
    
    
    # for each duplicated Entrez gene ID separately, we gather all rows with the corresponding gene expression data
    # and then extract the (rounded) mean expression value of all rows 
    for(i in 1:length(dupl_entrez)){
      
      #go through each ENTREZ IDs which occurs multiple times
      #determine all rows whose ENTREZ IDs correspond to currently considered ENTREZ ID 
      counts_dupl<-counts_expression_data_filterByExpr[counts_expression_data_filterByExpr$ENTREZID %in% unique(dupl_entrez)[i],]
      
      #compute the mean expression value of all rows that contain to the given Entrez gene ID 
      dupl_id<-round(colMeans(counts_dupl[,c(2:(ncol(counts_expression_data_filterByExpr)-1))]))
      #store rounded mean expression value in matrix 
      # this matrix is extended by a single row of gene expression data which corresponds to the 
      # (rounded) mean expression data that corresponds to the given Entrez gene ID 
      mean_entrez<-rbind(mean_entrez,dupl_id)
      
      
    }
    
    # after completing the for-loop, mean_entrez contains the mean expression measurements of each Entrez gene ID
    # which contains duplicates resulting from gene ID conversion
    mean_entrez
    
    #set corresponding ENTREZ gene IDs as rownames
    rownames(mean_entrez)<-unique(dupl_entrez)
    
    
    # test whether the number of rows in mean_entrez corresponds to the number ENTREZ IDs
    # that occur more than once 
    # result should be TRUE 
    nrow(mean_entrez)==length(dupl_entrez)
    
    # remove all rows from the expression data whose ENTREZ ID has at least one duplicate
    # intuition: we have just dealt with the corresponding rows and do not want them to be considered
    # in the second step (which deals with case 2
    
    exprdat_filterByExpr_dupl2<-counts_expression_data_filterByExpr[!counts_expression_data_filterByExpr$ENTREZID %in% dupl_entrez,]
    
    # test whether number of rows in resulting data set equals nrow of inital data set 
    # minus number of genes with at least one duplicate
    nrow(exprdat_filterByExpr_dupl2)==nrow(counts_expression_data_filterByExpr)-nrow(duplicated_conversion_entrez)
    dim(exprdat_filterByExpr_dupl2)
    
 
    
    #2. remove duplicated ENSEMBL IDs
    #caution: single ENSEMBL IDs that are mapped to multiple ENTREZ ID naturally generate
    #identical count data for all corresponding ENTREZ IDs
    #->pointless to compute mean expression values
    #verifiable by looking at data set only containing those ENSEMBL IDs that are
    #mapped by multiple ENTREZ IDs:
    #test_dupl_ensembl<-expression_data_filterByExpr[expression_data_filterByExpr$Row.names %in% dupl_ensembl,]
    #View(test_dupl_ensembl)
    
    #therefore: proceed as in option 1 and use ENTREZ ID that occurs first, remove the rest
    exprdat_filterByExpr_dupl2<-exprdat_filterByExpr_dupl2[!duplicated(exprdat_filterByExpr_dupl2$Row.names),]
    dim(exprdat_filterByExpr_dupl2)
    #set ENTREZ ID as rownames
    rownames(exprdat_filterByExpr_dupl2)<-exprdat_filterByExpr_dupl2$ENTREZID
    #remove any columns containing IDs
    exprdat_filterByExpr_dupl2<-subset(exprdat_filterByExpr_dupl2,select= -c(Row.names,ENTREZID))
    #add rows to data set that contain mean expression values of duplicate ENTREZ IDs
    exprdat_filterByExpr_dupl2<-rbind(exprdat_filterByExpr_dupl2,mean_entrez)
    #dimension of remaining expression data set:
    #dim(exprdat_dupl)
    
    
    #############################################
    ### final converted gene expression data set: 
    #############################################
    
    exprdat_filterByExpr_dupl2
    
    
    # note: exprdat_filterByExpr_dupl1 is NOT in the form that is required by edgeR and limma: 
    class(exprdat_filterByExpr_dupl2)
    # wherease edgeR and limma require a DGEList object 
    
    # We will therefore generate a DGEList object in the R script Instructions_Differential_Expression_Analysis.R
    
    #############################################
    ############################################
    
    
    
    
    ############
    ### option 3: among duplicates, keep row with highest overall expression values (i.e highest counts across all samples)
    ############
    

      #intuition: row with highest counts values has  highest power of detecting
      #differential expression 
      #as in option 2, this applies only to duplicates that result from multiple ENSEMBL IDs
      #that are mapped to the same ENTREZ ID
      
      
      #case 2: (case 1 below) multiple ENSEMBL IDs that are converted to the same single ENTREZ ID
      
      # generate matrix to later contain row with highest count values among ID duplicates
      # this data set is to be filled gradually and with each iteration of the follwing for-loop
      highest_count_entrez<-matrix(, nrow=0, ncol=ncol(expression_data_filterByExpr))
    
    
    
    
      # for each duplicated Entrez gene ID separately, we gather all rows with the corresponding gene expression data
      # and then extract the row with the highest overall magnitude of counts 
    
      for(i in 1:length(dupl_entrez)){
        
        # go through each ENTREZ IDs which occurs multiple times
        # determine all rows whose ENTREZ IDs correspond to currently considered ENTREZ ID 
        counts_dupl<-counts_expression_data_filterByExpr[counts_expression_data_filterByExpr$ENTREZID %in% unique(dupl_entrez)[i],]
        
        # order rows in decreasing manner by their number of read counts across all samples 
        order_rowsums<-order(rowSums(counts_dupl[,2:(ncol(counts_dupl)-1)]),decreasing=TRUE)
        #detect row with highest number of read counts across all samples (i.e. row with rank 1)
        dupl_id<-counts_dupl[order_rowsums==1,]
        #store corresponding expression 
        highest_count_entrez<-rbind(highest_count_entrez,dupl_id)
        #View(highest_count_entrez)
        #remove rows in counts_dupl from count data set successively
      }
      
    
      #Remove all initial values with ENTREZ duplicates from the dataset initial gene expression data set 
      exprdat_filterByExpr_dupl3<-counts_expression_data_filterByExpr[! counts_expression_data_filterByExpr$ENTREZID %in% unique(dupl_entrez),]
      
      
      #case 1: single ENSEMBL ID that is mapped to multiple ENTREZ gene IDs 
      #as in option 2, pointless to detect row with highest count values as all rows
      #corresponding to the same ENSEMBL ID naturally contain identical count data
      #therefore: remove duplicate ENSEMBL ID that occurs first 
      exprdat_filterByExpr_dupl3<-exprdat_filterByExpr_dupl3[!duplicated(exprdat_filterByExpr_dupl3$Row.names),]
      
      # add gene expression rows of all Entrez gene IDs that were initially duplicated
      exprdat_filterByExpr_dupl3<-rbind(exprdat_filterByExpr_dupl3,highest_count_entrez )
      
      #Set ENTREZ IDs as rownames 
      rownames(exprdat_filterByExpr_dupl3)<-exprdat_filterByExpr_dupl3$ENTREZID
      #Remove any column that information on contains gene IDs
      exprdat_filterByExpr_dupl3<-subset(exprdat_filterByExpr_dupl3, select=-c(Row.names,ENTREZID))
      dim(exprdat_filterByExpr_dupl3)
      # we see that the sample size (number of columns) is now back at its initial number 
      
      
      #############################################
      ### final converted gene expression data set: 
      #############################################
      
      exprdat_filterByExpr_dupl3
      
      # note: exprdat_filterByExpr_dupl1 is NOT in the form that is required by edgeR and limma: 
      class(exprdat_filterByExpr_dupl3)
      # wherease edgeR and limma require a DGEList object 
      
      # We will therefore generate a DGEList object in the R script Instructions_Differential_Expression_Analysis.R

      
      
      
      
      
      
      
      
      
      
      
      
      
      ################################################################################
      # (i) expression_data_filterDESeq2 #############################################
      ################################################################################
      
      # inspect data set: 
      # View(expression_data_filterDESeq2)
      dim(expression_data_filterDESeq2)
      
      # ->> this expression data set contains gene expression measurements of 52580 genes and 69 samples 
      
      # ->> goal: convert Ensembl gene IDs to Entrez gene IDs 
      
      # here, we work with functionalities provided by clusterProfiler 
    
      
      
      
      #############################
      # (I) Conversion of gene IDs
      #############################

      
      
      #######################################################################################
      ### step 1: indicate which organism the gene expression measurements originate from ###
      #######################################################################################
      
      
      
      # assumption 1: gene expression measurements originate from humans 
      # -> this corresponds to argument OrgDb = org.Hs.eg.db
      
      # IMPORTANT: library "org.Hs.eg.db", which contains the corresponding genome wide annotation,  
      # must be loaded if the given organism is human: 
      library(org.Hs.eg.db)
      
      # note: if, for instance, the organism is mouse, then the library "org.Mm.eg.db" must be loaded via 
      # library(org.Mm.eg.db)
      # and the argument OrgDb must be set to OrgDb = org.Mm.eg.db
      
      # assumption 2: genes originate from Ensembl gene ID format 
      # -> this corresponds to the argument fromType = "ENSEMBL"
      
      # assumption 3: gene set analysis tool requires genes to be in Entrez gene ID format 
      # -> this corresponds to argument toType = "ENTREZID" 
      
      
      # NOTE: arguments "fromType" and "toType" must be set as one of the following, depending on the given and the 
      # required gene ID format
      # ACCNUM, ALIAS, ENSEMBL, ENSEMBLPROT, ENSEMBLTRANS, ENTREZID, ENZYME, EVIDENCE, 
      # EVIDENCEALL, GENENAME, GENETYPE, GO, GOALL, IPI, MAP, OMIM, ONTOLOGY, ONTOLOGYALL, PATH, PFAM, PMID, PROSITE, 
      # REFSEQ, SYMBOL, UCSCKG, UNIPROT.
      
      
      
      
      ##############################################################
      ### step 2: Gene ID conversion via clusterProfiler::bitr() ###
      ##############################################################
      
      
      
      # obtain mapping of initial (Ensembl) gene IDs to required (Entrez) gene IDs: 
      bitr_EnsToEntr <- bitr(rownames(expression_data_filterDESeq2), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
      #note: not all ENSEMBL IDs could be converted to a corresponding ENTREZ Gene ID: 
      dim(bitr_EnsToEntr)
      
      #results: 
      #- not all Ensembl gene IDs can be mapped to a corresponding Entrez gene ID
      #- some individual Ensembl gene IDs were mapped to multiple distinct Entrez gene IDs
      #- some distinct Ensembl gene IDs were mapped to an identical Entrez gene ID 
      
      ### merge
      
      # concatenate initial gene expression data set with the mapping from initial (Ensembl) to required (Entrez) gene ID format
      # merge by row names of expression data set and ENSEMBL ID of gene ID mapping
      counts_expression_data_filterDESeq2 <- merge(expression_data_filterDESeq2, bitr_EnsToEntr, by.x=0, by.y="ENSEMBL", all.y=TRUE, sort=TRUE)
      
      # note that we will work with counts_expression_data_filterDESeq2 at this point and will set the content of this to 
      # expression_data_filterDESeq2 in the end
      
      
      # note on function arguments: 
      # by.x and by.y specify by which columns expression_data_filterDESeq2 and bitr_EnsToEntr are concatenated: 
      # by.x = 0: use row names of expression_data_filterDESeq2 (the row names contain the gene IDs)
      # by.y = "ENSEMBL": use the column that contains the Ensembl gene IDs 
      
      
      # take a look at dimension of resulting data set
      dim(counts_expression_data_filterDESeq2)
      # - number of genes in counts_expression_data_filterByExpr has decreased to 9537
      # -> this is directly caused by the circumstance that typically, not all initial (Ensembl) gene IDs can be 
      # converted to the required (Entrez) gene ID 
      # - the expression data set has been extended by two columns:
      # (1.) the very last column (here called "ENTREZID" containing the converted (Entrez) gene ID
      # (2.) the very first column (here called Row.names) with the initial (Ensembl) gene IDs that were originally stored in the row names
      
      # note: in the resulting data set, the row names now correspond to the row number     
      ##############################################
      ### step 3: take closer look at duplicates ###
      ##############################################
      
      ########
      # CASE 1: single ENSEMBL IDs are mapped to multiple ENTREZ IDs
      ########
      
      # obtain number of cases in which an ENSEMBL gene ID was converted to several ENTREZ IDs, i.e. the number of 
      # times an Ensembl ID appears repeatedly in the mapping
      sum(duplicated(bitr_EnsToEntr$ENSEMBL)) 
      # determine all duplicated ENSEMBL gene IDS(i.e. all ENSEMBL gene IDs that were mapped to multiple distinct Entrez gene IDs):
      dupl_ensembl<-unique(bitr_EnsToEntr$ENSEMBL[duplicated(bitr_EnsToEntr$ENSEMBL)])
      # number of ENSEMBL IDs that have at least one duplicate
      length(dupl_ensembl)
      
      
      
      # in the following, we can inspect the conversion pattern for each Ensembl ID that was mapped to 
      # multiple distinct Entrez IDs 
      duplicated_conversion_ens<-bitr_EnsToEntr[bitr_EnsToEntr$ENSEMBL %in% dupl_ensembl,]
      # take a look at conversion pattern for duplicated Ensembl gene IDs:
      duplicated_conversion_ens
      # -> we can see perfectly that these Ensembl ID are each mapped to multiple individual Entrez IDs 
      
      
      ########
      # CASE 2: multiple ENSEMBL IDs are mapped to single entrez ID
      ########
      
      # obtain number of cases in which multiple distinct Ensembl IDs are converted to the same Entrez gene ID 
      # in these cases, the corresponding Entrez IDs appear repeatedly in the mapping:
      sum(duplicated(bitr_EnsToEntr$ENTREZ)) 
      # determine all Entrez gene IDs affected by this duplication
      dupl_entrez<-unique(bitr_EnsToEntr$ENTREZID[duplicated(bitr_EnsToEntr$ENTREZID)])
      # number of ENTREZ IDs that have at least one duplicate
      length(dupl_entrez)
      # display of conversion pattern of duplicated ENTREZ IDs
      duplicated_conversion_entrez<-bitr_EnsToEntr[bitr_EnsToEntr$ENTREZID %in% dupl_entrez,]
      # take a look at conversion pattern: 
      duplicated_conversion_entrez
      # if ordered by Entrez IDs, we can see perfectly that for each of these Entrez gene IDs, there are multiple
      # corresponding Ensembl gene IDs 
      dim(duplicated_conversion_entrez)
      
      
      
      ##############################################
      ### step 4: remove duplicated gene IDs  ######
      ##############################################
      
      
      # for subsequent analysis (whether differential expression analysis of gene set analysis), we need to deal 
      # with the duplicated gene IDs and remove them in a suitable way such that only one unique gene expression 
      # measurement among the duplicates remains
      
      # there is no recommended way to proceed, i.e. no common approach presented in an official scientific puplication, 
      # but instead several approaches suggested by users in corresponding user platforms, 
      # generally, the manner of duplicate gene ID removal is performed at the discretion of the user 
      
      # we therefore present three approaches to duplicate gene ID removal 
      
      
      
      ############
      ### option 1: keep first subscript among duplicates 
      ############
      
      # this is the simplest approach among the three 
      
      
      #1. remove duplicated ENTREZ gene IDs
      counts_filterDESeq2_dupl1<-counts_expression_data_filterDESeq2[!duplicated(counts_expression_data_filterDESeq2$ENTREZID),]
      
      
      #2. remove duplicated ENSEMBL gene IDs
      counts_filterDESeq2_dupl1<-counts_filterDESeq2_dupl1[!duplicated(counts_filterDESeq2_dupl1$Row.names),]
      dim(counts_filterDESeq2_dupl1)
      
      #3. set Entrez gene IDs as rownames 
      rownames(counts_filterDESeq2_dupl1)<-counts_filterDESeq2_dupl1$ENTREZID
      #Remove columns containing ENSEMBL and ENTREZ IDs
      counts_filterDESeq2_dupl1<-subset(counts_filterDESeq2_dupl1, select=-c(Row.names,ENTREZID))
      dim(counts_filterDESeq2_dupl1)
      # we particularly see that the number of columns is now at its initial number 
      
      
      #############################################
      ### final converted gene expression data set: 
      #############################################
      
      
      counts_filterDESeq2_dupl1
      
      # note that the initial gene expression data set expression_data_filterDESeq2 (that 
      # we obtain from Instructions_PreFiltering.R is in the form of a DESeqDataSet
      
      # in objects of the class DESeqDataSet, the count data, the dimension of the count 
      # data, the rownames (i.e. gene IDs) etc. are stored
      
      # while DESeqDataSet is able to adapt the dimension of the count data (which is 
      # the consequence of e.g. pre-filtering), it gives an error if the rownames or
      # colnames of the supplied assay is not identical to those of the receiving
      # DESeqDataSet object. 
      # -> in our case, a conversion of the gene IDs leads to an alteration of the 
      # rownames, so that we CANNOT simply provide the new count data to the DESeqDataSet
      # with the following code line:
      
      # expression_data_filterDESeq2 <- as.matrix(counts_filterDESeq2_dupl1)
      
      # THEREFORE: we need to create a new DESeqDataSet in the following scripts 
      # (such as Instructions_Differential_Expression_Analysis.R)
      
      
      #############################################
      #############################################
      
      
      
      
      ############
      ### option 2: keep (rounded) mean expression value of all duplicated gene IDs
      ############
      
      
      # here, we switch order and remove duplicated Entrez gene IDs (case 2) before removing 
      # duplicated Ensembl gene IDs (case 1)
      # the reason for this is elaborated below 
      
      #1.remove duplicated ENTREZ gene IDs (case 2)
      #i.e. multiple different ENSEMBL IDs that are mapped to the same single ENTREZ ID
      
      #generate matrix to contain (rounded) mean expression values of all rows that 
      #have same ENTREZ gene ID
      #ncol=ncol(expression_data_filterDESeq2)-2 since data set contains 2 columns with IDs at this point
      mean_entrez<-matrix(, nrow=0, ncol=ncol(counts_expression_data_filterDESeq2)-2)
      
      
      # for each duplicated Entrez gene ID separately, we gather all rows with the corresponding gene expression data
      # and then extract the (rounded) mean expression value of all rows 
      for(i in 1:length(dupl_entrez)){
        
        #go through each ENTREZ IDs which occurs multiple times
        #determine all rows whose ENTREZ IDs correspond to currently considered ENTREZ ID 
        counts_dupl <- counts_expression_data_filterDESeq2[counts_expression_data_filterDESeq2$ENTREZID %in% unique(dupl_entrez)[i],]
        
        #compute the mean expression value of all rows that contain to the given Entrez gene ID 
        dupl_id<-round(colMeans(counts_dupl[,c(2:(ncol(counts_expression_data_filterDESeq2)-1))]))
        #store rounded mean expression value in matrix 
        # this matrix is extended by a single row of gene expression data which corresponds to the 
        # (rounded) mean expression data that corresponds to the given Entrez gene ID 
        mean_entrez<-rbind(mean_entrez,dupl_id)
        
        
      }
      
      # after completing the for-loop, mean_entrez contains the mean expression measurements of each Entrez gene ID
      # which contains duplicates resulting from gene ID conversion
      mean_entrez
      
      #set corresponding ENTREZ gene IDs as rownames
      rownames(mean_entrez)<-unique(dupl_entrez)
      
      
      # test whether the number of rows in mean_entrez corresponds to the number ENTREZ IDs
      # that occur more than once 
      # result should be TRUE 
      nrow(mean_entrez)==length(dupl_entrez)
      
      # remove all rows from the expression data whose ENTREZ ID has at least one duplicate
      # intuition: we have just dealt with the corresponding rows and do not want them to be considered
      # in the second step (which deals with case 2
      
      counts_filterDESeq2_dupl2<-counts_expression_data_filterDESeq2[!counts_expression_data_filterDESeq2$ENTREZID %in% dupl_entrez,]
      
      # test whether number of rows in resulting data set equals nrow of inital data set 
      # minus number of genes with at least one duplicate
      nrow(counts_filterDESeq2_dupl2)==nrow(counts_expression_data_filterDESeq2)-nrow(duplicated_conversion_entrez)
      dim(counts_filterDESeq2_dupl2)
      
      
      
      #2. remove duplicated ENSEMBL IDs
      #caution: single ENSEMBL IDs that are mapped to multiple ENTREZ ID naturally generate
      #identical count data for all corresponding ENTREZ IDs
      #->pointless to compute mean expression values
      #verifiable by looking at data set only containing those ENSEMBL IDs that are
      #mapped by multiple ENTREZ IDs:
      #test_dupl_ensembl<-expression_data_filterDESeq2[expression_data_filterDESeq2$Row.names %in% dupl_ensembl,]
      #View(test_dupl_ensembl)
      
      #therefore: proceed as in option 1 and use ENTREZ ID that occurs first, remove the rest
      counts_filterDESeq2_dupl2<-counts_filterDESeq2_dupl2[!duplicated(counts_filterDESeq2_dupl2$Row.names),]
      dim(counts_filterDESeq2_dupl2)
      #set ENTREZ ID as rownames
      rownames(counts_filterDESeq2_dupl2)<-counts_filterDESeq2_dupl2$ENTREZID
      #remove any columns containing IDs
      counts_filterDESeq2_dupl2<-subset(counts_filterDESeq2_dupl2,select= -c(Row.names,ENTREZID))
      #add rows to data set that contain mean expression values of duplicate ENTREZ IDs
      counts_filterDESeq2_dupl2<-rbind(counts_filterDESeq2_dupl2,mean_entrez)
      #dimension of remaining expression data set:
      #dim(exprdat_dupl)
      
      
      #############################################
      ### final converted gene expression data set: 
      #############################################
      
      counts_filterDESeq2_dupl2
      
      
      # note that the initial gene expression data set expression_data_filterDESeq2 (that 
      # we obtain from Instructions_PreFiltering.R is in the form of a DESeqDataSet
      
      # in objects of the class DESeqDataSet, the count data, the dimension of the count 
      # data, the rownames (i.e. gene IDs) etc. are stored
      
      # while DESeqDataSet is able to adapt the dimension of the count data (which is 
      # the consequence of e.g. pre-filtering), it gives an error if the rownames or
      # colnames of the supplied assay is not identical to those of the receiving
      # DESeqDataSet object. 
      # -> in our case, a conversion of the gene IDs leads to an alteration of the 
      # rownames, so that we CANNOT simply provide the new count data to the DESeqDataSet
      # with the following code line:
      
      # expression_data_filterDESeq2 <- as.matrix(counts_filterDESeq2_dupl1)
      
      # THEREFORE: we need to create a new DESeqDataSet in the following scripts 
      # (such as Instructions_Differential_Expression_Analysis.R)
      
    
      
      #############################################
      ############################################
      
      
      
      
      ############
      ### option 3: among duplicates, keep row with highest overall expression values (i.e highest counts across all samples)
      ############
      
      
      #intuition: row with highest counts values has  highest power of detecting
      #differential expression 
      #as in option 2, this applies only to duplicates that result from multiple ENSEMBL IDs
      #that are mapped to the same ENTREZ ID
      
      
      #case 2: (case 1 below) multiple ENSEMBL IDs that are converted to the same single ENTREZ ID
      
      # generate matrix to later contain row with highest count values among ID duplicates
      # this data set is to be filled gradually and with each iteration of the follwing for-loop
      highest_count_entrez<-matrix(, nrow=0, ncol=ncol(counts_expression_data_filterDESeq2))
      
      
      
      
      # for each duplicated Entrez gene ID separately, we gather all rows with the corresponding gene expression data
      # and then extract the row with the highest overall magnitude of counts 
      
      for(i in 1:length(dupl_entrez)){
        
        # go through each ENTREZ IDs which occurs multiple times
        # determine all rows whose ENTREZ IDs correspond to currently considered ENTREZ ID 
        counts_dupl<-counts_expression_data_filterDESeq2[counts_expression_data_filterDESeq2$ENTREZID %in% unique(dupl_entrez)[i],]
        
        # order rows in decreasing manner by their number of read counts across all samples 
        order_rowsums<-order(rowSums(counts_dupl[,2:(ncol(counts_dupl)-1)]),decreasing=TRUE)
        #detect row with highest number of read counts across all samples (i.e. row with rank 1)
        dupl_id<-counts_dupl[order_rowsums==1,]
        #store corresponding expression 
        highest_count_entrez<-rbind(highest_count_entrez,dupl_id)
        #View(highest_count_entrez)
        #remove rows in counts_dupl from count data set successively
      }
      
      
      #Remove all initial values with ENTREZ duplicates from the dataset initial gene expression data set 
      counts_filterDESeq2_dupl3<-counts_expression_data_filterDESeq2[!counts_expression_data_filterDESeq2$ENTREZID %in% unique(dupl_entrez),]
      
      
      #case 1: single ENSEMBL ID that is mapped to multiple ENTREZ gene IDs 
      #as in option 2, pointless to detect row with highest count values as all rows
      #corresponding to the same ENSEMBL ID naturally contain identical count data
      #therefore: remove duplicate ENSEMBL ID that occurs first 
      counts_filterDESeq2_dupl3<-counts_filterDESeq2_dupl3[!duplicated(counts_filterDESeq2_dupl3$Row.names),]
      
      # add gene expression rows of all Entrez gene IDs that were initially duplicated
      counts_filterDESeq2_dupl3<-rbind(counts_filterDESeq2_dupl3,highest_count_entrez )
      
      #Set ENTREZ IDs as rownames 
      rownames(counts_filterDESeq2_dupl3)<-counts_filterDESeq2_dupl3$ENTREZID
      #Remove any column that information on contains gene IDs
      counts_filterDESeq2_dupl3<-subset(counts_filterDESeq2_dupl3, select=-c(Row.names,ENTREZID))
      dim(counts_filterDESeq2_dupl3)
      # we see that the sample size (number of columns) is now back at its initial number 
      
      
      #############################################
      ### final converted gene expression data set: 
      #############################################
      
      counts_filterDESeq2_dupl3
      
      
      
      # note that the initial gene expression data set expression_data_filterDESeq2 (that 
      # we obtain from Instructions_PreFiltering.R is in the form of a DESeqDataSet
      
      # in objects of the class DESeqDataSet, the count data, the dimension of the count 
      # data, the rownames (i.e. gene IDs) etc. are stored
      
      # while DESeqDataSet is able to adapt the dimension of the count data (which is 
      # the consequence of e.g. pre-filtering), it gives an error if the rownames or
      # colnames of the supplied assay is not identical to those of the receiving
      # DESeqDataSet object. 
      # -> in our case, a conversion of the gene IDs leads to an alteration of the 
      # rownames, so that we CANNOT simply provide the new count data to the DESeqDataSet
      # with the following code line:
      
      # expression_data_filterDESeq2 <- as.matrix(counts_filterDESeq2_dupl1)
      
      # THEREFORE: we need to create a new DESeqDataSet in the following scripts 
      # (such as Instructions_Differential_Expression_Analysis.R)
      
      #############################################
      ############################################      
      
      
      
      
#######################
### step 5: Save ######
#######################
      
# for the purpose of illustration of the entire GSA workflow, we will proceed with a single
# pre-filtered and converted gene expression data set using filterByExpr() (for edgeR and limma) 
# and for DESeq2 
      
# to keep it simple, we will proceed with the respective first data set for edgeR and DESeq2
# however: you can witch around at your discretion. 
      
# specify pre-filtered and converted gene expression data set for edgeR and limma
exprdat_filter_conv_filterByExpr <- exprdat_filterByExpr_dupl1
# specify pre-filtered and converted gene expression data for DESeq2 
exprdat_filter_conv_DESeq2 <- counts_filterDESeq2_dupl1

save(exprdat_filter_conv_filterByExpr, file = "./Results_GeneID_Conversion_DuplicatesRemoval/exprdat_filter_conv_filterByExpr.Rdata")
save(exprdat_filter_conv_DESeq2, file = "./Results_GeneID_Conversion_DuplicatesRemoval/exprdat_filter_conv_DESeq2.Rdata")
      



  
  
  