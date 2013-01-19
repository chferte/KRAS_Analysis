# Charles Ferté
# Sage Bionetworks
# 14 Sept 2012


# train a model of differential gene expression per KRAS codon in Lung Adenocarcinoma
# for G12C,  in TCGA LUAD

#load the different packages
options(stringsAsFactors=FALSE)

library(affy)
library(corpcor)
library(lattice)
library(limma)
library(caret)
library(glmnet)
library(snm)
source("/home/cferte/FELLOW/cferte/KRAS_Project/JUSTIN_PREDICT_CCLE/code/lung_analysis_functions.R")
synapseLogin("charles.ferte@sagebase.org","charles")

###############################################################
# load the CCLE data (snm normalized) 
###############################################################


liste <- list.celfiles(path="/home/bmecham/ccle/celFiles",full.names=T)
map <- read.delim("/home/cferte/FELLOW/cferte/KRAS_Project/CCLE/CCLE_map.txt",header=T)


map <- map[ map$Site.Primary =="lung",]

map <- map[ map$Hist.Subtype1 !="small_cell_carcinoma",]
liste <- liste[match(paste(map$ID,".CEL",sep=""),substr(liste,29,nchar(liste)))]
liste <- liste[!is.na(liste)]
ccle_EXP <- ReadAffy(filenames=liste)
sampleNames(ccle_EXP) <- map$CCLE.name[which(!is.na(liste))]
sampleNames(ccle_EXP) <-sapply(strsplit(split="_",x=sampleNames(ccle_EXP)),function(x){x[[1]]})




##################################################################################
# rma
##################################################################################

ccle_rma <- rma(ccle_EXP)

# retrieve the feature names 

library(hgu133plus2.db)
tmp <- unlist(mget(x=featureNames(ccle_rma),hgu133plus2SYMBOL,ifnotfound=NA))

combine_probes_2_gene <- function(expr, genes, method="svd"){
  
  if(is.list(genes)) genes <- unlist(genes)
  
  stopifnot(dim(expr)[1] ==  length(genes))
  ugenes <- unique(genes)
  ugenes <- sort(ugenes[!is.na(ugenes)])
  M <- matrix(NaN, ncol=dim(expr)[2],nrow=length(ugenes),
              dimnames=list(ugenes, colnames(expr)))
  
  for(gene in ugenes){
    sub.expr <- as.matrix(expr[which(genes == gene),])
    if(dim(sub.expr)[2] == 1){
      M[gene,] <- sub.expr
    }else{
      tmp <- svd(sub.expr - rowMeans(sub.expr))$v[,1]
      tmp.c <- mean(cor(tmp, t(sub.expr)))
      #cat(gene," ", tmp.c, "\n")
      multiplier <- ifelse(tmp.c < 0, -1, 1)
      M[gene,] <- tmp * multiplier
    }
  }
  M
}

ccle_rma <- combine_probes_2_gene(expr=ccle_rma,genes=tmp)
colnames(ccle_rma) <- sampleNames(ccle_EXP)
CCLE_RMA <- ccle_rma
rm(ccle_rma,ccle_eset,map,mapCdfName)


# 
# # save ccle rma
# ccle_rma <- Data(list(name = "CCLE_RMA", parentId = 'syn1589868'))
# ccle_rma <- createEntity(ccle_rma)
# 
# # add object into the data entity
# ccle_rma <- addObject(ccle_rma,CCLE_RMA)
# 
# # push the raw data into this entity
# ccle_rma <- storeEntity(entity=ccle_rma)
# Load the data 
CCLE_RMA <- loadEntity('syn1589873')
CCLE_RMA <- CCLE_RMA$objects$CCLE_RMA

# load the ccle drug information and make it coherent with ccle_EXP
ccle_EXP <- CCLE_RMA
rm(CCLE_RMA)
ccle_drug <- loadEntity('syn1354656')
ccle_drug <- ccle_drug$objects$ccle_drug
ccle_drug <- ccle_drug@data
rownames(ccle_drug) <- sapply(strsplit(split="_",x=rownames(ccle_drug)),function(x){x[[1]]})
ccle_drug <- ccle_drug[colnames(ccle_EXP),]
