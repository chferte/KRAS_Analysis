# Charles Ferté
# Sage Bionetworks
# 5th January 2013


# sanger data input

#load the different packages
library(affy)
library(corpcor)
library(lattice)
library(limma)
library(caret)
library(glmnet)
library(snm)
library(synapseClient)
options(stringsAsFactors=FALSE)
synapseLogin("charles.ferte@sagebase.org","charles")

#####################################################################################################################
# sanger gene expression
#####################################################################################################################

sanger_exp <- loadEntity("syn1417725")
sanger_exp <- sanger_exp$objects$eset

library(org.Hs.eg.db)
tmp <- unlist(mget(x=sub(pattern="_mt",replacement="",x=featureNames(sanger_exp)),org.Hs.egSYMBOL,ifnotfound=NA))

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

sanger_exp1 <- combine_probes_2_gene(expr=sanger_exp,genes=tmp)
colnames(sanger_exp1) <- sampleNames(sanger_exp)
sanger_exp <- sanger_exp1
rm(sanger_exp1,tmp)

#####################################################################################################################
# input sanger mutations data
#####################################################################################################################
sanger_mut <- read.csv("/home/cferte/cell_line_data/Sanger_gdsc_mutation_w3.csv",header=TRUE)
rownames(sanger_mut) <- sanger_mut$Cell.Line
sanger_mut <- sanger_mut[, -which(colnames(sanger_mut) %in% c("Cell.Line","Cancer.Type","Cosmic_ID","Tissue"))]
sanger_mut <- as.matrix(sanger_mut)
sanger_mut[grep(pattern="wt",x=sanger_mut)] <- "0"
sanger_mut[grep(pattern="na",x=sanger_mut)] <- "NA"
sanger_mut[sanger_mut==""] <- "NA"
sanger_mut[!sanger_mut %in% c("0","NA")] <- "1"
colnames(sanger_mut) <- toupper(colnames(sanger_mut))
rownames(sanger_mut) <- toupper(rownames(sanger_mut))
rownames(sanger_mut) <- gsub(pattern="-",x=rownames(sanger_mut),replacement="")

#####################################################################################################################
# input the sanger drug response (IC 50 only)
#####################################################################################################################
sanger_drug1 <- read.delim("/home/cferte/cell_line_data/sanger_drugs_pharma.txt",header=TRUE,na.strings = c("","NA"))
rownames(sanger_drug) <- sanger_drug$Cell.Line
rownames(sanger_drug) <- gsub(pattern="-",replacement="",x=toupper(rownames(sanger_drug)))
sanger_drug <- sanger_drug[,grep(pattern="IC_50",x=colnames(sanger_drug))]
sanger_drug <- sanger_drug[,which(substr(x=colnames(sanger_drug),start=nchar(colnames(sanger_drug))-4,stop=nchar(colnames(sanger_drug)))=="IC_50")]


#####################################################################################################################
# input the sanger cnv
#####################################################################################################################
sanger_cnv <- loadEntity("syn1417763")
sanger_cnv <- sanger_cnv$objects$eset
tmp <- unlist(mget(x=sub(pattern="_eg",replacement="",x=featureNames(sanger_cnv)),org.Hs.egSYMBOL,ifnotfound=NA))
sanger_cnv1 <- combine_probes_2_gene(expr=sanger_cnv,genes=tmp)
colnames(sanger_cnv1) <- sampleNames(sanger_cnv)
sanger_cnv <- sanger_cnv1
rm(sanger_cnv1,tmp)

#####################################################################################################################
## input the sanger info
#####################################################################################################################
 
sanger_info <- read.delim("/home/cferte/cell_line_data/sanger_drugs_pharma.txt",header=TRUE,na.strings = c("","NA"))
rownames(sanger_info) <- sanger_info$Cell.Line
sanger_info <- sanger_info[,1:4]
sanger_info$cell_line_id <- toupper(gsub(pattern="-",replacement="",x=paste(sanger_info$Cell.Line,"_",sanger_info$Tissue,sep="")))

colnames(sanger_info)

#####################################################################################################################
# input the sanger drug info
#####################################################################################################################
sanger_drug_info <- read.delim("/home/cferte/cell_line_data/Sanger_drug_info.txt",header=TRUE,na.strings = c("","NA"))


#####################################################################################################################
# make the sampleNames coherent between mut exp and mut
#####################################################################################################################
tmp1 <- sapply(strsplit(x=colnames(sanger_exp),split="_"),function(x){x[[1]]})
tmp <- match(rownames(sanger_mut),tmp1)
rownames(sanger_mut)  <- colnames(sanger_exp)[tmp]


#####################################################################################################################
# make a list of all objects and save it into synapse
#####################################################################################################################

sanger_data <- list(sanger_info=sanger_info,sanger_exp=sanger_exp,sanger_cnv=sanger_cnv,sanger_mut=sanger_mut,sanger_drug=sanger_drug,sanger_drug_info=sanger_drug_info)

sanger_all <- Data(list(name = "sanger_all", parentId = 'syn1670945'))
sanger_all <- createEntity(sanger_all)

# add object into the data entity
sanger_all <- addObject(sanger_all,sanger_data)

# push the raw data into this entity
sanger_all <- storeEntity(entity=sanger_all)

