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

## sanger gene expression
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

# input sanger mutations data
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

## input the sanger drug response
sanger_drug <- read.csv("/home/cferte/cell_line_data/gdsc_manova_output_w2.csv",header=TRUE)

drug <- matrix(NA,nrow=length(unique(sanger_drug$sanger_Name)),ncol=length(unique(sanger_drug$Compound)))
colnames(drug) <- unique(sanger_drug$Compound)
rownames(drug) <- unique(sanger_drug$sanger_Name)
drug <- drug[-which(rownames(drug)=="b"),]

for(i in rownames(drug)){
  for(k in colnames(drug)){
    tmp <- sanger_drug$ActArea_.raw.[which(sanger_drug$sanger_Name==i & sanger_drug$Compound==k)]
    if(length(tmp)!=0){drug[i,k] <-tmp} 
  }}

sanger_drug <- drug
rm(drug)

## input the sanger cnv
sanger_cnv <- loadEntity("syn1417763")
sanger_cnv <- sanger_cnv$objects$eset
tmp <- unlist(mget(x=sub(pattern="_eg",replacement="",x=featureNames(sanger_cnv)),org.Hs.egSYMBOL,ifnotfound=NA))
sanger_cnv1 <- combine_probes_2_gene(expr=sanger_cnv,genes=tmp)
colnames(sanger_cnv1) <- sampleNames(sanger_cnv)
sanger_cnv <- sanger_cnv1
rm(sanger_cnv1,tmp)

## input the sanger info
sanger_info <- read.delim("/home/cferte/cell_line_data/sanger_sample_info_file_2012-04-06.txt")

# make a list of all objects and save it into synapse
sanger_data <- list(sanger_info=sanger_info,sanger_exp=sanger_exp,sanger_cnv=sanger_cnv,sanger_mut=sanger_mut,sanger_drug=sanger_drug)

#sanger_all <- Data(list(name = "sanger_all", parentId = 'syn1670945'))
#sanger_all <- createEntity(sanger_all)

# add object into the data entity
sanger_all <- addObject(sanger_all,sanger_data)

# push the raw data into this entity
sanger_all <- storeEntity(entity=sanger_all)



# # make the sampleNames coherent between mut exp and cnv
# tmp <- intersect(colnames(sanger_exp),colnames(sanger_cnv))
# tmp <- intersect(tmp,colnames(sanger_mut))
# tmp <- intersect(tmp,rownames(sanger_drug))
# sanger_cnv <- sanger_cnv[,tmp]
# sanger_exp <- sanger_exp[,tmp]
# sanger_mut <- sanger_mut[,tmp]
# 
