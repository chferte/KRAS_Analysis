# Charles Ferté
# Sage Bionetworks
# 5th January 2013


# ccle data input

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

## ccle gene expression
ccle_exp <- loadEntity("syn1417729")
ccle_exp <- ccle_exp$objects$eset

library(org.Hs.eg.db)
tmp <- unlist(mget(x=sub(pattern="_mt",replacement="",x=featureNames(ccle_exp)),org.Hs.egSYMBOL,ifnotfound=NA))
 
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

ccle_exp1 <- combine_probes_2_gene(expr=ccle_exp,genes=tmp)
colnames(ccle_exp1) <- sampleNames(ccle_exp)
ccle_exp <- ccle_exp1
rm(ccle_exp1)

## input ccle hybrid capture
ccle_mut <- read.delim("/home/cferte/cell_line_data/CCLE_hybrid_capture1650_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07.maf",header=TRUE)
ccle_mut <- ccle_mut[which(ccle_mut$Hugo_Symbol!="Unknown"),]
ccle_mut <- ccle_mut[,c("Hugo_Symbol","Tumor_Sample_Barcode","Protein_Change", "Genome_Change")]
ccle_mut$Protein_Change[ccle_mut$Protein_Change==""] <- ccle_mut$Genome_Change[ ccle_mut$Protein_Change==""]
MATMUT<-matrix(0,nrow=length(unique(ccle_mut$Hugo_Symbol)),ncol=length(unique(ccle_mut$Tumor_Sample_Barcode)))
colnames(MATMUT) <- unique(ccle_mut$Tumor_Sample_Barcode)
rownames(MATMUT) <- unique(ccle_mut$Hugo_Symbol)
for(i in rownames(MATMUT)){
  MATMUT[i,c(ccle_mut$Tumor_Sample_Barcode[which(ccle_mut$Hugo_Symbol==i)])] <- c(ccle_mut$Protein_Change[which(ccle_mut$Hugo_Symbol==i)])  
}
ccle_mut <- MATMUT
rm(MATMUT)


## input the ccle drug response as ActArea Norm
ccle_drug <- read.delim2("/home/cferte/cell_line_data/ccle_drug_response_v2.txt",header=TRUE,as.is=TRUE)
drug <- matrix(NA,nrow=length(unique(ccle_drug$CCLE_Name)),ncol=length(unique(ccle_drug$Compound)))
colnames(drug) <- unique(ccle_drug$Compound)
rownames(drug) <- unique(ccle_drug$CCLE_Name)
drug <- drug[-which(rownames(drug)=="b"),]
for(i in rownames(drug)){
  for(k in colnames(drug)){
    tmp <- ccle_drug$ActArea_.norm.[which(ccle_drug$CCLE_Name==i & ccle_drug$Compound==k)]
  if(length(tmp)!=0){drug[i,k] <-tmp} 
}}
ccle_drug_ActAreaNorm <- drug
rm(drug)

## input the ccle drug response as IC50 Norm
ccle_drug <- read.delim2("/home/cferte/cell_line_data/ccle_drug_response_v2.txt",header=TRUE,as.is=TRUE)
drug <- matrix(NA,nrow=length(unique(ccle_drug$CCLE_Name)),ncol=length(unique(ccle_drug$Compound)))
colnames(drug) <- unique(ccle_drug$Compound)
rownames(drug) <- unique(ccle_drug$CCLE_Name)
drug <- drug[-which(rownames(drug)=="b"),]
for(i in rownames(drug)){
  for(k in colnames(drug)){
    tmp <- ccle_drug$IC50_.microM...norm.[which(ccle_drug$CCLE_Name==i & ccle_drug$Compound==k)]
    if(length(tmp)!=0){drug[i,k] <-tmp} 
  }}
ccle_drug_IC50Norm <- drug
rm(drug)


ccle_drug <- list(ccle_drug_IC50Norm=ccle_drug_IC50Norm,ccle_drug_ActAreaNorm=ccle_drug_ActAreaNorm)


## input the ccle cnv
ccle_cnv <- loadEntity("syn1417763")
ccle_cnv <- ccle_cnv$objects$eset
tmp <- unlist(mget(x=sub(pattern="_eg",replacement="",x=featureNames(ccle_cnv)),org.Hs.egSYMBOL,ifnotfound=NA))
ccle_cnv1 <- combine_probes_2_gene(expr=ccle_cnv,genes=tmp)
colnames(ccle_cnv1) <- sampleNames(ccle_cnv)
ccle_cnv <- ccle_cnv1
rm(ccle_cnv1,tmp)

## input the ccle info
ccle_info <- read.delim("/home/cferte/cell_line_data/CCLE_sample_info_file_2012-04-06.txt")
rownames(ccle_info) <- ccle_info$CCLE.name

## input the ccle_drugs_info
ccle_drugs_info <- read.delim(file="/home/cferte/FELLOW/cferte/KRAS_Analysis/drug_sensitivity_inference/CCLE_drugs.txt",header=T,skip=2)

# make a list of all objects 
ccle_data <- list(ccle_info=ccle_info,ccle_exp=ccle_exp,ccle_cnv=ccle_cnv,ccle_mut=ccle_mut,ccle_drug=ccle_drug,ccle_drugs_info=ccle_drugs_info)

# ...and save it into synapse
ccle_all <- Data(list(name = "ccle_all", parentId = 'syn1670945'))
ccle_all <- createEntity(ccle_all)

# add object into the data entity
ccle_all <- addObject(ccle_all,ccle_data)

# push the raw data into this entity
ccle_all <- storeEntity(entity=ccle_all)