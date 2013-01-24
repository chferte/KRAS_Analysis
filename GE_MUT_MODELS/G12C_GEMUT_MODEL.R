# Charles Ferté
# Sage Bionetworks
# 24 Jan 2013


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
# run the analysis for the G12C in the TCGA LUAD  & CHEMORES (G12C vs any other KRAS)
###############################################################

###############################################################
# load the LUAD data  
###############################################################

source("/home/cferte/FELLOW/cferte/KRAS_Analysis/data_input/TCGA_LUAD_input.R")
par(mfrow=c(1,1))

# see if here is some signal with KRAS status
tmp <- ifelse(KRAS_LUAD!="WT",1,0)
fit <- eBayes(lmFit(LUAD_EXP,model.matrix(~tmp)))
hist(fit$p.value[,2],breaks=100, main="TCGA exp ~ KRAS")
rm(tmp)

###############################################################
# load the CCLE data (snm normalized) 
###############################################################

source("/home/cferte/FELLOW/cferte/KRAS_Analysis/data_input/ccle_load_data.R")
names(KRAS_CCLE) <- sub(pattern="_LUNG",replacement="",x=names(KRAS_CCLE))
colnames(CCLE_EXP) <- sub(pattern="_LUNG",replacement="",x=colnames(CCLE_EXP))

# see if here is some signal with KRAS status
tmp <- ifelse(KRAS_CCLE!="WT",1,0)
fit <- eBayes(lmFit(CCLE_EXP,model.matrix(~tmp)))
hist(fit$p.value[,2],breaks=100, main="CCLE exp ~ KRAS")

#####################################################################################
# make all datasets comparable 
#####################################################################################

# # first make all the db have coherent features
tmp1 <- intersect(rownames(CCLE_EXP),rownames(LUAD_EXP))
LUAD_EXP <- LUAD_EXP[tmp1,]
CCLE_EXP <- CCLE_EXP[tmp1,]
rm(tmp1)

# rescale the data (chemores) to have the same mean and variance than the LUAD
# Justin's function to rescale the VS to get the same mean/var than the TS
normalize_to_X <- function(mean.x, sd.x, Y){
  m.y <- rowMeans(Y)
  sd.y <- apply(Y, 1, sd)
  Y.adj <- (Y - m.y) * sd.x / sd.y  + mean.x 
  Y.adj[sd.y == 0] <- mean.x[sd.y==0]
  Y.adj
}

CCLE_EXP <- normalize_to_X(rowMeans(LUAD_EXP),apply(LUAD_EXP,1,sd),CCLE_EXP)

# get rid of the probes that are the less variant
tmp <- apply(LUAD_EXP,1,sd)
tmp1 <- which(tmp>quantile(tmp,probs=.2))
LUAD_EXP <- LUAD_EXP[tmp1,]
CCLE_EXP <- CCLE_EXP[tmp1,]
rm(tmp1,tmp)

#############################################################################
# input the LUAD TCGA mutations data
#############################################################################
load(file="/home/cferte/FELLOW/cferte/KRAS_Analysis/MATMUT_GENE_LUAD.RData")
colnames(MATMUT_GENE_LUAD) <- substr(x=colnames(MATMUT_GENE_LUAD),1,16)
tmp <- intersect(colnames(MATMUT_GENE_LUAD),colnames(LUAD_EXP))
mutations.luad <- MATMUT_GENE_LUAD[,tmp]
rm(MATMUT_GENE_LUAD)

# ################################################################################################################
# load the ccle data
# ################################################################################################################
MAF1 <- read.delim("/home/cferte/cell_line_data/CCLE_hybrid_capture1650_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07.maf",header=TRUE)
MAF1 <- MAF1[which(MAF1$Hugo_Symbol!="Unknown"),]
MAF1 <- MAF1[,c("Hugo_Symbol","Tumor_Sample_Barcode","Protein_Change", "Genome_Change")]
MAF1 <- MAF1[(grep(pattern="LUNG",x=MAF1$Tumor_Sample_Barcode)),]
MAF1$Tumor_Sample_Barcode <- sapply(strsplit(split="_",x=MAF1$Tumor_Sample_Barcode),function(x){x[[1]]})
MAF1$Protein_Change[MAF1$Protein_Change==""] <- MAF1$Genome_Change[ MAF1$Protein_Change==""]

# create a binary matrix with colnames equal to the sampel names and rownames equal the MUTIDs
MATMUT<-matrix(0,nrow=length(unique(MAF1$Hugo_Symbol)),ncol=length(unique(MAF1$Tumor_Sample_Barcode)))
colnames(MATMUT) <- unique(MAF1$Tumor_Sample_Barcode)
rownames(MATMUT) <- unique(MAF1$Hugo_Symbol)
for(i in rownames(MATMUT)){
  MATMUT[i,c(MAF1$Tumor_Sample_Barcode[which(MAF1$Hugo_Symbol==i)])] <- c(MAF1$Protein_Change[which(MAF1$Hugo_Symbol==i)])  
}
mutations.ccle <- MATMUT
rm(MAF1,MATMUT)
mutations.ccle[mutations.ccle!=0] <-1 
tmp <- rownames(mutations.ccle)
mutations.ccle <- apply(mutations.ccle,2,as.numeric)
rownames(mutations.ccle) <- tmp

# make the names of the samples coherents for ccle
tmp <- intersect(colnames(CCLE_EXP),colnames(mutations.ccle))
mutations.ccle <- mutations.ccle[,tmp]
CCLE_EXP <- CCLE_EXP[,tmp]
KRAS_CCLE <- KRAS_CCLE[tmp]

#############################################################################
# focus on G12C only
#############################################################################

tmp <- names(KRAS_LUAD)[which(KRAS_LUAD %in% c("G12C","WT"))]
KRAS_LUAD <- KRAS_LUAD[tmp]
LUAD_EXP <- LUAD_EXP[,tmp]
mutations.luad <- mutations.luad[,tmp]
rm(tmp)

tmp <- names(KRAS_CCLE)[which(KRAS_CCLE %in% c("G12C","WT"))]
KRAS_CCLE <- KRAS_CCLE[tmp]
CCLE_EXP <- CCLE_EXP[,tmp]
mutations.ccle <- mutations.ccle[,tmp]
rm(tmp)

#########################################################################
# make coherent the genes of mutations.ccle and .luad 

tmp <- intersect(rownames(mutations.ccle),rownames(mutations.luad))
mutations.ccle <- mutations.ccle[tmp,]
mutations.luad <- mutations.luad[tmp,]
tmp <- paste(tmp,"_mut",sep="")
identical(rownames(mutations.ccle),rownames(mutations.luad))
rownames(mutations.ccle) <- tmp
rownames(mutations.luad) <- tmp

###################################################################################################################
# train our predictive model of G12C in TCGA LUAD
# using a penalized regression approach  
# with alpha=.1 (more ridge) and determine lambda using nfolds= 5
# the robustness of the model is increased by boostrapping (n=100)
# validate each model in BATTLE
###################################################################################################################
par(mfrow=c(1,1))
require(glmnet)
G12C_LUAD <- ifelse(KRAS_LUAD=="G12C",1,0)
table(G12C_LUAD)
N <- 100
fit <- c()
#features <- c()
selected <- c()
yhat_CCLE <- c()
models <- 0
i <- 0
while(models<N)
{
  
  j <- c(names(which(G12C_LUAD==1)), sample(names(which(G12C_LUAD==0)),replace=TRUE))
  cv.fit <- cv.glmnet(t(rbind(LUAD_EXP[,j],mutations.luad[,j])), factor(G12C_LUAD[j]), nfolds=3, alpha=.1, family="binomial")
  fit <- glmnet(x=t(rbind(LUAD_EXP[,j],mutations.luad[,j])),y=factor(G12C_LUAD[j]),family="binomial",alpha=.1,lambda=cv.fit$lambda.1se)
  if(length(which(abs(as.numeric(fit$beta))> 10^-5))>10)
  {
    i=i+1
    print(i)
    #features <- c(features,length(which(abs(as.numeric(fit$beta))> 10^-5)))
    selected <- cbind(selected , as.numeric(fit$beta))
    
    yhat_CCLE <- cbind(yhat_CCLE,predict(fit, t(rbind(CCLE_EXP,mutations.ccle)),type="response"))
    
    models <- dim(yhat_CCLE)[2]
  } }

rownames(selected) <- rownames(fit$beta)



#######################################################
#evaluate the performance of the model in ccle
#######################################################
tmp <- ifelse(KRAS_CCLE=="G12C",1,0)
names(tmp) <- names(KRAS_CCLE)
G12C_CCLE <- tmp
rm(tmp)


boxplot(yhat_CCLE~G12C_CCLE,xlab=c("KRAS G12C mutational status"),ylab="model of G12C",main="predicting KRAS G12C in ccle
(modele trained in TCGA)")
stripchart(yhat_CCLE~G12C_CCLE,vertical=TRUE,pch=20,add=T,col="royalblue",cex=.7,method="jitter")

##################################################################
# plot a ROC curve to asses the performance of our model in ccle
##################################################################

AUC_CCLE <- c()
for (i in c(1:N)){
  Pred <- prediction(as.numeric(yhat_CCLE[,i]),as.numeric(G12C_CCLE))
  Perf <- performance(prediction.obj=Pred,"tpr","fpr")
  AUC <- performance(prediction.obj=Pred,"auc")
  AUC_CCLE <- c(AUC_CCLE,as.numeric(AUC@y.values))
}
boxplot(AUC_CCLE)


################
sort(apply(selected,1,sum),decreasing=TRUE)[1:50]


# ###############################################################################
# # infer drug correlations in CCLE and Sanger
# ###############################################################################
source("/home/cferte/FELLOW/cferte/KRAS_Analysis/drug_sensitivity_inference/G12C_drug_sens_correlations.R")

