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
# run the analysis for the G12V in the TCGA LUAD  & CHEMORES (G12V vs any other KRAS)
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

# see if here is some signal with KRAS status
tmp <- ifelse(KRAS_CCLE!="WT",1,0)
fit <- eBayes(lmFit(CCLE_EXP,model.matrix(~tmp)))
hist(fit$p.value[,2],breaks=100, main="CCLE exp ~ KRAS")


# see if the signal is different between G12C and G12V
tmp1 <- names(which(fit$p.value[,2]<.01))
tmp <- names(KRAS_CCLE)[which(KRAS_CCLE %in% c("G12C","G12V"))]
s <- svd(CCLE_EXP[tmp1,tmp]-rowMeans(CCLE_EXP[tmp1,tmp]))
boxplot(s$v[,1]~KRAS_CCLE[tmp],main="difference in G12C and G12V according to the Princ. Component 1 (CCLE)")
stripchart(s$v[,1]~KRAS_CCLE[tmp],cex=.7,col="royalblue",method="jitter",vertical=TRUE,pch=20,add=TRUE)
title(sub=paste("P=",format(wilcox.test(s$v[,1]~KRAS_CCLE[tmp])$p.value,digits=3),"(wilcoxon test)"))


##########################################
# load the Sanger data 
##########################################
source("~/FELLOW/cferte/KRAS_Analysis/data_input/sanger_load_data.R")

# see if here is some signal with KRAS status
tmp <- ifelse(KRAS_SANGER!="WT",1,0)
fit <- eBayes(lmFit(SANGER_EXP,model.matrix(~tmp)))
hist(fit$p.value[,2],breaks=100, main="SANGER exp ~ KRAS")

# see if the signal is different between G12C and G12V
tmp1 <- names(which(fit$p.value[,2]<.01))
tmp <- names(KRAS_SANGER)[which(KRAS_SANGER %in% c("G12C","G12V"))]
s <- svd(SANGER_EXP[tmp1,tmp]-rowMeans(SANGER_EXP[tmp1,tmp]))
boxplot(s$v[,1]~KRAS_SANGER[tmp],main="difference in G12C and G12V according to the Princ. Component 1 (SANGER)")
stripchart(s$v[,1]~KRAS_SANGER[tmp],cex=.7,col="royalblue",method="jitter",vertical=TRUE,pch=20,add=TRUE)
title(sub=paste("P=",format(wilcox.test(s$v[,1]~KRAS_SANGER[tmp])$p.value,digits=3),"(wilcoxon test)"))



#####################################################################################
# restrict to the probes that are highly correlated between sanger ccle
#####################################################################################
A <- SANGER_EXP
B <- CCLE_EXP
tmp <- intersect(rownames(SANGER_EXP),rownames(CCLE_EXP))
A <- A[tmp,]
B <- B[tmp,]
colnames(B) <- gsub(x=colnames(B),pattern="_LUNG",replacement="")

tmp <- intersect(colnames(A),colnames(B))
rat <- c()
raton <- c()
for(i in c(1:dim(A)[1]))
{
  raton  <- cor(A[i,tmp],B[i,tmp],method="spearman")
  rat <- c(rat,raton)
}
plot(density(rat))
tmp <- rownames(A)[which(rat>.5)]
SANGER_EXP <- SANGER_EXP[tmp,]
CCLE_EXP <- CCLE_EXP[tmp,]
rm(A,B,rat,tmp,raton)

#####################################################################################
# make all datasets comparable 
#####################################################################################

# # first make all the db have coherent features
tmp1 <- intersect(rownames(LUAD_EXP),intersect(rownames(BATTLE_EXP),intersect(rownames(SANGER_EXP),rownames(CHEMORES_EXP)))))
LUAD_EXP <- LUAD_EXP[tmp1,]
CHEMORES_EXP <- CHEMORES_EXP[tmp1,]
CCLE_EXP <- CCLE_EXP[tmp1,]
BATTLE_EXP <- BATTLE_EXP[tmp1,]
SANGER_EXP <- SANGER_EXP[tmp1,]
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

CHEMORES_EXP <- normalize_to_X(rowMeans(LUAD_EXP),apply(LUAD_EXP,1,sd),CHEMORES_EXP)
CCLE_EXP <- normalize_to_X(rowMeans(LUAD_EXP),apply(LUAD_EXP,1,sd),CCLE_EXP)
BATTLE_EXP <- normalize_to_X(rowMeans(LUAD_EXP),apply(LUAD_EXP,1,sd),BATTLE_EXP)
SANGER_EXP <- normalize_to_X(rowMeans(LUAD_EXP),apply(LUAD_EXP,1,sd),SANGER_EXP)

# get rid of the probes that are the less variant
tmp <- apply(LUAD_EXP,1,sd)
tmp1 <- which(tmp>quantile(tmp,probs=.2))
LUAD_EXP <- LUAD_EXP[tmp1,]
CHEMORES_EXP <- CHEMORES_EXP[tmp1,]
CCLE_EXP <- CCLE_EXP[tmp1,]
BATTLE_EXP <- BATTLE_EXP[tmp1,]
SANGER_EXP <- SANGER_EXP[tmp1,]
rm(tmp1,tmp)

#############################################################################
# focus on the drug sensitivity
#############################################################################

# read drugs & targets file
ccle.drugs <- read.delim(file="/home/cferte/FELLOW/cferte/KRAS_Analysis/drug_sensitivity_inference/CCLE_drugs.txt",header=T,skip=2)
sanger.drugs <- read.delim(file="/home/cferte/FELLOW/cferte/KRAS_Analysis/drug_sensitivity_inference/Sanger_drugs.txt",header=T)
sanger.drugs <- sanger.drugs[which(duplicated(sanger.drugs$DRUG.NAME)!=TRUE),c("DRUG.NAME","TARGET")]
sanger.drugs$DRUG.NAME <- sub(pattern="-",replacement=".",x=sanger.drugs$DRUG.NAME)
drug.names <-   ccle.drugs$Compound..code.or.generic.name.[grep(pattern="MEK",ccle.drugs$Target.s.)]



