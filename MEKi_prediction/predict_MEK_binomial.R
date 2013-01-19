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
CCLE_EXP <- ccle_EXP

##########################################
# load the Sanger data 
##########################################
source("~/FELLOW/cferte/KRAS_Analysis/data_input/sanger_load_data.R")

# make the gene expresson data coherent
tmp <- intersect(rownames(CCLE_EXP),rownames(SANGER_EXP))
SANGER_EXP <- SANGER_EXP[tmp,]
CCLE_EXP <- CCLE_EXP[tmp,]

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
par(mfrow=c(1,1))
hist(rat,breaks=100, main="correlation between the gene expression data from Sanger & CCLE",ylab="genes (N)")
tmp <- rownames(A)[which(rat>.7)]
SANGER_EXP <- SANGER_EXP[tmp,]
CCLE_EXP <- CCLE_EXP[tmp,]
rm(A,B,rat,tmp,raton)

#####################################################################################
# rescale the gene expression data so the ccle and sanger gene expression are comparable
#####################################################################################

# rescale the data (chemores) to have the same mean and variance than the LUAD
# Justin's function to rescale the VS to get the same mean/var than the TS
normalize_to_X <- function(mean.x, sd.x, Y){
  m.y <- rowMeans(Y)
  sd.y <- apply(Y, 1, sd)
  Y.adj <- (Y - m.y) * sd.x / sd.y  + mean.x 
  Y.adj[sd.y == 0] <- mean.x[sd.y==0]
  Y.adj
}

SANGER_EXP <- normalize_to_X(rowMeans(CCLE_EXP),apply(CCLE_EXP,1,sd),SANGER_EXP)

#####################################################################################
# get rid of the probes that are the less variant
#####################################################################################

tmp <- apply(CCLE_EXP,1,sd)
tmp1 <- which(tmp>quantile(tmp,probs=.1))
CCLE_EXP <- CCLE_EXP[tmp1,]
SANGER_EXP <- SANGER_EXP[tmp1,]
rm(tmp1,tmp)

#############################################################################
# focus on the MEK inhibitors
#############################################################################

# identify the names of the drugs targetting MEK
ccle.drugs <- read.delim(file="/home/cferte/FELLOW/cferte/KRAS_Analysis/drug_sensitivity_inference/CCLE_drugs.txt",header=T,skip=2)
sanger.drugs <- read.delim(file="/home/cferte/FELLOW/cferte/KRAS_Analysis/drug_sensitivity_inference/Sanger_drugs.txt",header=T)
sanger.drugs <- sanger.drugs[which(duplicated(sanger.drugs$DRUG.NAME)!=TRUE),c("DRUG.NAME","TARGET")]
sanger.drugs$DRUG.NAME <- sub(pattern="-",replacement=".",x=sanger.drugs$DRUG.NAME)

MEK.ccle <-   ccle.drugs$Compound..code.or.generic.name.[grep(pattern="MEK",ccle.drugs$Target.s.)]
MEK.sanger <- sanger.drugs$DRUG.NAME[grep(pattern="MEK",sanger.drugs$TARGET)]
rm(sanger.drugs,ccle.drugs)

# plot the distribution of the IC50 of the MEKi drugs
par(mfrow=c(2,3), oma=c(1,1,5,1))
plot(sort(SangerDrug[,MEK.sanger[1]]),main=paste(MEK.sanger[1]),ylab="IC50",pch=20)
abline(h=quantile(SangerDrug[,MEK.sanger[1]],probs=.2,na.rm=TRUE),col="red")
plot(sort(SangerDrug[,MEK.sanger[2]]),main=paste(MEK.sanger[2]),ylab="IC50",pch=20)
abline(h=quantile(SangerDrug[,MEK.sanger[2]],probs=.2,na.rm=TRUE),col="red")
plot(sort(SangerDrug[,MEK.sanger[3]]),main=paste(MEK.sanger[3]),ylab="IC50",pch=20)
abline(h=quantile(SangerDrug[,MEK.sanger[3]],probs=.2,na.rm=TRUE),col="red")
plot(sort(SangerDrug[,MEK.sanger[4]]),main=paste(MEK.sanger[4]),ylab="IC50",pch=20)
abline(h=quantile(SangerDrug[,MEK.sanger[4]],probs=.2,na.rm=TRUE),col="red")
plot(sort(ccle_drug[,MEK.ccle[1]]),main=paste(MEK.ccle[1]),ylab="IC50",pch=20)
abline(h=quantile(ccle_drug[,MEK.ccle[1]],probs=.2,na.rm=TRUE),col="red")
plot(sort(ccle_drug[,MEK.ccle[2]]),main=paste(MEK.ccle[2]),ylab="IC50",pch=20)
abline(h=quantile(ccle_drug[,MEK.ccle[2]],probs=.2,na.rm=TRUE),col="red")
title(main= "distribution of the IC50 across the MEK inhibitors in SANGER & CCLE db",outer=TRUE)

# identify the cells that are evaluated for MEK inhibitors in sanger
MEK.cells.sanger <- rownames(SangerDrug)[-unique(c(which(is.na(SangerDrug[,MEK.sanger[1]])), which(is.na(SangerDrug[,MEK.sanger[2]])), which(is.na(SangerDrug[,MEK.sanger[3]])), which(is.na(SangerDrug[,MEK.sanger[4]]))))]
M <- SangerDrug[MEK.cells.sanger,MEK.sanger]
MEK.cells.ccle <- rownames(ccle_drug)[-unique(c(which(is.na(ccle_drug[,MEK.ccle[1]])), which(is.na(ccle_drug[,MEK.ccle[2]]))))]
M2 <- ccle_drug[MEK.cells.ccle,MEK.ccle]

require(car)
scatterplotMatrix(M,main="correlations in the IC50 of the MEK inhibitors in Sanger")
scatterplotMatrix(normalizeCyclicLoess(M),main="correlations in the IC50 of the MEK inhibitors in Sanger (Loess normalized)")
cor(M,method="spearman",use="pairwise.complete.obs")
cor(normalizeCyclicLoess(M),method="spearman",use="pairwise.complete.obs")

scatterplotMatrix(normalizeCyclicLoess(M2),main="correlations in the IC50 of the MEK inhibitors in ccle (Loess normalized)")
cor(normalizeCyclicLoess(M2),method="spearman",use="pairwise.complete.obs")

par(mfrow=c(2,3))
fit <- eBayes(lmFit(SANGER_EXP[,MEK.cells.sanger],model.matrix(~M[,1])))
hist(fit$p.value[,2],breaks=100, main=paste("Sanger Expr ~",colnames(M)[1]))
table(fit$p.value[,2]<.05)
fit <- eBayes(lmFit(SANGER_EXP[,MEK.cells.sanger],model.matrix(~M[,2])))
hist(fit$p.value[,2],breaks=100, main=paste("Sanger Expr ~",colnames(M)[2]))
table(fit$p.value[,2]<.05)
fit <- eBayes(lmFit(SANGER_EXP[,MEK.cells.sanger],model.matrix(~M[,3])))
hist(fit$p.value[,2],breaks=100, main=paste("Sanger Expr ~",colnames(M)[3]))
table(fit$p.value[,2]<.05)
fit <- eBayes(lmFit(SANGER_EXP[,MEK.cells.sanger],model.matrix(~M[,4])))
hist(fit$p.value[,2],breaks=100, main=paste("Sanger Expr ~",colnames(M)[4]))
table(fit$p.value[,2]<.05)
fit <- eBayes(lmFit(CCLE_EXP[,MEK.cells.ccle],model.matrix(~M2[,1])))
hist(fit$p.value[,2],breaks=100, main=paste("ccle Expr ~",colnames(M2)[1]))
table(fit$p.value[,2]<.05)
fit <- eBayes(lmFit(CCLE_EXP[,MEK.cells.ccle],model.matrix(~M2[,2])))
hist(fit$p.value[,2],breaks=100, main=paste("ccle Expr ~",colnames(M2)[2]))
table(fit$p.value[,2]<.05)
title(main="univariate difefrential expression for sensitivity to MEKi across SANGER & CCLE",outer=TRUE)

# ###################################################################################################################
# # restrict to the common cell lines between ccle and sanger
# ###################################################################################################################
# 
# tmp <- intersect(rownames(M),rownames(M2))
# 
# 
# ################################################################################################################
# # load the mutations data
# ################################################################################################################
# mutations.ccle <- loadEntity("syn1528027")
# mutations.ccle <- mutations.ccle$objects$eSet_hybrid
# mutations.ccle <- exprs(mutations.ccle)
# tmp <- intersect(tmp,colnames(mutations.ccle))
# mutations.ccle <- mutations.ccle[,tmp]

# binomial transformation of the independant variables
M3 <- apply(M,2,function(y){ifelse(y<quantile(y,probs=.2),1,0)})
M4 <- apply(M2,2,function(y){ifelse(y<quantile(y,probs=.2),1,0)})

###################################################################################################################
# GOLD standard: correlations between IC50 of CCLE and Sanger
###################################################################################################################

cor(M[tmp,],M2[tmp,],method="spearman")
scatterplotMatrix(normalizeCyclicLoess(cbind(M[tmp,],M2[tmp,])))
title(main="correlations in the IC50 of the MEK inhibitors 
between CCLE & Sanger (Loess normalized)", outer=TRUE)

###################################################################################################################
# train our predictive model of MEK response in ccle
# using a penalized regression approach  
# with alpha=.1 (more ridge) and determine lambda using nfolds= 5
# the robustness of the model is increased by boostrapping (n=100)
# validate each model in BATTLE
###################################################################################################################
ccle.top <- ifelse(apply(M4,1,sum)==2,1,0)
names(ccle.top) <- rownames(M4)
# 
# IC50.MEK.ccle <- apply(normalizeCyclicLoess(M2[tmp,]),1,mean)
# IC50.MEK.sanger <- apply(normalizeCyclicLoess(M[tmp,]),1,mean)
# names(IC50.MEK.sanger) <- rownames(M[tmp,])
# names(IC50.MEK.ccle) <- rownames(M2[tmp,])

par(mfrow=c(1,1))
require(glmnet)

N <- 10
fit <- c()
selected <- c()
yhat <- c()
models <- 0
i <- 0
while(models<N)
{
j <- c(names(which(ccle.top==1)),sample(names(which(ccle.top==0)),replace=TRUE))
cv.fit <- cv.glmnet(t(CCLE_EXP[,j]), y=ccle.top[j], nfolds=3, alpha=.01)
fit <- glmnet(x=t(CCLE_EXP[,j]),y=ccle.top[j],alpha=.01,lambda=cv.fit$lambda.min,family="binomial")
if(length(which(abs(as.numeric(fit$beta))> 10^-5))>10)
{
i=i+1
print(i)
selected <- cbind(selected , as.numeric(fit$beta))
yhat <- cbind(yhat,predict(fit, t(SANGER_EXP[,rownames(M3)]),type="response"))
models <- dim(yhat)[2]
} }

AUC_RDEA119 <- c()
#par(mfrow=c(2,2))
require(ROCR)
for (i in c(1:N)){
Pred <- prediction(as.numeric(yhat[,i]),as.numeric(M3[,1]))
Perf <- performance(prediction.obj=Pred,"tpr","fpr")
AUC <- performance(prediction.obj=Pred,"auc")
AUC_RDEA119 <- c(AUC_RDEA119,as.numeric(AUC@y.values))
}
AUC_CI.1040 <- c()
#par(mfrow=c(2,2))
require(ROCR)
for (i in c(1:N)){
Pred <- prediction(as.numeric(yhat[,i]),as.numeric(M3[,2]))
Perf <- performance(prediction.obj=Pred,"tpr","fpr")
AUC <- performance(prediction.obj=Pred,"auc")
AUC_CI.1040 <- c(AUC_CI.1040,as.numeric(AUC@y.values))
}
AUC_PD.0325901 <- c()
#par(mfrow=c(2,2))
require(ROCR)
for (i in c(1:N)){
Pred <- prediction(as.numeric(yhat[,i]),as.numeric(M3[,3]))
Perf <- performance(prediction.obj=Pred,"tpr","fpr")
AUC <- performance(prediction.obj=Pred,"auc")
AUC_PD.0325901 <- c(AUC_PD.0325901,as.numeric(AUC@y.values))
}
AUC_AZD6244 <- c()
#par(mfrow=c(2,2))
require(ROCR)
for (i in c(1:N)){
Pred <- prediction(as.numeric(yhat[,i]),as.numeric(M3[,4]))
Perf <- performance(prediction.obj=Pred,"tpr","fpr")
AUC <- performance(prediction.obj=Pred,"auc")
AUC_AZD6244 <- c(AUC_AZD6244,as.numeric(AUC@y.values))
}
boxplot(cbind(AUC_RDEA119,AUC_CI.1040,AUC_PD.0325901,AUC_AZD6244),outline=FALSE,ylab="prediction of sensitivity (AUC)")
stripchart(list(RDEA119=AUC_RDEA119,CI.1040=AUC_CI.1040,PD.0325901=AUC_PD.0325901,AZD6244=AUC_AZD6244),add=TRUE,method="jitter",vertical=TRUE,col="royalblue",pch=20)
title( main=" bootstrapped elasticnet models of gene expression + hybrid capture sequencing
trained in 57 ccle treated by PD.0325901 or AZD6244
validation in the 71 cell lines processed by the Sanger group",outer=TRUE)
