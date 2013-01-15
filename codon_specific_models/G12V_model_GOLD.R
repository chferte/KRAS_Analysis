# Charles Ferté
# Sage Bionetworks
# 14 Sept 2012


# train a model of differential gene expression per KRAS codon in Lung Adenocarcinoma
# for G12V,  in TCGA LUAD

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

# 
# # see if the signal is different between G12C and G12V
# tmp1 <- names(which(fit$p.value[,2]<.05))
# tmp <- names(KRAS_LUAD)[which(KRAS_LUAD %in% c("G12C","G12V"))]
# 
# s <- svd(LUAD_EXP[tmp1,tmp]-rowMeans(LUAD_EXP[tmp1,tmp]))
# boxplot(s$v[,1]~KRAS_LUAD[tmp],main="difference in G12C and G12V according to the Princ. Component 1")
# stripchart(s$v[,1]~KRAS_LUAD[tmp],cex=.7,col="royalblue",method="jitter",vertical=TRUE,pch=20,add=TRUE)
# title(sub=paste("P=",format(wilcox.test(s$v[,1]~KRAS_LUAD[tmp])$p.value,digits=3),"(wilcoxon test)"))
# 
# fit <- eBayes(lmFit(LUAD_EXP[tmp1,tmp],model.matrix(~KRAS_LUAD[tmp])))
# hist(fit$p.value[,2],breaks=100, main="TCGA exp ~ KRAS G12C/V")
# 

##########################################
# load the Battle data 
##########################################

BATTLE_EXP  <- loadEntity('syn1467735')
BATTLE_EXP <- BATTLE_EXP$objects$Battle_snm

KRAS_BATTLE <- loadEntity('syn1465678')
KRAS_BATTLE <- KRAS_BATTLE$objects$KRAS_BATTLE

# see if here is some signal with KRAS status
tmp <- ifelse(KRAS_BATTLE != "WT",1,0)
fit <- eBayes(lmFit(BATTLE_EXP,model.matrix(~tmp)))
hist(fit$p.value[,2],breaks=100,main="BATTLE EXP ~ KRAS")

# see if the signal is different between G12C and G12V
tmp1 <- names(which(fit$p.value[,2]<.05))
tmp <- names(KRAS_BATTLE)[which(KRAS_BATTLE %in% c("G12C","G12V"))]
s <- svd(BATTLE_EXP[tmp1,tmp]-rowMeans(BATTLE_EXP[tmp1,tmp]))
boxplot(s$v[,1]~KRAS_BATTLE[tmp],main="difference in G12C and G12V according to the Princ. Component 1 (BATTLE)")
stripchart(s$v[,1]~KRAS_BATTLE[tmp],cex=.7,col="royalblue",method="jitter",vertical=TRUE,pch=20,add=TRUE)
title(sub=paste("P=",format(wilcox.test(s$v[,1]~KRAS_BATTLE[tmp])$p.value,digits=3),"(wilcoxon test)"))

###############################################################
# load the CHEMORES data  (without SCC and SCLC)
###############################################################

source("/home/cferte/FELLOW/cferte/KRAS_Analysis/data_input/Chemores_load_data.R")

# see if here is some signal with KRAS status
tmp <- ifelse(KRAS_CHEMORES!="WT",1,0)
fit <- eBayes(lmFit(CHEMORES_EXP,model.matrix(~tmp)))
hist(fit$p.value[,2],breaks=100, main="chemores exp ~ KRAS")

# see if the signal is different between G12C and G12V
tmp1 <- names(which(fit$p.value[,2]<.05))
tmp <- names(KRAS_CHEMORES)[which(KRAS_CHEMORES %in% c("G12C","G12V"))]
s <- svd(CHEMORES_EXP[tmp1,tmp]-rowMeans(CHEMORES_EXP[tmp1,tmp]))
boxplot(s$v[,1]~KRAS_CHEMORES[tmp],main="difference in G12C and G12V according to the Princ. Component 1 (CHEMORES)")
stripchart(s$v[,1]~KRAS_CHEMORES[tmp],cex=.7,col="royalblue",method="jitter",vertical=TRUE,pch=20,add=TRUE)
title(sub=paste("P=",format(wilcox.test(s$v[,1]~KRAS_CHEMORES[tmp])$p.value,digits=3),"(wilcoxon test)"))


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
tmp1 <- intersect(rownames(CCLE_EXP),intersect(rownames(LUAD_EXP),intersect(rownames(BATTLE_EXP),intersect(rownames(SANGER_EXP),rownames(CHEMORES_EXP)))))
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
# focus on G12V only
#############################################################################

tmp <- names(KRAS_LUAD)[which(KRAS_LUAD %in% c("G12V","WT"))]
KRAS_LUAD <- KRAS_LUAD[tmp]
LUAD_EXP <- LUAD_EXP[,tmp]
rm(tmp)

tmp <- names(KRAS_CHEMORES)[which(KRAS_CHEMORES %in% c("G12V","WT"))]
KRAS_CHEMORES <- KRAS_CHEMORES[tmp]
CHEMORES_EXP <- CHEMORES_EXP[,tmp]
CHEMORES_CLIN <- CHEMORES_CLIN[tmp,]
rm(tmp)

tmp <- names(KRAS_CCLE)[which(KRAS_CCLE %in% c("G12V","WT"))]
KRAS_CCLE <- KRAS_CCLE[tmp]
CCLE_EXP <- CCLE_EXP[,tmp]
rm(tmp)

tmp <- names(KRAS_BATTLE)[which(KRAS_BATTLE %in% c("G12V","WT"))]
KRAS_BATTLE <- KRAS_BATTLE[tmp]
BATTLE_EXP <- BATTLE_EXP[,tmp]
rm(tmp)

tmp <- names(KRAS_SANGER)[which(KRAS_SANGER %in% c("G12V","WT"))]
KRAS_SANGER <- KRAS_SANGER[tmp]
SANGER_EXP <- SANGER_EXP[,tmp]
rm(tmp)

##############################################################################################
# explore the biological meaning of the significant overlapping genes (fisher exact test) to see what is the biological signal 
##############################################################################################

tmp <- ifelse(KRAS_LUAD =="G12V",1,0)
fit <- eBayes(lmFit(LUAD_EXP[,names(tmp)],model.matrix(~tmp)))
hist(fit$p.value[,2],breaks=100)
table(fit$p.value[,2]<.05)

tmp <- ifelse(KRAS_CHEMORES=="G12V",1,0)
fit2 <- eBayes(lmFit(CHEMORES_EXP,model.matrix(~tmp)))
hist(fit2$p.value[,2],breaks=100)
table(fit2$p.value[,2]<.05)

tmp <- ifelse(KRAS_BATTLE=="G12V",1,0)
fit3 <- eBayes(lmFit(BATTLE_EXP,model.matrix(~tmp)))
hist(fit3$p.value[,2],breaks=100)
table(fit3$p.value[,2]<.05)

tmp <- ifelse(KRAS_CCLE=="G12V",1,0)
fit4 <- eBayes(lmFit(CCLE_EXP,model.matrix(~tmp)))
hist(fit4$p.value[,2],breaks=100)
table(fit4$p.value[,2]<.05)

tmp <- ifelse(KRAS_SANGER=="G12V",1,0)
fit5 <- eBayes(lmFit(SANGER_EXP,model.matrix(~tmp)))
hist(fit5$p.value[,2],breaks=100)
table(fit5$p.value[,2]<.05)

# biological meaning of the intersect
BIO_p <- intersect(intersect(names(which(fit$p.value[,2]<.1)), names(which(fit3$p.value[,2]<.1))),names(which(fit2$p.value[,2]<.1)))
BIO_p_cl <- intersect(names(which(fit4$p.value[,2]<.1)), names(which(fit5$p.value[,2]<.1)))
paste(BIO_p,collapse=" ")
paste(BIO_p_cl,collapse=" ")
rm(BIO_p,tmp)

###################################################################################################################
# train our predictive model of G12V in TCGA LUAD
# using a penalized regression approach  
# with alpha=.1 (more ridge) and determine lambda using nfolds= 5
# the robustness of the model is increased by boostrapping (n=100)
# validate each model in BATTLE
###################################################################################################################
par(mfrow=c(1,1))
require(glmnet)
G12V_LUAD <- ifelse(KRAS_LUAD=="G12V",1,0)
table(G12V_LUAD)
N <- 100
fit <- c()
#features <- c()
yhat_BATTLE <- c()
yhat_CHEMORES <- c()
yhat_CCLE <- c()
yhat_SANGER <- c()
selected <- c()
models <- 0
i <- 0
while(models<N)
{
  
  j <- c(names(which(G12V_LUAD==1)), sample(names(which(G12V_LUAD==0)),replace=TRUE))
  cv.fit <- cv.glmnet(t(LUAD_EXP[,j]), factor(G12V_LUAD[j]), nfolds=3, alpha=.1, family="binomial")
  fit <- glmnet(x=t(LUAD_EXP[,j]),y=factor(G12V_LUAD[j]),family="binomial",alpha=.1,lambda=cv.fit$lambda.1se)
    if(length(which(abs(as.numeric(fit$beta))> 10^-5))>10)
    {
      i=i+1
      print(i)
      #features <- c(features,length(which(abs(as.numeric(fit$beta))> 10^-5)))
      selected <- cbind(selected , as.numeric(fit$beta))
      yhat_BATTLE <- cbind(yhat_BATTLE,predict(fit, t(BATTLE_EXP),type="response"))
  yhat_CHEMORES <- cbind(yhat_CHEMORES,predict(fit, t(CHEMORES_EXP),type="response"))
  yhat_CCLE <- cbind(yhat_CCLE,predict(fit, t(CCLE_EXP),type="response"))
  yhat_SANGER <- cbind(yhat_SANGER,predict(fit, t(SANGER_EXP),type="response"))
  models <- dim(yhat_SANGER)[2]
  } }
rownames(selected) <- rownames(fit$beta)

###############################################################################
# assess the model performance in terms of AUC
###############################################################################
source("/home/cferte/FELLOW/cferte/KRAS_Analysis/model_performance/G12V_WT_performance.R")

###############################################################################
# infer drug correlations in CCLE and Sanger
###############################################################################
source("/home/cferte/FELLOW/cferte/KRAS_Analysis/drug_sensitivity_inference/G12V_drug_sens_correlations.R")
