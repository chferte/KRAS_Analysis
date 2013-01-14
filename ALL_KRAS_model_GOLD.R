# Charles Ferté
# Sage Bionetworks
# 14 Sept 2012


# train a model of differential gene expression per KRAS codon in Lung Adenocarcinoma
# both for G12C and G12V, combining CHEMORES and in TCGA LUAD

#load the different packages
options(stringsAsFactors=FALSE)

library(affy)
library(corpcor)
library(lattice)
library(limma)
library(caret)
library(glmnet)
source("/home/cferte/FELLOW/cferte/KRAS_Project/JUSTIN_PREDICT_CCLE/code/lung_analysis_functions.R")
synapseLogin("charles.ferte@sagebase.org","charles")

###############################################################
# run the analysis for All KRAS in the TCGA LUAD  
###############################################################
load("/home/cferte/FELLOW/cferte/KRAS_Project/OBJECTS/LUAD_EXP.RData")
load("/home/cferte/FELLOW/cferte/KRAS_Project/OBJECTS/KRAS_LUAD_STATUS.RData")
table(KRAS_LUAD)

# keep only the patients with available KRAS WT or G12C 
tmp <- intersect(names(KRAS_LUAD),colnames(LUAD_EXP))
LUAD_EXP <- LUAD_EXP[,tmp]
KRAS_LUAD <- KRAS_LUAD[tmp]
KRAS_LUAD <- ifelse(KRAS_LUAD=="WT",0,1)
table(KRAS_LUAD)

# check if there is any latent structure in the exp data
s <- svd(LUAD_EXP)
plot(s$v[,1],s$v[,2], col=ifelse(KRAS_LUAD==1,"royalblue","orange"),pch=20,xlab="PC1",ylab="PC2",main="LUAD EXP n=122")

# train the model on LUAD
fit <- eBayes(lmFit(LUAD_EXP,model.matrix(~KRAS_LUAD)))
hist(fit$p.value[,2],breaks=100,main="unadjusted p values (G12C vs. rest TCGA)")
abline(v=.05,col="red",lty=2)

###############################################################
# load the CHEMORES data  
###############################################################

load("/home/cferte/FELLOW/cferte/KRAS_Project/OBJECTS/CHEMORES_EXP.RData")
load("/home/cferte/FELLOW/cferte/KRAS_Project/OBJECTS/CHEMORES_CLIN.RData")
load("/home/cferte/FELLOW/cferte/KRAS_Project/OBJECTS/KRAS_CHEMORES_STATUS.RData")

# retain only the patients that are evaluated for KRAS
KRAS_CHEMORES <- KRAS_CHEMORES[which(!is.na(KRAS_CHEMORES))]
CHEMORES_EXP <- CHEMORES_EXP[,names(KRAS_CHEMORES)]
CHEMORES_CLIN <- CHEMORES_CLIN[names(KRAS_CHEMORES),]

# retain only the ADENOS only from CHEMORES 
CHEMORES_CLIN <- CHEMORES_CLIN[ CHEMORES_CLIN$Histology=="AC",]
CHEMORES_EXP <- CHEMORES_EXP[,rownames(CHEMORES_CLIN)]
KRAS_CHEMORES <- KRAS_CHEMORES[rownames(CHEMORES_CLIN)]

# transform KRAS_CHEMORES into a numeric variable
KRAS_CHEMORES <- ifelse(KRAS_CHEMORES=="WT",0,1)
table(KRAS_CHEMORES)

# check if there is any latent structure in the exp data
s <- svd(CHEMORES_EXP)
plot(s$v[,1],s$v[,2], col=ifelse(KRAS_CHEMORES==1,"royalblue","orange"),pch=20,xlab="PC1",ylab="PC2",main="CHEMORES EXP n=43")

# train the model on CHEMORES
fit <- eBayes(lmFit(CHEMORES_EXP,model.matrix(~KRAS_CHEMORES)))
hist(fit$p.value[,2],breaks=100,main="unadjusted p values 
  (G12C vs. rest CHEMORES)")
abline(v=.05,col="red",lty=2)
#hist(p.adjust(fit$p.value[,2],method="BH"),breaks=200)

################################################################
# merge the two dataset into the GOLD one to generate the model
################################################################


################################################################
# load the ccle data from synapse
################################################################
ccle_eset <- loadEntity('syn1354643')
ccle_drug <- loadEntity('syn1354656')
ccle_kras <- loadEntity('syn1443160')

# make all the db have the same features

tmp <- intersect(featureNames(ccle_eset$objects$ccle_eset),intersect(rownames(LUAD_EXP),rownames(CHEMORES_EXP)))
LUAD_EXP <- LUAD_EXP[tmp,]
CHEMORES_EXP <- CHEMORES_EXP[tmp,]
ccle_eset <- ccle_eset$objects$ccle_eset[tmp,]

# second quantile Normalization
require(caret)
CHEMO_EXP2 <- normalize2Reference(CHEMORES_EXP,rowMeans(LUAD_EXP))

# rescale the data (chemores) to have the same mean and variance than the LUAD
# Justin's function to rescale the VS to get the same mean/var than the TS
normalize_to_X <- function(mean.x, sd.x, Y){
  m.y <- rowMeans(Y)
  sd.y <- apply(Y, 1, sd)
  Y.adj <- (Y - m.y) * sd.x / sd.y  + mean.x 
  Y.adj[sd.y == 0] <- mean.x[sd.y==0]
  Y.adj
}

CHEMO_EXP2 <- normalize_to_X(rowMeans(LUAD_EXP),apply(LUAD_EXP,1,sd),CHEMO_EXP2)

# attempt to draw a qq plot
#qqplot(CHEMO_EXP2,LUAD_EXP)

# combine the two datasets
GOLD_EXP <- cbind(LUAD_EXP,CHEMO_EXP2)
DB_GOLD <- c(rep(1,times=137),rep(0,times=51))
KRAS_GOLD <- c(KRAS_LUAD,KRAS_CHEMORES)

# check if there is any latent structure in the exp data
s <- svd(GOLD_EXP)
plot(s$v[,1],s$v[,2], col=ifelse(KRAS_GOLD==1,"royalblue","orange"),pch=20,xlab="PC1",ylab="PC2",main="GOLD EXP n=188")
plot(s$v[,1],s$v[,2], col=ifelse(DB_GOLD==1,"black","red"),pch=20,xlab="PC1",ylab="PC2",main="GOLD EXP n=188")

# SNM to  remove the batch and adjustment variables effect
require(snm)
bio.var <- model.matrix(~ as.factor(KRAS_GOLD))
adj.var <- model.matrix(~ DB_GOLD)
snm.fit <- snm(GOLD_EXP, 
               bio.var=bio.var, 
               adj.var=adj.var, 
               rm.adj=TRUE)

new.dat <- snm.fit$norm.dat

# svd plots again
s <- svd(new.dat)
plot(s$v[,1],s$v[,2], col=ifelse(KRAS_GOLD==1,"royalblue","orange"),pch=20,xlab="PC1",ylab="PC2",main="GOLD EXP n=43")
plot(s$v[,1],s$v[,2], col=ifelse(DB_GOLD==1,"black","red"),pch=20,xlab="PC1",ylab="PC2",main="GOLD EXP n=43")

# train the model on GOLD1
GOLD1 <- new.dat
colnames(GOLD1) <- colnames(GOLD_EXP)
rownames(GOLD1) <- rownames(GOLD_EXP)
fit <- eBayes(lmFit(GOLD1,model.matrix(~KRAS_GOLD)))
hist(fit$p.value[,2],breaks=100,main="unadjusted p values 
  (G12C vs. all GOLD)")
abline(v=.05,col="red",lty=2)
GOLD1_names <- fit$p.value[,2]
table(fit$p.value[,2]<.05)

#################################################################################################################
# apply the models in the CCLE database and correlate them with drug sensitivity
#################################################################################################################

ccle_kras <- ccle_kras$objects$KRAS_CCLE
ccle_kras <- ifelse(ccle_kras=="0",0,1)
table(ccle_kras)

ccle_drug <- ccle_drug$objects$ccle_drug
ccle_drug <- ccle_drug[names(ccle_kras),]
drug <- ccle_drug@data

ccle_eset <- ccle_eset[,names(ccle_kras)]

# load the cell lines metadata and exclude the SCLC and the undifferenciated carcinoma
mData <- loadEntity('syn1443797')
mData <- mData$objects$CCLE_metadata
rownames(mData) <- mData$CCLE.name
mData <- mData[sampleNames(ccle_eset),]
table(mData$Hist.Subtype1)
#mData <- mData[mData$Hist.Subtype1!= "small_cell_carcinoma",]
#mData <- mData[mData$Hist.Subtype1!= "undifferentiated_carcinoma",]
#mData <- mData[mData$Hist.Subtype1!= "squamous_cell_carcinoma",]
ccle_eset <- ccle_eset[,rownames(mData)]
ccle_kras <- ccle_kras[rownames(mData)]
drug <- drug[rownames(mData),]
table(ccle_kras)

# perform quantile Normalization on ccle_eset with ref =GOLD1
ccle_eset1 <- normalize2Reference(exprs(ccle_eset),rowMeans(GOLD1))

# rescale the ccle_eset1 to have the same mean and variance than the GOLD1
# Justin's function to rescale the VS to get the same mean/var than the TS
normalize_to_X <- function(mean.x, sd.x, Y){
  m.y <- rowMeans(Y)
  sd.y <- apply(Y, 1, sd)
  Y.adj <- (Y - m.y) * sd.x / sd.y  + mean.x 
  Y.adj[sd.y == 0] <- mean.x[sd.y==0]
  Y.adj
}

ccle_eset1 <- normalize_to_X(rowMeans(GOLD1),apply(GOLD1,1,sd),ccle_eset1)

# check if there is any latent structure in the data
s <- svd(ccle_eset1)
plot(s$v[,1],s$v[,2], col=ifelse(KRAS_GOLD==1,"royalblue","orange"),pch=20,xlab="PC1",ylab="PC2",main="ccle_eset n=82")

# check if there is any signal in the ccle data

fit <- eBayes(lmFit(ccle_eset1,model.matrix(~ccle_kras+s$v[,4]+s$v[,5])))
hist(fit$p.value[,2],breaks=100,main="unadjusted p values 
  (G12C vs. G12V ccle)")
abline(v=.05,col="red",lty=2)
ccle_names <- fit$p.value[,2]
table(fit$p.value[,2]<.05)

# remove the 4th, 5th and 6th SVDw from ccle_eset since this is improving 
# remove the first svd of the battle
s$d[4] <- 0
s$d[5] <- 0
ccle_eset2 <- s$u %*% diag(s$d) %*% t(s$v)
rownames(ccle_eset2) <- rownames(ccle_eset1)
colnames(ccle_eset2) <- colnames(ccle_eset1)
ccle_eset1 <- ccle_eset2
rm(ccle_eset2)

# check if there is any latent structure in the data
s <- svd(ccle_eset1)
plot(s$v[,1],s$v[,2], col=ifelse(KRAS_GOLD==1,"royalblue","orange"),pch=20,xlab="PC1",ylab="PC2",main="ccle_eset n=71")

## plot the PVAL across GOLD1 and ccle
require(geneplotter)
require("RColorBrewer")

dim <- c(0,1)
par(mfrow=c(1,1),oma=c(0,0,3,0))
tmp <- intersect(names(GOLD1_names),names(ccle_names))
x <- GOLD1_names[tmp]
y <- ccle_names[tmp]
colors  <- densCols(x,y)
plot(x,y, col=colors, pch=20,xlab="GOLD1 (TCGA+CHEMORES) p values",ylab="CCLE p values",xlim=dim,ylim=dim)
abline(h=.05,col="red",lty=2)
abline(v=.05,col="red",lty=2)

table(GOLD1_names<.05)
table(ccle_names<.05)

#################################################################
# perform a feature selection on the overlapping probes only
#################################################################
tmp <- intersect(names(GOLD1_names)[which(GOLD1_names<.2)],names(ccle_names)[which(ccle_names<.2)])
length(tmp)

###################################################################################################################
# perform the prediction using the glmnet package with alpha=.1 (more ridge) and determine lambda using nfolds= 10
###################################################################################################################

require(glmnet)
set.seed(1234567)
cv.fit <- cv.glmnet(t(GOLD1[tmp,]), factor(KRAS_GOLD), nfolds=10, alpha=.1, family="binomial")
plot(cv.fit)
fit <- glmnet(x=t(GOLD1[tmp,]),y=KRAS_GOLD,family="binomial",alpha=.1,lambda=cv.fit$lambda.1se)
table(as.numeric(fit$beta)!=0)    



####################################
# apply this model in the ccle data
####################################
y_hat <- predict(fit, t(ccle_eset1[tmp,]),type="response")

########################################################
#evaluate the performance of the model
########################################################

boxplot(y_hat~ccle_kras,xlab=c("KRAS G12C mutational status"),ylab="model of G12C",main="predicting KRAS in ccle 
(modele trained in tcga + chemores)")
stripchart(y_hat~ccle_kras,vertical=TRUE,pch=20,add=T,col="royalblue",cex=.7 )

########################################################
# plot a ROC curve to asses the performance of our model
########################################################

require(ROCR)
Pred <- prediction(as.numeric(y_hat),as.numeric(ccle_kras))
Perf <- performance(prediction.obj=Pred,"tpr","fpr")
AUC <- performance(prediction.obj=Pred,"auc")
plot(Perf, col="royalblue",main="predicting KRAS G12C in ccle 
(modele trained in tcga + chemores)")
text(x=.7,y=.4,labels=paste("AUC=",format(x=AUC@y.values,digits=2)),col="royalblue")

###############################################################################
# boxplot the model across the different KRAS codon mutations within ccle
###############################################################################
codon_kras <- loadEntity('syn1443160')
codon_kras_ccle <- codon_kras$objects$KRAS_CCLE[colnames(ccle_eset1)]
codon_kras_ccle[codon_kras_ccle=="0"] <- "WT"
boxplot(y_hat~codon_kras_ccle,xlab=c("KRAS G12C mutational status"),ylab="KRAS lung model",main="predicting KRAS in ccle 
(modele trained in tcga + chemores)")
stripchart(y_hat~codon_kras_ccle,vertical=TRUE,pch=20,add=T,col="royalblue",cex=.7 )

###################################################################################################################
# improve the model by a dynamic penalization using the G12V/G12C diff p val vector
###################################################################################################################
load("/home/cferte/FELLOW/cferte/KRAS_Project/OBJECTS/KRAS_LUAD_STATUS.RData")
load("/home/cferte/FELLOW/cferte/KRAS_Project/OBJECTS/KRAS_CHEMORES_STATUS.RData")
KCV <- c(KRAS_LUAD,KRAS_CHEMORES)
KCV <- KCV[colnames(GOLD1)]
KCV <- KCV[which(KCV %in% c("G12C","G12V"))]
GOLDCV <- GOLD1[,names(KCV)]

fit <- eBayes(lmFit(GOLDCV,model.matrix(~KCV)))
hist(fit$p.value[,2],breaks=100,main="unadjusted p values (G12C vs.G12V)")
abline(v=.05,col="red",lty=2)
table(fit$p.value[,2]<.05)
pen <- fit$p.value[,2]
pen <- ifelse(pen<.05,.5,1)

set.seed(1234567)
cv.fit <- cv.glmnet(t(GOLD1), factor(KRAS_GOLD), nfolds=10, alpha=.1, family="binomial",penalty.factor=pen)
plot(cv.fit)
fit <- glmnet(x=t(GOLD1[tmp,]),y=KRAS_GOLD,family="binomial",alpha=.1,lambda=cv.fit$lambda.1se,penalty.factor=pen)
table(as.numeric(fit$beta)!=0)    
y_hat2 <- predict(fit, t(ccle_eset1[tmp,]),type="response")

###############################################################################
# boxplot the model across the different KRAS codon mutations within ccle
###############################################################################
codon_kras <- loadEntity('syn1443160')
codon_kras_ccle <- codon_kras$objects$KRAS_CCLE[colnames(ccle_eset1)]
codon_kras_ccle[codon_kras_ccle=="0"] <- "WT"
boxplot(y_hat2~codon_kras_ccle,xlab=c("KRAS G12C mutational status"),ylab="KRAS lung model",main="predicting KRAS in ccle 
(modele trained in tcga + chemores)")
stripchart(y_hat2~codon_kras_ccle,vertical=TRUE,pch=20,add=T,col="royalblue",cex=.7 )

table(codon_kras_ccle)



########################################################
# correlate the y_hat2 with the drug response
########################################################
drug1 <- drug[rownames(y_hat2),]
dim(drug1)
p <- c()
D <- c()
yhat2 <- sapply(colnames(drug1),function(x){
  p <- cor.test(drug[,x],y_hat2,method="spearman",use="pairwise.complete.obs") 
  D <- rbind(D,c(p$estimate,p$p.value))
})
rownames(yhat2) <- c("rho","p.value")
class(yhat2)
yhat2 <- as.data.frame(t(yhat2))

yhat2 <- yhat2[sort(yhat2$rho,decreasing=TRUE,index.return=TRUE)$ix,]

yhat2

yhat2$color <- ifelse(yhat2$p.value<.05,"red",ifelse(yhat2$p.value<.1,"orange","black"))

# we prefer to get rid of irinotecan given there are too many missing values
yhat2 <- yhat2[-which(rownames(yhat2)=="Irinotecan"),]

dotchart(yhat2$rho,labels=row.names(yhat2),pch=20,xlab="correlation (spearman rho)",main="correlation of drug sensitivity 
with G12C model in lung ccle",col=yhat2$color)
abline(v=0,lty=3,lwd=.3)

# do it in the G12C only ones...

# redo your model for KRAS only
# train in TCGA KRAS mutant only validate it in CHEMORES kras only - try to predict G12C
# if you get similar AUC in ccle then infer drug sensitivity information




