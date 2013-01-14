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
# run the analysis for the G12C in the TCGA LUAD  
###############################################################
load("/home/cferte/FELLOW/cferte/KRAS_Project/OBJECTS/LUAD_EXP.RData")
load("/home/cferte/FELLOW/cferte/KRAS_Project/OBJECTS/KRAS_LUAD_STATUS.RData")
table(KRAS_LUAD)

# keep only the patients with available KRAS WT or G12C 
#KRAS_LUAD <- KRAS_LUAD[which(KRAS_LUAD %in% c("WT","G12C"))]
tmp <- intersect(names(KRAS_LUAD),colnames(LUAD_EXP))
LUAD_EXP <- LUAD_EXP[,tmp]
KRAS_LUAD <- KRAS_LUAD[tmp]
table(KRAS_LUAD)
KRAS_LUAD <- ifelse(KRAS_LUAD=="G12C",1,ifelse(KRAS_LUAD=="G12V",2,ifelse(KRAS_LUAD=="WT",3,NA)))
table(KRAS_LUAD)

# get rid of the NAs in KRAS_LUAD and LUAD_EXP
KRAS_LUAD <- KRAS_LUAD[which(!is.na(KRAS_LUAD))]
LUAD_EXP <- LUAD_EXP[,names(KRAS_LUAD)]

# check if there is any latent structure in the exp data
s <- svd(LUAD_EXP)
plot(s$v[,1],s$v[,2], col=c("royalblue","aquamarine3","orange")[KRAS_LUAD],pch=20,xlab="PC1",ylab="PC2",main="LUAD EXP n=122")

# train the model on LUAD
fit <- eBayes(lmFit(LUAD_EXP,model.matrix(~KRAS_LUAD+s$v[,2])))
hist(fit$p.value[,2],breaks=100,main="unadjusted p values (G12C vs. rest TCGA)")
abline(v=.05,col="red",lty=2)
#hist(p.adjust(fit$p.value[,2],method="BH"),breaks=200)

# run the permutation test with n=100,
# with alpha =.001
# compute the FDR
#PERMUT <- replicate(100, eBayes(lmFit(LUAD_EXP,
#                                model.matrix(~KRAS_LUAD[sample(length(KRAS_LUAD))])))$p.value[,2])
#alpha <- .001
#N <- mean(apply(PERMUT,2,function(x){sum(x<alpha)}))
# M <- sum(fit$p.value[,2]<alpha)
# FDR_LUAD <- N/M
# paste("FDR LUAD=",FDR_LUAD, "(with alpha=",alpha,")")
# PERMUT_LUAD <- PERMUT
# rm(PERMUT)

###############################################################
# load the CHEMORES data  
###############################################################

load("/home/cferte/FELLOW/cferte/KRAS_Project/OBJECTS/CHEMORES_EXP.RData")
load("/home/cferte/FELLOW/cferte/KRAS_Project/OBJECTS/CHEMORES_CLIN.RData")
load("/home/cferte/FELLOW/cferte/KRAS_Project/OBJECTS/KRAS_CHEMORES_STATUS.RData")

# transform KRAS_CHEMORES into a numeric variable 1 G12C, 2 G12V, 3 WT
KRAS_CHEMORES <- ifelse(KRAS_CHEMORES=="G12C",1,ifelse(KRAS_CHEMORES=="G12V",2,ifelse(KRAS_CHEMORES=="WT",3,NA)))
table(KRAS_CHEMORES)

# retain only the patients that are evaluated for KRAS
KRAS_CHEMORES <- KRAS_CHEMORES[which(!is.na(KRAS_CHEMORES))]
CHEMORES_EXP <- CHEMORES_EXP[,names(KRAS_CHEMORES)]
CHEMORES_CLIN <- CHEMORES_CLIN[names(KRAS_CHEMORES),]

# retain only the ADENOS only from CHEMORES 
CHEMORES_CLIN <- CHEMORES_CLIN[ CHEMORES_CLIN$Histology=="AC",]
CHEMORES_EXP <- CHEMORES_EXP[,rownames(CHEMORES_CLIN)]
KRAS_CHEMORES <- KRAS_CHEMORES[rownames(CHEMORES_CLIN)]
table(KRAS_CHEMORES)

# check if there is any latent structure in the exp data
s <- svd(CHEMORES_EXP)
plot(s$v[,1],s$v[,2], col=c("royalblue","aquamarine3","orange")[KRAS_CHEMORES],pch=20,xlab="PC1",ylab="PC2",main="CHEMORES EXP n=43")

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
DB_GOLD <- c(rep(1,times=127),rep(0,times=47))
KRAS_GOLD <- c(KRAS_LUAD,KRAS_CHEMORES)

# check if there is any latent structure in the exp data
s <- svd(GOLD_EXP)
par(mfrow=c(2,2))
plot(s$v[,1],s$v[,2], col=c("royalblue","aquamarine3","orange")[KRAS_GOLD],pch=20,xlab="PC1",ylab="PC2",main="GOLD EXP n=188")
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
colnames(new.dat) <- colnames(GOLD_EXP)

# svd plots again
s <- svd(new.dat)
plot(s$v[,1],s$v[,2], col=c("royalblue","aquamarine3","orange")[KRAS_GOLD],pch=20,xlab="PC1",ylab="PC2",main="GOLD EXP n=43")
plot(s$v[,1],s$v[,2], col=ifelse(DB_GOLD==1,"black","red"),pch=20,xlab="PC1",ylab="PC2",main="GOLD EXP n=43")

# train the model on GOLD1
par(mfrow=c(1,1))
GOLD1 <- new.dat
fit <- eBayes(lmFit(GOLD1,model.matrix(~KRAS_GOLD +s$v[,2])))
hist(fit$p.value[,2],breaks=100,main="unadjusted p values 
  (G12C vs. all GOLD)")
abline(v=.05,col="red",lty=2)
GOLD1_names <- fit$p.value[,2]
table(fit$p.value[,2]<.05)

# remove the second SVD from GOLD1 since this is improving 
# remove the first svd of the battle
s$d[2] <- 0
GOLD2 <- GOLD1
GOLD2 <- s$u %*% diag(s$d) %*% t(s$v)
rownames(GOLD2) <- rownames(GOLD1)
colnames(GOLD2) <- colnames(GOLD1)
GOLD1 <- GOLD2
rm(GOLD2)    

#################################################################################################################
# apply the models in the CCLE database and correlate them with drug sensitivity
#################################################################################################################

ccle_kras <- ccle_kras$objects$KRAS_CCLE
table(ccle_kras)

ccle_kras <- ifelse(ccle_kras=="p.G12C",1,ifelse(ccle_kras=="p.G12V",2,ifelse(ccle_kras=="0",3,NA)))
ccle_kras <- ccle_kras[which(!is.na(ccle_kras))]
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
plot(s$v[,1],s$v[,2], col=c("royalblue","aquamarine3","orange")[ccle_kras],pch=20,xlab="PC1",ylab="PC2",main="ccle_eset n=70")

# check if there is any signal in the ccle data
fit <- eBayes(lmFit(ccle_eset1,model.matrix(~ccle_kras+s$v[,4])))
hist(fit$p.value[,2],breaks=100,main="unadjusted p values 
  (G12C vs. G12V ccle)")
abline(v=.05,col="red",lty=2)
ccle_names <- fit$p.value[,2]
table(fit$p.value[,2]<.05)

# remove the 4th SVDw from ccle_eset since this is improving 
# remove the first svd of the battle
s$d[4] <- 0
ccle_eset2 <- s$u %*% diag(s$d) %*% t(s$v)
colnames(ccle_eset2) <- colnames(ccle_eset1)
rownames(ccle_eset2) <- rownames(ccle_eset1)
ccle_eset1 <- ccle_eset2
rm(ccle_eset2)

# check if there is any latent structure in the data
s <- svd(ccle_eset1)
plot(s$v[,1],s$v[,2], col=c("royalblue","aquamarine3","orange")[ccle_kras],pch=20,xlab="PC1",ylab="PC2",main="ccle_eset n=70")

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
tmp <- intersect(names(GOLD1_names)[which(GOLD1_names<.1)],names(ccle_names)[which(ccle_names<.1)])
length(tmp)

###################################################################################################################
# perform the prediction using the glmnet package with alpha=.1 (more ridge) and determine lambda using nfolds= 10
###################################################################################################################

require(glmnet)
set.seed(1234567)
cv.fit <- cv.glmnet(t(GOLD1[tmp,]), factor(KRAS_GOLD), nfolds=10, alpha=.1, family="multinomial")
plot(cv.fit)
fit <- glmnet(x=t(GOLD1[tmp,]),y=KRAS_GOLD,family="multinomial",alpha=.1,lambda=cv.fit$lambda.1se)
 


##############################################################################################################################################
# create a better model based on the non penalization of features that are significantly diffenretially expressed between G12C and G12V alone
##############################################################################################################################################

tmp1 <- names(KRAS_GOLD)[which(KRAS_GOLD!=3)]
train_set <- GOLD1[,tmp1]
train_kras <- KRAS_GOLD[tmp1]
table(train_kras)
fit1 <- eBayes(lmFit(train_set,model.matrix(~train_kras )))
hist(fit1$p.value[,2])


tmp2 <-names(ccle_kras)[which(ccle_kras!=3)] 
test_set <- ccle_eset1[,tmp2]
test_kras <- ccle_kras[tmp2]
table(test_kras)
fit2 <- eBayes(lmFit(test_set,model.matrix(~test_kras)))
hist(fit2$p.value[,2])

pen <- fit1$p.value[,2]
set.seed(1234567)
cv.fit <- cv.glmnet(t(GOLD1[tmp,]), factor(KRAS_GOLD), nfolds=10, alpha=.1, family="multinomial",penalty.factor=pen)
plot(cv.fit)
fit <- glmnet(x=t(GOLD1[tmp,]),y=KRAS_GOLD,family="multinomial",alpha=.1,lambda=cv.fit$lambda.1se,penalty.factor=pen)


####################################
# apply this model in the ccle data
####################################
y_hat <- predict(fit, t(ccle_eset1[tmp,]),type="response")

######################################
#evaluate the performance of the model
######################################
y_hat <- as.data.frame(y_hat)
boxplot(y_hat)
boxplot(y_hat~ccle_kras,xlab=c("KRAS mutational status (1=G12C, 2= G12V, 3=WT)"),ylab="multinomial model",main="predicting KRAS status in ccle 
using a multinomial modele trained in tcga + chemores")
stripchart(y_hat~ccle_kras,vertical=TRUE,pch=20,add=T,col=c("royalblue","aquamarine3","orange")[ccle_kras],cex=.7 )

########################################################
# assess the performance of the multinomial model
########################################################



require(ROCR)
Pred <- prediction(as.numeric(d_G12C),as.numeric(ccle_kras))
Perf <- performance(prediction.obj=Pred,"tpr","fpr")
AUC <- performance(prediction.obj=Pred,"auc")
plot(Perf, col="royalblue",main="predicting KRAS G12C in ccle 
(modele trained in tcga + chemores)")
text(x=.7,y=.4,labels=paste("AUC=",format(x=AUC@y.values,digits=2)),col="royalblue")





########################################################
# correlate the G12C with the drug response
########################################################

p <- c()
D <- c()
G12C_SENS <- sapply(colnames(drug),function(x){
  p <- cor.test(drug[,x],d_G12C,method="spearman",use="pairwise.complete.obs") 
  D <- rbind(D,c(p$estimate,p$p.value))
})
rownames(G12C_SENS) <- c("rho","p.value")
G12C_SENS <- as.data.frame(t(G12C_SENS))
G12C_SENS <- G12C_SENS[sort(G12C_SENS$rho,decreasing=TRUE,index.return=TRUE)$ix,]

rownames(G12C_SENS)

G12C_SENS$color <- ifelse(G12C_SENS$p.value<.05,"red",ifelse(G12C_SENS$p.value<.1,"orange","black"))

# we prefer to get rid of irinotecan given there are too many missing values
G12C_SENS <- G12C_SENS[-which(rownames(G12C_SENS)=="Irinotecan"),]

dotchart(G12C_SENS$rho,labels=row.names(G12C_SENS),pch=20,xlab="correlation (spearman rho)",main="correlation of drug sensitivity 
with G12C model in lung ccle",col=G12C_SENS$color)
abline(v=0,lty=3,lwd=.3)




