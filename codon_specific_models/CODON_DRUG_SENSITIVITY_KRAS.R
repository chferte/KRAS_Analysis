## charles ferté,MD
## Sage Bionetworks 
## "2012-07-01"

## KRAS project LUAD

## mutational profile in LUAD & CRC
require(randomForest)
require(Biobase)
require(synapseClient)
require(HDclassif)
require(mclust)
require(Biobase)
require(ggplot2)
require(corpcor)
require(survival)
library(affy)
library(corpcor)
library(lattice)
library(snm)
library(WGCNA)
library(glmnet)

## synapse Login
synapseLogin("charles.ferte@sagebase.org","charles")
source("/home/cferte/FELLOW/cferte/KRAS_Project/JUSTIN_PREDICT_CCLE/code/lung_analysis_functions.R")


### magic option 
options(stringsAsFactors=FALSE)

ccle <- getCCLE()
ccle_eset <- ccle[[1]]
ccle_drug <- ccle[[2]]
idxs <- match(sampleNames(ccle_eset), sampleNames(ccle_drug))
ccle_eset <- ccle_eset[,!is.na(idxs)]
ccle_response <- ccle_drug[na.omit(idxs),]

#restrict the eset to the only 82 LUNG cell lines
lung.mask <- grepl("LUNG", sampleNames(ccle_eset))
ccle_eset <- ccle_eset[,lung.mask]

# create a vector with the KRAS status for these 82 cell lines
KRs <- ccle_eset@phenoData@data$KRAS
names(KRs) <- rownames(ccle_eset@phenoData@data)
table(KRs)

# retrieve the names of the cell lines associated with their mutations
mafCCLE <- loadEntity('syn427004')
MAF <- read.delim(file.path(mafCCLE$cacheDir,mafCCLE$files))
rm(mafCCLE)
 
# restrict to the sample names that are intersecting in the the MAF file & the eset
tmp <- intersect(unique(MAF$Tumor_Sample_Barcode),sampleNames(ccle_eset))
ccle_eset <- ccle_eset[,tmp]
MAF <- MAF[MAF$Tumor_Sample_Barcode %in% tmp,]
 
#restrict the MAF file & the eset to the 47 coherent LUNG cell lines
MAF <- MAF[grep(pattern="LUNG",x=MAF$Tumor_Sample_Barcode),]
ccle_eset <- ccle_eset[,grep(pattern="LUNG",sampleNames(ccle_eset))]
 
# restrict the MAF file to the KRAS mutations
MAF <- MAF[ MAF$Hugo_Symbol=="KRAS",c("Protein_Change","Tumor_Sample_Barcode")]
MAF <- MAF[-which(duplicated(MAF)),]

# call again the 82 ccle lung since we want to input the KRAS mutations protein change information with each cell line
ccle <- getCCLE()
ccle_eset <- ccle[[1]]
ccle_drug <- ccle[[2]]
idxs <- match(sampleNames(ccle_eset), sampleNames(ccle_drug))
ccle_eset <- ccle_eset[,!is.na(idxs)]
ccle_response <- ccle_drug[na.omit(idxs),]
lung.mask <- grepl("LUNG", sampleNames(ccle_eset))
ccle_eset <- ccle_eset[,lung.mask]

# input a vector with the KRAS protein change status in these 82 cell lines
KRs <- ccle_eset@phenoData@data$KRAS
names(KRs) <- rownames(ccle_eset@phenoData@data)
KRs[ MAF$Tumor_Sample_Barcode] <- MAF$Protein_Change

KRAS_CCLE <- KRs
table(KRAS_CCLE)


# store KRAS_CCLE in Synapse
KRAS_ccle <- Data(list(name = "KRAS CCLE", parentId = 'syn1337457'))
KRAS_ccle <- createEntity(KRAS_ccle)

# add files into the data entity
KRAS_ccle <- addObject(KRAS_ccle,KRAS_CCLE)

# push the raw data into this entity
KRAS_ccle <- storeEntity(entity=KRAS_ccle)


# see what are the differences in terms of drug sensitivty per codon:
DRUG <- ccle_drug@data[names(KRs),]
DRUG$KRs <- KRs

# look at the sensitivity for G12C
m <- c()
G12C_DRUG <- DRUG[DRUG$KRs=="p.G12C",-25]
for(i in c(1:24)){
  m <- c(m,mean(G12C_DRUG[,i],na.rm=T))
}
names(m) <- colnames(G12C_DRUG)

n <- c()
G12V_DRUG <- DRUG[DRUG$KRs=="p.G12V",-25]
for(i in c(1:24)){
  n <- c(n,mean(G12V_DRUG[,i],na.rm=T))
}
names(n) <- colnames(G12V_DRUG)

K61 <- c()
KRAS61_DRUG <- DRUG[ grep(DRUG$KRs,pattern="61"),-25]
for(i in c(1:24)){
  K61 <- c(K61,mean(KRAS61_DRUG[,i],na.rm=T))
}
names(K61) <- colnames(KRAS61_DRUG)

K12 <- c()
KRAS12_DRUG <- DRUG[ grep(DRUG$KRs,pattern="12"),-25]
for(i in c(1:24)){
  K12 <- c(K12,mean(KRAS12_DRUG[,i],na.rm=T))
}
names(K12) <- colnames(KRAS12_DRUG)



par(mfrow=c(1,1))
boxplot(G12C_DRUG,ylab="IC50",las=2,main="IC50 in 11 cell lines with KRAS G12C mutation")
boxplot(G12V_DRUG,ylab="IC50",las=2,main="IC50 in 6 cell lines with KRAS G12V mutation")
boxplot(KRAS61_DRUG,ylab="IC50",las=2,main="IC50 in 4 cell lines with KRAS G12V mutation")
boxplot(KRAS12_DRUG,ylab="IC50",las=2,main="IC50 in 21 cell lines with KRAS G12V mutation")

DRUG1 <- DRUG[,-c(25,6)]
require(gplots)
table(KRs)
tmp <- sort(as.numeric(as.factor(KRs)),index.return=TRUE)$ix
DRUG1 <- DRUG1[tmp,]
KRs1 <- KRs[tmp]
heatmap.2(t(DRUG1),trace="none",col=greenred(10),ColSideColors=rainbow(11)[as.numeric(as.factor(KRs1))],Colv=NULL,scale="none")


# only the KRAS mutant
DRUG2 <- DRUG[DRUG$KRs != "0",]
KRS2 <- DRUG2$KRs
DRUG2 <- DRUG2[,-c(25,6)]
require(gplots)
table(KRS2)
tmp <- sort(as.numeric(as.factor(KRS2)),index.return=TRUE)$ix
DRUG2 <- DRUG2[tmp,]
KRS22 <- KRS2[tmp]
heatmap.2(t(DRUG2),trace="none",col=greenred(10),ColSideColors=rainbow(11)[as.numeric(as.factor(KRS22))],Colv=NULL,scale="none")
KRS22

# KRAS G12C and G12V only
DRUG2 <- DRUG[DRUG$KRs %in% c("p.G12C","p.G12V"),]
KRS2 <- DRUG2$KRs
DRUG2 <- DRUG2[,-c(25,6)]
require(gplots)
table(KRS2)
tmp <- sort(as.numeric(as.factor(KRS2)),index.return=TRUE)$ix
DRUG2 <- DRUG2[tmp,]
KRS22 <- KRS2[tmp]
heatmap.2(t(DRUG2),trace="none",col=greenred(10),ColSideColors=rainbow(11)[as.numeric(as.factor(KRS22))],scale="none")
KRS22



# create a vector of KRAS status names KRs
#KRs <- rep(0,times=length(unique(MAF$Tumor_Sample_Barcode)))
#names(KRs) <- unique(MAF$Tumor_Sample_Barcode)
#KRs[unique(MAF$Tumor_Sample_Barcode[grep(pattern="KRAS",x=MAF$Hugo_Symbol)])] <- 1
#table(KRs)

# this is the training set
TS_EXP <- exprs(ccle_eset)
TS_y <- KRs

# load the four lung datasets
WILK_EXP <- loadEntity('syn1338139')
VS_EXP <- WILK_EXP$objects$WILK_EXP
WILK_CLIN  <-  loadEntity('syn1338095')
VS_y <- WILK_CLIN$objects$WILK_CLIN
VS_y <- VS_y$KRAS

LUAD_EXP <- loadEntity('syn1337634')
VS2_EXP <- LUAD_EXP$objects$LUAD_EXP
KRAS_LUAD <- loadEntity('syn1338272')
VS2_y <- ifelse(KRAS_LUAD$objects$KRAS_LUAD=="WT",0,1)
tmp <- intersect(colnames(VS2_EXP),names(VS2_y))
VS2_y <- VS2_y[tmp]

CHEMORES_EXP <- loadEntity('syn1337594')
VS3_EXP <- CHEMORES_EXP$objects$CHEMORES_EXP
CHEMORES_CLIN <- loadEntity('syn1337568')
VS3_y <- CHEMORES_CLIN$objects$CHEMORES_CLIN
VS3_y <- VS3_y$KRAS
names(VS3_y) <- rownames(CHEMORES_CLIN$objects$CHEMORES_CLIN)
VS3_y <- VS3_y[!is.na(VS3_y)]
VS3_EXP <- VS3_EXP[,names(VS3_y)]

BATTLE_EXP <- loadEntity('syn1337529')
VS4_EXP <- BATTLE_EXP$objects$BATTLE_EXP
colnames(VS4_EXP) <- substr(colnames(VS4_EXP),1,9)
BATTLE_CLIN <- loadEntity('syn1337513')
VS4_y <- BATTLE_CLIN$objects$BATTLE_CLIN
VS4_y <- VS4_y$KRAS
names(VS4_y) <- rownames(BATTLE_CLIN$objects$BATTLE_CLIN)
VS4_y <- VS4_y[!is.na(VS4_y)]
tmp <- names(VS4_y)
VS4_y <- ifelse(VS4_y=="WT",0,1)
names(VS4_y) <- tmp
VS4_EXP <- VS4_EXP[,names(VS4_y)]
battle.m <- rowMeans(VS4_EXP)
pc <- svd(VS4_EXP - battle.m)
plot(pc$v[,1],pc$v[,2])

#remove the first svd of the battle
pc$d[1] <- 0
EXP4 <- VS4_EXP
EXP4 <- pc$u %*% diag(pc$d) %*% t(pc$v) + battle.m 
rownames(EXP4) <- rownames(VS4_EXP)

# first make all the exp set with the same features
tmp <- intersect(rownames(TS_EXP),rownames(VS_EXP))
tmp <- intersect(tmp,rownames(VS2_EXP))
tmp <- intersect(tmp,rownames(VS3_EXP))
tmp <- intersect(tmp,rownames(EXP4))
length(tmp)

TS_EXP <- TS_EXP[tmp,]
VS_EXP <- VS_EXP[tmp,]
VS2_EXP <- VS2_EXP[tmp,]
VS3_EXP <- VS3_EXP[tmp,]
EXP4 <- EXP4[tmp,]
rm(tmp)

#optimisation (selection) of the number of features according to the top variant ones
 
# tmp <- apply(data.frame(apply(TS_EXP,1,var),apply(VS_EXP,1,var),apply(VS2_EXP,1,var),apply(VS3_EXP,1,var),apply(VS4_EXP,1,var)),1,mean)
# tmp <- names(which(tmp>quantile(tmp,probs=.60)))
# TS_EXP <- TS_EXP[tmp,]
# VS_EXP <- VS_EXP[tmp,]
# VS2_EXP <- VS2_EXP[tmp,]
# VS3_EXP <- VS3_EXP[tmp,]
# VS4_EXP <- VS4_EXP[tmp,]
# rm(tmp)

# get rid of the NAs in the KRAS vectors and the EXP of VS,VS3 and VS4
tmp <- which(!is.na(VS_y))
VS_y <- VS_y[tmp]
VS_EXP <- VS_EXP[,tmp]

tmp <- which(!is.na(VS3_y))
VS3_y <- VS3_y[tmp]
VS3_EXP <- VS3_EXP[,tmp]

tmp <- which(!is.na(VS4_y))
VS4_y <- VS4_y[tmp]
VS4_EXP <- EXP4[,tmp]

# predict the KRAS status fitted in the ccle (TS) into the 4 various datasets

# perform the prediction using the glmnet package with alpha=.1 (more ridge) and determine lambda using nfolds= 10
set.seed(1234567)
cv.fit <- cv.glmnet(t(TS_EXP), factor(TS_y), nfolds=20, alpha=.1, family="binomial")
plot(cv.fit)
fit <- glmnet(x=t(TS_EXP),y=TS_y,family="binomial",alpha=.1,lambda=cv.fit$lambda.1se)

# quantile normalization 
VSN <- normalize2Reference(VS_EXP, rowMeans(TS_EXP))
VS2N <- normalize2Reference(VS2_EXP, rowMeans(TS_EXP))
VS3N <- normalize2Reference(VS3_EXP, rowMeans(TS_EXP))
#VS4N <- normalize2Reference(VS4_EXP, rowMeans(TS_EXP))
VS4N <- normalize2Reference(EXP4, rowMeans(TS_EXP))

# Justin's function to rescale the VS to get the same mean/var than the TS
normalize_to_X <- function(mean.x, sd.x, Y){
  m.y <- rowMeans(Y)
  sd.y <- apply(Y, 1, sd)
  Y.adj <- (Y - m.y) * sd.x / sd.y  + mean.x 
  Y.adj[sd.y == 0] <- mean.x[sd.y==0]
  Y.adj
}


VSN <- normalize_to_X(rowMeans(TS_EXP),apply(TS_EXP,1,sd),VSN)
VS2N <- normalize_to_X(rowMeans(TS_EXP),apply(TS_EXP,1,sd),VS2N)
VS3N <- normalize_to_X(rowMeans(TS_EXP),apply(TS_EXP,1,sd),VS3N)
VS4N <- normalize_to_X(rowMeans(TS_EXP),apply(TS_EXP,1,sd),VS4N)


# predict it in the VS normalized (VSN)
y_hat <- predict(fit, t(VSN),type="response")
y_hat2 <- predict(fit, t(VS2N),type="response")
y_hat3 <- predict(fit, t(VS3N),type="response")
y_hat4 <- predict(fit, t(VS4N),type="response")

# plot a ROC curve to asses the response
require(ROCR)
erPred <- prediction(as.numeric(y_hat),as.numeric(VS_y))
erPerf <- performance(prediction.obj=erPred,"tpr","fpr")
erAUC <- performance(prediction.obj=erPred,"auc")
plot(erPerf, col="blue")
text(x=.7,y=.4,labels=paste("WILKERSON","AUC=",format(x=erAUC@y.values,digits=2)),col="blue")

erPred <- prediction(as.numeric(y_hat2),as.numeric(VS2_y))
erPerf <- performance(prediction.obj=erPred,"tpr","fpr")
erAUC <- performance(prediction.obj=erPred,"auc")
plot(erPerf, col="green",add=T)
text(x=.7,y=.35,labels=paste("TCGA LUAD","AUC=",format(x=erAUC@y.values,digits=2)),col="green")

erPred <- prediction(as.numeric(y_hat3),as.numeric(VS3_y))
erPerf <- performance(prediction.obj=erPred,"tpr","fpr")
erAUC <- performance(prediction.obj=erPred,"auc")
plot(erPerf, col="red",add=T)
text(x=.7,y=.3,labels=paste("CHEMORES","AUC=",format(x=erAUC@y.values,digits=2)),col="red")

erPred <- prediction(as.numeric(y_hat4),as.numeric(VS4_y))
erPerf <- performance(prediction.obj=erPred,"tpr","fpr")
erAUC <- performance(prediction.obj=erPred,"auc")
plot(erPerf, col="black",add=T)
text(x=.7,y=.25,labels=paste("BATTLE","AUC=",format(x=erAUC@y.values,digits=2)),col="black")
title("Performance of the KRAS model (CCLE) in 4 Adenocarcinoma datasets")

#############################################################################################
# boxplot the KRAS CCLE model across the different codons (CHEMORES and TCGA)
#############################################################################################

KCC <- CHEMORES_CLIN$objects$CHEMORES_CLIN$KRAS..NM_004448.2..exons.2_3
names(table(KCC))
KCC[KCC=="c.183A>W  p.Gln61His     (Q61H)    rs17851045 exon 3"] <- "Q61H"
KCC[KCC=="c.35G>K    p.Gly12Val (G12V) "] <- "G12V"
KCC[KCC=="c.35G>S    p.Gly12Ala (G12A) "] <- "G12A"
KCC[KCC=="c.34G>K    p.Gly12Cys (G12C) "] <- "G12C"
KCC[KCC=="c.34G>K    p.Gly12Cys (G12D) "] <- "G12D"
table(KCC)
names(KCC) <- rownames(CHEMORES_CLIN$objects$CHEMORES_CLIN)
WT_C <- predict(fit, t(VS3N[,names(which(KCC=="WT"))]),type="response")
G12C_C <- predict(fit, t(VS3N[,names(which(KCC=="G12C"))]),type="response")
G12V_C <- predict(fit, t(VS3N[,names(which(KCC=="G12V"))]),type="response")
boxplot(list(G12C=as.numeric(G12C_C),G12V=as.numeric(G12V_C),WT=as.numeric(WT_C)),ylab="RIS")
title("CHEMORES")
s <- svd(VS3N)
plot(s$v[,1],s$v[,2],col=as.factor(KCC),pch=20)

KCT <- loadEntity('syn1338272')
KCT <- KCT$objects$KRAS_LUAD
KCT <- KCT[intersect(colnames(VS2N),names(KCT))]
WT_T <- predict(fit, t(VS2N[,names(which(KCT=="WT"))]),type="response")
G12C_T <- predict(fit, t(VS2N[,names(which(KCT=="G12C"))]),type="response")
G12V_T <- predict(fit, t(VS2N[,names(which(KCT=="G12V"))]),type="response")
boxplot(list(G12C=as.numeric(G12C_T),G12V=as.numeric(G12V_T),WT=as.numeric(WT_T)),ylab="RIS")
title("TCGA_LUNG")
s <- svd(VS2N)
plot(s$v[,1],s$v[,2],col=as.factor(KCT),pch=20)

#############################################################################################
# boxplot the KRAS CCLE model across the different codons (CHEMORES and TCGA)
#############################################################################################

KCC <- CHEMORES_CLIN$objects$CHEMORES_CLIN$KRAS..NM_004448.2..exons.2_3
names(table(KCC))
KCC[KCC=="c.183A>W  p.Gln61His     (Q61H)    rs17851045 exon 3"] <- "Q61H"
KCC[KCC=="c.35G>K    p.Gly12Val (G12V) "] <- "G12V"
KCC[KCC=="c.35G>S    p.Gly12Ala (G12A) "] <- "G12A"
KCC[KCC=="c.34G>K    p.Gly12Cys (G12C) "] <- "G12C"
KCC[KCC=="c.34G>K    p.Gly12Cys (G12D) "] <- "G12D"
table(KCC)
names(KCC) <- rownames(CHEMORES_CLIN$objects$CHEMORES_CLIN)
WT_C <- predict(fit, t(VS3N[,names(which(KCC=="WT"))]),type="response")
G12C_C <- predict(fit, t(VS3N[,names(which(KCC=="G12C"))]),type="response")
G12V_C <- predict(fit, t(VS3N[,names(which(KCC=="G12V"))]),type="response")
boxplot(list(G12C=as.numeric(G12C_C),G12V=as.numeric(G12V_C),WT=as.numeric(WT_C)),ylab="RIS")
title("CHEMORES")
s <- svd(VS3N)
plot(s$v[,1],s$v[,2],col=as.factor(KCC),pch=20)

KCT <- loadEntity('syn1338272')
KCT <- KCT$objects$KRAS_LUAD
KCT <- KCT[intersect(colnames(VS2N),names(KCT))]
WT_T <- predict(fit, t(VS2N[,names(which(KCT=="WT"))]),type="response")
G12C_T <- predict(fit, t(VS2N[,names(which(KCT=="G12C"))]),type="response")
G12V_T <- predict(fit, t(VS2N[,names(which(KCT=="G12V"))]),type="response")
boxplot(list(G12C=as.numeric(G12C_T),G12V=as.numeric(G12V_T),WT=as.numeric(WT_T)),ylab="RIS")
title("TCGA_LUNG")
s <- svd(VS2N)
plot(s$v[,1],s$v[,2],col=as.factor(KCT),pch=20)

#############################################################################################
# boxplot the KRAS CCLE model across the different codons (CHEMORES and TCGA)
#############################################################################################

KCC <- CHEMORES_CLIN$objects$CHEMORES_CLIN$KRAS..NM_004448.2..exons.2_3
names(table(KCC))
KCC[KCC=="c.183A>W  p.Gln61His     (Q61H)    rs17851045 exon 3"] <- "Q61H"
KCC[KCC=="c.35G>K    p.Gly12Val (G12V) "] <- "G12V"
KCC[KCC=="c.35G>S    p.Gly12Ala (G12A) "] <- "G12A"
KCC[KCC=="c.34G>K    p.Gly12Cys (G12C) "] <- "G12C"
KCC[KCC=="c.34G>K    p.Gly12Cys (G12D) "] <- "G12D"
table(KCC)
names(KCC) <- rownames(CHEMORES_CLIN$objects$CHEMORES_CLIN)
WT_C <- predict(fit, t(VS3N[,names(which(KCC=="WT"))]),type="response")
G12C_C <- predict(fit, t(VS3N[,names(which(KCC=="G12C"))]),type="response")
G12V_C <- predict(fit, t(VS3N[,names(which(KCC=="G12V"))]),type="response")
boxplot(list(G12C=as.numeric(G12C_C),G12V=as.numeric(G12V_C),WT=as.numeric(WT_C)),ylab="RIS")
title("CHEMORES")
s <- svd(VS3N)
plot(s$v[,1],s$v[,2],col=as.factor(KCC),pch=20)

KCT <- loadEntity('syn1338272')
KCT <- KCT$objects$KRAS_LUAD
KCT <- KCT[intersect(colnames(VS2N),names(KCT))]
WT_T <- predict(fit, t(VS2N[,names(which(KCT=="WT"))]),type="response")
G12C_T <- predict(fit, t(VS2N[,names(which(KCT=="G12C"))]),type="response")
G12V_T <- predict(fit, t(VS2N[,names(which(KCT=="G12V"))]),type="response")
boxplot(list(G12C=as.numeric(G12C_T),G12V=as.numeric(G12V_T),WT=as.numeric(WT_T)),ylab="RIS")
title("TCGA_LUNG")
s <- svd(VS2N)
plot(s$v[,1],s$v[,2],col=as.factor(KCT),pch=20)

#############################################################################################
# train the KRAS PATIENTS model 
#############################################################################################
TEST <- cbind(VSN,VS4N)
KTEST <- c(VS_y,VS4_y)
set.seed(1234567)
cv.fit <- cv.glmnet(t(TEST), factor(KTEST), nfolds=20, alpha=.1, family="binomial")
plot(cv.fit)
fit <- glmnet(x=t(TEST),y=KTEST,family="binomial",alpha=.1,lambda=cv.fit$lambda.1se)

#############################################################################################
# boxplot the KRAS PATIENTS model across the different codons (CHEMORES and TCGA)
#############################################################################################
KCC <- CHEMORES_CLIN$objects$CHEMORES_CLIN$KRAS..NM_004448.2..exons.2_3
names(table(KCC))
KCC[KCC=="c.183A>W  p.Gln61His     (Q61H)    rs17851045 exon 3"] <- "Q61H"
KCC[KCC=="c.35G>K    p.Gly12Val (G12V) "] <- "G12V"
KCC[KCC=="c.35G>S    p.Gly12Ala (G12A) "] <- "G12A"
KCC[KCC=="c.34G>K    p.Gly12Cys (G12C) "] <- "G12C"
KCC[KCC=="c.35G>R    p.Gly12Asp (G12D) "] <- "G12D"
table(KCC)
names(KCC) <- rownames(CHEMORES_CLIN$objects$CHEMORES_CLIN)
WT_C <- predict(fit, t(VS3N[,names(which(KCC=="WT"))]),type="response")
G12C_C <- predict(fit, t(VS3N[,names(which(KCC=="G12C"))]),type="response")
G12V_C <- predict(fit, t(VS3N[,names(which(KCC=="G12V"))]),type="response")
par(mfrow=c(1,2))
boxplot(list(G12C=as.numeric(G12C_C),G12V=as.numeric(G12V_C),WT=as.numeric(WT_C)),ylab="RIS")
title("CHEMORES")
s <- svd(VS3N)
#plot(s$v[,1],s$v[,2],col=as.factor(KCC),pch=20)

KCT <- loadEntity('syn1338272')
KCT <- KCT$objects$KRAS_LUAD
KCT <- KCT[intersect(colnames(VS2N),names(KCT))]
WT_T <- predict(fit, t(VS2N[,names(which(KCT=="WT"))]),type="response")
G12C_T <- predict(fit, t(VS2N[,names(which(KCT=="G12C"))]),type="response")
G12V_T <- predict(fit, t(VS2N[,names(which(KCT=="G12V"))]),type="response")
boxplot(list(G12C=as.numeric(G12C_T),G12V=as.numeric(G12V_T),WT=as.numeric(WT_T)),ylab="RIS")
title("TCGA_LUNG")
s <- svd(VS2N)
#plot(s$v[,1],s$v[,2],col=as.factor(KCT),pch=20)

#####################################################################################################
### train G12C and G12V models from chemores +TCGA
#####################################################################################################
TEST <- cbind(VS2N,VS3N)
KCT <- loadEntity('syn1338272')
KCT <- KCT$objects$KRAS_LUAD
KCC <- CHEMORES_CLIN$objects$CHEMORES_CLIN$KRAS..NM_004448.2..exons.2_3
names(table(KCC))
KCC[KCC=="c.183A>W  p.Gln61His     (Q61H)    rs17851045 exon 3"] <- "Q61H"
KCC[KCC=="c.35G>K    p.Gly12Val (G12V) "] <- "G12V"
KCC[KCC=="c.35G>S    p.Gly12Ala (G12A) "] <- "G12A"
KCC[KCC=="c.34G>K    p.Gly12Cys (G12C) "] <- "G12C"
KCC[KCC=="c.35G>R    p.Gly12Asp (G12D) "] <- "G12D"
table(KCC)
names(KCC) <- rownames(CHEMORES_CLIN$objects$CHEMORES_CLIN)
KTEST <- c(KCT,KCC)
KTEST <- KTEST[colnames(TEST)]
table(KTEST)

KTEST_G12C <- ifelse(KTEST=="G12C",1,0)
table(KTEST_G12C)
cv.fit <- cv.glmnet(x=t(TEST),y=KTEST_G12C,nfolds=20,alpha=.1, family="binomial")
plot(cv.fit)
fit_G12C <- glmnet(x=t(TEST),y=KTEST_G12C,family="binomial",alpha=.1,lambda=cv.fit$lambda.min)

KTEST_G12V <- ifelse(KTEST=="G12V",1,0)
table(KTEST_G12V)
cv.fit <- cv.glmnet(x=t(TEST),y=KTEST_G12V,nfolds=20,alpha=.1, family="binomial")
plot(cv.fit)
fit_G12V <- glmnet(x=t(TEST),y=KTEST_G12V,family="binomial",alpha=.1,lambda=cv.fit$lambda.min)

#################################################################################################################
# apply the models in the CCLE database and correlate them with drug sensitivity
#################################################################################################################


ccle_response@data <- ccle_response@data[colnames(TS_EXP),]
drug <- ccle_response@data
sapply(colnames(drug),function(x){table(is.na(drug[,x]))})


d_G12C <- WT_T <- predict(fit_G12C, t(TS_EXP),type="response")
d_G12C

d_G12V <- WT_T <- predict(fit_G12V, t(TS_EXP),type="response")
d_G12V

# correlate the G12C and G12V with the drug response
p <- c()
D <- c()
G12C_SENS <- sapply(colnames(drug),function(x){
p <- cor.test(drug[,x],d_G12C,method="spearman",use="pairwise.complete.obs") 
D <- rbind(D,c(p$estimate,p$p.value))
})
rownames(G12C_SENS) <- c("rho","p.value")
G12C_SENS <- as.data.frame(t(G12C_SENS))
G12C_SENS <- G12C_SENS[sort(G12C_SENS$rho,index.return=TRUE)$ix,]

p <- c()
D <- c()
G12V_SENS <- sapply(colnames(drug),function(x){
  p <- cor.test(drug[,x],d_G12V,method="spearman",use="pairwise.complete.obs") 
  D <- rbind(D,p$estimate,c(p$p.value))
})
rownames(G12V_SENS) <- c("rho","p.value")
G12V_SENS <- as.data.frame(t(G12V_SENS))
G12V_SENS <- G12V_SENS[sort(G12V_SENS$rho,index.return=TRUE)$ix,]
