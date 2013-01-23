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

# load the KRAS_CCLE
ccle_kras <- loadEntity('syn1443160')
KRAS_CCLE <- ccle_kras$objects$KRAS_CCLE

KRAS_CCLE <- substr((KRAS_CCLE),3,nchar(KRAS_CCLE))
KRAS_CCLE[KRAS_CCLE==""] <- "WT"
table(KRAS_CCLE)
tmp <- ifelse(KRAS_CCLE %in% c("G12A","G12S","G13D","G13C","Q61H","Q61K","Q61L"),"rare",KRAS_CCLE)
names(tmp) <- names(KRAS_CCLE)
KRAS_CCLE <- tmp
table(KRAS_CCLE)

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
tmp <- rownames(A)[which(rat>.6)]
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
tmp1 <- which(tmp>quantile(tmp,probs=.2))
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
x <- log(sort(SangerDrug[,MEK.sanger[1]]))
plot(x,main=paste(MEK.sanger[1]),ylab="IC50",pch=20,col="gray60")
abline(h=quantile(x,probs=.2,na.rm=TRUE),col="red")
x <- log(sort(SangerDrug[,MEK.sanger[2]]))
plot(x,main=paste(MEK.sanger[2]),ylab="IC50",pch=20,col="gray60")
abline(h=quantile(x,probs=.2,na.rm=TRUE),col="red")
x <- log(sort(SangerDrug[,MEK.sanger[3]]))
plot(x,main=paste(MEK.sanger[3]),ylab="IC50",pch=20,col="gray60")
abline(h=quantile(x,probs=.2,na.rm=TRUE),col="red")
x <- log(sort(SangerDrug[,MEK.sanger[4]]))
plot(x,main=paste(MEK.sanger[4]),ylab="IC50",pch=20,col="gray60")
abline(h=quantile(x,probs=.2,na.rm=TRUE),col="red")
x <- log(sort(ccle_drug[,MEK.ccle[1]]))
plot(x,main=paste(MEK.ccle[1]),ylab="IC50",pch=20,col="gray60")
abline(h=quantile(x,probs=.2,na.rm=TRUE),col="red")
x <- log(sort(ccle_drug[,MEK.ccle[2]]))
plot(x,main=paste(MEK.ccle[2]),ylab="IC50",pch=20,col="gray60")
abline(h=quantile(x,probs=.2,na.rm=TRUE),col="red")
title(main= "distribution of the IC50 across the MEK inhibitors in SANGER & CCLE db",outer=TRUE)

# identify the cells that are evaluated for MEK inhibitors in sanger
MEK.cells.sanger <- rownames(SangerDrug)[-unique(c(which(is.na(SangerDrug[,MEK.sanger[1]])), which(is.na(SangerDrug[,MEK.sanger[2]])), which(is.na(SangerDrug[,MEK.sanger[3]])), which(is.na(SangerDrug[,MEK.sanger[4]]))))]
M <- SangerDrug[MEK.cells.sanger,MEK.sanger]
MEK.cells.ccle <- rownames(ccle_drug)[-unique(c(which(is.na(ccle_drug[,MEK.ccle[1]])), which(is.na(ccle_drug[,MEK.ccle[2]]))))]
M2 <- ccle_drug[MEK.cells.ccle,MEK.ccle]

rownames(M)
rownames(M2)

require(car)
scatterplotMatrix(M,main="correlations in the IC50 of the MEK inhibitors in Sanger")
scatterplotMatrix(normalizeCyclicLoess(M),main="correlations in the IC50 of the MEK inhibitors in Sanger (Loess normalized)")
cor(M,method="spearman",use="pairwise.complete.obs")
cor(normalizeCyclicLoess(M),method="spearman",use="pairwise.complete.obs")

scatterplotMatrix(normalizeCyclicLoess(M2),main="correlations in the IC50 of the MEK inhibitors in ccle (Loess normalized)")
cor(normalizeCyclicLoess(M2),method="spearman",use="pairwise.complete.obs")


# binomial transformation of the independant variables
M3 <- apply(M,2,function(y){ifelse(y<quantile(y,probs=.2),1,0)})
M4 <- apply(M2,2,function(y){ifelse(y<quantile(y,probs=.2),1,0)})


# is there any signal in the differential expression
par(mfrow=c(2,3))
fit <- eBayes(lmFit(SANGER_EXP[,MEK.cells.sanger],model.matrix(~M3[,1])))
hist(fit$p.value[,2],breaks=100, main=paste("Sanger Expr ~",colnames(M3)[1]))
table(fit$p.value[,2]<.05)
fit <- eBayes(lmFit(SANGER_EXP[,MEK.cells.sanger],model.matrix(~M3[,2])))
hist(fit$p.value[,2],breaks=100, main=paste("Sanger Expr ~",colnames(M3)[2]))
table(fit$p.value[,2]<.05)
fit <- eBayes(lmFit(SANGER_EXP[,MEK.cells.sanger],model.matrix(~M3[,3])))
hist(fit$p.value[,2],breaks=100, main=paste("Sanger Expr ~",colnames(M3)[3]))
table(fit$p.value[,2]<.05)
fit <- eBayes(lmFit(SANGER_EXP[,MEK.cells.sanger],model.matrix(~M3[,4])))
hist(fit$p.value[,2],breaks=100, main=paste("Sanger Expr ~",colnames(M3)[4]))
table(fit$p.value[,2]<.05)
fit <- eBayes(lmFit(CCLE_EXP[,MEK.cells.ccle],model.matrix(~M4[,1])))
hist(fit$p.value[,2],breaks=100, main=paste("ccle Expr ~",colnames(M4)[1]))
table(fit$p.value[,2]<.05)
fit <- eBayes(lmFit(CCLE_EXP[,MEK.cells.ccle],model.matrix(~M4[,2])))
hist(fit$p.value[,2],breaks=100, main=paste("ccle Expr ~",colnames(M4)[2]))
table(fit$p.value[,2]<.05)
title(main="univariate difefrential expression for sensitivity to MEKi across SANGER & CCLE",outer=TRUE)


# ################################################################################################################
# # load the  mutations data
# ################################################################################################################

# ################################################################################################################
# first load the ccle data
# ################################################################################################################
MAF1 <- read.delim("/home/cferte/CCLE_hybrid_capture1650_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07.maf",header=TRUE)
MAF1 <- MAF1[which(MAF1$Hugo_Symbol!="Unknown"),]
MAF1 <- MAF1[,c("Hugo_Symbol","Tumor_Sample_Barcode","Protein_Change", "Genome_Change")]
MAF1 <- MAF1[(grep(pattern="LUNG",x=MAF1$Tumor_Sample_Barcode)),]
MAF1$Tumor_Sample_Barcode <- sapply(strsplit(split="_",x=MAF1$Tumor_Sample_Barcode),function(x){x[[1]]})
MAF1$Protein_Change[MAF1$Protein_Change==""] <- MAF1$Genome_Change[ MAF1$Protein_Change==""]

# create a bianry matrix full of zeros with colnames equal to the sampel names and rownames equal the MUTIDs
MATMUT<-matrix(0,nrow=length(unique(MAF1$Hugo_Symbol)),ncol=length(unique(MAF1$Tumor_Sample_Barcode)))
colnames(MATMUT) <- unique(MAF1$Tumor_Sample_Barcode)
rownames(MATMUT) <- unique(MAF1$Hugo_Symbol)

# assign the protein or genome change information to any sample mutated for any MUTID
for(i in rownames(MATMUT)){
  MATMUT[i,c(MAF1$Tumor_Sample_Barcode[which(MAF1$Hugo_Symbol==i)])] <- c(MAF1$Protein_Change[which(MAF1$Hugo_Symbol==i)])  
}
mutations.ccle <- MATMUT
rm(MAF1,MATMUT)

# make the names of the samples coherents for ccle
tmp <- intersect(rownames(M4),colnames(mutations.ccle))
mutations.ccle <- mutations.ccle[,tmp]
M2 <- M2[tmp,]
M4 <- M4[tmp,]

# second, load the mutation data from Sanger
mutations.sanger <- read.csv("/home/cferte/FELLOW/cferte/Sanger_gdsc_mutation_w3.csv",header=TRUE)
rownames(mutations.sanger) <- toupper(gsub(pattern="-",replacement="",mutations.sanger$Cell.Line))
mutations.sanger <- t(mutations.sanger[ grep(pattern="NSCLC",x=mutations.sanger$Tissue),])
mutations.sanger[grep(pattern="na",x=mutations.sanger)] <- NA
mutations.sanger[grep(pattern="wt",x=mutations.sanger)] <- "0"

# make the names of the genes coherents between sanger and ccle
tmp <- intersect(rownames(mutations.sanger),rownames(mutations.ccle))
mutations.ccle <- mutations.ccle[tmp,]
mutations.sanger <- mutations.sanger[tmp,]

# make the names of the samples coherents for sanger
tmp <- intersect(rownames(M3),colnames(mutations.sanger))
mutations.sanger <- mutations.sanger[,tmp]
M3 <- M3[tmp,]
M <- M[tmp,]


# get rid of the NA's in Sanger mutations matrix
KRAS_SANGER <- mutations.sanger["KRAS",]
KRAS_SANGER[KRAS_SANGER!=0] <- sapply(strsplit(KRAS_SANGER[KRAS_SANGER!=0],split="::"),function(x){x[[1]]})
mutations.sanger[grep(pattern="p.",x=mutations.sanger)] <- "1"

# transform mutations.sanger into numeric and get rid of NAs
tmp <- apply(mutations.sanger,2,as.numeric)
rownames(tmp) <- rownames(mutations.sanger)
tmp <- tmp[names(which(!is.na(tmp[,1]))),]
mutations.sanger <- tmp
rm(tmp)

# restrict (again) to the common genes analyzed for both sanger and ccle
tmp <- intersect(rownames(mutations.ccle),rownames(mutations.sanger))
mutations.sanger <- mutations.sanger[tmp,]
mutations.ccle <- mutations.ccle[tmp,]
rm(tmp)

# create the KRAS_CCLE
KRAS_CCLE <- as.character(mutations.ccle["KRAS",])
names(KRAS_CCLE) <- colnames(mutations.ccle)
table(KRAS_CCLE)
table(KRAS_SANGER)

# set KRAS_CCLE and KRAS_SANGER to be a factor with same levels across ccle and sanger
# KRAS_CCLE[KRAS_CCLE %in% c("p.Q61H","p.Q61K","p.Q61L")] <- "p.Q61" 
# KRAS_SANGER[KRAS_SANGER %in% c("p.Q61H","p.Q61K","p.Q61L")] <- "p.Q61" 
# KRAS_CCLE[KRAS_CCLE %in% c("p.G13C","p.G13D")] <- "p.G13" 
# KRAS_SANGER[KRAS_SANGER %in% c("p.G13C","p.G13D")] <- "p.G13" 
# KRAS_SANGER[KRAS_SANGER %in% c("p.G12S","p.G12F","p.G12D","p.G12A")] <- "rare" 
# KRAS_CCLE[KRAS_CCLE %in% c("p.G12S","p.G12F","p.G12D","p.G12A")] <- "rare" 

theseLevels  <- unique(c(KRAS_CCLE, KRAS_SANGER))
table(KRAS_CCLE)
table(KRAS_SANGER)
KRAS_CCLE <- factor(KRAS_CCLE, levels=theseLevels)
KRAS_SANGER <- factor(KRAS_SANGER,levels=theseLevels)

kras.info.ccle <- t(model.matrix(~ -1 + KRAS_CCLE))
rownames(kras.info.ccle) <- sub(pattern="_CCLE",replacement= "", x=rownames(kras.info.ccle), fixed=T)
colnames(kras.info.ccle) <- names(KRAS_CCLE)
kras.info.sanger <- t(model.matrix(~ -1 + KRAS_SANGER))
rownames(kras.info.sanger) <- sub(pattern="_SANGER",replacement= "", x=rownames(kras.info.sanger), fixed=T)
colnames(kras.info.sanger) <- names(KRAS_SANGER)

# transform then mutations.ccle into a binary matrix
mutations.ccle[mutations.ccle!="0"] <- "1"
tmp <- apply(mutations.ccle,2,as.numeric)
rownames(tmp) <- paste(rownames(mutations.ccle),"_mut",sep="")
colnames(tmp) <- colnames(mutations.ccle)
mutations.ccle <- tmp
rm(tmp)
rownames(mutations.sanger) <- paste(rownames(mutations.sanger),"_mut",sep="")

###################################################################################################################
# imput the cnv data
###################################################################################################################
cnv.ccle <- read.delim("/home/cferte/CCLE_copynumber_byGene_2012-09-29.txt",header=TRUE)
CNV1 <- cnv.ccle
cnv.ccle <- cnv.ccle[,grep(pattern="LUNG",colnames(cnv.ccle))]
rownames(cnv.ccle) <- CNV1$geneName
colnames(cnv.ccle) <- sub(pattern="_LUNG",replacement="",x=colnames(cnv.ccle))
cnv.ccle <- cnv.ccle[,rownames(M4)]

cnv.sanger <- read.delim("/external-data/DAT_032__Sanger_Cell_Lines/CNV/Sanger800_copy_number_by_gene.txt",header=TRUE,na.strings = c("NA"," ","NaN"))
cnv.sanger <- loadEntity('syn464292')

cnv.sanger <- cnv.sanger$objects$eset
sampleNames(cnv.sanger)
featureNames(cnv.sanger)
rownames(cnv.sanger) <- cnv.sanger$Gene.Symb
tmp <- sapply(strsplit(x=colnames(cnv.sanger),split="_"),function(x){x[[1]]})
tmp <- toupper(x=tmp)
tmp <- gsub(pattern=".",replacement="",x=tmp,fixed = TRUE)
colnames(cnv.sanger) <- tmp
tmp <- intersect(rownames(M3),colnames(cnv.sanger))
cnv.sanger <- cnv.sanger[,tmp]
M3 <- M3[tmp,]

# get rid of the Nas in cnv.sanger
# pblm here
x <- 1
tmp <- apply(cnv.sanger,2,function(x){rownames(cnv.sanger)[which(is.na(x))]})
tmp <- unique(unlist(tmp))
cnv.sanger <- cnv.sanger[tmp,]
table(is.na(cnv.sanger))

# rescale the cnv data so the ccle and sanger cnv are comparable
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


cnv.sanger <- normalize_to_X(rowMeans(cnv.ccle),apply(cnv.ccle,1,sd),cnv.sanger)


tmp <- intersect(rownames(cnv.sanger),rownames(cnv.ccle))




################################
# restrict to KRAS mutant samples only 
################################

KC <- colnames(kras.info.ccle)
KS <- colnames(kras.info.sanger)
#KC <- names(which(kras.info.ccle["KRAS0",]==1))
#KS <- names(which(kras.info.sanger["KRAS0",]==1))

###################################################################################################################
# GOLD standard: correlations between IC50 of CCLE and Sanger
###################################################################################################################
tmp <- intersect(rownames(M3),rownames(M4))
cor(M3[tmp,-1],M4[tmp,],method="spearman")
scatterplotMatrix(normalizeCyclicLoess(cbind(M[tmp,-1],M2[tmp,])))
title(main="correlations in the IC50 of the MEK inhibitors 
between CCLE & Sanger (Loess normalized)", outer=TRUE)

par(mfrow=c(1,2))
###################################################################################################################
# train our predictive model of MEK response in ccle
# using a penalized regression approach  
# with alpha=.1 (more ridge) and determine lambda using nfolds= 5
# the robustness of the model is increased by boostrapping (n=100)
# validate each model in BATTLE
###################################################################################################################



require(glmnet)
#ccle.top <- M4[,2]
ccle.top <- ifelse(apply(M4[KC,],1,sum)>0,1,0)
table(ccle.top)
names(ccle.top) <- rownames(M4[KC,])

N <- 20
fit <- c()
selected <- c()
yhat <- c()
models <- 0
i <- 0
bal <- c()
trainex <- c()
validex <- c()
while(models<N)
{
j <- c(names(which(ccle.top==1)),sample(names(which(ccle.top==0)),replace=TRUE))
trainex <- rbind(rbind(CCLE_EXP[,j],mutations.ccle[,j]),kras.info.ccle[,j])
cv.fit <- cv.glmnet(t(trainex), y=ccle.top[j], nfolds=3, alpha=.1)
fit <- glmnet(x=t(trainex),y=ccle.top[j],alpha=.1,lambda=cv.fit$lambda.1se,family="binomial")
if(length(which(abs(as.numeric(fit$beta))> 10^-5))>10)
{
i=i+1
print(i)
selected <- cbind(selected , as.numeric(fit$beta))

dev <- which(rownames(M3[KS,]) %in% intersect(j,rownames(M3[KS,])))
val <- rownames(M3[KS,])[-dev]
validex <- rbind(rbind(SANGER_EXP[,val],mutations.sanger[,val]),kras.info.sanger[,val])
yhat <- c(yhat,list(predict(fit, t(validex),type="response")))
models <- length(yhat)
} }

# assess and plot the performance
require(ROCR)

AUC_RDEA119 <- c()
for (i in c(1:N)){
Pred <- prediction(as.numeric(yhat[[i]]),as.numeric(M3[rownames(yhat[[i]]),1]))
Perf <- performance(prediction.obj=Pred,"tpr","fpr")
AUC <- performance(prediction.obj=Pred,"auc")
AUC_RDEA119 <- c(AUC_RDEA119,as.numeric(AUC@y.values))
}

AUC_CI.1040 <- c()
for (i in c(1:N)){
Pred <- prediction(as.numeric(yhat[[i]]),as.numeric(M3[rownames(yhat[[i]]),2]))
Perf <- performance(prediction.obj=Pred,"tpr","fpr")
AUC <- performance(prediction.obj=Pred,"auc")
AUC_CI.1040 <- c(AUC_CI.1040,as.numeric(AUC@y.values))
}

AUC_PD.0325901 <- c()
for (i in c(1:N)){
Pred <- prediction(as.numeric(yhat[[i]]),as.numeric(M3[rownames(yhat[[i]]),3]))
Perf <- performance(prediction.obj=Pred,"tpr","fpr")
AUC <- performance(prediction.obj=Pred,"auc")
AUC_PD.0325901 <- c(AUC_PD.0325901,as.numeric(AUC@y.values))
}

AUC_AZD6244 <- c()
for (i in c(1:N)){
Pred <- prediction(as.numeric(yhat[[i]]),as.numeric(M3[rownames(yhat[[i]]),4]))
Perf <- performance(prediction.obj=Pred,"tpr","fpr")
AUC <- performance(prediction.obj=Pred,"auc")
AUC_AZD6244 <- c(AUC_AZD6244,as.numeric(AUC@y.values))
}


boxplot(cbind(AUC_RDEA119,AUC_CI.1040,AUC_PD.0325901,AUC_AZD6244),outline=FALSE,ylab="prediction of sensitivity (AUC)",ylim=c(0.3,.9))
stripchart(list(RDEA119=AUC_RDEA119,CI.1040=AUC_CI.1040,PD.0325901=AUC_PD.0325901,AZD6244=AUC_AZD6244),add=TRUE,method="jitter",vertical=TRUE,col="royalblue",pch=20)
title( main=" bootstrapped elasticnet models of gene expression + hybrid capture sequencing
trained in 68 ccle treated by PD.0325901 or AZD6244
validation in 35 cell lines processed by the Sanger group",outer=TRUE)

# extract the biological meaning
rownames(selected) <- rownames(fit$beta)
mus <- selected
mus[mus!=0] <- 1
#plot(density(apply(mus,2,sum)))
rownames(selected)[which(apply(mus,1,sum)>quantile(apply(mus,1,sum),probs=.99))]
sort(apply(mus,1,sum),decreasing=TRUE)[1:100]




###################################################################################################################
# train our predictive model of MEK response in sanger
# using a penalized regression approach  
# with alpha=.1 (more ridge) and determine lambda using nfolds= 5
# the robustness of the model is increased by boostrapping (n=100)
# validate each model in BATTLE
###################################################################################################################


require(glmnet)
sanger.top <- ifelse(apply(M3[KS,-1],1,sum)>2,1,0)
#ccle.top <- M3[,4]
names(sanger.top) <- rownames(M3[KS,])

N <- 30
fit <- c()
selected <- c()
yhat <- c()
models <- 0
i <- 0
bal <- c()

while(models<N)
{
  j <- c(names(which(sanger.top==1)),sample(names(which(sanger.top==0)),replace=TRUE))
  
  cv.fit <- cv.glmnet(t(rbind(rbind(SANGER_EXP[,j],mutations.sanger[,j]),kras.info.sanger[,j])), y=sanger.top[j], nfolds=3, alpha=.1)
  fit <- glmnet(x=t(rbind(rbind(SANGER_EXP[,j],mutations.sanger[,j]),kras.info.sanger[,j])),y=sanger.top[j],alpha=.1,lambda=cv.fit$lambda.1se,family="binomial")
  if(length(which(abs(as.numeric(fit$beta))> 10^-5))>10)
  {
    i=i+1
    print(i)
    selected <- cbind(selected , as.numeric(fit$beta))
    
    
    dev <- which(rownames(M4[KC,]) %in% intersect(j,rownames(M4[KC,])))
    val <- rownames(M4[KC,])[-dev]
    yhat <- c(yhat,list(predict(fit, t(rbind(rbind(CCLE_EXP[,val],mutations.ccle[,val]),kras.info.ccle[,val])),type="response")))
    models <- length(yhat)
  } }

# assess and plot the performance
require(ROCR)

AUC_PD.0325901 <- c()
for (i in c(1:N)){
  Pred <- prediction(as.numeric(yhat[[i]]),as.numeric(M4[rownames(yhat[[i]]),1]))
  Perf <- performance(prediction.obj=Pred,"tpr","fpr")
  AUC <- performance(prediction.obj=Pred,"auc")
  AUC_PD.0325901 <- c(AUC_PD.0325901,as.numeric(AUC@y.values))
}

AUC_AZD6244 <- c()
for (i in c(1:N)){
  Pred <- prediction(as.numeric(yhat[[i]]),as.numeric(M4[rownames(yhat[[i]]),2]))
  Perf <- performance(prediction.obj=Pred,"tpr","fpr")
  AUC <- performance(prediction.obj=Pred,"auc")
  AUC_AZD6244 <- c(AUC_AZD6244,as.numeric(AUC@y.values))
}

boxplot(cbind(AUC_PD.0325901,AUC_AZD6244),outline=FALSE,ylab="prediction of sensitivity (AUC)",ylim=c(0.3,.9))
stripchart(list(PD.0325901=AUC_PD.0325901,AZD6244=AUC_AZD6244),add=TRUE,method="jitter",vertical=TRUE,col="royalblue",pch=20)
title( main=" bootstrapped elasticnet models of gene expression + hybrid capture sequencing
trained in 68 ccle treated by PD.0325901 or AZD6244
       validation in 35 cell lines processed by the Sanger group",outer=TRUE)

# extract the biological meaning
rownames(selected) <- rownames(fit$beta)
mus <- selected
mus[mus!=0] <- 1
#plot(density(apply(mus,2,sum)))
rownames(selected)[which(apply(mus,1,sum)>quantile(apply(mus,1,sum),probs=.99))]
sort(apply(mus,1,sum),decreasing=TRUE)[1:50]

