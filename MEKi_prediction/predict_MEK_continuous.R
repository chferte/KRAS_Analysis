# Charles Ferté
# Sage Bionetworks
# 14 Sept 2012


# train a model of MEK sensibility (independant response = bionomial)

#load the different packages
options(stringsAsFactors=FALSE)

library(affy)
library(corpcor)
library(lattice)
library(limma)
library(caret)
library(glmnet)
library(snm)
library(synapseClient)
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
# get rid of the probes that are the less variant
#####################################################################################

tmp <- apply(CCLE_EXP,1,sd)
plot(density(tmp))
hist(tmp)
tmp1 <- which(tmp>quantile(tmp,probs=.1))
CCLE_EXP <- CCLE_EXP[tmp1,]
SANGER_EXP <- SANGER_EXP[tmp1,]
rm(tmp1,tmp)

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
#rm(sanger.drugs,ccle.drugs)

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

# identify the cell lines that are evaluated for MEK inhibitors in Sanger and in CCLE, and produce matrices of IC50
MEK.cells.sanger <- rownames(SangerDrug)[-unique(c(which(is.na(SangerDrug[,MEK.sanger[1]])), which(is.na(SangerDrug[,MEK.sanger[2]])), which(is.na(SangerDrug[,MEK.sanger[3]])), which(is.na(SangerDrug[,MEK.sanger[4]]))))]
mek.ic50.sanger <- SangerDrug[MEK.cells.sanger,MEK.sanger]
MEK.cells.ccle <- rownames(ccle_drug)[-unique(c(which(is.na(ccle_drug[,MEK.ccle[1]])), which(is.na(ccle_drug[,MEK.ccle[2]]))))]
mek.ic50.ccle <- ccle_drug[MEK.cells.ccle,MEK.ccle]

rownames(mek.ic50.sanger)
rownames(mek.ic50.ccle)


# correlation between mek inhib of sanger and ccle
require(car)
scatterplotMatrix(mek.ic50.sanger,main="correlations in the IC50 of the MEK inhibitors in Sanger")
scatterplotMatrix(normalizeCyclicLoess(mek.ic50.sanger),main="correlations in the IC50 of the MEK inhibitors in Sanger (Loess normalized)")
cor(mek.ic50.sanger,method="spearman",use="pairwise.complete.obs")
cor(normalizeCyclicLoess(mek.ic50.sanger),method="spearman",use="pairwise.complete.obs")

scatterplotMatrix(normalizeCyclicLoess(mek.ic50.ccle),main="correlations in the IC50 of the MEK inhibitors in ccle (Loess normalized)")
cor(normalizeCyclicLoess(mek.ic50.ccle),method="spearman",use="pairwise.complete.obs")

tmp <- intersect(rownames(mek.ic50.ccle),rownames(mek.ic50.sanger))
cor(mek.ic50.sanger[tmp,],mek.ic50.ccle[tmp,],method="spearman",use="pairwise.complete.obs")


# is there any signal in the differential expression
par(mfrow=c(2,3))
fit <- eBayes(lmFit(SANGER_EXP[,MEK.cells.sanger],model.matrix(~mek.ic50.sanger[,1])))
hist(fit$p.value[,2],breaks=100, main=paste("Sanger Expr ~",colnames(mek.ic50.sanger)[1]))
table(fit$p.value[,2]<.05)
fit <- eBayes(lmFit(SANGER_EXP[,MEK.cells.sanger],model.matrix(~mek.ic50.sanger[,2])))
hist(fit$p.value[,2],breaks=100, main=paste("Sanger Expr ~",colnames(mek.ic50.sanger)[2]))
table(fit$p.value[,2]<.05)
fit <- eBayes(lmFit(SANGER_EXP[,MEK.cells.sanger],model.matrix(~mek.ic50.sanger[,3])))
hist(fit$p.value[,2],breaks=100, main=paste("Sanger Expr ~",colnames(mek.ic50.sanger)[3]))
table(fit$p.value[,2]<.05)
fit <- eBayes(lmFit(SANGER_EXP[,MEK.cells.sanger],model.matrix(~mek.ic50.sanger[,4])))
hist(fit$p.value[,2],breaks=100, main=paste("Sanger Expr ~",colnames(mek.ic50.sanger)[4]))
table(fit$p.value[,2]<.05)
fit <- eBayes(lmFit(CCLE_EXP[,MEK.cells.ccle],model.matrix(~mek.ic50.ccle[,1])))
hist(fit$p.value[,2],breaks=100, main=paste("CCLE Expr ~",colnames(mek.ic50.ccle)[1]))
table(fit$p.value[,2]<.05)
fit <- eBayes(lmFit(CCLE_EXP[,MEK.cells.ccle],model.matrix(~mek.ic50.ccle[,2])))
hist(fit$p.value[,2],breaks=100, main=paste("CCLE Expr ~",colnames(mek.ic50.ccle)[2]))
table(fit$p.value[,2]<.05)
title(main="univariate difefrential expression for sensitivity to MEKi across SANGER & CCLE",outer=TRUE)


# ################################################################################################################
# # load the  mutations data
# ################################################################################################################

# ################################################################################################################
# first load the ccle data
# ################################################################################################################
MAF1 <- read.delim("/home/cferte/cell_line_data/CCLE_hybrid_capture1650_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07.maf",header=TRUE)
MAF1 <- MAF1[which(MAF1$Hugo_Symbol!="Unknown"),]
MAF1 <- MAF1[,c("Hugo_Symbol","Tumor_Sample_Barcode","Protein_Change", "Genome_Change")]
MAF1 <- MAF1[(grep(pattern="LUNG",x=MAF1$Tumor_Sample_Barcode)),]
MAF1$Tumor_Sample_Barcode <- sapply(strsplit(split="_",x=MAF1$Tumor_Sample_Barcode),function(x){x[[1]]})
MAF1$Protein_Change[MAF1$Protein_Change==""] <- MAF1$Genome_Change[ MAF1$Protein_Change==""]

# create a binary matrix full of zeros with colnames equal to the sampel names and rownames equal the MUTIDs
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
tmp <- intersect(rownames(mek.ic50.ccle),colnames(mutations.ccle))
mutations.ccle <- mutations.ccle[,tmp]
mek.ic50.ccle <- mek.ic50.ccle[tmp,]

# second, load the mutation data from Sanger
mutations.sanger <- read.csv("/home/cferte/cell_line_data/Sanger_gdsc_mutation_w3.csv",header=TRUE)
rownames(mutations.sanger) <- toupper(gsub(pattern="-",replacement="",mutations.sanger$Cell.Line))
mutations.sanger <- t(mutations.sanger[ grep(pattern="NSCLC",x=mutations.sanger$Tissue),])
mutations.sanger[grep(pattern="na",x=mutations.sanger)] <- NA
mutations.sanger[grep(pattern="wt",x=mutations.sanger)] <- "0"

# make the names of the genes coherents between sanger and ccle
tmp <- intersect(rownames(mutations.sanger),rownames(mutations.ccle))
mutations.ccle <- mutations.ccle[tmp,]
mutations.sanger <- mutations.sanger[tmp,]

# make the names of the samples coherents for sanger
tmp <- intersect(rownames(mek.ic50.sanger),colnames(mutations.sanger))
mutations.sanger <- mutations.sanger[,tmp]
mek.ic50.sanger <- mek.ic50.sanger[tmp,]

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

# set KRAS_CCLE and KRAS_SANGER to be a factor with same levels across ccle and sanger
theseLevels  <- unique(c(KRAS_CCLE, KRAS_SANGER))
table(KRAS_CCLE)
table(KRAS_SANGER)
KRAS_CCLE <- factor(KRAS_CCLE, levels=theseLevels)
KRAS_SANGER <- factor(KRAS_SANGER,levels=theseLevels)
table(KRAS_CCLE)
table(KRAS_SANGER)
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

################################
# restrict to KRAS mutant samples only 
################################

KC <- names(which(kras.info.ccle["KRAS0",]==1))
KS <- names(which(kras.info.sanger["KRAS0",]==1))

################################################################################################################
# is there any signal in the differential expression
################################################################################################################

par(mfrow=c(2,3))
fit <- eBayes(lmFit(SANGER_EXP[,rownames(mek.ic50.sanger)],model.matrix(~mek.ic50.sanger[,1])))
hist(fit$p.value[,2],breaks=100, main=paste("Sanger Expr ~",colnames(mek.ic50.sanger)[1]))
table(fit$p.value[,2]<.05)
fit <- eBayes(lmFit(SANGER_EXP[,rownames(mek.ic50.sanger)],model.matrix(~mek.ic50.sanger[,2])))
hist(fit$p.value[,2],breaks=100, main=paste("Sanger Expr ~",colnames(mek.ic50.sanger)[2]))
table(fit$p.value[,2]<.05)
fit <- eBayes(lmFit(SANGER_EXP[,rownames(mek.ic50.sanger)],model.matrix(~mek.ic50.sanger[,3])))
hist(fit$p.value[,2],breaks=100, main=paste("Sanger Expr ~",colnames(mek.ic50.sanger)[3]))
table(fit$p.value[,2]<.05)
fit <- eBayes(lmFit(SANGER_EXP[,rownames(mek.ic50.sanger)],model.matrix(~mek.ic50.sanger[,4])))
hist(fit$p.value[,2],breaks=100, main=paste("Sanger Expr ~",colnames(mek.ic50.sanger)[4]))
table(fit$p.value[,2]<.05)
fit <- eBayes(lmFit(CCLE_EXP[,rownames(mek.ic50.ccle)],model.matrix(~mek.ic50.ccle[,1])))
hist(fit$p.value[,2],breaks=100, main=paste("CCLE Expr ~",colnames(mek.ic50.ccle)[1]))
table(fit$p.value[,2]<.05)
fit <- eBayes(lmFit(CCLE_EXP[,rownames(mek.ic50.ccle)],model.matrix(~mek.ic50.ccle[,2])))
hist(fit$p.value[,2],breaks=100, main=paste("CCLE Expr ~",colnames(mek.ic50.ccle)[2]))
table(fit$p.value[,2]<.05)
title(main="univariate difefrential expression for sensitivity to MEKi across SANGER & CCLE",outer=TRUE)



################################################################################################################
# is there any signal in the mutation
################################################################################################################

par(mfrow=c(2,3))
fit <- eBayes(lmFit(mutations.sanger[,rownames(mek.ic50.sanger)],model.matrix(~mek.ic50.sanger[,1])))
hist(fit$p.value[,2], main=paste("Sanger Expr ~",colnames(mek.ic50.sanger)[1]))
table(fit$p.value[,2]<.05)
fit <- eBayes(lmFit(mutations.sanger[,rownames(mek.ic50.sanger)],model.matrix(~mek.ic50.sanger[,2])))
hist(fit$p.value[,2], main=paste("Sanger Expr ~",colnames(mek.ic50.sanger)[2]))
table(fit$p.value[,2]<.05)
fit <- eBayes(lmFit(mutations.sanger[,rownames(mek.ic50.sanger)],model.matrix(~mek.ic50.sanger[,3])))
hist(fit$p.value[,2], main=paste("Sanger Expr ~",colnames(mek.ic50.sanger)[3]))
table(fit$p.value[,2]<.05)
fit <- eBayes(lmFit(mutations.sanger[,rownames(mek.ic50.sanger)],model.matrix(~mek.ic50.sanger[,4])))
hist(fit$p.value[,2], main=paste("Sanger Expr ~",colnames(mek.ic50.sanger)[4]))
table(fit$p.value[,2]<.05)
fit <- eBayes(lmFit(mutations.ccle[,rownames(mek.ic50.ccle)],model.matrix(~mek.ic50.ccle[,1])))
hist(fit$p.value[,2], main=paste("CCLE Expr ~",colnames(mek.ic50.ccle)[1]))
table(fit$p.value[,2]<.05)
fit <- eBayes(lmFit(mutations.ccle[,rownames(mek.ic50.ccle)],model.matrix(~mek.ic50.ccle[,2])))
hist(fit$p.value[,2], main=paste("CCLE Expr ~",colnames(mek.ic50.ccle)[2]))
table(fit$p.value[,2]<.05)
title(main="univariate mutational association for sensitivity to MEKi across SANGER & CCLE",outer=TRUE)


###################################################################################################################
# GOLD standard: correlations between IC50 of CCLE and Sanger
###################################################################################################################
tmp <- intersect(rownames(mek.ic50.ccle),rownames(mek.ic50.sanger))
cor(mek.ic50.ccle[tmp,],mek.ic50.sanger[tmp,],method="spearman")
scatterplotMatrix(normalizeCyclicLoess(cbind(mek.ic50.ccle[tmp,],mek.ic50.sanger[tmp,])))

ic50.train <- mek.ic50.sanger
ic50.val <- mek.ic50.ccle

# is there any signal in the differential expression
par(mfrow=c(2,3))
fit <- eBayes(lmFit(SANGER_EXP[,rownames(ic50.train)],model.matrix(~ic50.train[,1])))
hist(fit$p.value[,2],breaks=100, main=paste("ccle Expr ~",colnames(ic50.train)[1]))
table(fit$p.value[,2]<.05)
fit <- eBayes(lmFit(SANGER_EXP[,rownames(ic50.train)],model.matrix(~ic50.train[,2])))
hist(fit$p.value[,2],breaks=100, main=paste("ccle Expr ~",colnames(ic50.train)[2]))
table(fit$p.value[,2]<.05)
fit <- eBayes(lmFit(SANGER_EXP[,rownames(ic50.train)],model.matrix(~ic50.train[,3])))
hist(fit$p.value[,2],breaks=100, main=paste("ccle Expr ~",colnames(ic50.train)[3]))
table(fit$p.value[,2]<.05)
fit <- eBayes(lmFit(SANGER_EXP[,rownames(ic50.train)],model.matrix(~ic50.train[,4])))
hist(fit$p.value[,2],breaks=100, main=paste("ccle Expr ~",colnames(ic50.train)[4]))
table(fit$p.value[,2]<.05)
fit <- eBayes(lmFit(CCLE_EXP[,rownames(ic50.val)],model.matrix(~ic50.val[,1])))
hist(fit$p.value[,2],breaks=100, main=paste("Sanger Expr ~",colnames(ic50.val)[1]))
table(fit$p.value[,2]<.05)
fit <- eBayes(lmFit(CCLE_EXP[,rownames(ic50.val)],model.matrix(~ic50.val[,2])))
hist(fit$p.value[,2],breaks=100, main=paste("Sanger Expr ~",colnames(ic50.val)[2]))
table(fit$p.value[,2]<.05)
title(main="univariate difefrential expression for sensitivity to MEKi across SANGER & CCLE",outer=TRUE)




###################################################################################################################
# train our predictive model of MEK response in ccle
# using a penalized regression approach  
# with alpha=.1 (more ridge) and determine lambda using nfolds= 5
# the robustness of the model is increased by boostrapping (n=100)
# train in sanger cells
# validate each model in ccle
###################################################################################################################

# # get rid of the non informative cell lines
# tmp <- sapply(colnames(ic50.train),function(x){rownames(ic50.train)[which(ic50.train[,x]<quantile(ic50.train[,x],probs=.6) & ic50.train[,x]>quantile(ic50.train[,x],probs=.4))]})
# tmp <- as.character(tmp)
# trash <- which(rownames(ic50.train) %in% unique(tmp[duplicated(tmp)]))


par(mfrow=c(1,1))
require(glmnet)
N <- 30
fit <- c()
selected <- c()
yhat <- c()
models <- 0
i <- 0
while(models<N)
{
  dev <- c()
  val <- c()
  j <- sample(rownames(ic50.train),replace=TRUE)
  table(j)
  #trainex <- SANGER_EXP[,j]
  trainex <- rbind(rbind(SANGER_EXP[,j],mutations.sanger[,j]),kras.info.sanger[,j])
  vec.train <- apply(normalizeCyclicLoess(ic50.train[j,]),1,mean)
  cv.fit <- cv.glmnet(t(trainex), y=vec.train,nfolds=5, alpha=.01)
  fit <- glmnet(x=t(trainex),y=vec.train,alpha=.01,lambda=cv.fit$lambda.1se)
  if(length(which(abs(as.numeric(fit$beta))> 10^-5))>10)
  {
    i=i+1
    print(i)
    selected <- cbind(selected , as.numeric(fit$beta))
    
    dev <- which(rownames(ic50.val) %in% intersect(j,rownames(ic50.val)))
    val <- rownames(ic50.val)
    #val <- rownames(ic50.val)[-dev]
    #intersect(j,val)
    #validex <- CCLE_EXP[,val]
    validex <- rbind(rbind(CCLE_EXP[,val],mutations.ccle[,val]),kras.info.ccle[,val])
    yhat <- c(yhat,list(predict(fit, t(validex))))
    models <- length(yhat)
  } }

rownames(selected) <- rownames(trainex)

y <- c()
for(i in c(1:length(yhat)))
{
  y <- cbind(y, yhat[[i]])
}
rownames(y) <- val

boxplot(list(PD0325901=cor(y,ic50.val[,1],method="spearman"), AZD6244=cor(y,ic50.val[,2],method="spearman")),outline=FALSE,ylim=c(0,1))
stripchart(list(PD0325901=cor(y,ic50.val[,1],method="spearman"), AZD6244=cor(y,ic50.val[,2],method="spearman")),method="jitter",vertical=TRUE,add=TRUE,col="royalblue",pch=20)
tmp <- intersect(rownames(ic50.train),rownames(ic50.val))
a <- cor(ic50.val[tmp,],apply(ic50.train[tmp,],1,mean),method="spearman")[1]
b <- cor(ic50.val[tmp,],apply(ic50.train[tmp,],1,mean),method="spearman")[2]
stripchart(list(PD0325901=a, AZD6244=b),pch=20,cex=3,col="orange",add=TRUE,vertical=TRUE)
abline(h=c(0,.2,.4,.6,.8,1),lty=2)

###################################################################################################################
# train our predictive model of MEK response in ccle
# using a penalized regression approach  
# with alpha=.1 (more ridge) and determine lambda using nfolds= 5
# the robustness of the model is increased by boostrapping (n=100)
# validate each model in BATTLE
###################################################################################################################

require(glmnet)
ccle.mean <- rank(M2[,1])
sanger.mean <- M[,2]
names(ccle.mean) <- rownames(M2)
names(sanger.mean) <- rownames(M)

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
j <- sample(rownames(M2),replace=TRUE)
trainex <- rbind(rbind(CCLE_EXP[,j],mutations.ccle[,j]),kras.info.ccle[,j])
cv.fit <- cv.glmnet(t(trainex), y=ccle.mean[j], nfolds=3, alpha=.1)
fit <- glmnet(x=t(trainex),y=ccle.mean[j],alpha=.1,lambda=cv.fit$lambda.1se,family="gaussian")
if(length(which(abs(as.numeric(fit$beta))> 10^-5))>10)
{
i=i+1
print(i)
selected <- cbind(selected , as.numeric(fit$beta))

dev <- which(rownames(M) %in% intersect(j,rownames(M)))

val <- rownames(M)[-dev]
validex <- rbind(rbind(SANGER_EXP[,val],mutations.sanger[,val]),kras.info.sanger[,val])
yhat <- c(yhat,list(predict(fit, t(validex),type="link")))
models <- length(yhat)
} }

# assess and plot the performance
length(yhat)
yhat[[4]]

colnames(M)
COR_AZD6244 <- c()
for (i in c(1:N)){
COR_AZD6244 <- c(COR_AZD6244,cor(as.numeric(yhat[[i]]),rank(as.numeric(M[rownames(yhat[[i]]),4])),method="spearman"))
}

# COR_CI.1040 <- c()
# for (i in c(1:N)){
#   COR_CI.1040 <- c(COR_CI.1040,cor(as.numeric(yhat[[i]]),as.numeric(M[rownames(yhat[[i]]),2]),method="spearman"))
# }

boxplot(COR_AZD6244)




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

