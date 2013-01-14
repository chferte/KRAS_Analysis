# charles Ferté
# 14 octobre 2012
# Sage Bionetworks

# look a the drug sensitivity per KRAS codon in Lung cell lines using ccle
# mutational profile in LUAD & CRC

#load the packages
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

#load the data from Synapse
ccle_eset <- loadEntity('syn1354643')
ccle_drug <- loadEntity('syn1354656')
ccle_kras <- loadEntity('syn1443160')

# create a DRUG object
ccle_kras <- ccle_kras$objects$KRAS_CCLE
DRUG <- ccle_drug$objects$ccle_drug@data[names(ccle_kras),]
DRUG$kras <- ccle_kras

table(ccle_kras)
# get rid of irinotecan since it comprise many missing data
DRUG <- DRUG[,-which(colnames(DRUG)=="Irinotecan")]

###########################################################################################
# describe the sensitivity to the drugs across G12C G12V
###########################################################################################

# look at the sensitivity for G12C and set m as a vector of median IC50 per drug
m <- c()
G12C_DRUG <- DRUG[DRUG$kras=="p.G12C",-24]
for(i in c(1:23)){
  m <- c(m,median(G12C_DRUG[,i],na.rm=T))
}
names(m) <- colnames(G12C_DRUG)

# look at the sensitivity for G12V and set n as a vector of median IC50 per drug
n <- c()
G12V_DRUG <- DRUG[DRUG$kras=="p.G12V",-24]
for(i in c(1:23)){
  n <- c(n,median(G12V_DRUG[,i],na.rm=T))
}
names(n) <- colnames(G12V_DRUG)

# look at the sensitivity for codon 61 and set K61 as a vector of median IC50 per drug
K61 <- c()
KRAS61_DRUG <- DRUG[ grep(DRUG$kras,pattern="61"),-24]
for(i in c(1:23)){
  K61 <- c(K61,median(KRAS61_DRUG[,i],na.rm=T))
}
names(K61) <- colnames(KRAS61_DRUG)

# look at the sensitivity for codon 12 and set K12 as a vector of median IC50 per drug
K12 <- c()
KRAS12_DRUG <- DRUG[ grep(DRUG$kras,pattern="12"),-24]
for(i in c(1:23)){
  K12 <- c(K12,mean(KRAS12_DRUG[,i],na.rm=T))
}
names(K12) <- colnames(KRAS12_DRUG)


par(mfrow=c(1,1),oma=c(5,2,0,2))
boxplot(G12C_DRUG[,sort(m,index.return=T)$ix],ylab="IC50",las=2,main="IC50 in 11 cell lines with KRAS G12C mutation")
boxplot(G12V_DRUG[,sort(n,index.return=T)$ix],ylab="IC50",las=2,main="IC50 in 6 cell lines with KRAS G12V mutation")
#boxplot(KRAS61_DRUG[,sort(K61,index.return=T)$ix],ylab="IC50",las=2,main="IC50 in 4 cell lines with KRAS G12V mutation")
#boxplot(KRAS12_DRUG[,sort(K12,index.return=T)$ix],ylab="IC50",las=2,main="IC50 in 21 cell lines with KRAS G12V mutation")

###########################################################################################
# draw a heatmap across different kras codons
###########################################################################################

# draw a heatmap with all the kras status
par(oma=c(2,2,2,2))
DRUG1 <- DRUG[,-c(24)]
require(gplots)
table(ccle_kras)
tmp <- sort(as.numeric(as.factor(ccle_kras)),index.return=TRUE)$ix
DRUG1 <- DRUG1[tmp,]
KRS1 <- ccle_kras[tmp]
heatmap.2(t(DRUG1),trace="none",col=greenred(100),ColSideColors=rainbow(11)[as.numeric(as.factor(KRS1))],Colv=NULL,scale="none")

# draw a heatmap with the KRAS mutant only
par(oma=c(2,2,2,2))
DRUG2 <- DRUG[DRUG$kras != "0",]
KRS2 <- DRUG2$kras
DRUG2 <- DRUG2[,-c(24)]
require(gplots)
table(KRS2)
tmp <- sort(as.numeric(as.factor(KRS2)),index.return=TRUE)$ix
DRUG2 <- DRUG2[tmp,]
KRS22 <- KRS2[tmp]
heatmap.2(t(DRUG2),trace="none",col=greenred(100),ColSideColors=rainbow(11)[as.numeric(as.factor(KRS22))],Colv=NULL,scale="none")
KRS22

# KRAS G12C and G12V only
par(oma=c(2,2,2,2))
DRUG2 <- DRUG[DRUG$kras %in% c("p.G12C","p.G12V"),]
KRS2 <- DRUG2$kras
DRUG2 <- DRUG2[,-c(24)]
require(gplots)
table(KRS2)
tmp <- sort(as.numeric(as.factor(KRS2)),index.return=TRUE)$ix
DRUG2 <- DRUG2[tmp,]
KRS22 <- KRS2[tmp]
heatmap.2(t(DRUG2),trace="none",col=greenred(100),ColSideColors=rainbow(11)[as.numeric(as.factor(KRS22))],scale="none",Colv=NULL)
heatmap.2(t(DRUG2),trace="none",col=greenred(100),ColSideColors=rainbow(11)[as.numeric(as.factor(KRS22))],scale="none")
KRS22

##########################################################################
# compute the difference of drug sensitivity between codon G12C and G12V:
##########################################################################



diff <- c()
DRUG1 <- DRUG[ DRUG$kras %in% c("p.G12C","p.G12V"),]
for(i in c(1:23)){
  diff <- c(diff,wilcox.test(DRUG1[,i]~as.factor(DRUG1$kras))$p.value)
}
names(diff) <- colnames(DRUG1)[1:23]
diff <- sort(diff)
which(diff<.05)

par(mfrow=c(2,2))
# boxplot the IC50 of the 4 drugs where there are significant differences between G12C and G12V
boxplot(DRUG1[,"RAF265"]~DRUG1$kras,main="RAF265 in ccle data",ylab="IC50")
stripchart(DRUG1[,"RAF265"]~DRUG1$kras,add=TRUE,vertical=TRUE,col="red",pch=20)

# boxplot the IC50 of the 4 drugs where there are significant differences between G12C and G12V
boxplot(DRUG1[,"PLX4720"]~DRUG1$kras,main="PLX4720 in ccle data",ylab="IC50")
stripchart(DRUG1[,"PLX4720"]~DRUG1$kras,add=TRUE,vertical=TRUE,col="red",pch=20)

# boxplot the IC50 of the 4 drugs where there are significant differences between G12C and G12V
boxplot(DRUG1[,"TKI258"]~DRUG1$kras,main="TKI258 in ccle data",ylab="IC50")
stripchart(DRUG1[,"TKI258"]~DRUG1$kras,add=TRUE,vertical=TRUE,col="red",pch=20)

# boxplot the IC50 of the 4 drugs where there are significant differences between G12C and G12V
boxplot(DRUG1[,"17-AAG"]~DRUG1$kras,main="17-AAG in ccle data",ylab="IC50")
stripchart(DRUG1[,"17-AAG"]~DRUG1$kras,add=TRUE,vertical=TRUE,col="red",pch=20)
