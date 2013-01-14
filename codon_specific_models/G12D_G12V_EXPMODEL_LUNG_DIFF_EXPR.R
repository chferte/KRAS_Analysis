# Charles Ferté
# Sage Bionetworks
# 14 Sept 2012


# train a model of differential gene expression per KRAS codon in Lung Adenocarcinoma
# both in G12C and G12V in CHEMORES and in TCGA LUAD
# chechk if there is some signal !

options(stringsAsFactors=FALSE)

library(affy)
library(corpcor)
library(lattice)
library(limma)


###############################################################
# run the analysis for the G12D G12V in the TCGA LUAD  
###############################################################

load("/home/cferte/FELLOW/cferte/KRAS_Project/OBJECTS/LUAD_EXP.RData")
load("/home/cferte/FELLOW/cferte/KRAS_Project/OBJECTS/KRAS_LUAD_STATUS.RData")
table(KRAS_LUAD)

# keep only the patients with available KRAS WT or G12C 
tmp <- intersect(names(KRAS_LUAD),colnames(LUAD_EXP))
LUAD_EXP <- LUAD_EXP[,tmp]
KRAS_LUAD <- KRAS_LUAD[tmp]
table(KRAS_LUAD)

KRAS_LUAD <- ifelse(KRAS_LUAD=="WT",0,1)
fit <- eBayes(lmFit(LUAD_EXP,model.matrix(~KRAS_LUAD)))
par(mfrow=c(1,2))
hist(fit$p.value[,2],breaks=100)
hist(p.adjust(fit$p.value[,2],method="BH"),breaks=200)

# run the permutation test with n=100,
# with alpha =.001
# compute the FDR
PERMUT <- replicate(100, eBayes(lmFit(LUAD_EXP,
      model.matrix(~KRAS_LUAD[sample(length(KRAS_LUAD))])))$p.value[,2])
alpha <- .001
N <- mean(apply(PERMUT,2,function(x){sum(x<alpha)}))
M <- sum(fit$p.value[,2]<alpha)
FDR <- N/M
paste("FDR=",FDR, "(with alpha=",alpha,")")


load("/home/cferte/FELLOW/cferte/KRAS_Project/OBJECTS/LUAD_EXP.RData")
load("/home/cferte/FELLOW/cferte/KRAS_Project/OBJECTS/KRAS_LUAD_STATUS.RData")
table(KRAS_LUAD)

# keep only the patients with available KRAS WT or G12C 
KRAS_LUAD <- KRAS_LUAD[which(KRAS_LUAD %in% c("G12V","G12C"))]
tmp <- intersect(names(KRAS_LUAD),colnames(LUAD_EXP))
LUAD_EXP <- LUAD_EXP[,tmp]
KRAS_LUAD <- KRAS_LUAD[tmp]
table(KRAS_LUAD)

fit <- eBayes(lmFit(LUAD_EXP,model.matrix(~KRAS_LUAD)))
hist(fit$p.value[,2],breaks=100)
hist(p.adjust(fit$p.value[,2],method="BH"),breaks=100)


G12C_LUAD_PVAL <- fit$p.value[,2]
G12C_LUAD_COEF <- fit$coefficients[,2]
rm(fit,tmp,KRAS_LUAD,LUAD_EXP)

###############################################################
# run the analysis for the G12V vs WT in the TCGA LUAD  
###############################################################

load("/home/cferte/FELLOW/cferte/KRAS_Project/OBJECTS/LUAD_EXP.RData")
load("/home/cferte/FELLOW/cferte/KRAS_Project/OBJECTS/KRAS_LUAD_STATUS.RData")

# keep only the patients with available KRAS WT or G12V (= 37 pts)
KRAS_LUAD <- KRAS_LUAD[which(KRAS_LUAD %in% c("WT","G12V"))]


tmp <- intersect(names(KRAS_LUAD),colnames(LUAD_EXP))
LUAD_EXP <- LUAD_EXP[,tmp]
KRAS_LUAD <- KRAS_LUAD[tmp]
table(KRAS_LUAD)

fit <- eBayes(lmFit(LUAD_EXP,model.matrix(~KRAS_LUAD)))

G12V_LUAD_PVAL <- fit$p.value[,2]
G12V_LUAD_COEF <- fit$coefficients[,2]
rm(fit,tmp,KRAS_LUAD,LUAD_EXP)


###############################################################
# load the CHEMORES data  
###############################################################

load("/home/cferte/FELLOW/cferte/KRAS_Project/OBJECTS/CHEMORES_EXP.RData")
load("/home/cferte/FELLOW/cferte/KRAS_Project/OBJECTS/CHEMORES_CLIN.RData")
load("/home/cferte/FELLOW/cferte/KRAS_Project/OBJECTS/KRAS_CHEMORES_STATUS.RData")

# extract the ADENOS only from CHEMORES 
CHEMORES_CLIN <- CHEMORES_CLIN[ CHEMORES_CLIN$Histology=="AC",]
CHEMORES_EXP <- CHEMORES_EXP[,rownames(CHEMORES_CLIN)]
KRAS_CHEMORES <- KRAS_CHEMORES[colnames(CHEMORES_EXP)]

# keep only the patients with available KRAS G12V or G12D (= 43 pts)
KRAS_CHEMORES <- KRAS_CHEMORES[which(KRAS_CHEMORES %in% c("G12V","G12D"))]
table(KRAS_CHEMORES)

CHEMORES_CLIN <- CHEMORES_CLIN[names(KRAS_CHEMORES),]
CHEMORES_EXP <- CHEMORES_EXP[,names(KRAS_CHEMORES)]

fit <- eBayes(lmFit(CHEMORES_EXP,model.matrix(~KRAS_CHEMORES)))

G12C_CHEMORES_PVAL <- fit$p.value[,2]
G12C_CHEMORES_COEF <- fit$coefficients[,2]
rm(CHEMORES_CLIN,CHEMORES_EXP,KRAS_CHEMORES,fit)



##############################################################
# get the interesect between the G12C
##############################################################
tmp <- names(G12C_LUAD_COEF)[intersect(which(G12C_CHEMORES_PVAL<.05),which(G12C_LUAD_PVAL<.05))]
paste(tmp,collapse=" ")

##############################################################
# get the interesect between the G12C
##############################################################
tmp1 <- names(G12V_LUAD_COEF)[intersect(which(G12V_CHEMORES_PVAL<.05),which(G12C_LUAD_PVAL<.05))]
paste(tmp1,collapse=" ")

save(G12C_CHEMORES_COEF,G12C_CHEMORES_PVAL,G12C_LUAD_COEF,G12C_LUAD_PVAL,file="/home/cferte/FELLOW/cferte/KRAS_Project/OBJECTS/G12C_CODON_LUNG_DIFF_EXP_PCVAL.RData")
save(G12V_CHEMORES_COEF,G12V_CHEMORES_PVAL,G12V_LUAD_COEF,G12V_LUAD_PVAL,file="/home/cferte/FELLOW/cferte/KRAS_Project/OBJECTS/G12V_CODON_LUNG_DIFF_EXP_PCVAL.RData")


