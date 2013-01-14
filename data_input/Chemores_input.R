# charles ferté
# sage bionetworks
# November 13th

#######################################################################################################################
# load packages and synapse
#######################################################################################################################

#load the different packages
options(stringsAsFactors=FALSE)
library(affy)
library(limma)
library(snm)
source("/home/cferte/FELLOW/cferte/KRAS_Project/JUSTIN_PREDICT_CCLE/code/lung_analysis_functions.R")
synapseLogin("charles.ferte@sagebase.org","charles")

#######################################################################################################################
# normalize the data
#######################################################################################################################
targets <- readTargets("/external-data/DAT_042__chemores/2012-05/GE/chemores_GE_target.txt")


#######################################################################################################################
# 3. Use LIMMA's read.maimages function to load your data into an RGList object. 
# Remember to set the path to the location where your Agilent FE text files are stored.
#######################################################################################################################
RG <- read.maimages(targets, path="/external-data/DAT_042__chemores/2012-05/GE/data/", source="agilent.median")


#######################################################################################################################
# 4. Subtract the background (optional, but recommended). 
# Setting an offset >1 ensures there will be no negative A values after log-transforming the data.
#######################################################################################################################
RG.bc <- backgroundCorrect(RG, method="normexp", offset=70)

#######################################################################################################################
#5. Normalize within arrays. Select a method appropriate for your data - see the user guide for details.
#######################################################################################################################
MA <- normalizeWithinArrays(RG.bc, method="loess",iterations=6)

#######################################################################################################################
# 6. Use the avereps function to average replicate spots
#######################################################################################################################
#load("/external-data/DAT_042__chemores/2012-05/GE/MA_loessNorm.RData")
MA.avg <- avereps(MA, ID=MA$genes$ProbeName)

#######################################################################################################################
# get rid of the dye bias using the dyes swap : (T/N-N/T)/2
#######################################################################################################################
MA.dcor <- (MA.avg$M[,(1:123)] -MA.avg$M[,(124:246)])/2  

#######################################################################################################################
# make all objects be cohenrent as regards to the samples
#######################################################################################################################
CHEMORES_EXP <- MA.dcor
targets <- readTargets("/external-data/DAT_042__chemores/2012-05/GE/chemores_GE_target.txt")
colnames(CHEMORES_EXP) <- targets$PatientID[1:123]
CHEMORES_CLIN <- loadEntity('syn1337568')
CHEMORES_CLIN <- CHEMORES_CLIN$objects$CHEMORES_CLIN
tmp <- intersect(colnames(CHEMORES_EXP),CHEMORES_CLIN$Patient_ID)
rownames(CHEMORES_CLIN) <- CHEMORES_CLIN$Patient_ID
CHEMORES_CLIN <- CHEMORES_CLIN[tmp,]
CHEMORES_EXP <- CHEMORES_EXP[,tmp]

#######################################################################################################################
# create KRAS_CHEMORES object 
#######################################################################################################################
KRAS_CHEMORES <- CHEMORES_CLIN$KRAS..NM_004448.2..exons.2_3
names(KRAS_CHEMORES) <- rownames(CHEMORES_CLIN)
KRAS_CHEMORES[grep(pattern="G12C",KRAS_CHEMORES)] <- "G12C"
KRAS_CHEMORES[grep(pattern="G12D",KRAS_CHEMORES)] <- "G12D"
KRAS_CHEMORES[grep(pattern="G12V",KRAS_CHEMORES)] <- "G12V"
KRAS_CHEMORES[grep(pattern="G12A",KRAS_CHEMORES)] <- "G12A"
KRAS_CHEMORES[grep(pattern="Q61H",KRAS_CHEMORES)] <- "Q61H"
KRAS_CHEMORES[grep(pattern="WT",KRAS_CHEMORES)] <- "WT"
table(KRAS_CHEMORES)


#######################################################################################################################
# get rid of the KRAS =NA
#######################################################################################################################
tmp <- names(KRAS_CHEMORES)[!is.na(KRAS_CHEMORES)]
KRAS_CHEMORES <- KRAS_CHEMORES[tmp]
CHEMORES_EXP <- CHEMORES_EXP[,tmp]
CHEMORES_CLIN <- CHEMORES_CLIN[tmp,]

#######################################################################################################################
# retrieve the annotations of CHEMORES_EXP
#######################################################################################################################
rat <- read.delim(file="/external-data/DAT_042__chemores/2012-05/GE/array_design/014850_D_AA_20070207.txt",header=T, as.is=T)
rownames(rat) <- rat$ProbeID
identical(rownames(CHEMORES_EXP)[match(rownames(rat),rownames(CHEMORES_EXP))],rownames(rat))
CHEMORES_EXP <- CHEMORES_EXP[match(rownames(rat),rownames(CHEMORES_EXP)),]
identical(rownames(rat),rownames(CHEMORES_EXP))
rownames(CHEMORES_EXP) <- rat$GeneSymbol
table(rownames(CHEMORES_EXP)=="")
CHEMORES_EXP <- CHEMORES_EXP[rownames(CHEMORES_EXP)!="",]


###################################################################################
# explore the signal
###################################################################################
tmp <- ifelse(KRAS_CHEMORES=="WT",1,0)
s <- svd(CHEMORES_EXP)
fit2 <- eBayes(lmFit(CHEMORES_EXP,model.matrix(~tmp)))
hist(fit2$p.value[,2],breaks=100)
table(fit2$p.value[,2]<.05)

Cancer.cell.rate <- ifelse(as.numeric(CHEMORES_CLIN$Cancer_Cell_Rate)<50,0,1)
Cancer.cell.rate[is.na(Cancer.cell.rate)] <- 0
names(Cancer.cell.rate) <- rownames(CHEMORES_CLIN)
fit2 <- eBayes(lmFit(CHEMORES_EXP,model.matrix(~tmp + Cancer.cell.rate)))
hist(fit2$p.value[,2],breaks=100)
table(fit2$p.value[,2]<.05)

tmp2 <- rownames(CHEMORES_CLIN)[which(CHEMORES_CLIN$Histology!="SCC")]
tmp <- ifelse(KRAS_CHEMORES[tmp2]=="WT",1,0)
Cancer.cell.rateACLCC <- Cancer.cell.rate[tmp2] 
fit2 <- eBayes(lmFit(CHEMORES_EXP[,tmp2],model.matrix(~tmp)))
hist(fit2$p.value[,2],breaks=100)
table(fit2$p.value[,2]<.05)

# 
# ###################################################################################
# # perform supervised normalization using snm 
# ###################################################################################
# 
# # the cancer cell rate (<50% or >= 50%)  is a  confounder of the expression
# Cancer.cell.rate <- ifelse(as.numeric(CHEMORES_CLIN$Cancer_Cell_Rate)<50,0,1)
# Cancer.cell.rate[is.na(Cancer.cell.rate)] <- 0
# s <- svd(CHEMORES_EXP)
# plot(s$v[,1],s$v[,2],col=rainbow(2)[as.factor(Cancer.cell.rate)],pch=20)
# adj.var <- model.matrix(~ Cancer.cell.rate )
# 
# # we know that KRAS and Histology are biological & study variables of interest 
# KRAS <- ifelse(KRAS_CHEMORES=="WT",0,1)
# plot(s$v[,1],s$v[,2],col=rainbow(6)[as.factor(KRAS_CHEMORES)],pch=20)
# plot(s$v[,1],s$v[,2],col=rainbow(4)[as.factor(CHEMORES_CLIN$Histology)],pch=20)
# 
# bio.var <- model.matrix(~ KRAS_CHEMORES + CHEMORES_CLIN$Histology)
# 
# # run snm
# snm.fit <- snm(CHEMORES_EXP, 
#                bio.var=bio.var, 
#                adj.var=adj.var, 
#                rm.adj=TRUE)
# 
# new.expr <- snm.fit$norm.dat
# rownames(new.expr) <- rownames(CHEMORES_EXP)
# colnames(new.expr) <- colnames(CHEMORES_EXP)
# 
# s <- svd(new.expr)
# plot(s$v[,1],s$v[,2],col=rainbow(2)[as.factor(Cancer.cell.rate)],pch=20)
# plot(s$v[,1],s$v[,2],col=rainbow(2)[as.factor(KRAS)],pch=20)
# plot(s$v[,1],s$v[,2],col=rainbow(6)[as.factor(KRAS_CHEMORES)],pch=20)
# plot(s$v[,1],s$v[,2],col=rainbow(4)[as.factor(CHEMORES_CLIN$Histology)],pch=20)
# 
# CHEMORES_EXP <- new.expr

###################################################################################
# save the CHEMORES objects 
###################################################################################

save(file= "/home/cferte/FELLOW/cferte/KRAS_Project/OBJECTS/CHEMORES_EXP.RData",CHEMORES_EXP)
save(file= "/home/cferte/FELLOW/cferte/KRAS_Project/OBJECTS/KRAS_CHEMORES_STATUS.RData",KRAS_CHEMORES)
save(file= "/home/cferte/FELLOW/cferte/KRAS_Project/OBJECTS/CHEMORES_CLIN.RData",CHEMORES_CLIN)


###################################################################################
# do not run 
###################################################################################

###################################################################################
# identification of flag proibes with agilent
###################################################################################

# 
# 
# 
# setwd("/external-data/DAT_042__chemores/2012-05/GE/data/")
# filenames <- dir()
# for(i in c(filenames[1:10])){
#   
# 
# FileName <- i
# # Step2: Load the second block of data to collect Intensity values and QC columns
# cat("\nStep2: Chargement des données mArray:")
# eset <- read.csv(paste(FileName, sep = ""), header = T, skip = 9, sep = "\t")
# 
# keepCol <- which(as.character(colnames(eset)) %in% 
#   c(	"ProbeName", "SystematicName",
#      "gMedianSignal", "rMedianSignal",
#      "gIsSaturated", "rIsSaturated",
#      "gIsFeatNonUnifOL", "rIsFeatNonUnifOL",
#      "gIsWellAboveBG", "rIsWellAboveBG"))
# 
# eset <- eset[which(substr(eset$ProbeName, 1, 1)== "A"), keepCol]  	# Select QC columns and rows containing A_xxx probes
# 
# # Step4: Flags suppression
# Flags = TRUE
# pflags <- NA
# if(Flags){
#   flags.info <- eset[,5:10]
#   flags <- which(	flags.info[,1]== 1 | 		    # 1 = gTsSaturated invalide
#                   flags.info[,2]== 1 | 				# 1 = rTsSaturated invalide
#                   flags.info[,3]== 1 | 				# 1 = gIsFeatureNonUnifOL invalide
#                   flags.info[,4]== 1  #|				# 1 = rIsFeatureNonUnifOL invalide
#                   #flags.info[,5]== 0 | 				# 0 = gIsWellAboveBG invalide
#                   #flags.info[,6]== 0
#                   )					# 0 = rIsWellAboveBG invalide
#   
#   eset$gMedianSignal[flags] = NA
#   eset$rMedianSignal[flags] = NA
#   #eset<- eset[-flags,]
#   pflags <- round(length(flags)/nrow(eset)*100, 2)
#   cat(pflags, "%\n")
#   }
# }
# 
# 
# 
# 
# for(i in c(filenames[1:10])){
#   FileName <- i
#   # Step2: Load the second block of data to collect Intensity values and QC columns
#   eset <- read.csv(paste(FileName, sep = ""), header = T, skip = 9, sep = "\t")
#   
#   keepCol <- which(as.character(colnames(eset)) %in% 
#     c(  "ProbeName", "SystematicName",
#         "gMedianSignal", "rMedianSignal",
#         "gIsSaturated", "rIsSaturated",
#         "gIsFeatNonUnifOL", "rIsFeatNonUnifOL",
#         "gIsWellAboveBG", "rIsWellAboveBG"))
#   
#     eset <- eset[which(substr(eset$ProbeName, 1, 1)== "A"), keepCol]  	# Select QC columns and rows containing A_xxx probes
#   
#   # Step4: Flags suppression
#     flags.info <- eset[,5:10]
#     cat(i, '\n')
#     cat('flag1: ', round(length(which(flags.info[,1]== 1))/nrow(eset)*100, 2), '\n')
#     cat('flag2: ', round(length(which(flags.info[,2]== 1))/nrow(eset)*100, 2), '\n')
#     cat('flag3: ', round(length(which(flags.info[,3]== 1))/nrow(eset)*100, 2), '\n')
#     cat('flag4: ', round(length(which(flags.info[,4]== 1))/nrow(eset)*100, 2), '\n')
#     cat('flag5: ', round(length(which(flags.info[,5]== 0))/nrow(eset)*100, 2), '\n')
#     cat('flag6: ', round(length(which(flags.info[,6]== 0))/nrow(eset)*100, 2), '\n')
# }