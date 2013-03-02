# Charles Ferté
# Sage Bionetworks
# 14 Feb 2012


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
synapseLogin("charles.ferte@sagebase.org","charles")

###############################################################
# load the sanger data (snm normalized) 
###############################################################
sanger_all <- loadEntity("syn1690438")
sanger_all <- sanger_all$objects$sanger_data
assign(x=names(sanger_all)[1],sanger_all[[1]])
assign(x=names(sanger_all)[2],sanger_all[[2]])
assign(x=names(sanger_all)[3],sanger_all[[3]])
assign(x=names(sanger_all)[4],sanger_all[[4]])
assign(x=names(sanger_all)[5],sanger_all[[5]])
assign(x=names(sanger_all)[6],sanger_all[[6]])

a <- rownames(sanger_drug)
b <- sapply(strsplit(x=colnames(sanger_exp),split="_"),function(x){x[[1]]})
tmp <- match(a,b)
sanger_drug$cell_new_id <- colnames(sanger_exp)[tmp]
sanger_drug <- sanger_drug[!is.na(sanger_drug$cell_new_id),]
rownames(sanger_drug) <- sanger_drug$cell_new_id
rm(a,b,tmp)
sanger_drug$cell_new_id <- NULL
colnames(sanger_drug) <- toupper(sub(pattern="_IC_50",replacement="",x=colnames(sanger_drug)))
colnames(sanger_drug) <- sub(pattern=".",replacement="-",x=colnames(sanger_drug),fixed=TRUE)
tmp <- rownames(sanger_drug)
sanger_drug <- apply(sanger_drug,2,function(x){as.numeric(sub(pattern=",",replacement=".",x=x,fixed=TRUE))})
rownames(sanger_drug) <- tmp
rm(tmp)

#########################################################################################################
## Make the data coherent between all the datasets
#########################################################################################################
tmp <- intersect(colnames(sanger_cnv),colnames(sanger_exp))
tmp <- intersect(rownames(sanger_mut),tmp)
tmp <- intersect(rownames(sanger_drug),tmp)
sanger_exp <- sanger_exp[,tmp]
sanger_cnv <- sanger_cnv[,tmp]
sanger_mut <- sanger_mut[tmp,]
sanger_drug <- sanger_drug[tmp,]
rm(tmp)

# identify the mek inhibitors
mek.inhib.sanger <-   sanger_drug_info$Name[grep(pattern="MEK",sanger_drug_info$Targets)][c(1,3)]

# identify the cells that have been consistently evaluated for with both MEK inhibitors
mek.sanger.cells <- rownames(sanger_drug)[-unique(c(which(is.na(sanger_drug[,mek.inhib.sanger[1]])), which(is.na(sanger_drug[,mek.inhib.sanger[2]]))))]

# identify the ActArea of the mek inhibs within the cells
mek.sanger.ic50 <- sanger_drug[mek.sanger.cells,mek.inhib.sanger]

# reduce the matrices to the mek.cells
sanger_exp <- sanger_exp[,mek.sanger.cells]
sanger_cnv <- sanger_cnv[,mek.sanger.cells]
sanger_mut <- sanger_mut[mek.sanger.cells,]
sanger_drug <- sanger_drug[mek.sanger.cells,]

# ####################################################################################################
# # define the NSCLC Breast Lung Melanoma Glioma & heMal (hematological malignacies) cells
# ####################################################################################################
# 
# carcinoma.mek.cells <-  intersect(mek.cells,sanger_info$sanger.name[sanger_info$Histology =="carcinoma"])
# nsclc.mek.cells <- carcinoma.mek.cells[grep(pattern="LUNG",x=carcinoma.mek.cells)]
# nsclc.mek.cells <- intersect(nsclc.mek.cells,sanger_info$sanger.name[ sanger_info$Hist.Subtype1 !="small_cell_carcinoma"])
# crc.mek.cells <- carcinoma.mek.cells[grep(pattern="LARGE_INTESTINE",x=carcinoma.mek.cells)]
# breast.mek.cells <- carcinoma.mek.cells[grep(pattern="BREAST",x=carcinoma.mek.cells)]
# melanoma.mek.cells <-  intersect(mek.cells,sanger_info$sanger.name[sanger_info$Histology =="malignant_melanoma"])
# glioma.mek.cells <-  intersect(mek.cells,sanger_info$sanger.name[sanger_info$Histology =="glioma"])
# hemal.mek.cells <- intersect(mek.cells,sanger_info$sanger.name[sanger_info$Histology %in% c("haematopoietic_neoplasm","lymphoid_neoplasm")])
# 
