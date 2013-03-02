# Charles Ferté
# Sage Bionetworks
# 14 Feb 2012


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
synapseLogin("charles.ferte@sagebase.org","charles")

###############################################################
# load the CCLE data (snm normalized) 
###############################################################
ccle_all <- loadEntity("syn1671195")
ccle_all <- ccle_all$objects$ccle_data
assign(x=names(ccle_all)[1],ccle_all[[1]])
assign(x=names(ccle_all)[2],ccle_all[[2]])
assign(x=names(ccle_all)[3],ccle_all[[3]])
assign(x=names(ccle_all)[4],ccle_all[[4]])
assign(x=names(ccle_all)[5],ccle_all[[5]])
assign(x=names(ccle_all)[6],ccle_all[[6]])
assign(x=names(ccle_drug)[1],ccle_drug[[1]])
assign(x=names(ccle_drug)[2],ccle_drug[[2]])

# let's use the ActArea as ccle_drug: "_IC50Norm" or "_ActAreaNorm"
ccle_drug <- ccle_drug_ActAreaNorm

#modify the ccle_mut into a binary matrix
ccle_mut[ccle_mut!="0"] <- 1
tmp <- rownames(ccle_mut)
ccle_mut <- apply(ccle_mut,2,as.numeric)
rownames(ccle_mut) <- tmp
rm(tmp)


#########################################################################################################################################
# make the data coherent between exp and cnv to ultimately compute the eigengenes vector
#########################################################################################################################################

tmp <- intersect(colnames(ccle_exp),ccle_info$CCLE.name)
tmp <- intersect(tmp,colnames(ccle_cnv))
ccle_exp <- ccle_exp[,tmp]
ccle_cnv <- ccle_cnv[,tmp]
ccle_info <- ccle_info[which(ccle_info$CCLE.name %in% intersect(tmp,ccle_info$CCLE.name)),]
rm(tmp)

global.matrix <- rbind(ccle_exp,ccle_cnv)
tissue.origin <- as.factor(ccle_info$Site.Primary)
# s <- fast.svd(global.matrix-rowMeans(global.matrix))
# par(mfrow=c(1,1))
# 
# # percentage variance explained
# plot(s$d^2/sum(s$d^2),pch=20)
# 
# # plot first and second svd
# plot(s$v[,1],s$v[,2],col=rainbow(24)[tissue.origin],pch=20,cex=1.5, xlab="PC1",ylab="PC2")
# plot.new()
# legend(0,1,legend=levels(tissue.origin),cex=.8,col=rainbow(24),pch=20,text.width=.4,text.font=2,
#        text.col=rainbow(24))
# 
# # assign values for the eigengenes since they discriminate well the tissue specificity
# eigengenes <- s$v
# rownames(eigengenes) <- colnames(ccle_exp)
# colnames(eigengenes) <- paste("PC_",seq(from=1,to=ncol(eigengenes),by=1),sep="")

#########################################################################################################
## feature selection for low variance
#########################################################################################################
tmp <- apply(ccle_exp,1,sd)
ccle_exp <- ccle_exp[which(tmp>quantile(x=tmp,probs=.25)),]
rm(tmp)

tmp <- apply(ccle_cnv,1,sd)
ccle_cnv <- ccle_cnv[which(tmp>quantile(x=tmp,probs=.25)),]
rm(tmp)

#########################################################################################################
## Make the data coherent between all the datasets
#########################################################################################################
tmp <- intersect(colnames(ccle_cnv),colnames(ccle_exp))
tmp <- intersect(colnames(ccle_mut),tmp)
tmp <- intersect(rownames(ccle_drug),tmp)
ccle_exp <- ccle_exp[,tmp]
ccle_cnv <- ccle_cnv[,tmp]
ccle_mut <- ccle_mut[,tmp]
ccle_drug <- ccle_drug[tmp,]
ccle_info <- ccle_info[which(ccle_info$CCLE.name %in% intersect(tmp,ccle_info$CCLE.name)),]
rm(tmp)

# identify the mek inhibitors
mek.inhib <-   ccle_drugs_info$Compound..code.or.generic.name.[grep(pattern="MEK",ccle_drugs_info$Target.s.)]

# identify the cells that have been consistently evaluated for with both MEK inhibitors
mek.cells <- rownames(ccle_drug)[-unique(c(which(is.na(ccle_drug[,mek.inhib[1]])), which(is.na(ccle_drug[,mek.inhib[2]]))))]

# identify the ActArea of the mek inhibs within the cells
mek.ActArea <- ccle_drug[mek.cells,mek.inhib]

# reduce the matrices to the mek.cells
ccle_exp <- ccle_exp[,mek.cells]
ccle_cnv <- ccle_cnv[,mek.cells]
ccle_mut <- ccle_mut[,mek.cells]
ccle_drug <- ccle_drug[mek.cells,]

####################################################################################################
# define the NSCLC Breast Lung Melanoma Glioma & heMal (hematological malignacies) cells
####################################################################################################

carcinoma.mek.cells <-  intersect(mek.cells,ccle_info$CCLE.name[ccle_info$Histology =="carcinoma"])
nsclc.mek.cells <- carcinoma.mek.cells[grep(pattern="LUNG",x=carcinoma.mek.cells)]
nsclc.mek.cells <- intersect(nsclc.mek.cells,ccle_info$CCLE.name[ ccle_info$Hist.Subtype1 !="small_cell_carcinoma"])
crc.mek.cells <- carcinoma.mek.cells[grep(pattern="LARGE_INTESTINE",x=carcinoma.mek.cells)]
breast.mek.cells <- carcinoma.mek.cells[grep(pattern="BREAST",x=carcinoma.mek.cells)]
melanoma.mek.cells <-  intersect(mek.cells,ccle_info$CCLE.name[ccle_info$Histology =="malignant_melanoma"])
glioma.mek.cells <-  intersect(mek.cells,ccle_info$CCLE.name[ccle_info$Histology =="glioma"])
hemal.mek.cells <- intersect(mek.cells,ccle_info$CCLE.name[ccle_info$Histology %in% c("haematopoietic_neoplasm","lymphoid_neoplasm")])

