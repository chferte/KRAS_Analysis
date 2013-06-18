# Charles Ferté
# Sage Bionetworks
# June 2, 2013

# load the lung virtual drug sensitivity values
load("/home/cferte/RESULTS/vds_lung_mut.Rda")
y_hat <- vds

#load the data
load("/home/jguinney/data/TCGA_ds.rda")

library(synapseClient)
synapseLogin(username="charles.ferte@sagebase.org",password="charles")
luad_all <- loadEntity("syn1676707")
luad_all <- luad_all$objects$luad_data
luad_mut <- luad_all[[2]]
head(luad_mut)
luad_mut <- luad_mut==1
colnames(luad_mut) <- gsub(colnames(luad_mut),pattern="-",replacement=".",fixed=TRUE)
colnames(luad_mut) <- substr(colnames(luad_mut),1,12)
luad[[4]]  <- luad_mut
rm(luad_all,luad_mut)


# replace the luad mut by the luad mut that i have (much more samples)
source("/home/cferte/FELLOW/cferte/KRAS_Analysis/MEKi_prediction/MEK_framework/cellline_2_tcga_pipeline.R")
fMat <- build_feature_matrix(luad, gene.dict="cosmic",with.rppa=FALSE)

drug_fmat <- find_drug_features(y_hat,fMat, rep(TRUE, nrow(fMat)), 
                                beta_threshold=10^-3,num.bootstraps=250)
drug <- "MEK inhibitors"
top=30
par(mfrow=c(1,1))
plot_features_2(drug_fmat, 
                paste(drug," in luad cancer n=",ncol(fMat),sep=""),
                top=top,text.cex=.7)

# plot the distribution of the top features selected in the model as a heatmap
M <- fMat[rownames(drug_fmat[[1]])[1:top],]
M <- ifelse(M==TRUE,1,0)
library(gplots)
par(oma=c(0,0,0,4))
heatmap.2(M,scale="none",trace="none",col=c("gray80","orangered"),key=FALSE)

# correlation matrix for the top sensitive ones
top.sens <- rownames(drug_fmat$df)[1:top][drug_fmat$df$posFreq[1:top]>.7]

foo <- fMat[top.sens,]
library(e1071)
H <- hamming.distance(foo)

#COR <- cor(foo,method="spearman")
heatmap.2(H,trace="none")

# correlation matrix for the top resistant ones
top.res <- rownames(drug_fmat$df)[1:top][drug_fmat$df$posFreq[1:top]<.3]
foo <- t(fMat[top.res,])
COR <- cor(foo,method="spearman")
heatmap.2(COR,trace="none")

# draw a more elaborated heatmap on the top features
par(oma=c(0,0,0,5))
top.sens <- which(drug_fmat$df$posFreq[1:top]>.7)
top.res <- which(drug_fmat$df$posFreq[1:top]<.3)
order <- c(top.sens,top.res)
M <- fMat[rownames(drug_fmat[[1]])[order],]
M <- ifelse(M==TRUE,1,0)
vec <- sort(y_hat[colnames(M)])
M <- M[,names(vec)]
heatmap.2(M,Rowv=FALSE,Colv="Rowv",dendrogram="none",
          scale="none",trace="none",
          keysize = 1.25, density.info = 'none',key=FALSE,
          ColSideColors=greenred(ncol(M))[c(ncol(M):1)],
          col = colorpanel(2, 'gray70', 'orangered'),
          main=paste("TCGA crc cancer \nn=", ncol(M)))





# decision tree
library(rpart)
par(oma=c(0,0,0,0),mar=c(1,1,1,1))
top_genes <- as.character(drug_fmat$df$genes[1:50])
idxs <- groupMatch(names(y_hat), colnames(fMat))
y_hat.m <- y_hat[idxs[[1]]]
sub_fmat <- fMat[top_genes, idxs[[2]]]

M <- t(sub_fmat)

fit <- rpart(y_hat.m ~ M,method="anova",
             control=rpart.control(cp=.01,
          xval=50,minsplit=10))
plot(fit)
text(fit,use.n=TRUE,cex=.7)

