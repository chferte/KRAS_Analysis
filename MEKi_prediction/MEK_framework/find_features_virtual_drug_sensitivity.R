# Charles Ferté
# Sage Bionetworks
# June 2, 2013

# load the lung virtual drug sensitivity values
load("/home/cferte/RESULTS/vds_lung_mut.Rda")
y_hat <- vds
summary(vds)

#load the data
load("/home/jguinney/data/TCGA_ds_ver2.rda")

# # replace the luad mut by the luad mut that i have (much more samples)
# library(synapseClient)
# synapseLogin(username="charles.ferte@sagebase.org",password="charles")
# luad_all <- loadEntity("syn1676707")
# luad_all <- luad_all$objects$luad_data
# luad_mut <- luad_all[[2]]
# luad_mut <- luad_mut==1
# colnames(luad_mut) <- gsub(colnames(luad_mut),pattern="-",replacement=".",fixed=TRUE)
# colnames(luad_mut) <- substr(colnames(luad_mut),1,12)
# luad[[4]]  <- luad_mut
# rm(luad_all,luad_mut)


# load the functions to find the features
source("/home/cferte/FELLOW/cferte/KRAS_Analysis/MEKi_prediction/MEK_framework/cellline_2_tcga_pipeline.R")
#fMat <- build_feature_matrix(crc, gene.dict="cosmic",with.rppa=FALSE)
fMat <- build_feature_matrix(luad, with.rppa=FALSE)

drug_fmat <- find_drug_features(y_hat,fMat, rep(TRUE, nrow(fMat)), 
                                beta_threshold=10^-3,num.bootstraps=200)
drug <- "MEK inhibitors"
top=30
par(mfrow=c(1,1))
plot_features_2(drug_fmat, 
                paste(drug," in luad cancer n=",ncol(fMat),sep=""),
                top=top,text.cex=.6)

# plot the distribution of the top features selected in the model as a heatmap
top=15
M <- fMat[rownames(drug_fmat[[1]])[1:top],]
M <- ifelse(M==TRUE,1,0)
library(gplots)
par(oma=c(0,0,0,4))
tmp <- intersect(names(y_hat),colnames(M))
y_hat <- y_hat[tmp]
M <- M[,tmp]
vec <- (y_hat-min(y_hat))/max(y_hat-min(y_hat))
vec <- round(200*vec+1,digits=0)
#vec <- order(y_hat)
col.vec <- redgreen(ncol(M))[vec ]
heatmap.2(M,scale="none",trace="none",col=c("gray80","orangered"),key=FALSE, ColSideColors=col.vec)

# draw a more elaborated heatmap on the top features
top=15
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
          main=paste("TCGA luad cancer \nn=", ncol(M)))

# decision tree
library(rpart)
par(oma=c(0,0,0,0),mar=c(1,1,1,1))
top_genes <- as.character(drug_fmat$df$genes[1:top])
grep(pattern="mut_",x=top_genes,value=FALSE)
idxs <- groupMatch(names(y_hat), colnames(fMat))
y_hat.m <- y_hat[idxs[[1]]]
sub_fmat <- fMat[top_genes, idxs[[2]]]

M <- t(sub_fmat)

fit <- rpart(y_hat.m ~ M,method="anova",
             control=rpart.control(cp=.01,
          xval=50,minsplit=10))
plot(fit)
text(fit,use.n=TRUE,cex=.7)


# apply the model in cell lines 
source("/home/cferte/FELLOW/cferte/KRAS_Analysis/MEKi_prediction/MEK_framework/load_mek_ccle.R")
foo <- sapply(strsplit(colnames(M),split="_"),function(x)  x[[length(x)]]  )
foo <- unique(intersect(rownames(ccle_mut),foo))


# heatmap with clustering
ccle.matrix <- ccle_mut[foo,nsclc.mek.cells]
vec <- apply(mek.ActArea[nsclc.mek.cells,],1,sum)
vec <- order(vec)
col.vec <- redgreen(ncol(ccle.matrix))[vec]
heatmap.2(ccle.matrix,scale="none",trace="none",col=c("gray80","orangered"),key=FALSE, ColSideColors=col.vec)

# heatmap without clustering

vec <- apply(mek.ActArea[nsclc.mek.cells,],1,sum)
vec <- sort(vec)
ccle.matrix <- ccle_mut[foo,names(vec)]
heatmap.2(ccle.matrix,Rowv=FALSE,Colv="Rowv",dendrogram="none",
          scale="none",trace="none",
          keysize = 1.25, density.info = 'none',key=FALSE,
          ColSideColors=greenred(ncol(ccle.matrix))[c(ncol(ccle.matrix):1)],
          col = colorpanel(2, 'gray70', 'orangered'),
          main=paste("TCGA lung cancer \nn=", ncol(ccle.matrix)))

