# charles ferté
# sage bionetworks
# June 2, 2013

# load the lung virtual rug sensitivity values
load("/home/cferte/RESULTS/vds_lung_mut.Rda")
y_hat <- vds

#load the data
load("/home/jguinney/data/TCGA_ds.rda")

library(synapseClient)
synapseLogin(username="charles.ferte@sagebase.org",password="charles")
luad_all <- loadEntity("syn1676707")
luad_all <- luad_all$objects$luad_data
luad_mut <- luad_all[[2]]
luad_mut <- luad_mut==1
colnames(luad_mut) <- gsub(colnames(luad_mut),pattern="-",replacement=".",fixed=TRUE)
colnames(luad_mut) <- substr(colnames(luad_mut),1,12)
luad[[4]]  <- luad_mut
rm(luad_all,luad_mut)


# replace the luad mut by the luad mut that i have (much more samples)
source("/home/cferte/FELLOW/cferte/KRAS_Analysis/MEKi_prediction/MEK_framework/cellline_2_tcga_pipeline.R")
fMat <- build_feature_matrix(luad, gene.dict="cosmic",with.rppa=FALSE)
drug_fmat <- find_drug_features(y_hat,fMat, rep(TRUE, nrow(fMat)), 
                                beta_threshold=10^-3,num.bootstraps=200)
drug <- "MEK inhibitors"
plot_features_2(drug_fmat, 
                paste(drug," in lung adenocarcinoma (n=327)",sep=""),
                top=30,text.cex=.75)

# here are the top features selected in the model
head(drug_fmat[[1]])

library(rpart)
top_genes <- as.character(drug_fmat$df$genes[1:50])
idxs <- groupMatch(names(y_hat), colnames(fMat))
y_hat.m <- y_hat[idxs[[1]]]
sub_fmat <- fMat[top_genes, idxs[[2]]]

M <- t(sub_fmat)

fit <- rpart(y_hat.m ~ M,method="anova",
             control=rpart.control(cp=.01,xval=50,minsplit=10))
plot(fit)
text(fit,use.n=TRUE,cex=.7)
