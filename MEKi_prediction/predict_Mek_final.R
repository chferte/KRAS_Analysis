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

# let's use the ActArea as ccle_drug
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
s <- fast.svd(global.matrix-rowMeans(global.matrix))
par(mfrow=c(1,1))

# percentage variance explained
plot(s$d^2/sum(s$d^2),pch=20)

# plot first and second svd
plot(s$v[,1],s$v[,2],col=rainbow(24)[tissue.origin],pch=20,cex=1.5, xlab="PC1",ylab="PC2")
plot.new()
legend(0,1,legend=levels(tissue.origin),cex=.8,col=rainbow(24),pch=20,text.width=.4,text.font=2,
       text.col=rainbow(24))

# assign values for the eigengenes since they discriminate well the tissue specificity
eigengenes <- s$v
rownames(eigengenes) <- colnames(ccle_exp)
colnames(eigengenes) <- paste("PC_",seq(from=1,to=ncol(eigengenes),by=1),sep="")

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

#######################################################
# predictive modeling
# we predict the ic50 
# training all cells global matrix without eigengenes and with eigengenes in parallel
#######################################################

# define globalmatrix (without eigengenes)
global.matrix <- rbind(ccle_exp,ccle_cnv,ccle_mut)
rownames(global.matrix) <- c(paste(rownames(ccle_exp),"_exp",sep=""),paste(rownames(ccle_cnv),"_cnv",sep=""),paste(rownames(ccle_mut),"_mut",sep=""))

# define globalmatrix2 (with eigengenes)
eigengenes <- t(eigengenes[colnames(global.matrix),c(1:30)])
global.matrix2 <- rbind(global.matrix,eigengenes)
rownames(global.matrix2) <- c(rownames(global.matrix),paste("PC",c(1:30),sep=""))


N=200
models <- 0
i <- 0

yhat.all <- c()
yhat.breast <- c()
yhat.nsclc <- c()
yhat.crc <- c()
yhat.glioma <- c()
yhat.melanoma <- c()
yhat.hemal  <- c()
selected <- c()

yhat.all2 <- c()
yhat.breast2 <- c()
yhat.nsclc2 <- c()
yhat.crc2 <- c()
yhat.glioma2 <- c()
yhat.melanoma2 <- c()
yhat.hemal2  <- c()
selected2 <- c()

while(models<N)
{
  par(mfrow=c(1,1))
  train <- sample(mek.cells,replace=TRUE)
  val <-mek.cells[-which(mek.cells %in% train)]
  vec.train <-apply(ccle_drug[train,mek.inhib],1,mean)
  
  cv.fit <- cv.glmnet(t(global.matrix[,train]), y=vec.train,nfolds=3, alpha=.1)
  fit <- glmnet(x=t(global.matrix[,train]),y=vec.train,alpha=.1,lambda=cv.fit$lambda.1se)
  selected <- c(selected,list(fit$beta))
  yhat.all <- c(yhat.all,list(predict(fit, t(global.matrix[,val]))))
  yhat.breast <- c(yhat.breast,list(predict(fit,t(global.matrix[,breast.mek.cells[-which(breast.mek.cells %in% train)]]))))
  yhat.nsclc <- c(yhat.nsclc,list(predict(fit,t(global.matrix[,nsclc.mek.cells[-which(nsclc.mek.cells %in% train)]]))))
  yhat.crc <- c(yhat.crc,list(predict(fit,t(global.matrix[,crc.mek.cells[-which(crc.mek.cells %in% train)]]))))
  yhat.glioma <- c(yhat.glioma,list(predict(fit,t(global.matrix[,glioma.mek.cells[-which(glioma.mek.cells %in% train)]]))))
  yhat.melanoma <- c(yhat.melanoma,list(predict(fit,t(global.matrix[,melanoma.mek.cells[-which(melanoma.mek.cells %in% train)]]))))
  yhat.hemal <- c(yhat.hemal,list(predict(fit,t(global.matrix[,hemal.mek.cells[-which(hemal.mek.cells %in% train)]]))))
 
  cv.fit2 <- cv.glmnet(t(global.matrix2[,train]), y=vec.train,nfolds=3, alpha=.1)
  fit2 <- glmnet(x=t(global.matrix2[,train]),y=vec.train,alpha=.1,lambda=cv.fit$lambda.1se)
  selected2 <- c(selected2,list(fit2$beta))
  yhat.all2 <- c(yhat.all2,list(predict(fit2, t(global.matrix2[,val]))))
  yhat.breast2 <- c(yhat.breast2,list(predict(fit2,t(global.matrix2[,breast.mek.cells[-which(breast.mek.cells %in% train)]]))))
  yhat.nsclc2 <- c(yhat.nsclc2,list(predict(fit2,t(global.matrix2[,nsclc.mek.cells[-which(nsclc.mek.cells %in% train)]]))))
  yhat.crc2 <- c(yhat.crc2,list(predict(fit2,t(global.matrix2[,crc.mek.cells[-which(crc.mek.cells %in% train)]]))))
  yhat.glioma2 <- c(yhat.glioma2,list(predict(fit2,t(global.matrix2[,glioma.mek.cells[-which(glioma.mek.cells %in% train)]]))))
  yhat.melanoma2 <- c(yhat.melanoma2,list(predict(fit2,t(global.matrix2[,melanoma.mek.cells[-which(melanoma.mek.cells %in% train)]]))))
  yhat.hemal2 <- c(yhat.hemal2,list(predict(fit2,t(global.matrix2[,hemal.mek.cells[-which(hemal.mek.cells %in% train)]]))))
  
  
  i=1+i
  print(i)
  models <- length(yhat.all2)
  }  

#################################################################################################################
# explore what is retained in the models
#################################################################################################################

# abc <- matrix(NA,ncol=length(selected),nrow=nrow(global.matrix))
# rownames(abc) <- rownames(global.matrix)
# for(i in c(1:length(selected)))
# {  abc[,i]<- as.numeric(selected[[i]])  
# }

# fit.betas.mek.models <- selected
# 
# # store it into synapse
# fits <- Data(list(name = "fitbetas_mek_models", parentId = 'syn1670945'))
# fits<- createEntity(fits)
# 
# # add object into the data entity
# fits <- addObject(fits,fit.betas.mek.models)
# 
# # push the raw data into this entity
# fits <- storeEntity(entity=fits)
# 
# par(mfrow=c(1,1))
# abc[abc!=0] <-1 
# hist(log10(rowSums(abs(abc))),breaks=50,col="red")
# h <- rowSums(abs(abc))
# sort(h,decreasing=TRUE)[1:100]
# which(h>quantile(h,probs=.99))

###### SPEARMAN ###########
par(mfrow=c(2,4),oma=c(0,0,6,0))
method.cor <- "spearman"
cex=1.5
source("/home/cferte/FELLOW/cferte/KRAS_Analysis/MEKi_prediction/Mek_performance.R")
title(main="Performance of 200 bootstrapped models\npredicting sensitivity to MEK inhibitors\ntraining in all CCLE cell lines",
      sub="sensitivity was assessed by ActArea",outer=TRUE)

###### PEARSON ###########
par(mfrow=c(2,4),oma=c(0,0,6,0))
method.cor <- "pearson"
cex=1.5
source("/home/cferte/FELLOW/cferte/KRAS_Analysis/MEKi_prediction/Mek_performance.R")
title(main="Performance of 200 bootstrapped models\npredicting sensitivity to MEK inhibitors\ntraining in CCLE all cell lines",
      sub="sensitivity was assessed by ActArea",outer=TRUE)

#############################
# plot the distribution of MEK ActArea
############################

par(mfrow=c(2,4),oma=c(0,0,6,0))
plot(density(apply(ccle_drug[mek.cells,mek.inhib],1,mean)),xlim=c(-2,8.5),main="ALL CELLS",ylim=c(0,.8))
plot(density(apply(ccle_drug[nsclc.mek.cells,mek.inhib],1,mean)),xlim=c(-2,8.5),main="NSCLC",ylim=c(0,.8))
plot(density(apply(ccle_drug[breast.mek.cells,mek.inhib],1,mean)),xlim=c(-2,8.5),main="BREAST",ylim=c(0,.8))
plot(density(apply(ccle_drug[crc.mek.cells,mek.inhib],1,mean)),xlim=c(-2,8.5),main="COLORECTAL",ylim=c(0,.8))
plot(density(apply(ccle_drug[hemal.mek.cells,mek.inhib],1,mean)),xlim=c(-2,8.5),main="Hematologic\nMalignancies",ylim=c(0,.8))
plot(density(apply(ccle_drug[glioma.mek.cells,mek.inhib],1,mean)),xlim=c(-2,8.5),main="GLIOMA",ylim=c(0,.8))
plot(density(apply(ccle_drug[melanoma.mek.cells,mek.inhib],1,mean)),xlim=c(-2,8.5),main="MELANOMA",ylim=c(0,.8))
title(main="Distribution of the sensitivity of cell lines to MEK inhibitors according to tissue type \n(sensitivity assessed by ActArea)",
      sub="sensitivity was assessed by ActArea",outer=TRUE)
#############################
# 
############################
par(mfrow=c(1,1),oma=c(1,1,1,1))
boxplot(apply(ccle_drug[mek.cells,mek.inhib],1,mean),main="ALL MEK CELLS")
par(mfrow=c(2,3))
boxplot(apply(ccle_drug[nsclc.mek.cells,mek.inhib],1,mean),main="NSCLC",ylim=c(0,8))
boxplot(apply(ccle_drug[breast.mek.cells,mek.inhib],1,mean),main="BREAST",ylim=c(0,8))
boxplot(apply(ccle_drug[crc.mek.cells,mek.inhib],1,mean),main="COLORECTAL",ylim=c(0,8))
boxplot(apply(ccle_drug[hemal.mek.cells,mek.inhib],1,mean),main="Hematologic\nMalignancies",ylim=c(0,8))
boxplot(apply(ccle_drug[glioma.mek.cells,mek.inhib],1,mean),main="GLIOMA",ylim=c(0,8))
boxplot(apply(ccle_drug[melanoma.mek.cells,mek.inhib],1,mean),main="MELANOMA",ylim=c(0,8))


############################
# plot the confidence intervall for each IC50 prediction
############################

# first plot the density of the IC50
k <- list(mek.cells,nsclc.mek.cells,breast.mek.cells,crc.mek.cells,hemal.mek.cells,glioma.mek.cells,melanoma.mek.cells)
yhats <- list(yhat.all,yhat.nsclc,yhat.breast,yhat.crc,yhat.hemal,yhat.glioma,yhat.melanoma)
names(k) <- c("ALL CELLS","NSCLC","BREAST","COLORECTAL","Hematologic\nMalignancies","GLIOMA","MELANOMA")
par(mfrow=c(2,4))
for(i in c(1:length(k)))
  {

tissue.ic50 <- apply(ccle_drug[k[[i]],mek.inhib],1,mean)
plot(density(tissue.ic50),xlim=c(-2,6),main=paste(names(k)[i]),ylim=c(0,3.5))

# then create a matrix abc containing all the yhat per cells
abc <- matrix(NA,ncol=N,nrow=length(k[[i]]))
rownames(abc) <- k[[i]]
for(j in c(1:length(yhats[[i]]))) { abc[rownames(yhats[[i]][[j]]),j]<- yhats[[i]][[j]]}

# and compute the rmse for each cell
RMSE <- c()
for(n in c(1:length(tissue.ic50))){

  tmp <- sapply(c(1:ncol(abc)),function(x){(tissue.ic50[n]-abc[n,x])^2})
  RMSE <- c(RMSE,sqrt(mean(tmp,na.rm=TRUE)))
  }
names(RMSE) <- k[[i]]

# plot the variance of the predictions accodring to their actual IC50
# to show that the prediction performance is dependant on the distribution of the model features in our model
tmp <- names(sort(apply(ccle_drug[k[[i]],mek.inhib],1,mean)))
lines(tissue.ic50[tmp],RMSE[tmp],col="red",type="p",pch=20)
q <- loess(RMSE[tmp]~tissue.ic50[tmp])
lines(q$x,q$fitted,col="blue",lwd=3)
}

title(main="Distribution of the sensitivity of cell lines to MEK inhibitors according to tissue type \n(sensitivity assessed by ActArea)",,outer=TRUE)
