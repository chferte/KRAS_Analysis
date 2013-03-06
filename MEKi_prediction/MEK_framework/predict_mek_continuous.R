# Charles Ferté
# Sage Bionetworks
# 14 Feb 2012

#############################
# load the data
#############################
source("/home/cferte/FELLOW/cferte/KRAS_Analysis/MEKi_prediction/MEK_framework/load_mek_ccle.R")
require(mixtools)

#######################################################
# predictive modeling
# we predict the ic50 
# training all cells global matrix without eigengenes and with eigengenes in parallel
#######################################################
# define globalmatrix (without eigengenes)
global.matrix <- rbind(ccle_exp,ccle_cnv,ccle_mut)
rownames(global.matrix) <- c(paste(rownames(ccle_exp),"_exp",sep=""),paste(rownames(ccle_cnv),"_cnv",sep=""),paste(rownames(ccle_mut),"_mut",sep=""))

# define globalmatrix2 (with eigengenes)
# eigengenes <- t(eigengenes[colnames(global.matrix),c(1:30)])
# global.matrix2 <- rbind(global.matrix,eigengenes)
# rownames(global.matrix2) <- c(rownames(global.matrix),paste("PC",c(1:30),sep=""))

par(mfrow=c(2,4),oma=c(0,0,6,0))
cells <- list(mek.cells,nsclc.mek.cells,breast.mek.cells,crc.mek.cells,hemal.mek.cells,glioma.mek.cells,melanoma.mek.cells)
cell.names <- list("ALL CELLS","NSCLC","BREAST","CRC","Hematologic\nMalignancies","GLIOMA","MELANOMA")
for(i in c(1:length(cells))){
  plot(density(apply(ccle_drug[cells[[i]],mek.inhib],1,mean)),xlim=c(-2,8.5),main=paste(cell.names[[i]]),ylim=c(0,.8),lwd=3)
}
title(main="Distribution of the sensitivity of cell lines to MEK inhibitors according to tissue type \n(sensitivity assessed by ActArea)",
      sub="sensitivity was assessed by ActArea",outer=TRUE)


blah <- apply(ccle_drug[crc.mek.cells,mek.inhib],1,mean)
abc <- normalmixEM(blah)
plot(abc,density=T)

# assign weights 
# mean.mek.sens <- apply(ccle_drug[mek.cells,mek.inhib],1,mean)
# density.id <- sapply(c(1:length(mean.mek.sens)),function(x){min(which(density(mean.mek.sens)$x>mean.mek.sens[x]))})
# density.weights <- log(1/(density(mean.mek.sens)$y[density.id]))
# names(density.weights) <- names(mean.mek.sens)
# rm(density.id,mean.mek.sens)
# density.weights  <- rep(1,times=length(density.weights))                    
# 
# summary(density.weights)
# # boxplot(density.weights)
# density.weights <- rep(0.5,times=length(mek.cells))
# names(density.weights) <- mek.cells
# density.weights[names(new.weight)] <- new.weight

require(multicore)

N=100
#models <- 0
i <- 0

# 
# # set up the Q1:Q4 for running balanced models
# q25 <- quantile(apply(ccle_drug[mek.cells,mek.inhib],1,mean),probs=.35) 
# q50 <- quantile(apply(ccle_drug[mek.cells,mek.inhib],1,mean),probs=.5)
# q75 <- quantile(apply(ccle_drug[mek.cells,mek.inhib],1,mean),probs=.65)
# Q1 <- names(which(apply(ccle_drug[mek.cells,mek.inhib],1,mean) < q25))
# Q2 <- names(which(apply(ccle_drug[mek.cells,mek.inhib],1,mean) > q25 & apply(ccle_drug[mek.cells,mek.inhib],1,mean) < q50 ))
# Q3 <- names(which(apply(ccle_drug[mek.cells,mek.inhib],1,mean) > q50 & apply(ccle_drug[mek.cells,mek.inhib],1,mean) < q75 ))
# Q4 <- names(which(apply(ccle_drug[mek.cells,mek.inhib],1,mean) > q75))
# rm(q25,q50,q75)


PARAL <- mclapply(X=1:N,FUN=function(x){
  print(i)
  i <- 1+1
#while(models<N)
#{
  
  #standard sampling
  train <- sample(mek.cells,replace=TRUE)
  
  # balanced model
  #train <- c(sample(x=c(Q1),replace=TRUE,size=220),sample(Q4,replace=TRUE,size=220))
  
  # weighted model on the distribution
  #train <- sample(mek.cells,replace=TRUE)
  
  
  # tissue specific models
  #train <- sample(hemal.mek.cells,replace=TRUE)
  
  # weighted lung models
  #train <- c(sample(mek.cells,replace=TRUE,size=144),sample(nsclc.mek.cells,replace=TRUE))
  
  # mixed balanced & tissue specific model
  #train  <- c(sample(nsclc.mek.cells,replace=TRUE),sample(x=c(Q1,Q2,Q3),replace=TRUE,size=36),sample(Q4,replace=TRUE,size=36))
  
  vec.train <-apply(ccle_drug[train,mek.inhib],1,mean)
  
  #cv.fit <- cv.glmnet(t(global.matrix[,train]), y=vec.train,nfolds=3, alpha=.1,weights=density.weights[train])
  #fit <- glmnet(x=t(global.matrix[,train]),y=vec.train,alpha=.1,lambda=cv.fit$lambda.1se,weights=density.weights[train])
  cv.fit <- cv.glmnet(t(global.matrix[,train]), y=vec.train,nfolds=3, alpha=.1)
  fit <- glmnet(x=t(global.matrix[,train]),y=vec.train,alpha=.1,lambda=cv.fit$lambda.1se)
  return(list(fit,train)) },mc.set.seed=TRUE,mc.cores=6)
  


yhat.all <- c()
yhat.breast <- c()
yhat.nsclc <- c()
yhat.crc <- c()
yhat.glioma <- c()
yhat.melanoma <- c()
yhat.hemal  <- c()
selected <- c()

for(i in c(1:N)){
  train <- PARAL[[i]][[2]]
  fit <- PARAL[[i]][[1]]
  val <-mek.cells[-which(mek.cells %in% train)]
  selected <- c(selected,list(fit$beta))
  yhat.all <- c(yhat.all,list(predict(fit, t(global.matrix[,val]))))
  yhat.breast <- c(yhat.breast,list(predict(fit,t(global.matrix[,breast.mek.cells[-which(breast.mek.cells %in% train)]]))))
  yhat.nsclc <- c(yhat.nsclc,list(predict(fit,t(global.matrix[,nsclc.mek.cells[-which(nsclc.mek.cells %in% train)]]))))
  yhat.crc <- c(yhat.crc,list(predict(fit,t(global.matrix[,crc.mek.cells[-which(crc.mek.cells %in% train)]]))))
  yhat.glioma <- c(yhat.glioma,list(predict(fit,t(global.matrix[,glioma.mek.cells[-which(glioma.mek.cells %in% train)]]))))
  yhat.melanoma <- c(yhat.melanoma,list(predict(fit,t(global.matrix[,melanoma.mek.cells[-which(melanoma.mek.cells %in% train)]]))))
  yhat.hemal <- c(yhat.hemal,list(predict(fit,t(global.matrix[,hemal.mek.cells[-which(hemal.mek.cells %in% train)]]))))
  print(i)}  


# 
# 
# ###### SPEARMAN ###########
# par(mfrow=c(2,4),oma=c(0,0,6,0))
# method.cor <- "spearman"
# cex=1.5
# source("/home/cferte/FELLOW/cferte/KRAS_Analysis/MEKi_prediction/Mek_performance.R")
# title(main="Performance of 100 bootstrapped models\npredicting sensitivity to MEK inhibitors\ntraining in all CCLE cell lines",
#       sub="sensitivity was assessed by ActArea",outer=TRUE)
# 
# ###### PEARSON ###########
# par(mfrow=c(2,4),oma=c(0,0,6,0))
# method.cor <- "pearson"
# cex=1.5
# source("/home/cferte/FELLOW/cferte/KRAS_Analysis/MEKi_prediction/Mek_performance.R")
# title(main="Performance of 100 bootstrapped models\npredicting sensitivity to MEK inhibitors\ntraining in CCLE all cell lines",
#       sub="sensitivity was assessed by ActArea",outer=TRUE)
# 
