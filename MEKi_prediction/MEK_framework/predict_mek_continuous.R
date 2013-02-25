# Charles Ferté
# Sage Bionetworks
# 14 Feb 2012

#############################
# load the data
#############################
source("/home/cferte/FELLOW/cferte/KRAS_Analysis/MEKi_prediction/MEK_framework/load_mek_data.R")

#######################################################
# predictive modeling
# we predict the ic50 
# training all cells global matrix without eigengenes and with eigengenes in parallel
#######################################################

# define globalmatrix (without eigengenes)
global.matrix <- rbind(ccle_exp,ccle_cnv,ccle_mut)
rownames(global.matrix) <- c(paste(rownames(ccle_exp),"_exp",sep=""),paste(rownames(ccle_cnv),"_cnv",sep=""),paste(rownames(ccle_mut),"_mut",sep=""))
# 
# # define globalmatrix2 (with eigengenes)
# eigengenes <- t(eigengenes[colnames(global.matrix),c(1:30)])
# global.matrix2 <- rbind(global.matrix,eigengenes)
# rownames(global.matrix2) <- c(rownames(global.matrix),paste("PC",c(1:30),sep=""))
# 

N=100
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

# set up the Q1:Q4 for running balanced models
q25 <- quantile(apply(ccle_drug[mek.cells,mek.inhib],1,mean),probs=.25) 
q50 <- quantile(apply(ccle_drug[mek.cells,mek.inhib],1,mean),probs=.5)
q75 <- quantile(apply(ccle_drug[mek.cells,mek.inhib],1,mean),probs=.75)
Q1 <- names(which(apply(ccle_drug[mek.cells,mek.inhib],1,mean) < q25))
Q2 <- names(which(apply(ccle_drug[mek.cells,mek.inhib],1,mean) > q25 & apply(ccle_drug[mek.cells,mek.inhib],1,mean) < q50 ))
Q3 <- names(which(apply(ccle_drug[mek.cells,mek.inhib],1,mean) > q50 & apply(ccle_drug[mek.cells,mek.inhib],1,mean) < q75 ))
Q4 <- names(which(apply(ccle_drug[mek.cells,mek.inhib],1,mean) > q75))
rm(q25,q50,q75)

while(models<N)
{
  par(mfrow=c(1,1))
  #train <- sample(mek.cells,replace=TRUE)
  train <- c(sample(Q1,replace=TRUE),sample(Q2,replace=TRUE),sample(Q3,replace=TRUE),sample(Q4,replace=TRUE))
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
   
  i=1+i
  print(i)
  models <- length(yhat.all)
  }  

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
