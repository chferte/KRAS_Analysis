# charles ferte
# feb 21 th 2012
# sage bionetworks

#############################
# load the data
#############################

source("/home/cferte/FELLOW/cferte/KRAS_Analysis/MEKi_prediction/MEK_framework/load_mek_data.R")

#############################
# first,plot the distribution of MEK ActArea
############################


par(mfrow=c(1,2))
plot(density(ccle_drug[mek.cells,mek.inhib[1]]),xlim=c(-2,8.5),main=paste("Sensitivity to",mek.inhib[1], "(ActArea)","\nin all cells"),ylim=c(0,.8))
abline(v=quantile(ccle_drug[mek.cells,mek.inhib[1]],probs=.8),col="red")
abline(v=quantile(ccle_drug[mek.cells,mek.inhib[1]],probs=.2),col="red")
plot(density(ccle_drug[mek.cells,mek.inhib[2]]),xlim=c(-2,8.5),main=paste("Sensitivity to",mek.inhib[2], "(ActArea)","\nin all cells"),ylim=c(0,.8))
abline(v=quantile(ccle_drug[mek.cells,mek.inhib[2]],probs=.8),col="green")
abline(v=quantile(ccle_drug[mek.cells,mek.inhib[2]],probs=.2),col="green")


#############################
# discretise the ActArea
#############################

a <- ifelse(ccle_drug[mek.cells,mek.inhib[1]]>quantile(ccle_drug[mek.cells,mek.inhib[1]],probs=.7),1,0)
b <- ifelse(ccle_drug[mek.cells,mek.inhib[2]]>quantile(ccle_drug[mek.cells,mek.inhib[2]],probs=.7),1,0)
discrete.mek.response <- ifelse(a==1 & b==1,1,0)

table(discrete.mek.response)

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


N=50
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
  train <- c(sample(mek.cells[which(discrete.mek.response==1)],replace=TRUE), sample(mek.cells[which(discrete.mek.response==0)],size=112,replace=TRUE))
  trainval <-mek.cells[-which(mek.cells %in% train)]
  vec.train <-discrete.mek.response[train]
  
  cv.fit <- cv.glmnet(t(global.matrix[,train]), y=vec.train,nfolds=3, alpha=.1,family="binomial")
  #plot(cv.fit)
  fit <- glmnet(x=t(global.matrix[,train]),y=vec.train,alpha=.1,lambda=cv.fit$lambda.1se,family="binomial")
  selected <- c(selected,list(fit$beta))
  yhat.all <- c(yhat.all,list(predict(fit, t(global.matrix[,val]),type="response")))
  yhat.breast <- c(yhat.breast,list(predict(fit,t(global.matrix[,breast.mek.cells[-which(breast.mek.cells %in% train)]]),type="response")))
  yhat.nsclc <- c(yhat.nsclc,list(predict(fit,t(global.matrix[,nsclc.mek.cells[-which(nsclc.mek.cells %in% train)]]),type="response")))
  yhat.crc <- c(yhat.crc,list(predict(fit,t(global.matrix[,crc.mek.cells[-which(crc.mek.cells %in% train)]]),type="response")))
  yhat.glioma <- c(yhat.glioma,list(predict(fit,t(global.matrix[,glioma.mek.cells[-which(glioma.mek.cells %in% train)]]),type="response")))
  yhat.melanoma <- c(yhat.melanoma,list(predict(fit,t(global.matrix[,melanoma.mek.cells[-which(melanoma.mek.cells %in% train)]]),type="response")))
  yhat.hemal <- c(yhat.hemal,list(predict(fit,t(global.matrix[,hemal.mek.cells[-which(hemal.mek.cells %in% train)]]),type="response")))

  
  i=1+i
  print(i)
  models <- length(yhat.all)
}  

# assess and plot the performance
require(ROCR)
# AUC.all <- c()
# for (i in c(1:N)){
#   Pred <- prediction(as.numeric(yhat.all[[i]]),discrete.mek.response[rownames(yhat.all[[i]])])
#   Perf <- performance(prediction.obj=Pred,"tpr","fpr")
#   AUC <- performance(prediction.obj=Pred,"auc")
#   AUC.all <- c(AUC.all,as.numeric(AUC@y.values)) }
# 
# AUC.breast <- c()
# for (i in c(1:N)){
#   Pred <- prediction(as.numeric(yhat.breast[[i]]),discrete.mek.response[rownames(yhat.breast[[i]])])
#   Perf <- performance(prediction.obj=Pred,"tpr","fpr")
#   AUC <- performance(prediction.obj=Pred,"auc")
#   AUC.breast <- c(AUC.breast,as.numeric(AUC@y.values)) }


abc <- matrix(NA,nrow=length(nsclc.mek.cells),ncol=N)
rownames(abc) <- nsclc.mek.cells
colnames(abc) <- c(1:N)
for(i in c(1:N)){abc[rownames(yhat.nsclc[[i]]),i] <- yhat.nsclc[[i]]}
tmp <- apply(abc,1,median,na.rm=TRUE)
Pred <- prediction(tmp,discrete.mek.response[nsclc.mek.cells])
Perf <- performance(prediction.obj=Pred,"tpr","fpr")
AUC <- performance(prediction.obj=Pred,"auc")
AUC.nsclc <- c(AUC.nsclc,as.numeric(AUC@y.values)) }

# AUC.crc <- c()
# for (i in c(1:N)){
#   Pred <- prediction(as.numeric(yhat.crc[[i]]),discrete.mek.response[rownames(yhat.crc[[i]])])
#   Perf <- performance(prediction.obj=Pred,"tpr","fpr")
#   AUC <- performance(prediction.obj=Pred,"auc")
#   AUC.crc <- c(AUC.crc,as.numeric(AUC@y.values)) }
# 
# AUC.glioma <- c()
# for (i in c(1:N)){
#   Pred <- prediction(as.numeric(yhat.glioma[[i]]),discrete.mek.response[rownames(yhat.glioma[[i]])])
#   Perf <- performance(prediction.obj=Pred,"tpr","fpr")
#   AUC <- performance(prediction.obj=Pred,"auc")
#   AUC.glioma <- c(AUC.glioma,as.numeric(AUC@y.values)) }
# 
# AUC.melanoma <- c()
# for (i in c(1:N)){
#   Pred <- prediction(as.numeric(yhat.melanoma[[i]]),discrete.mek.response[rownames(yhat.melanoma[[i]])])
#   Perf <- performance(prediction.obj=Pred,"tpr","fpr")
#   AUC <- performance(prediction.obj=Pred,"auc")
#   AUC.melanoma <- c(AUC.melanoma,as.numeric(AUC@y.values)) }
# 
# AUC.hemal <- c()
# for (i in c(1:N)){
#   Pred <- prediction(as.numeric(yhat.hemal[[i]]),discrete.mek.response[rownames(yhat.hemal[[i]])])
#   Perf <- performance(prediction.obj=Pred,"tpr","fpr")
#   AUC <- performance(prediction.obj=Pred,"auc")
#   AUC.hemal <- c(AUC.hemal,as.numeric(AUC@y.values)) }
# 
# 
# M <- list(ALL=AUC.all, BREAST=AUC.breast, NSCLC=AUC.nsclc, CRC=AUC.crc, MELANOMA=AUC.melanoma, GLIOMA=AUC.glioma, HEMAL=AUC.hemal)
# boxplot(M,outline=FALSE)
# stripchart(M,add=TRUE,vertical=TRUE,method="jitter",col="red")
