# charles ferté
# feb 28
# sage bionetworks

# script to plot the concordance between the IC50 and ActArea from Sanger and CCLE for the 2 MEK inhibs.

###########################################################
# load the data
###########################################################
source("/home/cferte/FELLOW/cferte/KRAS_Analysis/MEKi_prediction/MEK_framework/load_mek_ccle.R")
source("/home/cferte/FELLOW/cferte/KRAS_Analysis/MEKi_prediction/MEK_framework/load_mek_sanger.R")

###################################################################################################
# plot the concordant cell lines between ccle and sanger
###################################################################################################
par(mfrow=c(1,2))

plot(density(apply(mek.sanger.ic50[tmp,],1,mean)),main="sanger ic50")
abline(v=c(apply(mek.sanger.ic50[bad.concordant,],1,mean)),col="red")
abline(v=c(apply(mek.sanger.ic50[true.concordant,],1,mean)),col="blue")

plot(density(-1*apply(mek.ActArea[tmp,],1,mean)),main="ccle AA")
abline(v=c(-1*apply(mek.ActArea[bad.concordant,],1,mean)),col="red")
abline(v=c(-1*apply(mek.ActArea[true.concordant,],1,mean)),col="blue")
tmp <- intersect(rownames(mek.sanger.ic50),rownames(mek.ActArea))
x <- rank(apply(mek.sanger.ic50[tmp,],1,mean))
y <- rank(-1*apply(mek.ActArea[tmp,],1,mean))

par(mfrow=c(1,1), oma=c(0,0,0,0))
fit <- rlm(y~x)
plot(x,y,pch=20,xlab="sanger ranks",ylab= "ccle ranks",col="gray50")
abline(fit,lty=2)

good.bad <- rep(0,times=length(tmp))
names(good.bad) <- tmp
good.bad[true.concordant] <-  1
plot(apply(mek.sanger.ic50[tmp,],1,mean),apply(mek.ActArea[tmp,],1,mean),col=c("blue","red")[as.factor(good.bad)],pch=20)

# # asesss if the distribution of sensitivity is affected by the tissue type
# tissue.color <- as.numeric(factor(gsub("^.*?_(.*)","\\1",tmp)))
# plot(x,y,pch=20,col=rainbow(23)[as.numeric(tissue.color)])

# set the true concordant cells as the ones where the difference in ranks is in the first 33th percentile
true.concordant <- names(which(abs(x-y)<quantile(abs(x-y),probs=.33)))
bad.concordant <- tmp[-which(tmp%in%true.concordant)]
points(x[true.concordant],y[true.concordant],pch=19)
fit <- rlm(y[true.concordant]~x[true.concordant])
abline(fit,col="red",lwd=3)

# Map each sanger rank result to have the ActArea of the corresponding rank in ccle 
predicted.sanger.ActArea <- sort(apply(mek.ActArea[names(sort(y)),],1,mean))
names(predicted.sanger.ActArea) <- names(sort(x))

# Map each ccle rank result to have the ic50 of the corresponding rank in sanger 
predicted.ccle.ic50 <- apply(mek.sanger.ic50[names(sort(x)),],1,mean)
names(predicted.ccle.ic50) <- names(sort(y))

# # Map each sanger rank result to have the ActArea of the corresponding rank in ccle 
# # for each sanger rank , assign the predicted ActArea based on the true concordant fit
# predicted.sanger.ActArea <- c()
# predicted.rank <- as.numeric(format(x*fit$coefficients[2] +fit$coefficients[1],digits=1))
# names(predicted.rank) <- names(x)
# for(i in 1:length(predicted.rank)){
#   predicted.sanger.ActArea <-c(predicted.sanger.ActArea, mean(mek.ActArea[names(which(y==predicted.rank[i])),]))
# }
# names(predicted.sanger.ActArea) <- names(predicted.rank)

global.matrix2 <- global.matrix

#train model
N <- 20
PARAL1 <- mclapply(X=1:N,FUN=function(x){
train <- sample(true.concordant,replace=TRUE)
vec.train1 <-apply(ccle_drug[train,mek.inhib],1,mean)
cv.fit1 <- cv.glmnet(t(global.matrix[,train]), y=vec.train1,nfolds=5, alpha=.1)
fit1 <- glmnet(x=t(global.matrix[,train]),y=vec.train1,alpha=.1,lambda=cv.fit1$lambda.1se)
return(list(fit1,train))},mc.set.seed=TRUE,mc.cores=6)

# now investigate whether we do better or not if we train in true.concordant rather than in all
PARAL2 <- mclapply(X=1:N,FUN=function(x){
train <- sample(true.concordant,replace=TRUE)
vec.train2 <-apply(sanger_drug[train,mek.inhib],1,mean)
cv.fit2 <- cv.glmnet(t(global.matrix2[,train]), y=vec.train2,nfolds=5, alpha=.1)
fit2 <- glmnet(x=t(global.matrix2[,train]),y=vec.train2,alpha=.1,lambda=cv.fit2$lambda.1se)
return(list(fit2,train))},mc.set.seed=TRUE,mc.cores=6)


# compute the yhats of ccle and sanger
yhat.ccle <- c()
for(i in c(1:N)){
  train <- PARAL1[[i]][[2]]
  fit1 <- PARAL1[[i]][[1]]
  val <-true.concordant[-which(true.concordant %in% train)]
  #val <- bad.concordant
  yhat.ccle <- c(yhat.ccle,list(predict(fit1, t(global.matrix[,val]))))
}

yhat.sanger <- c()
for(i in c(1:N)){
  train <- PARAL2[[i]][[2]]
  fit2 <- PARAL2[[i]][[1]]
  val <-true.concordant[-which(true.concordant %in% train)]
  #val <- bad.concordant
  yhat.sanger <- c(yhat.sanger,list(predict(fit2, t(global.matrix2[,val]))))
}

# compute the RMSE for our strategy yhat - the actual actArea from ccle or from the  sanger ActArea infered from a true concordant model, 
# then, compute the RMSE using the same fit using the sanger IC50 and the ccle IC50 infered from sanger 
par(mfrow=c(1,2))
par(oma=c(5,0,5,0))                                 
RMSE1 <- sapply(1:N,function(z){sqrt(mean((yhat.ccle[[z]] - apply(ccle_drug[rownames(yhat.ccle[[z]]),mek.inhib],1,mean))^2))})
RMSE2 <- sapply(1:N,function(z){sqrt(mean((yhat.ccle[[z]] - predicted.sanger.ActArea[rownames(yhat.ccle[[z]])])^2))})
boxplot(list(ccle.ActArea=RMSE1,sanger.rank.infered.ActArea=RMSE2),las=2,cex.axis=.8,ylim=c(0,3))
stripchart(list(ccle.ActArea=RMSE1,sanger.rank.infered.ActArea=RMSE2),method="jitter",cex.axis=.8,vertical=TRUE,pch=20,col="red",add=TRUE,las=2,cex=.8,ylim=c(0,3))
                               
RMSE3 <- sapply(1:N,function(z){sqrt(mean((yhat.sanger[[z]] - apply(sanger_drug[rownames(yhat.sanger[[z]]),mek.inhib],1,mean))^2))})
RMSE4 <- sapply(1:N,function(z){sqrt(mean((yhat.sanger[[z]] - predicted.ccle.ic50[rownames(yhat.sanger[[z]])])^2))})
boxplot(list(sanger.ic50=RMSE3,ccle.rank.infered.ic50=RMSE4),las=2,cex.axis=.8,ylim=c(0,3))
stripchart(list(sanger.ic50=RMSE3,ccle.rank.infered.ic50=RMSE4),method="jitter",cex.axis=.8,vertical=TRUE,pch=20,col="blue",add=TRUE,las=2,cex=.8,ylim=c(0,3))

title("validation in the true concordant",outer=TRUE)
                                 