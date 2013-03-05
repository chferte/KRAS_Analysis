# charles ferté
# feb 28
# sage bionetworks

# script to plot the concordance between the IC50 and ActArea from Sanger and CCLE for the 2 MEK inhibs.

###########################################################
# load the data
###########################################################
source("/home/cferte/FELLOW/cferte/KRAS_Analysis/MEKi_prediction/MEK_framework/load_mek_sanger.R")
source("/home/cferte/FELLOW/cferte/KRAS_Analysis/MEKi_prediction/MEK_framework/load_mek_ccle.R")

###################################################################################################
# plot the concordant cell lines between ccle and sanger
###################################################################################################
tmp <- intersect(rownames(mek.sanger.ic50),rownames(mek.ActArea))
x <- rank(apply(mek.sanger.ic50[tmp,],1,mean))
y <- rank(-1*apply(mek.ActArea[tmp,],1,mean))

par(mfrow=c(1,2))
fit <- rlm(y~x)
#fit <- loess(y~x)
plot(x,y,pch=19,xlab="sanger ranks",ylab= "ccle ranks")
abline(fit)
#abline(a=fit$coefficients[1],b=fit$coefficients[2],col="orangered",lty=2,lwd=4)

# abline(v=quantile(abs(x-y),probs=c(.75)),col="red")


# # asesss if the distribution of sensitivity is affected by the tissue type
# tissue.color <- as.numeric(factor(gsub("^.*?_(.*)","\\1",tmp)))
# plot(x,y,pch=20,col=rainbow(23)[as.numeric(tissue.color)])



# plot the true concordant ones
plot(density(abs(x-y)))
abline(v=quantile(abs(x-y),probs=.4),col="red")
true.concordant <- names(which(abs(x-y)<quantile(abs(x-y),probs=.4)))
length(true.concordant)
fit <- rlm(y[true.concordant]~x[true.concordant])
abline(fit,col="red")
#true.concordant <- names(which(abs(fit$residuals)<quantile(abs(fit$residuals),probs=.5)))
plot(x[true.concordant],y[true.concordant],pch=19)

# plot the Actarea and IC50 of these cells
#plot(apply(mek.sanger.ic50[true.concordant,],1,mean),apply(mek.ActArea[true.concordant,],1,mean),xlab="sanger ic50",ylab="ccle ActArea",pch=20)
#fit <- rlm(apply(mek.ActArea[true.concordant,],1,mean)~apply(mek.sanger.ic50[true.concordant,],1,mean))
#new.fit <- rlm(y[true.concordant]~x[true.concordant])
#abline(new.fit)
#plot(x,y,pch=19)
#abline(a=new.fit$coefficients[1],b=new.fit$coefficients[2],lty=2,col="orangered",lwd=3)

# for each sanger rank result, infer the ActArea from our model:
predicted.sanger.ActArea <- c()
predicted.rank <- as.numeric(format(x*fit$coefficients[2] +fit$coefficients[1],digits=1))
names(predicted.rank) <- names(x)
for(i in 1:length(predicted.rank)){
  predicted.sanger.ActArea <-c(predicted.sanger.ActArea, mean(mek.ActArea[names(which(y==predicted.rank[i])),]))
}
names(predicted.sanger.ActArea) <- names(predicted.rank)
plot(x,predicted.sanger.ActArea)

#train model
N <- 10
PARAL1 <- mclapply(X=1:N,FUN=function(x){
train <- sample(true.concordant,replace=TRUE)
vec.train1 <-apply(ccle_drug[train,mek.inhib],1,mean)

cv.fit1 <- cv.glmnet(t(global.matrix[,train]), y=vec.train1,nfolds=5, alpha=.1)
fit1 <- glmnet(x=t(global.matrix[,train]),y=vec.train1,alpha=.1,lambda=cv.fit1$lambda.1se)
  
  return(list(fit1,train)) },
                  mc.set.seed=TRUE,mc.cores=6)

yhat.ccle <- c()
for(i in c(1:N)){
  train <- PARAL1[[i]][[2]]
  fit1 <- PARAL1[[i]][[1]]
  val <-tmp[-which(tmp %in% true.concordant)]
  yhat.ccle <- c(yhat.ccle,list(predict(fit1, t(global.matrix[,val]))))
}

# compute the RMSE for our strategy yhat - y ccle
RMSE1 <- sapply(1:N,function(z){sqrt(mean((yhat.ccle[[z]] - apply(ccle_drug[rownames(yhat.ccle[[z]]),mek.inhib],1,mean))^2))})
RMSE2 <- sapply(1:N,function(z){sqrt(mean((yhat.ccle[[z]] - predicted.sanger.ActArea[rownames(yhat.ccle[[z]])])^2))})


boxplot(list(RMSE1=RMSE1,RMSE2=RMSE2))
stripchart(list(RMSE1=RMSE1,RMSE2=RMSE2),method="jitter",vertical=TRUE,pch=19,col="red",add=TRUE)

