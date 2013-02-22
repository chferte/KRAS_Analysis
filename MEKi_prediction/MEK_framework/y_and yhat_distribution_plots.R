# charles ferté
# sage bionetworks
# Feb 22th, 2013

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
# same but with boxplots 
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
# plot the RMSE for the prediction of each cell lines
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
