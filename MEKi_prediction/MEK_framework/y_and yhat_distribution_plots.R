# charles ferté
# sage bionetworks
# Feb 22th, 2013


#############################
# load the data
#############################
source("/home/cferte/FELLOW/cferte/KRAS_Analysis/MEKi_prediction/MEK_framework/load_mek_ccle.R")


################################################################################################################
# compute mixture models to assign sensitive and resistant cells
################################################################################################################

require(flexmix)

Nclust <- 2

cells <- list(mek.cells,nsclc.mek.cells,breast.mek.cells,crc.mek.cells,hemal.mek.cells,glioma.mek.cells,melanoma.mek.cells)
cell.names <- list("ALL CELLS","NSCLC","BREAST","CRC","Hematologic\nMalignancies","GLIOMA","MELANOMA")

#################
all.prob <- c()
all.prob.sens <- c()
all.prob.res <- c()
total.cell.status <-c()
for(i in c(1:length(cells))){
  
  blah <- apply(ccle_drug[cells[[i]],mek.inhib],1,mean)
  blah <- blah[names(sort(blah))]
  rank(blah)
  Nclust <- 2
  
  tissue.post.prob <- c()
  K <- 30
  while(length(tissue.post.prob)<K) {
    sample.out <- as.numeric(format(length(names(blah))*.07,digits=1))
    iter.cells <- sample(names(blah),replace=FALSE,size=length(names(blah))-sample.out)
    #iter.cells <- sample(rownames(blah),replace=FALSE,size=length(rownames(blah))-2)
    #tryclust <- flexmix(blah[iter.cells,]~1,k=Nclust, model=FLXMCmvnorm(diag=FALSE))
    tryclust <- flexmix(blah[iter.cells]~1,k=Nclust)
    clust <- tryclust@cluster
    table(clust)
    if(length(table(clust))<Nclust){} else{
    names(clust) <- iter.cells
    
    post.prob <- tryclust@posterior$scaled
    if(mean(blah[names(which(clust==1))]) > mean(blah[names(which(clust==2))])) {
      colnames(post.prob) <- c("sens","res")
    } else { colnames(post.prob) <- c("res","sens")}
    post.prob <- post.prob[,c("sens","res")]
    rownames(post.prob) <- iter.cells
    tissue.post.prob  <- c(tissue.post.prob,list(post.prob))
    }
  }
  

final.prob.sens <- matrix(NA,ncol=length(tissue.post.prob),nrow=length(cells[[i]])) 
rownames(final.prob.sens) <- names(blah)
for(x in  1:length(tissue.post.prob)){final.prob.sens[rownames(tissue.post.prob[[x]]),x]<-(tissue.post.prob[[x]][,"sens"])}
  
  all.prob.sens <- c(all.prob.sens,list(final.prob.sens))
  
  final.prob.res <- matrix(NA,ncol=length(tissue.post.prob),nrow=length(cells[[i]])) 
  rownames(final.prob.res) <- names(blah)
  for(x in  1:length(tissue.post.prob)){final.prob.res[rownames(tissue.post.prob[[x]]),x]<-(tissue.post.prob[[x]][,"res"])}
  all.prob.res <- c(all.prob.res,list(final.prob.res))
  
final.prob <- cbind(apply(final.prob.sens,1,function(x){exp(mean(log(x[!is.na(x)])))}),apply(final.prob.res,1,function(x){exp(mean(log(x[!is.na(x)])))}))
  colnames(final.prob) <- c("sens","res")
  
  all.prob <- c(all.prob,list(final.prob))

  }


###################
# plot it 
##################
par(mfrow=c(2,4),oma=c(0,0,0,0))
for(i in 1:length(all.prob))
{
#col.vec <- (apply(all.prob[[i]],1,function(x){names(which(x==max(x)))}))
#col.vec <- ifelse(col.vec=="sens",1,2)
col.vec <- factor(-1*(all.prob[[i]][,1]+1))
plot(density(apply(ccle_drug[cells[[i]],mek.inhib],1,mean)),lwd=3,main=paste(cell.names[[i]]),xlim=c(-1.5,8),ylim=c(0,1))
ypos <- abs(rnorm(sd=.04,mean=.8,n=length(rownames(all.prob[[i]]))))
points(apply(ccle_drug[rownames(all.prob[[i]]),mek.inhib],1,mean),ypos,col=greenred(length(col.vec))[col.vec],pch=20,cex=1)
}


# cell status: sensitive are the ones >.6, res are <.4, inter are the rest !
cell.status <- sapply(1:length(cell.names),function(x){ifelse(all.prob[[x]][,1]>.6,"sens",ifelse(all.prob[[x]][,1]<.4,"res","inter"))})


###########################################################################################################################
tmp <- sapply(c(1:N),function(x){return(PARAL[[x]][[2]])})
plot(density(apply(ccle_drug[tmp,mek.inhib],1,mean)),xlim=c(-2,8.5),main="training set",ylim=c(0,.8),lwd=3,col="red")
title(main="Distribution of the sensitivity of cell lines to MEK inhibitors according to tissue type \n(sensitivity assessed by ActArea)",
      sub="sensitivity was assessed by ActArea",outer=TRUE)

################################################################################################################
# plot the distribution of yhats of MEK ActArea (ie: plot the y)
################################################################################################################

par(mfrow=c(2,4),oma=c(0,0,4,0))

cells <- list(mek.cells,nsclc.mek.cells,breast.mek.cells,crc.mek.cells,hemal.mek.cells,glioma.mek.cells,melanoma.mek.cells)
cell.names <- list("ALL","NSCLC","BREAST","CRC","Hematologic\nMalignancies","GLIOMA","MELANOMA")
yhats <- list(yhat.all,yhat.nsclc,yhat.breast,yhat.crc,yhat.hemal,yhat.glioma,yhat.melanoma)
threshold <- .8
for(j in c(1:length(cells)))
  {
abc <- matrix(NA,nrow=length(cells[[j]]),ncol=N)
rownames(abc) <- cells[[j]]
colnames(abc) <- c(1:N)
for(i in c(1:N)){abc[rownames(yhats[[j]][[i]]),i] <- yhats[[j]][[i]]}
abc <- apply(abc,1,median,na.rm=TRUE)
plot(density(abc),xlim=c(-2,5.5),main=paste(cell.names[[j]],"yhats"),ylim=c(0,1.4),lwd=2)
abline(v=quantile(abc,probs=c(.2,.8)),col="red",lty=3)
top.predicted.cells <- ifelse(abc>quantile(abc,probs=threshold),1,0)
def <- apply(ccle_drug[cells[[j]],mek.inhib],1,mean)
top.true.cells <- ifelse(def>quantile(def,probs=threshold),1,0)
cat("predictive performance when fixing a threshold at the 80th quantile\n",paste("for ",cell.names[[j]]),"cell lines")
PPV <- length(which(top.predicted.cells==1 & top.true.cells==1))/length(which(top.predicted.cells==1))
NPV <- length(which(top.predicted.cells==0 & top.true.cells==0))/length(which(top.predicted.cells==0))
Accuracy <- (length(which(top.predicted.cells==1 & top.true.cells==1)) + length(which(top.predicted.cells==0 & top.true.cells==0)))/length(top.predicted.cells)
print(table(PREDICTED=top.predicted.cells,REALITY=top.true.cells))
print(paste("PPV=",PPV))
print(paste("NPV=",NPV))
print(paste("Accuracy=",Accuracy))
}

title(main="Distribution of the yhats of the MEK inhibitors according to tissue type \n(sensitivity assessed by ActArea)",outer=TRUE)


#################################################################################################################
# plot the distribution of yhats of MEK ActArea (ie: plot the y)
################################################################################################################

par(mfrow=c(2,3),oma=c(0,0,0,0))
cells <- list(mek.cells,nsclc.mek.cells,breast.mek.cells,crc.mek.cells,hemal.mek.cells,glioma.mek.cells,melanoma.mek.cells)
cell.names <- list("ALL CELLS","NSCLC","BREAST","CRC","Hematologic\nMalignancies","GLIOMA","MELANOMA")
yhats <- list(yhat.all,yhat.nsclc,yhat.breast,yhat.crc,yhat.hemal,yhat.glioma,yhat.melanoma)

for(j in c(1:length(cells)))
{
  PPV <- c()
  NPV <- c()
  Accuracy <- c()
  threshold <- seq(from=0,to=1,by=.005)
  for(k in threshold){
   abc <- matrix(NA,nrow=length(cells[[j]]),ncol=N)
  rownames(abc) <- cells[[j]]
  colnames(abc) <- c(1:N)
  for(i in c(1:N)){abc[rownames(yhats[[j]][[i]]),i] <- yhats[[j]][[i]]}
  
   abc <- apply(abc,1,median,na.rm=TRUE)
    top.predicted.cells <- ifelse(abc>=quantile(abc,probs=k),1,0)
  def <- apply(ccle_drug[cells[[j]],mek.inhib],1,mean)
  top.true.cells <- ifelse(def>=quantile(def,probs=k),1,0)
  PPV <- c(PPV,length(which(top.predicted.cells==1 & top.true.cells==1))/length(which(top.predicted.cells==1)))
  NPV <- c(NPV,length(which(top.predicted.cells==0 & top.true.cells==0))/length(which(top.predicted.cells==0)))
  Accuracy <- c(Accuracy,(length(which(top.predicted.cells==1 & top.true.cells==1)) + length(which(top.predicted.cells==0 & top.true.cells==0)))/length(top.predicted.cells))
 }
  plot(threshold,PPV,main=paste(cell.names[[j]]),ylim=c(0,1),type="l",lwd=3, cex.axis=.9)
  abline(h=c(.2,.4,.6,.8),lty=2,cex=.8,col="gray60")  
  abline(v=c(.2,.4,.6,.8),lty=2,cex=.8,col="gray60")
  plot(threshold,NPV,main=paste(cell.names[[j]]),ylim=c(0,1),type="l",lwd=3, cex.axis=.9)
  abline(h=c(.2,.4,.6,.8),lty=2,cex=.8,col="gray60")  
  abline(v=c(.2,.4,.6,.8),lty=2,cex=.8,col="gray60")
  plot(threshold,Accuracy,main=paste(cell.names[[j]]),ylim=c(0,1),type="l",lwd=3, cex.axis=.9)  
  abline(h=c(.2,.4,.6,.8),lty=2,cex=.8,col="gray60")  
  abline(v=c(.2,.4,.6,.8),lty=2,cex=.8,col="gray60")
}

################################################################################################################
# plot the RMSE for the prediction of each cell lines
################################################################################################################
par(oma=c(0,0,6,0))
# first plot the density of the IC50
k <- list(mek.cells,nsclc.mek.cells,breast.mek.cells,crc.mek.cells,hemal.mek.cells,glioma.mek.cells,melanoma.mek.cells)
yhats <- list(yhat.all,yhat.nsclc,yhat.breast,yhat.crc,yhat.hemal,yhat.glioma,yhat.melanoma)
names(k) <- c("ALL CELLS","NSCLC","BREAST","COLORECTAL","Hematologic\nMalignancies","GLIOMA","MELANOMA")
par(mfrow=c(2,4))
for(i in c(1:length(k)))
{
 
  tissue.ic50 <- apply(ccle_drug[k[[i]],mek.inhib],1,mean)
  plot(density(tissue.ic50),main=paste(names(k)[i]),ylim=c(0,3.5),xlim=c(-2,8.5),lwd=1) 
  
  #  create a matrix abc containing all the yhat per cells
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
  lines(tissue.ic50[tmp],RMSE[tmp],col="coral1",type="p",pch=20,cex=.7)
  axis(4, ylim=c(0,max(RMSE[tmp])+1),,col="red",col.axis="red")
  q <- loess(RMSE[tmp]~tissue.ic50[tmp])
  lines(q$x,q$fitted,col="red",lwd=3)
  #lines(density(tissue.ic50),main=paste(names(k)[i]),ylim=c(0,3.5),xlim=c(-2,8.5),lwd=2)
 lines(density(tissue.ic50),main=paste(names(k)[i]),lwd=3)
}

title(main="Distribution of the sensitivity of cell lines to MEK inhibitors according to tissue type \n(sensitivity assessed by ActArea)",,outer=TRUE)


