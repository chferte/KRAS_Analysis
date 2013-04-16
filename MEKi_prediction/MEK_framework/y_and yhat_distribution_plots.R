# charles ferté
# sage bionetworks
# Feb 22th, 2013

library(ROCR)

#############################
# load the data
#############################
source("/home/cferte/FELLOW/cferte/KRAS_Analysis/MEKi_prediction/MEK_framework/load_mek_ccle.R")
ccle_probs_status <- loadEntity("syn1709732")
ccle_probs_status <- ccle_probs_status$objects$ccle_probs_status
all.prob <- ccle_probs_status[[1]]
cell.status <- ccle_probs_status[[2]]



ccle_probs_status <- loadEntity("syn1709732")
ccle_probs_status <- ccle_probs_status$objects$ccle_probs_status
all.prob <- ccle_probs_status[[1]]
cell.status <- ccle_probs_status[[2]]

cells <- list(mek.cells,nsclc.mek.cells,breast.mek.cells,crc.mek.cells,hemal.mek.cells,melanoma.mek.cells,pancreas.mek.cells,ovary.mek.cells)
cell.names <- list("ALL CELLS","NSCLC","BREAST","CRC","Hematologic\nMalignancies","MELANOMA","PANCREAS","OVARY")

load("/home/cferte/RESULTS/MEKi/GLOBAL_MODEL/ROBJECTS/global_model_yhats.Rda")
yhats <- global_model_yhats


#load("/home/cferte/RESULTS/MEKi/GLOBAL_MODEL/ROBJECTS/lkb1_model_yhats.Rda")
#yhats <- lkb1_model_yhats

#yhats <- list(yhat.all,yhat.nsclc,yhat.breast,yhat.crc,yhat.hemal,yhat.glioma,yhat.melanoma)
#yhats <- pure_tissue_models_yhats

#yhats <- global_model_tissue_adjusted_yhats

#load("/home/cferte/RESULTS/MEKi/balanced_model_yhats.Rda")
#yhats <- balanced_model_yhats

# cell status: sensitive are the ones >.6, res are <.4, inter are the rest !
cell.status <- c()
limit.sens <- .55
limit.res <- .35
for(i in 1:length(cell.names)){
foo <- names(which(all.prob[[i]][,1]>=limit.sens))
min.sens.act.area <- min(apply(ccle_drug[foo,mek.inhib],1,mean))
foo <- names(which(all.prob[[i]][,1]<=limit.res))
max.res.act.area <- max(apply(ccle_drug[foo,mek.inhib],1,mean))
cell.status <- c(cell.status,
                 list(ifelse(apply(ccle_drug[cells[[i]],mek.inhib],1,mean) >= min.sens.act.area,"sens",
                             ifelse(apply(ccle_drug[cells[[i]],mek.inhib],1,mean) <= max.res.act.area,
                                    "res","inter"))))
}
names(cell.status) <- cell.names

# # try with the top 25 % percentile and lowest 25% percentile as sensitive and resistant, respectively
# gold <- cell.status
# for(i in 1:length(cell.names)){
# 
# foo <- apply(ccle_drug[cells[[i]],mek.inhib],1,mean)
# foo2 <- quantile(foo,probs=c(.25,.75))
# gold[[i]][names(which(foo>=foo2[2]))] <- "sens"
# gold[[i]][names(which(foo<=foo2[1]))] <- "res"
# gold[[i]][names(which(foo>foo2[1] & foo<foo2[2]))] <- "inter"
# 
# print(cell.names[[i]])
# print(table(gold[[i]]))
# }
# 
# cell.status <- gold
# rm(foo,foo2,gold)


#############################
# plot the ActArea & the status
#############################
par(mfrow=c(2,4),oma=c(0,0,4,0))
for(i in 1:length(all.prob))
{
  print(cell.names[[i]])
  print(table(cell.status[[i]]))
  #col.vec <- (apply(all.prob[[i]],1,function(x){names(which(x==max(x)))}))
  #col.vec <- ifelse(col.vec=="sens",1,2)
  #col.vec <- factor(-1*(all.prob[[i]][,1]+1))
  col.vec <- ifelse(cell.status[[i]]=="sens","green",ifelse(cell.status[[i]]=="res","red","grey60"))
  plot(density(apply(ccle_drug[cells[[i]],mek.inhib],1,mean)),lwd=3,main=paste(cell.names[[i]]),xlim=c(-1.5,8),ylim=c(0,1))
  ypos <- abs(rnorm(sd=.04,mean=.8,n=length(rownames(all.prob[[i]]))))
  points(apply(ccle_drug[cells[[i]],mek.inhib],1,mean),ypos,col=col.vec,pch=20,cex=1.1)
  #points(apply(ccle_drug[rownames(all.prob[[i]]),mek.inhib],1,mean),ypos,col=greenred(length(col.vec))[col.vec],pch=20,cex=1.1)
}

title(main="Distribution of the sensitivity of cell lines to MEK inhibitors according to tissue type \n(sensitivity assessed by ActArea, red & green colors stand for resistant & sensitive cells)",
      ,outer=TRUE)

################################################################################################################
# plot the distribution of yhats of MEK ActArea (ie: plot the y)
################################################################################################################
N <- length(yhats[[1]])
par(mfrow=c(2,4),oma=c(0,0,4,0))
threshold <- .8
for(j in c(1:length(cells)))
  {
abc <- c()
  abc <- as.data.frame(matrix(NA,nrow=length(cells[[j]]),ncol=N))
rownames(abc) <- cells[[j]]
colnames(abc) <- c(1:N)
for(l in c(1:N)){abc[rownames(yhats[[j]][[l]]),l] <- yhats[[j]][[l]]}
abc <- apply(abc,1,median,na.rm=TRUE)
plot(density(abc,na.rm=TRUE),xlim=c(-2,5.5),main=paste(cell.names[[j]],"yhats"),ylim=c(0,1.4),lwd=2)
abline(v=quantile(abc,probs=c(.2,.8),na.rm=TRUE),col="red",lty=3)
top.predicted.cells <- ifelse(abc>=quantile(abc,probs=threshold,na.rm=TRUE),1,0)
#def <- apply(ccle_drug[cells[[j]],mek.inhib],1,mean)
#top.true.cells <- ifelse(def>quantile(def,probs=threshold),1,0)
top.true.cells <- ifelse(cell.status[[j]]=="sens",1,0)
tmp <- intersect(names(top.true.cells),names(top.predicted.cells))
top.predicted.cells <- top.predicted.cells[tmp]
top.true.cells <- top.true.cells[tmp]
cat("predictive performance when fixing a threshold at the 80th quantile\n",paste("for ",cell.names[[j]]),"cell lines")
PPV <- length(which(top.predicted.cells==1 & top.true.cells==1))/length(which(top.predicted.cells==1))
NPV <- length(which(top.predicted.cells==0 & top.true.cells==0))/length(which(top.predicted.cells==0))
Accuracy <- (length(which(top.predicted.cells==1 & top.true.cells==1)) + length(which(top.predicted.cells==0 & top.true.cells==0)))/length(top.predicted.cells)


print(table(PREDICTED=top.predicted.cells,REALITY=top.true.cells))
print(paste("Accuracy=",format(Accuracy,digits=2)))
print(paste("NPV=",format(NPV,digits=2)))
print(paste("PPV=",format(PPV,digits=2)))
}

title(main="Distribution of the yhats of the MEK inhibitors according to tissue type \n(sensitivity assessed by ActArea)",outer=TRUE)

################################################################################################################
# plot the yhats (x axis) and MEK ActArea (y axis)
################################################################################################################
par(mfrow=c(2,4))
for(i in 1:length(cell.names)){
  
  y <- apply(ccle_drug[cells[[i]],mek.inhib],1,mean)
  abc <- c()
    abc <- as.data.frame(matrix(NA,nrow=length(cells[[i]]),ncol=N))
    rownames(abc) <- cells[[i]]
    colnames(abc) <- c(1:N)
    for(l in c(1:N)){abc[rownames(yhats[[i]][[l]]),l] <- yhats[[i]][[l]]}
    abc <- apply(abc,1,median,na.rm=TRUE)
  x <- abc
  rm(abc)
plot(x,y,pch=20,main=paste(cell.names[[i]]),
     xlab="yhats",ylab="ActArea",ylim=c(0,8))
abline(col="red",lw=2,lty=2,h=c(min(y[names(which(cell.status[[i]]=="sens"))]),max(y[names(which(cell.status[[i]]=="res"))])))
  
#abline(h=quantile(y,probs=c(.25,.75)),col="red",lw=2,lty=2)
g <- loess(formula=y~x)
lines(lwd=3,g$x[sort(g$x,index.return=TRUE)$ix],
      g$fitted[sort(g$x,index.return=TRUE)$ix],col="royalblue1")
print(cell.names[[i]])
  print(cor.test(x,y,method="spearman")$p.value)
  print(cor.test(x,y,method="spearman")$estimate)
  abline(v=quantile(x,probs=c(.2,.8),na.rm=TRUE),col="green",lwd=3)
  
}

title(main="Distribution of the sensitivity to the MEK inhibitors according their preedictied value (yhats) \n(sensitivity assessed by ActArea)",outer=TRUE)



#################################################################################################################
# plot the distribution of yhats of MEK ActArea (ie: plot the y)
################################################################################################################

par(mfrow=c(4,5),oma=c(0,0,0,0))
#cells <- list(mek.cells,nsclc.mek.cells,breast.mek.cells,crc.mek.cells,hemal.mek.cells,glioma.mek.cells,melanoma.mek.cells)
#cell.names <- list("ALL CELLS","NSCLC","BREAST","CRC","Hematologic\nMalignancies","GLIOMA","MELANOMA")
#yhats <- list(yhat.all,yhat.nsclc,yhat.breast,yhat.crc,yhat.hemal,yhat.glioma,yhat.melanoma)
#yhats <- pure_tissue_models_yhats
for(j in c(1:length(cells)))
{
  PPV <- c()
  NPV <- c()
  Accuracy <- c()
  threshold <- seq(from=0,to=1,by=.01)
  for(k in threshold){
  
  #assemble the yhats in a median vector abc  
    abc <- matrix(NA,nrow=length(cells[[j]]),ncol=N)
  rownames(abc) <- cells[[j]]
  colnames(abc) <- c(1:N)
  for(i in c(1:N)){abc[rownames(yhats[[j]][[i]]),i] <- yhats[[j]][[i]]}
   predicted.sensitivity <- apply(abc,1,median,na.rm=TRUE)
  top.predicted.cells <- ifelse(predicted.sensitivity>=quantile(predicted.sensitivity,probs=k),1,0)
  bottom.predicted.cells <- ifelse(predicted.sensitivity<=quantile(predicted.sensitivity,probs=k),1,0)
    #def <- apply(ccle_drug[cells[[j]],mek.inhib],1,mean)
  #top.true.cells <- ifelse(def>=quantile(def,probs=k),1,0)
   top.true.cells <- ifelse(cell.status[[j]]=="sens",1,0)
    bottom.true.cells <- ifelse(cell.status[[j]]=="res",1,0)
    tmp <- rownames(abc)
    top.predicted.cells <- top.predicted.cells[tmp]
    top.true.cells <- top.true.cells[tmp]
    bottom.predicted.cells <- bottom.predicted.cells[tmp]
    bottom.true.cells <- bottom.true.cells[tmp]
   PPV <- c(PPV,length(which(top.predicted.cells==1 & top.true.cells==1))/length(which(top.predicted.cells==1)))
  NPV <- c(NPV,length(which(bottom.predicted.cells==1 & bottom.true.cells==1))/length(which(bottom.predicted.cells==1)))
  Accuracy <- c(Accuracy,(length(which(top.predicted.cells==1 & top.true.cells==1)) + length(which(top.predicted.cells==0 & top.true.cells==0)))/length(top.predicted.cells))
   pred <- prediction(predicted.sensitivity, top.true.cells)
   perf <- performance(pred, measure = "tpr", x.measure = "fpr")
  auc <- performance(pred,"auc")
    auc <- format(unlist(slot(auc, "y.values")),digits=2)
  cor <- cor(abc,all.prob[[j]][rownames(abc),1],method="spearman",use="pairwise.complete.obs")
  
  }
  title.name <- paste(cell.names[[j]],"n=",length(cells[[j]]))
  boxplot(cor,main=paste(title.name,"\nCorrelation"),outline=FALSE,ylab="spearman rho",ylim=c(-0.4,1))
  abline(h=c(-.4,-.2,0,.2,.4,.6,.8,1),lty=2,cex=.8,col="gray60")
  stripchart(cor,col="coral",add=TRUE,vertical=TRUE,method="jitter",pch=19)
  print(paste("correlation", cell.names[[j]],format(median(cor,na.rm=TRUE),digits=2)))
  
  plot(perf, lwd=3,main=paste(title.name,"\nAUC"),col="coral",xlab="False positive rate (1-Specificity)",ylab="True positive rate (sensitivity)")
  text(x=.6,y=.2,labels=paste("AUC=",auc),cex=.8,font=2)
  abline(b=1,a=0,lty=2,cex=.8,col="gray60")
  print(paste("AUC=", cell.names[[j]],auc))

  plot(threshold,Accuracy,main=paste(title.name,"\nAccuracy"),ylim=c(0,1),type="l",lwd=3, col="coral",cex.axis=.9)  
  abline(v=c(.2,.4,.6,.8),lty=2,cex=.8,col="gray60")  
  abline(h=c(.2,.4,.6,.8),lty=2,cex=.8,col="gray60")
  
  plot(threshold,NPV,main=paste(title.name,"\nNegative Predicted value"),ylim=c(0,1),col="red",type="l",lwd=3, cex.axis=.9,xlim=c(0,1))
  abline(v=c(0,.2,.4,.6,.8,1),lty=2,cex=.8,col="gray60")  
  abline(h=c(.2,.4,.6,.8),lty=2,cex=.8,col="gray60")

  plot(threshold,PPV,main=paste(title.name,"\nPositive Predicted Value"),col="green",ylim=c(0,1),type="l",lwd=3, cex.axis=.9,xlim=c(0,1))
  abline(v=c(0,.2,.4,.6,.8,1),lty=2,cex=.8,col="gray60")  
  abline(h=c(.2,.4,.6,.8),lty=2,cex=.8,col="gray60")
}


################################################################################################################
# plot the RMSE for the prediction of each cell lines
################################################################################################################
par(oma=c(0,0,6,0))
# first plot the density of the IC50
k <- cells
names(k) <- cell.names
par(mfrow=c(2,4))
for(i in c(1:length(k)))
{
 
  tissue.ic50 <- apply(ccle_drug[k[[i]],mek.inhib],1,mean)
  plot(density(tissue.ic50),main=paste(names(k)[i]),ylim=c(0,3.5),xlim=c(-2,8.5),lwd=1,ylab="density of the ActArea") 
  
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
  lines(tissue.ic50[tmp],RMSE[tmp],col="red",type="p",pch=19,cex=.7)
  axis(4, ylim=c(0,max(RMSE[tmp])+1),,col="coral",col.axis="red",ylab="RMSE")
  q <- loess(RMSE[tmp]~tissue.ic50[tmp])
  lines(q$x,q$fitted,col="royalblue",lwd=3)
  #lines(density(tissue.ic50),main=paste(names(k)[i]),ylim=c(0,3.5),xlim=c(-2,8.5),lwd=2)
 lines(density(tissue.ic50),main=paste(names(k)[i]),lwd=3)
}

title(main="Distribution of the sensitivity of cell lines to MEK inhibitors according to tissue type \n(sensitivity assessed by ActArea)",,outer=TRUE)


