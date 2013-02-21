# charles ferté
# Feb 14th 2013
# sage bionetworks

# script to compute the correlations and display them



#############################################################################################################################################
# compute the correlations in All cell lines and display it!
#############################################################################################################################################

abc <- matrix(NA,ncol=N,nrow=length(mek.cells))
rownames(abc) <- mek.cells
for(i in c(1:length(yhat.all)))
{
  abc[rownames(yhat.all[[i]]),i]<- yhat.all[[i]]  
}


All1 <- cor(abc,apply(mek.ActArea[rownames(abc),],1,mean),method=method.cor,use="pairwise.complete.obs")

boxplot(All1, 
        main="All cell lines",
        ylab=paste(method.cor, "correlation"),ylim=c(0,1),outline=FALSE)
stripchart(All1,vertical=TRUE,method="jitter",add=TRUE,col="darkgreen",pch=20,cex=cex)
abline(h=seq(0,1,.1),lty=2)



#############################################################################################################################################
# compute the correlations in nsclc cell lines and display it!
#############################################################################################################################################
abc <- matrix(NA,ncol=N,nrow=length(nsclc.mek.cells))
rownames(abc) <- nsclc.mek.cells
for(i in c(1:length(yhat.nsclc)))
{
  abc[rownames(yhat.nsclc[[i]]),i]<- yhat.nsclc[[i]]  
}

nsclc1 <- cor(abc,apply(mek.ActArea[rownames(abc),],1,mean),method=method.cor,use="pairwise.complete.obs")

boxplot(nsclc1, 
        main="NSCLC",
        ylab=paste(method.cor, "correlation"),ylim=c(0,1),outline=FALSE)
stripchart(nsclc1,vertical=TRUE,method="jitter",add=TRUE,col="darkgreen",pch=20,cex=cex)
abline(h=seq(0,1,.1),lty=2)


#############################################################################################################################################
# compute the correlations in Breast cell lines and display it!
#############################################################################################################################################
abc <- matrix(NA,ncol=N,nrow=length(breast.mek.cells))
rownames(abc) <- breast.mek.cells
for(i in c(1:length(yhat.breast)))
{
  abc[rownames(yhat.breast[[i]]),i]<- yhat.breast[[i]]  
}


breast1 <- cor(abc,apply(mek.ActArea[rownames(abc),],1,mean),method=method.cor,use="pairwise.complete.obs")

boxplot(breast1, 
        main="Breast",
        ylab=paste(method.cor, "correlation"),ylim=c(0,1),outline=FALSE)
stripchart(breast1,vertical=TRUE,method="jitter",add=TRUE,col="darkgreen",pch=20,cex=cex)
abline(h=seq(0,1,.1),lty=2)

#############################################################################################################################################
# compute the correlations in crc cell lines and display it!
#############################################################################################################################################
abc <- matrix(NA,ncol=N,nrow=length(crc.mek.cells))
rownames(abc) <- crc.mek.cells
for(i in c(1:length(yhat.crc)))
{
  abc[rownames(yhat.crc[[i]]),i]<- yhat.crc[[i]]  
}

crc1 <- cor(abc,apply(mek.ActArea[rownames(abc),],1,mean),method=method.cor,use="pairwise.complete.obs")

boxplot(crc1, 
        main="COLORECTAL",
        ylab=paste(method.cor, "correlation"),ylim=c(0,1),outline=FALSE)
stripchart(crc1,vertical=TRUE,method="jitter",add=TRUE,col="darkgreen",pch=20,cex=cex)
abline(h=seq(0,1,.1),lty=2)


#############################################################################################################################################
# compute the correlations in hemal cell lines and display it!
#############################################################################################################################################
abc <- matrix(NA,ncol=N,nrow=length(hemal.mek.cells))
rownames(abc) <- hemal.mek.cells
for(i in c(1:length(yhat.hemal)))
{
  abc[rownames(yhat.hemal[[i]]),i]<- yhat.hemal[[i]]  
}


hemal1 <- cor(abc,apply(mek.ActArea[rownames(abc),],1,mean),method=method.cor,use="pairwise.complete.obs")

boxplot(hemal1, 
        main="Hematological\nMalignancies",
        ylab=paste(method.cor, "correlation"),ylim=c(0,1),outline=FALSE)
stripchart(hemal1,vertical=TRUE,method="jitter",add=TRUE,col="darkgreen",pch=20,cex=cex)
abline(h=seq(0,1,.1),lty=2)

#############################################################################################################################################
# compute the correlations in glioma cell lines and display it!
#############################################################################################################################################
abc <- matrix(NA,ncol=N,nrow=length(glioma.mek.cells))
rownames(abc) <- glioma.mek.cells
for(i in c(1:length(yhat.glioma)))
{
  abc[rownames(yhat.glioma[[i]]),i]<- yhat.glioma[[i]]  
}

glioma1 <- cor(abc,apply(mek.ActArea[rownames(abc),],1,mean),method=method.cor,use="pairwise.complete.obs")

boxplot(glioma1, 
        main="GLIOMA",
        ylab=paste(method.cor, "correlation"),ylim=c(0,1),outline=FALSE)
stripchart(glioma1,vertical=TRUE,method="jitter",add=TRUE,col="darkgreen",pch=20,cex=cex)
abline(h=seq(0,1,.1),lty=2)


#############################################################################################################################################
# compute the correlations in melanoma cell lines and display it!
#############################################################################################################################################
abc <- matrix(NA,ncol=N,nrow=length(melanoma.mek.cells))
rownames(abc) <- melanoma.mek.cells
for(i in c(1:length(yhat.melanoma)))
{
  abc[rownames(yhat.melanoma[[i]]),i]<- yhat.melanoma[[i]]  
}


melanoma1 <- cor(abc,apply(mek.ActArea[rownames(abc),],1,mean),method=method.cor,use="pairwise.complete.obs")

boxplot(melanoma1, 
        main="MELANOMA",
        ylab=paste(method.cor, "correlation"),ylim=c(0,1),outline=FALSE)
stripchart(melanoma1,vertical=TRUE,method="jitter",add=TRUE,col="darkgreen",pch=20,cex=cex)
abline(h=seq(0,1,.1),lty=2)

