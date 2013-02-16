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


All1 <- cor(abc,mek.ActArea[rownames(abc),1],method=method.cor,use="pairwise.complete.obs")
All2 <- cor(abc,mek.ActArea[rownames(abc),2],method=method.cor,use="pairwise.complete.obs")

boxplot(list(PD0325901=All1,AZD6244=All2), 
        main="All cell lines",
        ylab=paste(method.cor, "correlation"),ylim=c(0,1),outline=FALSE)
stripchart(list(PD0325901=All1,AZD6244=All2),vertical=TRUE,method="jitter",add=TRUE,col="aquamarine4",pch=20)
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


breast1 <- cor(abc,mek.ActArea[rownames(abc),1],method=method.cor,use="pairwise.complete.obs")
breast2 <- cor(abc,mek.ActArea[rownames(abc),2],method=method.cor,use="pairwise.complete.obs")

boxplot(list(PD0325901=breast1,AZD6244=breast2), 
        main="Breast",
        ylab=paste(method.cor, "correlation"),ylim=c(0,1),outline=FALSE)
stripchart(list(PD0325901=breast1,AZD6244=breast2),vertical=TRUE,method="jitter",add=TRUE,col="aquamarine4",pch=20)
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

nsclc1 <- cor(abc,mek.ActArea[rownames(abc),1],method=method.cor,use="pairwise.complete.obs")
nsclc2 <- cor(abc,mek.ActArea[rownames(abc),2],method=method.cor,use="pairwise.complete.obs")

boxplot(list(PD0325901=nsclc1,AZD6244=nsclc2), 
        main="NSCLC",
        ylab=paste(method.cor, "correlation"),ylim=c(0,1),outline=FALSE)
stripchart(list(PD0325901=nsclc1,AZD6244=nsclc2),vertical=TRUE,method="jitter",add=TRUE,col="aquamarine4",pch=20)
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

crc1 <- cor(abc,mek.ActArea[rownames(abc),1],method=method.cor,use="pairwise.complete.obs")
crc2 <- cor(abc,mek.ActArea[rownames(abc),2],method=method.cor,use="pairwise.complete.obs")

boxplot(list(PD0325901=crc1,AZD6244=crc2), 
        main="Colorectal",
        ylab=paste(method.cor, "correlation"),ylim=c(0,1),outline=FALSE)
stripchart(list(PD0325901=crc1,AZD6244=crc2),vertical=TRUE,method="jitter",add=TRUE,col="aquamarine4",pch=20)
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


melanoma1 <- cor(abc,mek.ActArea[rownames(abc),1],method=method.cor,use="pairwise.complete.obs")
melanoma2 <- cor(abc,mek.ActArea[rownames(abc),2],method=method.cor,use="pairwise.complete.obs")

boxplot(list(PD0325901=melanoma1,AZD6244=melanoma2), 
        main="Malignant melanoma",
        ylab=paste(method.cor, "correlation"),ylim=c(0,1),outline=FALSE)
stripchart(list(PD0325901=melanoma1,AZD6244=melanoma2),vertical=TRUE,method="jitter",add=TRUE,col="aquamarine4",pch=20)
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

glioma1 <- cor(abc,mek.ActArea[rownames(abc),1],method=method.cor,use="pairwise.complete.obs")
glioma2 <- cor(abc,mek.ActArea[rownames(abc),2],method=method.cor,use="pairwise.complete.obs")

boxplot(list(PD0325901=glioma1,AZD6244=glioma2), 
        main="Glioma",
        ylab=paste(method.cor, "correlation"),ylim=c(0,1),outline=FALSE)
stripchart(list(PD0325901=glioma1,AZD6244=glioma2),vertical=TRUE,method="jitter",add=TRUE,col="aquamarine4",pch=20)
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


hemal1 <- cor(abc,mek.ActArea[rownames(abc),1],method=method.cor,use="pairwise.complete.obs")
hemal2 <- cor(abc,mek.ActArea[rownames(abc),2],method=method.cor,use="pairwise.complete.obs")

boxplot(list(PD0325901=hemal1,AZD6244=hemal2), 
        main="Hematological\nMalignancies",
        ylab=paste(method.cor, "correlation"),ylim=c(0,1),outline=FALSE)
stripchart(list(PD0325901=hemal1,AZD6244=hemal2),vertical=TRUE,method="jitter",add=TRUE,col="aquamarine4",pch=20)
abline(h=seq(0,1,.1),lty=2)
