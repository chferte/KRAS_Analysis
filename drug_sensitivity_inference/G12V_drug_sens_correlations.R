# Charles Ferté
# Sage Bionetworks
# Dec 14th 2012


###################################################################################
# identify drug correlations between model G12V WT and CCLE and Sanger drug panel
###################################################################################

###############################################################################
# correlate the G12V with the drug response using the entire ccle lung dataset
###############################################################################
drug <- ccle_drug[names(KRAS_CCLE),]
drug.cor <- c()
drug.sig <- c()
p <- c()
D <- c()
par(mfrow=c(1,1),oma=c(4,2,2,2))
G12V_SENS <- sapply(colnames(drug),function(x){
  p <- cor(drug[,x],yhat_CCLE,method="spearman",use="pairwise.complete.obs") 
  D <- rbind(D,p)
})
colnames(G12V_SENS) <- colnames(drug)
drug.cor <- apply(G12V_SENS,2,function(x){median(x,na.rm=TRUE)})
tmp <- names(drug.cor)[sort(drug.cor,decreasing=TRUE,index.return=TRUE)$ix]
G12V_SENS <- G12V_SENS[,tmp]


##################################################################################
# compute the empirial p value for the CCLE corr
##################################################################################
res <- c()
emp.pval <- c()
drug <- ccle_drug[names(KRAS_CCLE),]
par(mfrow=c(1,1),oma=c(4,2,2,2))

for(k in colnames(drug)){
  drug.cor <- c()
  p <- c()
  D <- c()
  median.cor <- c()
  
  for(j in c(1:1000)){
    # print(j)     
    median.cor <- c(median.cor,median(cor(drug[sample(x=rownames(drug),replace=FALSE),k],yhat_CCLE,method="spearman",use="pairwise.complete.obs")))
  }
  print(k)
  emp.pval <- length(which(median.cor<median(G12V_SENS[,k])))/1000
  print(emp.pval)
  res <- c(res,emp.pval)
}

names(res) <- colnames(drug)
res <- res[colnames(G12V_SENS)]

G12V_SENS_CCLE_COR <- G12V_SENS
G12V_SENS_CCLE_PVAL <- res


###############################################################################
# plot the volcano plot for CCLE
###############################################################################

P <- .05
sig.drugs <- -log10(G12V_SENS_CCLE_PVAL) > -log10(P)
palette <- topo.colors(length(which(sig.drugs==TRUE)))
cols <- rep("gray60",length(G12V_SENS_CCLE_PVAL))
cols[sig.drugs] <- palette
pch <- rep(16,length(G12V_SENS_CCLE_PVAL))
pch[sig.drugs] <- 17
cex <- rep(.6, length(G12V_SENS_CCLE_PVAL))
cex[sig.drugs] <- 1
plot(apply(G12V_SENS_CCLE_COR,2,mean)*-1, -log10(G12V_SENS_CCLE_PVAL),type="p",pch=pch,cex=cex,col=cols,
     ylab="-log10(p)",xlab="spearman rho",xlim=c(-.5,.5))
abline(h=-log10(P),lty=2,col="gray60")
abline(v=c(-.5,.5),lty=2,col="gray60")
legend(-0.2,1,legend=names(G12V_SENS_CCLE_PVAL)[sig.drugs],
       pch=17,col=palette,cex=.8,xjust=.5)


###############################################################################
# correlate the G12V with the drug response using the entire SANGER lung dataset
###############################################################################
drug <- SangerDrug[names(KRAS_SANGER),]
drug.cor <- c()
drug.sig <- c()
p <- c()
D <- c()
par(mfrow=c(1,1),oma=c(4,2,2,2))
G12V_SENS <- sapply(colnames(drug),function(x){
  p <- cor(drug[,x],yhat_SANGER,method="spearman",use="pairwise.complete.obs") 
  D <- rbind(D,p)
})
colnames(G12V_SENS) <- colnames(drug)
drug.cor <- apply(G12V_SENS,2,function(x){median(x,na.rm=TRUE)})
tmp <- names(drug.cor)[sort(drug.cor,decreasing=TRUE,index.return=TRUE)$ix]
G12V_SENS <- G12V_SENS[,tmp]


##################################################################################
# compute the empirial p value for the SANGER corr
##################################################################################
res <- c()
emp.pval <- c()
drug <- SangerDrug[names(KRAS_SANGER),]
par(mfrow=c(1,1),oma=c(4,2,2,2))

for(k in colnames(drug)){
  drug.cor <- c()
  p <- c()
  D <- c()
  median.cor <- c()
  
  for(j in c(1:1000)){
    # print(j)     
    median.cor <- c(median.cor,median(cor(drug[sample(x=rownames(drug),replace=FALSE),k],yhat_SANGER,method="spearman",use="pairwise.complete.obs"),na.rm=TRUE))
  }
  print(k)
  emp.pval <- length(which(median.cor<median(G12V_SENS[,k])))/1000
  print(emp.pval)
  res <- c(res,emp.pval)
}

names(res) <- colnames(drug)
res <- res[colnames(G12V_SENS)]

G12V_SENS_SANGER_COR <- G12V_SENS
G12V_SENS_SANGER_PVAL <- res


###############################################################################
# plot the volcano plot for SANGER
###############################################################################
P <- .05
sig.drugs <- -log10(G12V_SENS_SANGER_PVAL) > -log10(P)
palette <- topo.colors(length(which(sig.drugs==TRUE)))
cols <- rep("gray60",length(G12V_SENS_SANGER_PVAL))
cols[sig.drugs] <- palette
pch <- rep(16,length(G12V_SENS_SANGER_PVAL))
pch[sig.drugs] <- 17
cex <- rep(.6, length(G12V_SENS_SANGER_PVAL))
cex[sig.drugs] <- 1
plot(apply(G12V_SENS_SANGER_COR,2,mean)*-1, -log10(G12V_SENS_SANGER_PVAL),type="p",pch=pch,cex=cex,col=cols,
     ylab="-log10(p)",xlab="spearman rho",xlim=c(-.5,.5))
abline(h=-log10(P),lty=2,col="gray60")
abline(v=c(-.5,.5),lty=2,col="gray60")
legend(-0.2,1.5,legend=names(G12V_SENS_SANGER_PVAL)[sig.drugs],
       pch=17,col=palette,cex=.8,xjust=.5)

###############################################################################
# save the correlation objects
###############################################################################

save(G12V_SENS_CCLE_COR,G12V_SENS_CCLE_PVAL,G12V_SENS_SANGER_COR,G12V_SENS_SANGER_PVAL,file="/home/cferte/FELLOW/cferte/KRAS_Analysis/results_drug_sens/G12V_WT_drug_correlations.rda")

