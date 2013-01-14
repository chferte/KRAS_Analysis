# Charles Ferté
# dec 12 2012

# read drugs & targets file
ccle.drugs <- read.delim(file="/home/cferte/FELLOW/cferte/KRAS_Analysis/drug_sensitivity_inference/CCLE_drugs.txt",header=T,skip=2)
sanger.drugs <- read.delim(file="/home/cferte/FELLOW/cferte/KRAS_Analysis/drug_sensitivity_inference/Sanger_drugs.txt",header=T)
sanger.drugs <- sanger.drugs[which(duplicated(sanger.drugs$DRUG.NAME)!=TRUE),c("DRUG.NAME","TARGET")]
sanger.drugs$DRUG.NAME <- sub(pattern="-",replacement=".",x=sanger.drugs$DRUG.NAME)
# compute the comparison between the correlatiosn of the G12C-WT and G12V-WT both with CCLE and Sanger drug sensitivity

load("/home/cferte/FELLOW/cferte/KRAS_Analysis/results_drug_sens/G12C_WT_drug_correlations.rda")
load("/home/cferte/FELLOW/cferte/KRAS_Analysis/results_drug_sens/G12V_WT_drug_correlations.rda")


# make sure the correlations are coherently sorted
tmp <- intersect(colnames(G12C_SENS_CCLE_COR),colnames(G12V_SENS_CCLE_COR))
G12C_SENS_CCLE_COR <- G12C_SENS_CCLE_COR[,tmp]
G12V_SENS_CCLE_COR <- G12V_SENS_CCLE_COR[,tmp]
rm(tmp)

tmp <- intersect(colnames(G12C_SENS_SANGER_COR),colnames(G12V_SENS_SANGER_COR))
G12C_SENS_SANGER_COR <- G12C_SENS_SANGER_COR[,tmp]
G12V_SENS_SANGER_COR <- G12V_SENS_SANGER_COR[,tmp]
rm(tmp)


# function to apply the r to z transformation  (two sided test)
diff.corr <- function( r1, n1, r2, n2 ){ 
  
  Z1 <- 0.5 * log( (1+r1)/(1-r1) ) 
  Z2 <- 0.5 * log( (1+r2)/(1-r2) ) 
  
  diff   <- Z1 - Z2 
  SEdiff <- sqrt( 1/(n1 - 3) + 1/(n2 - 3) ) 
  diff.Z  <- diff/SEdiff 
  
  p <- 2*pnorm( abs(diff.Z), lower=F) 
  p
} 


#  look at the differences in the correlations with drug sensitivity in CCLE between G12C and G12V
N_G12C_CCLE <- 64
N_G12V_CCLE <-- 59
R_G12C <- apply(G12C_SENS_CCLE_COR,2,median)
R_G12V <- apply(G12V_SENS_CCLE_COR,2,median)

M <- R_G12C-R_G12V
p.val <- c()
for(i in c(1:length(colnames(G12C_SENS_CCLE_COR))))
{ p.val <- c(p.val,diff.corr(r1=R_G12C[i],r2=R_G12V[i],n1=N_G12C_CCLE,n2=N_G12V_CCLE))}
sort(p.val)
# 
# #  volcano plot of the r to z transformation
# P <- .05
# sig.drugs <- -log10(p.val) > -log10(P)
# palette <- topo.colors(length(which(sig.drugs==TRUE)))
# cols <- rep("gray60",length(M))
# cols[sig.drugs] <- palette
# pch <- rep(16,length(M))
# pch[sig.drugs] <- 17
# cex <- rep(.6, length(M))
# cex[sig.drugs] <- 1
# plot(M, -log10(p.val),type="p",pch=pch,cex=cex,col=cols,
#      ylab="-log10(p)",xlab="difference in the correlations (r-to-z)",xlim=c(-.5,.5))
# abline(h=-log10(P),lty=2,col="gray60")
# abline(v=c(-.5,.5),lty=2,col="gray60")
# legend(-0,1.5,legend=names(p.val)[sig.drugs],
#        pch=17,col=palette,cex=.8,xjust=.5)

# plot the median spearman rho of G12C and G12V for the drugs that exhibit differential sensitivity between G12C and G12V
#drug.names <- names(which(sig.drugs==TRUE))
drug.names <-   ccle.drugs$Compound..code.or.generic.name.[grep(pattern="MEK",ccle.drugs$Target.s.)]


COR.G12C <- G12C_SENS_CCLE_COR[,drug.names]
COR.G12V <- G12V_SENS_CCLE_COR[,drug.names]

par(mfrow=c(1,length(drug.names)))
for(i in c(1:length(drug.names))) {
tmp <- cbind(COR.G12C[,i],COR.G12V[,i])
colnames(tmp) <- c("G12C","G12V")
boxplot(tmp,outline=FALSE, main=drug.names[i], ylab="Correlation with drug IC50 (Spearman Rho)")
stripchart(list(G12C=tmp[,1], G12V=tmp[,2]), add=T,vertical=TRUE,method="jitter",col="red",pch=20)
}
p.val[drug.names]

#  SANGER
#  look at the differences in the correlations with drug sensitivity in SANGER between G12C and G12V
N_G12C_SANGER <- 98
N_G12V_SANGER <- 96
R_G12C <- apply(G12C_SENS_SANGER_COR,2,median)
R_G12V <- apply(G12V_SENS_SANGER_COR,2,median)


M <- R_G12C-R_G12V
p.val <- c()
for(i in c(1:length(colnames(G12C_SENS_SANGER_COR))))
{ p.val <- c(p.val,diff.corr(r1=R_G12C[i],r2=R_G12V[i],n1=N_G12C_SANGER,n2=N_G12V_SANGER))}

# #  volcano plot of the r to z transformation
# par(mfrow=c(1,1))
# P <- .05
# sig.drugs <- -log10(p.val) > -log10(P)
# palette <- topo.colors(length(which(sig.drugs==TRUE)))
# cols <- rep("gray60",length(M))
# cols[sig.drugs] <- palette
# pch <- rep(16,length(M))
# pch[sig.drugs] <- 17
# cex <- rep(.6, length(M))
# cex[sig.drugs] <- 1
# plot(M, -log10(p.val),type="p",pch=pch,cex=cex,col=cols,
#      ylab="-log10(p)",xlab="difference in the correlations (r-to-z)",xlim=c(-.5,.5))
# abline(h=-log10(P),lty=2,col="gray60")
# abline(v=c(-.5,.5),lty=2,col="gray60")
# legend(-0,1.5,legend=names(p.val)[sig.drugs],
#        pch=17,col=palette,cex=.8,xjust=.5)



# plot the median spearman rho of G12C and G12V for the drugs that exhibit differential sensitivity between G12C and G12V
#drug.names <- names(which(sig.drugs==TRUE))
drug.names <- sanger.drugs$DRUG.NAME[grep(pattern="MEK",sanger.drugs$TARGET)]
COR.G12C <- G12C_SENS_SANGER_COR[,drug.names]
COR.G12V <- G12V_SENS_SANGER_COR[,drug.names]

par(mfrow=c(2,2))
for(i in c(1:length(drug.names))) {
  tmp <- cbind(COR.G12C[,i],COR.G12V[,i])
  colnames(tmp) <- c("G12C","G12V")
  boxplot(tmp,outline=FALSE, main=drug.names[i], ylab="Correlation with drug IC50 (Spearman Rho)")
  stripchart(list(G12C=tmp[,1], G12V=tmp[,2]), add=T,vertical=TRUE,method="jitter",col="red",pch=20)
}
p.val[drug.names]
