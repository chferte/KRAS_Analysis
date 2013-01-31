# Charles Ferté
# dec 12 2012

# read drugs & targets file
ccle.drugs <- read.delim(file="/home/cferte/FELLOW/cferte/KRAS_Analysis/drug_sensitivity_inference/CCLE_drugs.txt",header=T,skip=2)
sanger.drugs <- read.delim(file="/home/cferte/FELLOW/cferte/KRAS_Analysis/drug_sensitivity_inference/Sanger_drugs.txt",header=T)
sanger.drugs <- sanger.drugs[which(duplicated(sanger.drugs$DRUG.NAME)!=TRUE),c("DRUG.NAME","TARGET")]
sanger.drugs$DRUG.NAME <- sub(pattern="-",replacement=".",x=sanger.drugs$DRUG.NAME)

# compute the comparison between the correlations of the G12C-WT and G12V-WT both with CCLE and Sanger drug sensitivity
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

par(mfrow=c(2,2))
plot(apply(G12C_SENS_CCLE_COR,2,median),col="gray60",pch=20, main="G12C ~ CCLE", ylab="Spearman rho", ylim=c(-.4,.4))
abline(h=0,lty=2,lwd=.5)
plot(apply(G12V_SENS_CCLE_COR,2,median),col="gray60",pch=20, main="G12V ~ CCLE", ylab="Spearman rho", ylim=c(-.4,.4))
abline(h=0,lty=2,lwd=.5)
plot(apply(G12C_SENS_SANGER_COR,2,median),col="gray60",pch=20, main="G12C ~ SANGER", ylab="Spearman rho", ylim=c(-.6,.6))
abline(h=0,lty=2,lwd=.5)
plot(apply(G12V_SENS_SANGER_COR,2,median),col="gray60",pch=20, main="G12V ~ SANGER", ylab="Spearman rho", ylim=c(-.6,.6))
abline(h=0,lty=2,lwd=.5)

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

########################################################################################################################
#  CCLE
########################################################################################################################
#  look at the differences in the correlations with drug sensitivity in CCLE between G12C and G12V
N_G12C_CCLE <- 79
N_G12V_CCLE <- 79
identical(colnames(G12C_SENS_CCLE_COR),colnames(G12V_SENS_CCLE_COR))
R_G12C <- apply(G12C_SENS_CCLE_COR,2,median)
R_G12V <- apply(G12V_SENS_CCLE_COR,2,median)

M <- R_G12C-R_G12V
p.val <- c()
for(i in c(1:length(colnames(G12C_SENS_CCLE_COR))))
{ p.val <- c(p.val,diff.corr(r1=R_G12C[i],r2=R_G12V[i],n1=N_G12C_CCLE,n2=N_G12V_CCLE))}


# plot the median spearman rho of G12C and G12V for the drugs that exhibit differential sensitivity between G12C and G12V
par(mfrow=c(2,2))
for(i in names(which(p.val<.05))){
boxplot(data.frame(cbind(G12C=G12C_SENS_CCLE_COR[,i],G12V=G12V_SENS_CCLE_COR[,i])), main=paste(i), ylab="Spearman rho",sub=paste("P=",format(p.val[i],digits=2)),ylim=c(-.43,.5))
stripchart(list(G12C=G12C_SENS_CCLE_COR[,i], G12V=G12V_SENS_CCLE_COR[,i]), add=T,vertical=TRUE,method="jitter",col="red",pch=20)
abline(h=0,lty=2,lwd=.5)
}


# plot the median spearman rho of G12C and G12V for the MEK inhibitors that exhibit differential sensitivity between G12C and G12V
drug.names <- ccle.drugs$Compound..code.or.generic.name.[grep(pattern="MEK",ccle.drugs$Target.s.)]
par(mfrow=c(2,2))
for(i in drug.names){
  boxplot(data.frame(cbind(G12C=G12C_SENS_CCLE_COR[,i],G12V=G12V_SENS_CCLE_COR[,i])), main=paste(i), ylab="Spearman rho",sub=paste("P=",format(p.val[i],digits=2)))
  stripchart(list(G12C=G12C_SENS_CCLE_COR[,i], G12V=G12V_SENS_CCLE_COR[,i]), add=T,vertical=TRUE,method="jitter",col="red",pch=20)
  abline(h=0,lty=2,lwd=.5)
}

# plot the median spearman rho of G12C and G12V for the MET inhibitors that exhibit differential sensitivity between G12C and G12V
drug.names <- ccle.drugs$Compound..code.or.generic.name.[grep(pattern="MET",ccle.drugs$Target.s.)]
drug.names[drug.names=="PF-2341066"] <- "PF2341066"
par(mfrow=c(2,2))
for(i in drug.names){
  boxplot(data.frame(cbind(G12C=G12C_SENS_CCLE_COR[,i],G12V=G12V_SENS_CCLE_COR[,i])), main=paste(i), ylab="Spearman rho",sub=paste("P=",format(p.val[i],digits=2)))
  stripchart(list(G12C=G12C_SENS_CCLE_COR[,i], G12V=G12V_SENS_CCLE_COR[,i]), add=T,vertical=TRUE,method="jitter",col="red",pch=20)
  abline(h=0,lty=2,lwd=.5)
}

# plot the median spearman rho of G12C and G12V for the RAF inhibitors that exhibit differential sensitivity between G12C and G12V
drug.names <- ccle.drugs$Compound..code.or.generic.name.[grep(pattern="RAF",ccle.drugs$Target.s.)]
drug.names <- c(drug.names, ccle.drugs$Compound..code.or.generic.name.[grep(pattern="Raf",ccle.drugs$Target.s.)])
par(mfrow=c(2,2))
for(i in drug.names){
  boxplot(data.frame(cbind(G12C=G12C_SENS_CCLE_COR[,i],G12V=G12V_SENS_CCLE_COR[,i])), main=paste(i), ylab="Spearman rho",sub=paste("P=",format(p.val[i],digits=2)))
  stripchart(list(G12C=G12C_SENS_CCLE_COR[,i], G12V=G12V_SENS_CCLE_COR[,i]), add=T,vertical=TRUE,method="jitter",col="red",pch=20)
  abline(h=0,lty=2,lwd=.5)
}

# plot the median spearman rho of G12C and G12V for the FGFR inhibitors that exhibit differential sensitivity between G12C and G12V
drug.names <- ccle.drugs$Compound..code.or.generic.name.[grep(pattern="FGF",ccle.drugs$Target.s.)]
par(mfrow=c(2,2))
for(i in drug.names){
  boxplot(data.frame(cbind(G12C=G12C_SENS_CCLE_COR[,i],G12V=G12V_SENS_CCLE_COR[,i])), main=paste(i), ylab="Spearman rho",sub=paste("P=",format(p.val[i],digits=2)))
  stripchart(list(G12C=G12C_SENS_CCLE_COR[,i], G12V=G12V_SENS_CCLE_COR[,i]), add=T,vertical=TRUE,method="jitter",col="red",pch=20)
  abline(h=0,lty=2,lwd=.5)
}

# plot the median spearman rho of G12C and G12V for the Gamma Secretase inhibitors that exhibit differential sensitivity between G12C and G12V
drug.names <- ccle.drugs$Compound..code.or.generic.name.[grep(pattern="ecretase",ccle.drugs$Target.s.)]
par(mfrow=c(2,2))
for(i in drug.names){
  boxplot(data.frame(cbind(G12C=G12C_SENS_CCLE_COR[,i],G12V=G12V_SENS_CCLE_COR[,i])), main=paste(i), ylab="Spearman rho",sub=paste("P=",format(p.val[i],digits=2)))
  stripchart(list(G12C=G12C_SENS_CCLE_COR[,i], G12V=G12V_SENS_CCLE_COR[,i]), add=T,vertical=TRUE,method="jitter",col="red",pch=20)
  abline(h=0,lty=2,lwd=.5)
}

# plot the median spearman rho of G12C and G12V for the HSP90 that exhibit differential sensitivity between G12C and G12V
drug.names <- ccle.drugs$Compound..code.or.generic.name.[grep(pattern="HSP90",ccle.drugs$Target.s.)]
par(mfrow=c(2,2))
for(i in drug.names){
  boxplot(data.frame(cbind(G12C=G12C_SENS_CCLE_COR[,i],G12V=G12V_SENS_CCLE_COR[,i])), main=paste(i), ylab="Spearman rho",sub=paste("P=",format(p.val[i],digits=2)))
  stripchart(list(G12C=G12C_SENS_CCLE_COR[,i], G12V=G12V_SENS_CCLE_COR[,i]), add=T,vertical=TRUE,method="jitter",col="red",pch=20)
  abline(h=0,lty=2,lwd=.5)
}

# plot the median spearman rho of G12C and G12V for the IGF1R that exhibit differential sensitivity between G12C and G12V
drug.names <- ccle.drugs$Compound..code.or.generic.name.[grep(pattern="IGF",ccle.drugs$Target.s.)]
par(mfrow=c(2,2))
for(i in drug.names){
  boxplot(data.frame(cbind(G12C=G12C_SENS_CCLE_COR[,i],G12V=G12V_SENS_CCLE_COR[,i])), main=paste(i), ylab="Spearman rho",sub=paste("P=",format(p.val[i],digits=2)))
  stripchart(list(G12C=G12C_SENS_CCLE_COR[,i], G12V=G12V_SENS_CCLE_COR[,i]), add=T,vertical=TRUE,method="jitter",col="red",pch=20)
  abline(h=0,lty=2,lwd=.5)
}

# plot the median spearman rho of G12C and G12V for the VEGF that exhibit differential sensitivity between G12C and G12V
drug.names <- ccle.drugs$Compound..code.or.generic.name.[grep(pattern="VEGF",ccle.drugs$Target.s.)]
drug.names <- drug.names[-1]
par(mfrow=c(2,2))
for(i in drug.names){
  boxplot(data.frame(cbind(G12C=G12C_SENS_CCLE_COR[,i],G12V=G12V_SENS_CCLE_COR[,i])), main=paste(i), ylab="Spearman rho",sub=paste("P=",format(p.val[i],digits=2)))
  stripchart(list(G12C=G12C_SENS_CCLE_COR[,i], G12V=G12V_SENS_CCLE_COR[,i]), add=T,vertical=TRUE,method="jitter",col="red",pch=20)
  abline(h=0,lty=2,lwd=.5)
}


########################################################################################################################
#  SANGER
########################################################################################################################

#  look at the differences in the correlations with drug sensitivity in SANGER between G12C and G12V
N_G12C_SANGER <- 67
N_G12V_SANGER <- 67
identical(colnames(G12C_SENS_SANGER_COR),colnames(G12V_SENS_SANGER_COR))
R_G12C <- apply(G12C_SENS_SANGER_COR,2,median)
R_G12V <- apply(G12V_SENS_SANGER_COR,2,median)


M <- R_G12C-R_G12V
p.val <- c()
for(i in c(1:length(colnames(G12C_SENS_SANGER_COR))))
{ p.val <- c(p.val,diff.corr(r1=R_G12C[i],r2=R_G12V[i],n1=N_G12C_SANGER,n2=N_G12V_SANGER))}
p.val


# plot the median spearman rho of G12C and G12V for the MEK inhibitors that exhibit differential sensitivity between G12C and G12V
drug.names <- sanger.drugs$DRUG.NAME[grep(pattern="MEK",sanger.drugs$TARGET)]
par(mfrow=c(2,2))
for(i in drug.names){
  boxplot(data.frame(cbind(G12C=G12C_SENS_SANGER_COR[,i],G12V=G12V_SENS_SANGER_COR[,i])), main=paste(i), ylab="Spearman rho",sub=paste("P=",format(p.val[i],digits=2)))
  stripchart(list(G12C=G12C_SENS_SANGER_COR[,i], G12V=G12V_SENS_SANGER_COR[,i]), add=T,vertical=TRUE,method="jitter",col="red",pch=20)
  abline(h=0,lty=2,lwd=.5)
}

# plot the median spearman rho of G12C and G12V for the drugs that exhibit differential sensitivity between G12C and G12V
par(mfrow=c(2,2))
for(i in names(which(p.val<.05))){
  boxplot(data.frame(cbind(G12C=G12C_SENS_SANGER_COR[,i],G12V=G12V_SENS_SANGER_COR[,i])), main=paste(i), ylab="Spearman rho",sub=paste("P=",format(p.val[i],digits=2)))
  stripchart(list(G12C=G12C_SENS_SANGER_COR[,i], G12V=G12V_SENS_SANGER_COR[,i]), add=T,vertical=TRUE,method="jitter",col="red",pch=20)
  abline(h=0,lty=2,lwd=.5)
}


# plot the median spearman rho of G12C and G12V for the MET inhibitors that exhibit differential sensitivity between G12C and G12V
drug.names <- sanger.drugs$DRUG.NAME[grep(pattern="MET",sanger.drugs$TARGET)]
par(mfrow=c(2,2))
for(i in drug.names){
  boxplot(data.frame(cbind(G12C=G12C_SENS_SANGER_COR[,i],G12V=G12V_SENS_SANGER_COR[,i])), main=paste(i), ylab="Spearman rho",sub=paste("P=",format(p.val[i],digits=2)))
  stripchart(list(G12C=G12C_SENS_SANGER_COR[,i], G12V=G12V_SENS_SANGER_COR[,i]), add=T,vertical=TRUE,method="jitter",col="red",pch=20)
  abline(h=0,lty=2,lwd=.5)
}

# plot the median spearman rho of G12C and G12V for the RAF inhibitors that exhibit differential sensitivity between G12C and G12V
drug.names <- sanger.drugs$DRUG.NAME[grep(pattern="RAF",sanger.drugs$TARGET)]
par(mfrow=c(2,2))
for(i in drug.names){
  boxplot(data.frame(cbind(G12C=G12C_SENS_SANGER_COR[,i],G12V=G12V_SENS_SANGER_COR[,i])), main=paste(i), ylab="Spearman rho",sub=paste("P=",format(p.val[i],digits=2)))
  stripchart(list(G12C=G12C_SENS_SANGER_COR[,i], G12V=G12V_SENS_SANGER_COR[,i]), add=T,vertical=TRUE,method="jitter",col="red",pch=20)
  abline(h=0,lty=2,lwd=.5)
}


# plot the median spearman rho of G12C and G12V for the FGFR inhibitors that exhibit differential sensitivity between G12C and G12V
drug.names <- sanger.drugs$DRUG.NAME[grep(pattern="FGF",sanger.drugs$TARGET)]
drug.names
par(mfrow=c(2,2))
for(i in drug.names){
  boxplot(data.frame(cbind(G12C=G12C_SENS_SANGER_COR[,i],G12V=G12V_SENS_SANGER_COR[,i])), main=paste(i), ylab="Spearman rho",sub=paste("P=",format(p.val[i],digits=2)))
  stripchart(list(G12C=G12C_SENS_SANGER_COR[,i], G12V=G12V_SENS_SANGER_COR[,i]), add=T,vertical=TRUE,method="jitter",col="red",pch=20)
  abline(h=0,lty=2,lwd=.5)
}

# plot the median spearman rho of G12C and G12V for the HSP90 inhibitors that exhibit differential sensitivity between G12C and G12V
drug.names <- sanger.drugs$DRUG.NAME[grep(pattern="HSP90",sanger.drugs$TARGET)]
drug.names <- c(drug.names,sanger.drugs$DRUG.NAME[grep(pattern="Hsp90",sanger.drugs$TARGET)])
drug.names[drug.names == "17.AAG"] <- "X17.AAG" 
par(mfrow=c(2,2))
for(i in drug.names){
  boxplot(data.frame(cbind(G12C=G12C_SENS_SANGER_COR[,i],G12V=G12V_SENS_SANGER_COR[,i])), main=paste(i), ylab="Spearman rho",sub=paste("P=",format(p.val[i],digits=2)))
  stripchart(list(G12C=G12C_SENS_SANGER_COR[,i], G12V=G12V_SENS_SANGER_COR[,i]), add=T,vertical=TRUE,method="jitter",col="red",pch=20)
  abline(h=0,lty=2,lwd=.5)
}

# plot the median spearman rho of G12C and G12V for the gamma secretase inhibitors that exhibit differential sensitivity between G12C and G12V
drug.names <- sanger.drugs$DRUG.NAME[grep(pattern="ecretase",sanger.drugs$TARGET)]
drug.names[drug.names=="Z.LLNle-CHO"] <- "Z.LLNle.CHO" 
drug.names <- drug.names[-1]
par(mfrow=c(1,1))
for(i in drug.names){
  boxplot(data.frame(cbind(G12C=G12C_SENS_SANGER_COR[,i],G12V=G12V_SENS_SANGER_COR[,i])), main=paste(i), ylab="Spearman rho",sub=paste("P=",format(p.val[i],digits=2)),ylim=c(-.6,.2))
  stripchart(list(G12C=G12C_SENS_SANGER_COR[,i], G12V=G12V_SENS_SANGER_COR[,i]), add=T,vertical=TRUE,method="jitter",col="red",pch=20)
  abline(h=0,lty=2,lwd=.5)
}

# plot the median spearman rho of G12C and G12V for the IGF-1R inhibitors that exhibit differential sensitivity between G12C and G12V
drug.names <- sanger.drugs$DRUG.NAME[grep(pattern="IGF",sanger.drugs$TARGET)]
drug.names
par(mfrow=c(2,2))
for(i in drug.names){
  boxplot(data.frame(cbind(G12C=G12C_SENS_SANGER_COR[,i],G12V=G12V_SENS_SANGER_COR[,i])), main=paste(i), ylab="Spearman rho",sub=paste("P=",format(p.val[i],digits=2)))
  stripchart(list(G12C=G12C_SENS_SANGER_COR[,i], G12V=G12V_SENS_SANGER_COR[,i]), add=T,vertical=TRUE,method="jitter",col="red",pch=20)
  abline(h=0,lty=2,lwd=.5)
}


# plot the median spearman rho of G12C and G12V for the VEGFR inhibitors that exhibit differential sensitivity between G12C and G12V
drug.names <- sanger.drugs$DRUG.NAME[grep(pattern="VEGF",sanger.drugs$TARGET)]
drug.names
par(mfrow=c(2,2))
for(i in drug.names){
  boxplot(data.frame(cbind(G12C=G12C_SENS_SANGER_COR[,i],G12V=G12V_SENS_SANGER_COR[,i])), main=paste(i), ylab="Spearman rho",sub=paste("P=",format(p.val[i],digits=2)))
  stripchart(list(G12C=G12C_SENS_SANGER_COR[,i], G12V=G12V_SENS_SANGER_COR[,i]), add=T,vertical=TRUE,method="jitter",col="red",pch=20)
  abline(h=0,lty=2,lwd=.5)
}