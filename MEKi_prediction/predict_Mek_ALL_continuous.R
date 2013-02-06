# Charles Ferté
# Sage Bionetworks
# 14 Sept 2012


# train a model of MEK sensibility (independant response = bionomial)

#load the different packages
options(stringsAsFactors=FALSE)

library(affy)
library(corpcor)
library(lattice)
library(limma)
library(caret)
library(glmnet)
library(snm)
library(synapseClient)
synapseLogin("charles.ferte@sagebase.org","charles")

###############################################################
# load the CCLE data (snm normalized) 
###############################################################
ccle_all <- loadEntity("syn1671180")
ccle_all <- ccle_all$objects$ccle_data
assign(x=names(ccle_all)[1],ccle_all[[1]])
assign(x=names(ccle_all)[2],ccle_all[[2]])
assign(x=names(ccle_all)[3],ccle_all[[3]])
assign(x=names(ccle_all)[4],ccle_all[[4]])
assign(x=names(ccle_all)[5],ccle_all[[5]])
assign(x=names(ccle_all)[6],ccle_all[[6]])

# make the data coherent between all the datasets
tmp <- intersect(colnames(ccle_cnv),colnames(ccle_exp))
tmp <- intersect(colnames(ccle_mut),tmp)
tmp <- intersect(rownames(ccle_drug),tmp)
ccle_exp <- ccle_exp[,tmp]
ccle_cnv <- ccle_cnv[,tmp]
ccle_mut <- ccle_mut[,tmp]
ccle_drug <- ccle_drug[tmp,]
ccle_info <- ccle_info[which(ccle_info$CCLE.name %in% intersect(tmp,ccle_info$CCLE.name)),]
rm(tmp)

#modify the ccle_mut into a binary matrix
ccle_mut[ccle_mut!="0"] <- 1
ccle_mut <- apply(ccle_mut,2,as.numeric)

# identify the mek inhibitors
mek.inhib <-   ccle_drugs_info$Compound..code.or.generic.name.[grep(pattern="MEK",ccle_drugs_info$Target.s.)]
# identify the cells treated consistently with both MEK inhibitors
mek.cells <- rownames(ccle_drug)[-unique(c(which(is.na(ccle_drug[,mek.inhib[1]])), which(is.na(ccle_drug[,mek.inhib[2]]))))]
# identify the ic50 of the mek inhibs within the cells
mek.ic50 <- ccle_drug[mek.cells,mek.inhib]
  
# assess if there is any signal in the differential expression
par(mfrow=c(2,3))
fit <- eBayes(lmFit(ccle_exp[,mek.cells],model.matrix(~mek.ic50[,1])))
hist(fit$p.value[,2],breaks=100, main=paste("ccle expr ~",colnames(mek.ic50)[1]))
table(fit$p.value[,2]<.05)
fit <- eBayes(lmFit(ccle_cnv[,mek.cells],model.matrix(~mek.ic50[,1])))
hist(fit$p.value[,2],breaks=100, main=paste("ccle cnv ~",colnames(mek.ic50)[1]))
table(fit$p.value[,2]<.05)
fit <- eBayes(lmFit(ccle_mut[,mek.cells],model.matrix(~mek.ic50[,1])))
hist(fit$p.value[,2],breaks=100, main=paste("ccle mut ~",colnames(mek.ic50)[1]))
table(fit$p.value[,2]<.05)
fit <- eBayes(lmFit(ccle_exp[,mek.cells],model.matrix(~mek.ic50[,2])))
hist(fit$p.value[,2],breaks=100, main=paste("ccle expr ~",colnames(mek.ic50)[2]))
table(fit$p.value[,2]<.05)
fit <- eBayes(lmFit(ccle_cnv[,mek.cells],model.matrix(~mek.ic50[,2])))
hist(fit$p.value[,2],breaks=100, main=paste("ccle cnv ~",colnames(mek.ic50)[2]))
table(fit$p.value[,2]<.05)
fit <- eBayes(lmFit(ccle_mut[,mek.cells],model.matrix(~mek.ic50[,2])))
hist(fit$p.value[,2],breaks=100, main=paste("ccle mut ~",colnames(mek.ic50)[2]))
table(fit$p.value[,2]<.05)



title(main="univariate differential expression & cnv for sensitivity to MEKi in CCLE",outer=TRUE)

################################
# restrict to KRAS mutant samples only 
################################

KC <- names(which(kras.info.ccle["KRAS0",]==1))
KS <- names(which(kras.info.sanger["KRAS0",]==1))

################################################################################################################
# is there any signal in the differential expression
################################################################################################################

par(mfrow=c(2,3))
fit <- eBayes(lmFit(SANGER_EXP[,rownames(mek.ic50.sanger)],model.matrix(~mek.ic50.sanger[,1])))
hist(fit$p.value[,2],breaks=100, main=paste("Sanger Expr ~",colnames(mek.ic50.sanger)[1]))
table(fit$p.value[,2]<.05)
fit <- eBayes(lmFit(SANGER_EXP[,rownames(mek.ic50.sanger)],model.matrix(~mek.ic50.sanger[,2])))
hist(fit$p.value[,2],breaks=100, main=paste("Sanger Expr ~",colnames(mek.ic50.sanger)[2]))
table(fit$p.value[,2]<.05)
fit <- eBayes(lmFit(SANGER_EXP[,rownames(mek.ic50.sanger)],model.matrix(~mek.ic50.sanger[,3])))
hist(fit$p.value[,2],breaks=100, main=paste("Sanger Expr ~",colnames(mek.ic50.sanger)[3]))
table(fit$p.value[,2]<.05)
fit <- eBayes(lmFit(SANGER_EXP[,rownames(mek.ic50.sanger)],model.matrix(~mek.ic50.sanger[,4])))
hist(fit$p.value[,2],breaks=100, main=paste("Sanger Expr ~",colnames(mek.ic50.sanger)[4]))
table(fit$p.value[,2]<.05)
fit <- eBayes(lmFit(CCLE_EXP[,rownames(mek.ic50.ccle)],model.matrix(~mek.ic50.ccle[,1])))
hist(fit$p.value[,2],breaks=100, main=paste("CCLE Expr ~",colnames(mek.ic50.ccle)[1]))
table(fit$p.value[,2]<.05)
fit <- eBayes(lmFit(CCLE_EXP[,rownames(mek.ic50.ccle)],model.matrix(~mek.ic50.ccle[,2])))
hist(fit$p.value[,2],breaks=100, main=paste("CCLE Expr ~",colnames(mek.ic50.ccle)[2]))
table(fit$p.value[,2]<.05)
title(main="univariate difefrential expression for sensitivity to MEKi across SANGER & CCLE",outer=TRUE)



################################################################################################################
# is there any signal in the mutation
################################################################################################################

par(mfrow=c(2,3))
fit <- eBayes(lmFit(mutations.sanger[,rownames(mek.ic50.sanger)],model.matrix(~mek.ic50.sanger[,1])))
hist(fit$p.value[,2], main=paste("Sanger Expr ~",colnames(mek.ic50.sanger)[1]))
table(fit$p.value[,2]<.05)
fit <- eBayes(lmFit(mutations.sanger[,rownames(mek.ic50.sanger)],model.matrix(~mek.ic50.sanger[,2])))
hist(fit$p.value[,2], main=paste("Sanger Expr ~",colnames(mek.ic50.sanger)[2]))
table(fit$p.value[,2]<.05)
fit <- eBayes(lmFit(mutations.sanger[,rownames(mek.ic50.sanger)],model.matrix(~mek.ic50.sanger[,3])))
hist(fit$p.value[,2], main=paste("Sanger Expr ~",colnames(mek.ic50.sanger)[3]))
table(fit$p.value[,2]<.05)
fit <- eBayes(lmFit(mutations.sanger[,rownames(mek.ic50.sanger)],model.matrix(~mek.ic50.sanger[,4])))
hist(fit$p.value[,2], main=paste("Sanger Expr ~",colnames(mek.ic50.sanger)[4]))
table(fit$p.value[,2]<.05)
fit <- eBayes(lmFit(mutations.ccle[,rownames(mek.ic50.ccle)],model.matrix(~mek.ic50.ccle[,1])))
hist(fit$p.value[,2], main=paste("CCLE Expr ~",colnames(mek.ic50.ccle)[1]))
table(fit$p.value[,2]<.05)
fit <- eBayes(lmFit(mutations.ccle[,rownames(mek.ic50.ccle)],model.matrix(~mek.ic50.ccle[,2])))
hist(fit$p.value[,2], main=paste("CCLE Expr ~",colnames(mek.ic50.ccle)[2]))
table(fit$p.value[,2]<.05)
title(main="univariate mutational association for sensitivity to MEKi across SANGER & CCLE",outer=TRUE)


###################################################################################################################
# GOLD standard: correlations between IC50 of CCLE and Sanger
###################################################################################################################
tmp <- intersect(rownames(mek.ic50.ccle),rownames(mek.ic50.sanger))
cor(mek.ic50.ccle[tmp,],mek.ic50.sanger[tmp,],method="spearman")
scatterplotMatrix(normalizeCyclicLoess(cbind(mek.ic50.ccle[tmp,],mek.ic50.sanger[tmp,])))

ic50.train <- mek.ic50.sanger
ic50.val <- mek.ic50.ccle

# is there any signal in the differential expression
par(mfrow=c(2,3))
fit <- eBayes(lmFit(SANGER_EXP[,rownames(ic50.train)],model.matrix(~ic50.train[,1])))
hist(fit$p.value[,2],breaks=100, main=paste("ccle Expr ~",colnames(ic50.train)[1]))
table(fit$p.value[,2]<.05)
fit <- eBayes(lmFit(SANGER_EXP[,rownames(ic50.train)],model.matrix(~ic50.train[,2])))
hist(fit$p.value[,2],breaks=100, main=paste("ccle Expr ~",colnames(ic50.train)[2]))
table(fit$p.value[,2]<.05)
fit <- eBayes(lmFit(SANGER_EXP[,rownames(ic50.train)],model.matrix(~ic50.train[,3])))
hist(fit$p.value[,2],breaks=100, main=paste("ccle Expr ~",colnames(ic50.train)[3]))
table(fit$p.value[,2]<.05)
fit <- eBayes(lmFit(SANGER_EXP[,rownames(ic50.train)],model.matrix(~ic50.train[,4])))
hist(fit$p.value[,2],breaks=100, main=paste("ccle Expr ~",colnames(ic50.train)[4]))
table(fit$p.value[,2]<.05)
fit <- eBayes(lmFit(CCLE_EXP[,rownames(ic50.val)],model.matrix(~ic50.val[,1])))
hist(fit$p.value[,2],breaks=100, main=paste("Sanger Expr ~",colnames(ic50.val)[1]))
table(fit$p.value[,2]<.05)
fit <- eBayes(lmFit(CCLE_EXP[,rownames(ic50.val)],model.matrix(~ic50.val[,2])))
hist(fit$p.value[,2],breaks=100, main=paste("Sanger Expr ~",colnames(ic50.val)[2]))
table(fit$p.value[,2]<.05)
title(main="univariate difefrential expression for sensitivity to MEKi across SANGER & CCLE",outer=TRUE)

###################################################################################################################
# train our predictive model of MEK response in sanger
# using a penalized regression approach  
# with alpha=.1 (more ridge) and determine lambda using nfolds= 5
# the robustness of the model is increased by boostrapping (n=100)
# validate each model in ccle
###################################################################################################################
ic50.train <- mek.ic50.sanger
ic50.val <- mek.ic50.ccle
par(mfrow=c(1,2))
require(glmnet)
N <- 100
fit <- c()
selected <- c()
yhat <- c()
models <- 0
i <- 0
while(models<N)
{
  dev <- c()
  val <- c()
  j <- sample(rownames(ic50.train),replace=TRUE)
  table(j)
  #trainex <- SANGER_EXP[,j]
  trainex <- rbind(rbind(SANGER_EXP[,j],mutations.sanger[,j]),kras.info.sanger[,j])
  vec.train <- apply(ic50.train[j,c(3,4)],1,mean)
  cv.fit <- cv.glmnet(t(trainex), y=vec.train,nfolds=5, alpha=.5)
  fit <- glmnet(x=t(trainex),y=vec.train,alpha=.5,lambda=cv.fit$lambda.1se)
  if(length(which(abs(as.numeric(fit$beta))> 10^-5))>10)
  {
    i=i+1
    print(i)
    selected <- cbind(selected , as.numeric(fit$beta))
    
    dev <- which(rownames(ic50.val) %in% intersect(j,rownames(ic50.val)))
    val <- rownames(ic50.val)
    #val <- rownames(ic50.val)[-dev]
    #intersect(j,val)
    #validex <- CCLE_EXP[,val]
    validex <- rbind(rbind(CCLE_EXP[,val],mutations.ccle[,val]),kras.info.ccle[,val])
    yhat <- c(yhat,list(predict(fit, t(validex))))
    models <- length(yhat)
  } }

rownames(selected) <- rownames(trainex)
selected1 <- selected
y <- c()
for(i in c(1:length(yhat)))
{
  y <- cbind(y, yhat[[i]])
}
rownames(y) <- val

boxplot(list(PD0325901=cor(y,ic50.val[,1],method="spearman"), AZD6244=cor(y,ic50.val[,2],method="spearman")),outline=FALSE,ylim=c(0,1),cex.axis=.7)
stripchart(list(PD0325901=cor(y,ic50.val[,1],method="spearman"), AZD6244=cor(y,ic50.val[,2],method="spearman")),method="jitter",vertical=TRUE,add=TRUE,col="royalblue",pch=20)
tmp <- intersect(rownames(ic50.train),rownames(ic50.val))
a <- cor(ic50.val[tmp,],apply(ic50.train[tmp,],1,mean),method="spearman")[1]
b <- cor(ic50.val[tmp,],apply(ic50.train[tmp,],1,mean),method="spearman")[2]
stripchart(list(PD0325901=a, AZD6244=b),pch=20,cex=3,col="orange",add=TRUE,vertical=TRUE)
abline(h=c(0,.2,.4,.6,.8,1),lty=2)


###################################################################################################################
# train our predictive model of MEK response in ccle
# using a penalized regression approach  
# with alpha=.1 (more ridge) and determine lambda using nfolds= 5
# the robustness of the model is increased by boostrapping (n=100)
# validate each model in sanger
###################################################################################################################
ic50.train <- mek.ic50.ccle
ic50.val <- mek.ic50.sanger
require(glmnet)
N <- 100
fit <- c()
selected <- c()
yhat <- c()
models <- 0
i <- 0
while(models<N)
{
  dev <- c()
  val <- c()
  j <- sample(rownames(ic50.train),replace=TRUE)
  table(j)
  trainex <- rbind(rbind(CCLE_EXP[,j],mutations.ccle[,j]),kras.info.ccle[,j])
  vec.train <- apply(ic50.train[j,],1,mean)
  cv.fit <- cv.glmnet(t(trainex), y=vec.train,nfolds=5, alpha=.5)
  fit <- glmnet(x=t(trainex),y=vec.train,alpha=.5,lambda=cv.fit$lambda.1se)
  if(length(which(abs(as.numeric(fit$beta))> 10^-5))>10)
  {
    i=i+1
    print(i)
    selected <- cbind(selected , as.numeric(fit$beta))
    
    dev <- which(rownames(ic50.val) %in% intersect(j,rownames(ic50.val)))
    val <- rownames(ic50.val)
    #val <- rownames(ic50.val)[-dev]
    #intersect(j,val)
    #validex <- CCLE_EXP[,val]
    validex <- rbind(rbind(SANGER_EXP[,val],mutations.sanger[,val]),kras.info.sanger[,val])
    yhat <- c(yhat,list(predict(fit, t(validex))))
    models <- length(yhat)
  } }

rownames(selected) <- rownames(trainex)
selected2 <- selected
y <- c()
for(i in c(1:length(yhat)))
{
  y <- cbind(y, yhat[[i]])
}
rownames(y) <- val

boxplot(list(RDEA119=cor(y,ic50.val[,1],method="spearman"),CI.1040=cor(y,ic50.val[,2],method="spearman"), PD0325901=cor(y,ic50.val[,3],method="spearman"), AZD6244=cor(y,ic50.val[,4],method="spearman")),outline=FALSE,ylim=c(0,1),cex.axis=.7)
stripchart(list(RDEA119=cor(y,ic50.val[,1],method="spearman"),CI.1040=cor(y,ic50.val[,2],method="spearman"), PD0325901=cor(y,ic50.val[,3],method="spearman"), AZD6244=cor(y,ic50.val[,4],method="spearman")),method="jitter",vertical=TRUE,add=TRUE,col="royalblue",pch=20)
tmp <- intersect(rownames(ic50.train),rownames(ic50.val))
a <- cor(ic50.val[tmp,],apply(ic50.train[tmp,],1,mean),method="spearman")[1]
b <- cor(ic50.val[tmp,],apply(ic50.train[tmp,],1,mean),method="spearman")[2]
c <- cor(ic50.val[tmp,],apply(ic50.train[tmp,],1,mean),method="spearman")[3]
d <- cor(ic50.val[tmp,],apply(ic50.train[tmp,],1,mean),method="spearman")[4]
stripchart(list(RDEA119=a,CI.1040=b,PD0325901=c, AZD6244=d),pch=20,cex=3,col="orange",add=TRUE,vertical=TRUE)
abline(h=c(0,.2,.4,.6,.8,1),lty=2)

#######################################################################

# extract the biological meaning
mus1 <- selected1
mus1[mus1!=0] <- 1
x <- rownames(selected1)[which(apply(mus1,1,sum)>quantile(apply(mus1,1,sum),probs=.9))]
sort(apply(mus1,1,sum),decreasing=TRUE)[1:50]

mus2 <- selected2
mus2[mus2!=0] <- 1
y <- rownames(selected2)[which(apply(mus2,1,sum)>quantile(apply(mus2,1,sum),probs=.9))]
sort(apply(mus2,1,sum),decreasing=TRUE)[1:50]

intersect(x,y)
