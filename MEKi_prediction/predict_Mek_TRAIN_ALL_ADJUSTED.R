# Charles Ferté
# Sage Bionetworks
# 14 Feb 2012


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
ccle_all <- loadEntity("syn1671195")
ccle_all <- ccle_all$objects$ccle_data
assign(x=names(ccle_all)[1],ccle_all[[1]])
assign(x=names(ccle_all)[2],ccle_all[[2]])
assign(x=names(ccle_all)[3],ccle_all[[3]])
assign(x=names(ccle_all)[4],ccle_all[[4]])
assign(x=names(ccle_all)[5],ccle_all[[5]])
assign(x=names(ccle_all)[6],ccle_all[[6]])
assign(x=names(ccle_drug)[1],ccle_drug[[1]])
assign(x=names(ccle_drug)[2],ccle_drug[[2]])

# let's use the ActArea as ccle_drug
ccle_drug <- ccle_drug_ActAreaNorm

#modify the ccle_mut into a binary matrix
ccle_mut[ccle_mut!="0"] <- 1
tmp <- rownames(ccle_mut)
ccle_mut <- apply(ccle_mut,2,as.numeric)
rownames(ccle_mut) <- tmp
rm(tmp)


#########################################################################################################################################
# make the data coherent between exp and cnv to ultimately compute the eigengenes vector
#########################################################################################################################################

tmp <- intersect(colnames(ccle_exp),ccle_info$CCLE.name)
tmp <- intersect(tmp,colnames(ccle_cnv))
ccle_exp <- ccle_exp[,tmp]
ccle_cnv <- ccle_cnv[,tmp]
ccle_info <- ccle_info[which(ccle_info$CCLE.name %in% intersect(tmp,ccle_info$CCLE.name)),]
rm(tmp)

global.matrix <- rbind(ccle_exp,ccle_cnv)
tissue.origin <- as.factor(ccle_info$Site.Primary)
s <- fast.svd(global.matrix-rowMeans(global.matrix))
par(mfrow=c(1,1))

# percentage variance explained
plot(s$d^2/sum(s$d^2),pch=20)

# plot first and second svd
plot(s$v[,1],s$v[,2],col=rainbow(24)[tissue.origin],pch=20,cex=1.5, xlab="PC1",ylab="PC2")
plot.new()
legend(0,1,legend=levels(tissue.origin),cex=.8,col=rainbow(24),pch=20,text.width=.4,text.font=2,
       text.col=rainbow(24))

# assign values for the eigengenes since they discriminate well the tissue specificity
eigengenes <- s$v
rownames(eigengenes) <- colnames(ccle_exp)
colnames(eigengenes) <- paste("PC_",seq(from=1,to=ncol(eigengenes),by=1),sep="")

#########################################################################################################
## feature selection for low variance
#########################################################################################################
tmp <- apply(ccle_exp,1,sd)
ccle_exp <- ccle_exp[which(tmp>quantile(x=tmp,probs=.25)),]
rm(tmp)

tmp <- apply(ccle_cnv,1,sd)
ccle_cnv <- ccle_cnv[which(tmp>quantile(x=tmp,probs=.25)),]
rm(tmp)

#########################################################################################################
## Make the data coherent between all the datasets
#########################################################################################################
tmp <- intersect(colnames(ccle_cnv),colnames(ccle_exp))
tmp <- intersect(colnames(ccle_mut),tmp)
tmp <- intersect(rownames(ccle_drug),tmp)
ccle_exp <- ccle_exp[,tmp]
ccle_cnv <- ccle_cnv[,tmp]
ccle_mut <- ccle_mut[,tmp]
ccle_drug <- ccle_drug[tmp,]
ccle_info <- ccle_info[which(ccle_info$CCLE.name %in% intersect(tmp,ccle_info$CCLE.name)),]
rm(tmp)

# identify the mek inhibitors
mek.inhib <-   ccle_drugs_info$Compound..code.or.generic.name.[grep(pattern="MEK",ccle_drugs_info$Target.s.)]

# identify the cells that have been consistently evaluated for with both MEK inhibitors
mek.cells <- rownames(ccle_drug)[-unique(c(which(is.na(ccle_drug[,mek.inhib[1]])), which(is.na(ccle_drug[,mek.inhib[2]]))))]

# identify the ActArea of the mek inhibs within the cells
mek.ActArea <- ccle_drug[mek.cells,mek.inhib]

# reduce the matrices to the mek.cells
ccle_exp <- ccle_exp[,mek.cells]
ccle_cnv <- ccle_cnv[,mek.cells]
ccle_mut <- ccle_mut[,mek.cells]
ccle_drug <- ccle_drug[mek.cells,]

# assess if there is any signal in the differential expression
par(mfrow=c(2,2),oma=c(0,0,6,0))
fit <- eBayes(lmFit(ccle_exp[,mek.cells],model.matrix(~mek.ActArea[,1])))
hist(fit$p.value[,2],breaks=50, main=paste("gene expr ~",colnames(mek.ActArea)[1]),col="aquamarine4",xlab="p values")
abline(v=.05,col="red",lty=2,lwd=2)
table(fit$p.value[,2]<.05)
fit <- eBayes(lmFit(ccle_cnv[,mek.cells],model.matrix(~mek.ActArea[,1])))
hist(fit$p.value[,2],breaks=50, main=paste(" cnv ~",colnames(mek.ActArea)[1]),col="aquamarine4",xlab="p values")
abline(v=.05,col="red",lty=2,lwd=2)
table(fit$p.value[,2]<.05)
fit <- eBayes(lmFit(ccle_exp[,mek.cells],model.matrix(~mek.ActArea[,2])))
hist(fit$p.value[,2],breaks=50, main=paste("gene expr ~",colnames(mek.ActArea)[2]),col="aquamarine4",xlab="p values")
abline(v=.05,col="red",lty=2,lwd=2)
table(fit$p.value[,2]<.05)
fit <- eBayes(lmFit(ccle_cnv[,mek.cells],model.matrix(~mek.ActArea[,2])))
hist(fit$p.value[,2],breaks=50, main=paste(" cnv ~",colnames(mek.ActArea)[2]),col="aquamarine4",xlab="p values")
abline(v=.05,col="red",lty=2,lwd=2)
table(fit$p.value[,2]<.05)

title(main=paste("univariate differential gene expression & differential CNV \nfor sensitivity to MEK inhibitors (AZ6244 and PD0325901) \nin the Cancer Cell Line Encyclopedia (n=",length(mek.cells),")",sep=""),outer=TRUE)

####################################################################################################
# define the NSCLC Breast Lung Melanoma Glioma & heMal (hematological malignacies) cells
####################################################################################################

carcinoma.mek.cells <-  intersect(mek.cells,ccle_info$CCLE.name[ccle_info$Histology =="carcinoma"])
nsclc.mek.cells <- carcinoma.mek.cells[grep(pattern="LUNG",x=carcinoma.mek.cells)]
nsclc.mek.cells <- intersect(nsclc.mek.cells,ccle_info$CCLE.name[ ccle_info$Hist.Subtype1 !="small_cell_carcinoma"])
crc.mek.cells <- carcinoma.mek.cells[grep(pattern="LARGE_INTESTINE",x=carcinoma.mek.cells)]
breast.mek.cells <- carcinoma.mek.cells[grep(pattern="BREAST",x=carcinoma.mek.cells)]
melanoma.mek.cells <-  intersect(mek.cells,ccle_info$CCLE.name[ccle_info$Histology =="malignant_melanoma"])
glioma.mek.cells <-  intersect(mek.cells,ccle_info$CCLE.name[ccle_info$Histology =="glioma"])
hemal.mek.cells <- intersect(mek.cells,ccle_info$CCLE.name[ccle_info$Histology %in% c("haematopoietic_neoplasm","lymphoid_neoplasm")])

#############################
# predictive modeling
############################

# first define the global matrix with the eigen genes
global.matrix <- rbind(ccle_exp,ccle_cnv,ccle_mut)

# add non penalization on the eigengenes 
eigengenes <- t(eigengenes[colnames(global.matrix),c(1:10)])

global.matrix <- rbind(global.matrix,eigengenes)

rownames(global.matrix) <- c(paste(rownames(ccle_exp),"_exp",sep=""),paste(rownames(ccle_cnv),"_cnv",sep=""),paste(rownames(ccle_mut),"_mut",sep=""),paste("PC",c(1:10),sep=""))


N=100
models <- 0
yhat.all <- c()
yhat.breast <- c()
yhat.nsclc <- c()
yhat.crc <- c()
yhat.glioma <- c()
yhat.melanoma <- c()
yhat.hemal  <- c()
selected <- c()
i <- 0
while(models<N)
{
  par(mfrow=c(1,1))
  train <- sample(mek.cells,replace=TRUE)
  val <-mek.cells[-which(mek.cells %in% train)]
  vec.train <-apply(ccle_drug[train,mek.inhib],1,mean)
  cv.fit <- cv.glmnet(t(global.matrix[,train]), y=vec.train,nfolds=3, alpha=.1)
  fit <- glmnet(x=t(global.matrix[,train]),y=vec.train,alpha=.1,lambda=cv.fit$lambda.1se)
  selected <- c(selected,list(fit$beta))
  yhat.all <- c(yhat.all,list(predict(fit, t(global.matrix[,val]))))
  yhat.breast <- c(yhat.breast,list(predict(fit,t(global.matrix[,breast.mek.cells[-which(breast.mek.cells %in% train)]]))))
  yhat.nsclc <- c(yhat.nsclc,list(predict(fit,t(global.matrix[,nsclc.mek.cells[-which(nsclc.mek.cells %in% train)]]))))
  yhat.crc <- c(yhat.crc,list(predict(fit,t(global.matrix[,crc.mek.cells[-which(crc.mek.cells %in% train)]]))))
  yhat.glioma <- c(yhat.glioma,list(predict(fit,t(global.matrix[,glioma.mek.cells[-which(glioma.mek.cells %in% train)]]))))
  yhat.melanoma <- c(yhat.melanoma,list(predict(fit,t(global.matrix[,melanoma.mek.cells[-which(melanoma.mek.cells %in% train)]]))))
  yhat.hemal <- c(yhat.hemal,list(predict(fit,t(global.matrix[,hemal.mek.cells[-which(hemal.mek.cells %in% train)]]))))
 
  i=1+i
  print(i)
  models <- length(yhat.all)
  }  

## see what is retained in the models
abc <- matrix(NA,ncol=length(selected),nrow=nrow(global.matrix))
rownames(abc) <- rownames(global.matrix)
for(i in c(1:length(selected)))
{  abc[,i]<- as.numeric(selected[[i]])  
}

par(mfrow=c(1,1))
abc[abc!=0] <-1 
hist(log10(rowSums(abs(abc))),breaks=50,col="red")
h <- rowSums(abs(abc))
sort(h,decreasing=TRUE)[1:100]
which(h>quantile(h,probs=.99))
###### SPEARMAN ###########
par(mfrow=c(2,4),oma=c(0,0,6,0))
method.cor <- "spearman"
source("/home/cferte/FELLOW/cferte/KRAS_Analysis/MEKi_prediction/Mek_performance.R")
title("performance of 100 bootstrapped models\npredicting sensitivity to MEK inhibitors (as assessed by ActArea)\nTraining in All cell lines, ajusted on the 1st Princ. Component",outer=TRUE)

###### PEARSON ###########
par(mfrow=c(2,4),oma=c(0,0,6,0))
method.cor <- "pearson"
source("/home/cferte/FELLOW/cferte/KRAS_Analysis/MEKi_prediction/Mek_performance.R")
title("performance of 100 bootstrapped models\npredicting sensitivity to MEK inhibitors (as assessed by ActArea)\nTraining in All cell lines,ajusted on the 1st Princ. Component",outer=TRUE)
#

