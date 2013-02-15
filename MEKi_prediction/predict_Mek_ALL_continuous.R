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

# first define the global matrix
global.matrix <- rbind(ccle_exp,ccle_cnv,ccle_mut)
#eigengenes <- t(eigengenes[colnames(global.matrix),c(1:10)])
#global.matrix <- rbind(global.matrix,eigengenes)
rownames(global.matrix) <- c(paste(rownames(ccle_exp),"_exp",sep=""),paste(rownames(ccle_cnv),"_cnv",sep=""),paste(rownames(ccle_mut),"_mut",sep=""))

N=50
#selected <- c()
models <- 0


yhat.all <- c()
yhat.breast <- c()
yhat.nsclc <- c()
yhat.crc <- c()
yhat.glioma <- c()
yhat.melanoma <- c()
yhat.hemal  <- c()

i <- 0
while(models<N)
{
  par(mfrow=c(1,1))
  train <- sample(mek.cells,replace=TRUE)
  val <-mek.cells[-which(mek.cells %in% train)]

#pen <- c(rep(1,times=dim(trainex)[1]-1),0)
  vec.train <-apply(ccle_drug[train,mek.inhib],1,mean)
  cv.fit <- cv.glmnet(t(global.matrix[,train]), y=vec.train,nfolds=3, alpha=.1)
  fit <- glmnet(x=t(global.matrix[,train]),y=vec.train,alpha=.1,lambda=cv.fit$lambda.1se)

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






  #selected <- cbind(selected,as.numeric(fit$beta))
All1 <- c(All1,cor(yhat.all,mek.ActArea[rownames(yhat.all),1],method="pearson",use="pairwise.complete.obs"))
All2 <- c(All2,cor(yhat.all,mek.ActArea[rownames(yhat.all),2],method="pearson",use="pairwise.complete.obs"))

  breast1 <- c(breast1,cor(yhat.breast,mek.ActArea[rownames(yhat.breast),1],method="pearson",use="pairwise.complete.obs"))
  breast2 <- c(breast2,cor(yhat.breast,mek.ActArea[rownames(yhat.breast),2],method="pearson",use="pairwise.complete.obs"))
  
  nsclc1 <- c(nsclc1,cor(yhat.nsclc,mek.ActArea[rownames(yhat.nsclc),1],method="pearson",use="pairwise.complete.obs"))
  nsclc2 <- c(nsclc2,cor(yhat.nsclc,mek.ActArea[rownames(yhat.nsclc),2],method="pearson",use="pairwise.complete.obs"))
  
  crc1 <- c(crc1,cor(yhat.crc,mek.ActArea[rownames(yhat.crc),1],method="pearson",use="pairwise.complete.obs"))
  crc2 <- c(crc2,cor(yhat.crc,mek.ActArea[rownames(yhat.crc),2],method="pearson",use="pairwise.complete.obs"))
  
  glioma1 <- c(glioma1,cor(yhat.glioma,mek.ActArea[rownames(yhat.glioma),1],method="pearson",use="pairwise.complete.obs"))
  glioma2 <- c(glioma2,cor(yhat.glioma,mek.ActArea[rownames(yhat.glioma),2],method="pearson",use="pairwise.complete.obs"))
  
  hemal1 <- c(hemal1,cor(yhat.hemal,mek.ActArea[rownames(yhat.hemal),1],method="pearson",use="pairwise.complete.obs"))
  hemal2 <- c(hemal2,cor(yhat.hemal,mek.ActArea[rownames(yhat.hemal),2],method="pearson",use="pairwise.complete.obs"))
  
  melanoma1 <- c(melanoma1,cor(yhat.melanoma,mek.ActArea[rownames(yhat.melanoma),1],method="pearson",use="pairwise.complete.obs"))
  melanoma2 <- c(melanoma2,cor(yhat.melanoma,mek.ActArea[rownames(yhat.melanoma),2],method="pearson",use="pairwise.complete.obs"))
  

#}


All1 <- c()
All2 <- c()

breast1 <- c()
nsclc1 <- c()
crc1 <- c()
glioma1 <- c()
hemal1 <- c()
melanoma1 <- c()

breast2 <- c()
nsclc2 <- c()
crc2 <- c()
glioma2 <- c()
hemal2 <- c()
melanoma2 <- c()

# display TALL VALL
boxplot(list(PD0325901=All1,AZD6244=All2), 
        main="Performance of  bootstrapped models predicting \nthe sensitivity to Mek inhibitors (assessed by ActArea) \ntraining in all cell lines, validated in all cell lines not used for training",
        ylab="Pearson r",ylim=c(0,1),outline=FALSE)
stripchart(list(PD0325901=All1,AZD6244=All2),vertical=TRUE,method="jitter",add=TRUE,col="aquamarine4",pch=20)
abline(h=seq(0,1,.1),lty=2)

# display TALL Breast
boxplot(list(PD0325901=All1,AZD6244=All2), 
        main="Performance of  bootstrapped models predicting \nthe sensitivity to Mek inhibitors (assessed by ActArea) \ntraining in all cell lines, validated in all cell lines not used for training",
        ylab="Pearson r",ylim=c(0,1),outline=FALSE)
stripchart(list(PD0325901=All1,AZD6244=All2),vertical=TRUE,method="jitter",add=TRUE,col="aquamarine4",pch=20)
abline(h=seq(0,1,.1),lty=2)


###################################################################################################################
# train our predictive model of MEK response in the carcinoma ccle but half of the lung
# using a penalized regression approach  
# with alpha=.1 (more ridge) and determine lambda using nfolds= 5
# the robustness of the model is increased by boostrapping (n=100)
# validate each model in the rest of the lung
###################################################################################################################

# restrict to the carcinoma only
carcinoma.mek.cells <-  intersect(mek.cells,ccle_info$CCLE.name[ccle_info$Histology =="carcinoma"])
lung.mek.cells <- carcinoma.mek.cells[grep(pattern="LUNG",x=carcinoma.mek.cells)]
nsclc.mek.cells <- intersect(lung.mek.cells,ccle_info$CCLE.name[ ccle_info$Hist.Subtype1 !="small_cell_carcinoma"])


par(mfrow=c(1,1))
require(glmnet)
N <- 50
fit <- c()
selected <- c()
yhat <- c()
models <- 0
i <- 0
while(models<N)
{
  
  val <- sample(x=nsclc.mek.cells,size=length(nsclc.mek.cells)/2,replace=FALSE)
  train <- carcinoma.mek.cells[-which(carcinoma.mek.cells %in% val)]
  trainex <- rbind(ccle_exp[,train],ccle_cnv[,train])
  vec.train <- apply(ccle_drug[train,mek.inhib],1,mean)
  cv.fit <- cv.glmnet(t(trainex), y=vec.train,nfolds=5, alpha=.6)
  fit <- glmnet(x=t(trainex),y=vec.train,alpha=.6,lambda=cv.fit$lambda.1se)
  if(length(which(abs(as.numeric(fit$beta))> 10^-5))>10)
  {
    i=i+1
    print(i)
    selected <- cbind(selected , as.numeric(fit$beta))
    validex <- rbind(ccle_exp[,val],ccle_cnv[,val])
    yhat <- c(yhat,list(predict(fit, t(validex))))
    models <- length(yhat)
  } }

rownames(selected) <- rownames(trainex)
selected1 <- selected

y <- c()
for(i in c(1:length(yhat)))
{  y <- unique(c(y, rownames(yhat[[i]])))}

Y <- matrix(NA,nrow=length(unique(y)),ncol=length(yhat))
rownames(Y) <- y
colnames(Y) <- c(1:length(yhat))
for(i in c(1:length(yhat))){
  Y[rownames(yhat[[i]]),i] <- yhat[[i]]
}

par(mfrow=c(1,1))
boxplot(list(PD0325901=cor(Y,mek.ActArea[rownames(Y),1],method="spearman",use="pairwise.complete.obs"), AZD6244=cor(Y,mek.ActArea[rownames(Y),2],method="spearman",use="pairwise.complete.obs")),outline=FALSE,ylim=c(0,1),cex.axis=.7)
stripchart(list(PD0325901=cor(Y,mek.ActArea[rownames(Y),1],method="spearman",use="pairwise.complete.obs"), AZD6244=cor(Y,mek.ActArea[rownames(Y),2],method="spearman",use="pairwise.complete.obs")),method="jitter",vertical=TRUE,add=TRUE,col="aquamarine5",pch=20)
abline(h=c(0,.2,.4,.6,.8,1),lty=2)

# #######################################################################
# # extract the biological meaning
# #######################################################################
# 
# mus1 <- selected1
# mus1[mus1!=0] <- 1
# x <- rownames(selected1)[which(apply(mus1,1,sum)>quantile(apply(mus1,1,sum),probs=.9))]
# sort(apply(mus1,1,sum),decreasing=TRUE)[1:50]
# 
# mus2 <- selected2
# mus2[mus2!=0] <- 1
# y <- rownames(selected2)[which(apply(mus2,1,sum)>quantile(apply(mus2,1,sum),probs=.9))]
# sort(apply(mus2,1,sum),decreasing=TRUE)[1:50]
# 
# intersect(x,y)
