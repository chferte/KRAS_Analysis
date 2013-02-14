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
heMal.mek.cells <- intersect(mek.cells,ccle_info$CCLE.name[ccle_info$Histology %in% c("haematopoietic_neoplasm","lymphoid_neoplasm")])

#############################
# predictive modeling
############################

# first define the global matrix
global.matrix <- rbind(ccle_exp,ccle_cnv,ccle_mut)
#eigengenes <- t(eigengenes[colnames(global.matrix),c(1:10)])
#global.matrix <- rbind(global.matrix,eigengenes)
rownames(global.matrix) <- c(paste(rownames(ccle_exp),"_exp",sep=""),paste(rownames(ccle_cnv),"_cnv",sep=""),paste(rownames(ccle_mut),"_mut",sep=""))

selected <- c()
k <- c()
j <- c()
i <- 0
for(i in c(1:50))
{
  par(mfrow=c(1,1))
yhat <- c()
train <- sample(mek.cells,replace=TRUE)
  val <-mek.cells[-which(mek.cells %in% train)]
trainex <- global.matrix[,train]
#pen <- c(rep(1,times=dim(trainex)[1]-1),0)
vec.train <-apply(ccle_drug[train,mek.inhib],1,mean)
cv.fit <- cv.glmnet(t(trainex), y=vec.train,nfolds=3, alpha=.1)
fit <- glmnet(x=t(trainex),y=vec.train,alpha=1,lambda=cv.fit$lambda.1se)
validex <- global.matrix[,val]
  rownames(validex) <- rownames(trainex)
yhat <- predict(fit, t(validex))
  selected <- cbind(selected,as.numeric(fit$beta))
j <- c(j,cor(yhat,mek.ActArea[rownames(yhat),1],method="spearman",use="pairwise.complete.obs"))
k <- c(k,cor(yhat,mek.ActArea[rownames(yhat),2],method="spearman",use="pairwise.complete.obs"))
print(i)
}
boxplot(list(PD0325901=j,AZD6244=k), main="Correlation between models and true ActArea of Mek inhibitors \ training in all cell lines (n=444) and validated in the cell lines not used for training",ylab="Spearman rho",ylim=c(0,1))
stripchart(list(PD0325901=j,AZD6244=k),vertical=TRUE,method="jitter",add=TRUE,col="aquamarine4",pch=20)




rownames(selected) <- rownames(fit$beta)


s <- selected
s <- apply(abs(s),1,sum)
names(s) <- rownames(selected)
plot(sort(s))
sort(s,decreasing=TRUE)[1:20]
top1 <- names(sort(s,decreasing=TRUE)[1:3])

####################################################################################################
# tests RF
####################################################################################################

library(randomForest)
# let's find the optimal value of mtry
a <- tuneRF(x=t(trainex), y=vec.train, stepFactor=1.5, doBest=TRUE, trace=TRUE, ntreeTry=200,type="regression")

# run the model with this mtry value
fitRF <- randomForest(x=t(trainex), y=vec.train, mtry=a$mtry,do.trace=10, ntree=200, importance=TRUE,type="regression")
yhatRF <- predict(object=fitRF, newdata=t(validex),type="response")

cor(yhatRF,mek.ActArea[names(yhatRF),1],method="spearman",use="pairwise.complete.obs")
cor(yhatRF,mek.ActArea[names(yhatRF),2],method="spearman",use="pairwise.complete.obs")

varImpPlot(fitRF,n.var=10,sort=TRUE)

top2 <- names(sort(fitRF$importance[,1],decreasing=TRUE)[1:3])

####################################################################################################
# linear model based on the top genes
####################################################################################################

top <- c(top1[1:2],top2[1:2])
top <- top1[1:3]

a <- c()
b <- c()
for(i in c(1:50))
{
train1 <- sample(nsclc.mek.cells,size=36,replace=FALSE)
trainex1 <- rbind(ccle_exp[,train1],ccle_cnv[,train1],ccle_mut[,train1],cell.type.vec[train1])
rownames(trainex1) <- c(paste(rownames(ccle_exp),"_exp",sep=""),paste(rownames(ccle_cnv),"_cnv",sep=""),paste(rownames(ccle_mut),"_mut",sep=""),"cell.type.vec")
val1 <- nsclc.mek.cells[-which(nsclc.mek.cells %in% train1)]
validex1 <- rbind(ccle_exp[,val1],ccle_cnv[,val1],ccle_mut[,val1],cell.type.vec[val1])
rownames(validex1) <- rownames(trainex1)
dat1 <- trainex1[top,]
dat1 <- rbind(dat1,vec.train[colnames(dat1)])
dat1 <- as.data.frame(t(dat1))
colnames(dat1) <- c(top,"ActArea.train")
fit1 <- lm(ActArea.train ~ LRAT_exp + OSTalpha_exp + PAQR5_exp,data=dat1)
yhat1 <- predict(fit1, as.data.frame(t(validex1[top,])),type="response")
a <- c(a,cor(yhat1,mek.ActArea[names(yhat1),1],method="spearman",use="pairwise.complete.obs"))
b <- c(b,cor(yhat1,mek.ActArea[names(yhat1),2],method="spearman",use="pairwise.complete.obs"))
print(i)
}
boxplot(list(PD0325901=a,AZD6244=b), main="meki prediction in nsclc using top genes",ylab="spearman correlation with ActArea",ylim=c(0,1))
stripchart(list(PD0325901=a,AZD6244=b),vertical=TRUE,method="jitter",add=TRUE,col="aquamarine4",pch=20)
abline(h=seq(from=0,to=1,by=.1),lty=2)
fit1$effects
summary(fit1)

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
N <- 5
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
