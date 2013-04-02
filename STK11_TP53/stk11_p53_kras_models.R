# Charles Ferté
# Sage Bionetworks
# 24 Jan 2013


# train a model of differential gene expression in Lung Adenocarcinoma


#load the different packages
options(stringsAsFactors=FALSE)

library(affy)
library(corpcor)
library(lattice)
library(limma)
library(caret)
library(glmnet)
library(snm)
#source("/home/cferte/FELLOW/cferte/KRAS_Project/JUSTIN_PREDICT_CCLE/code/lung_analysis_functions.R")
synapseLogin("charles.ferte@sagebase.org","charles")


###############################################################
# load the LUAD data  
###############################################################
load("/home/jguinney/projects/AZRasPaper/data~/luad/luad_rnaseq_v3.1.6.rda")
load("/home/cferte/FELLOW/cferte/KRAS_Analysis/mutations_LUAD.RData")
load("/home/cferte/FELLOW/cferte/KRAS_Analysis/KRAS_LUAD.RData")
luad_exp <- RPKM.cqn
luad_mut <- MATMUT_LUAD
rm(Counts,RPKM.cqn,MATMUT_LUAD)
names(KRAS_LUAD) <- substr(names(KRAS_LUAD),start=1,stop=12)
# clean up the KRAS_LUAD variable
#KRAS_LUAD[!KRAS_LUAD %in% c("WT","G12C","G12A","G12V","G12D")] <- "rare"

# define the key variables
ALL_KRAS_LUAD <- luad_mut[which(rownames(luad_mut)=="KRAS"),]
STK11_LUAD <- luad_mut[which(rownames(luad_mut)=="STK11"),]
TP53_LUAD <- luad_mut[which(rownames(luad_mut)=="TP53"),]

# STK11 mutation is more frequent among KRAS mutant samples
table("KRAS"=ALL_KRAS_LUAD,"STK11"=STK11_LUAD)
prop.table(table("KRAS"=ALL_KRAS_LUAD,"STK11"=STK11_LUAD),1)
fisher.test(ALL_KRAS_LUAD,STK11_LUAD)

# Most of the STK11 mutations are concommitant with KRAS G12C mutations, among KRAS mutant samples.
table("KRAS"=KRAS_LUAD[KRAS_LUAD!="WT"],"STK11"=STK11_LUAD[KRAS_LUAD!="WT"])
prop.table(table("KRAS"=KRAS_LUAD[KRAS_LUAD!="WT"],"STK11"=STK11_LUAD[KRAS_LUAD!="WT"]),2)
fisher.test(KRAS_LUAD[KRAS_LUAD!="WT"],STK11_LUAD[KRAS_LUAD!="WT"])

# TP53 mutation is less frequent among KRAS mutant samples
table("KRAS"=ALL_KRAS_LUAD,"TP53"=TP53_LUAD)
prop.table(table("KRAS"=ALL_KRAS_LUAD,"TP53"=TP53_LUAD))
fisher.test(ALL_KRAS_LUAD,TP53_LUAD)

# TP53 mutation is distributed in G12C, G12V and in rare codons among KRAS mutant samples
table("KRAS"=KRAS_LUAD[KRAS_LUAD!="WT"],"TP53"=TP53_LUAD[KRAS_LUAD!="WT"])
prop.table(table("KRAS"=KRAS_LUAD[KRAS_LUAD!="WT"],"TP53"=TP53_LUAD[KRAS_LUAD!="WT"]),2)
fisher.test(KRAS_LUAD[KRAS_LUAD!="WT"],TP53_LUAD[KRAS_LUAD!="WT"])

# STK11 and TP53 are mutually exclusive
table("TP53"=TP53_LUAD,"STK11"=STK11_LUAD)
prop.table(table("TP53"=TP53_LUAD,"STK11"=STK11_LUAD),1)
fisher.test(STK11_LUAD,TP53_LUAD,alternative="less")

# STK11 and TP53 are also mutually exclusive within KRAS mutant samples
# 23% of KRAS mutant are mutated for STK11, 33% of KRAS mutant are mutated for TP53, 2% are comutated TP53 STK11
table("TP53"=TP53_LUAD[KRAS_LUAD!="WT"],"STK11"=STK11_LUAD[KRAS_LUAD!="WT"])
prop.table(table("TP53"=TP53_LUAD[KRAS_LUAD!="WT"],"STK11"=STK11_LUAD[KRAS_LUAD!="WT"]))
fisher.test(TP53_LUAD[KRAS_LUAD!="WT"],STK11_LUAD[KRAS_LUAD!="WT"],alternative="less")

# # find other mutually exclusive genes
# 
# # by identifying genes taht are mutually exclusive with STK11 and TP53
# kras.mat <- luad_mut[,KRAS_LUAD!="WT"]
# tmp <- names(which(apply(kras.mat,1,sum)>3 & apply(kras.mat,1,sum)<107))
# kras.mat <- kras.mat[tmp,]
# STK11_TP53 <- ifelse(STK11_LUAD[KRAS_LUAD!="WT"]==1 | TP53_LUAD[KRAS_LUAD!="WT"]==1,1,0)
# res <- apply(kras.mat,1,function(x){fisher.test(x,STK11_TP53,alternative="less")$p.value})
# names(res) <- rownames(kras.mat)
# gold <- res[which(res<.05)]
# print(gold)
# 
# # by identifying genes that are overlapping with  G12V
# kras.mat <- luad_mut[,KRAS_LUAD!="WT"]
# tmp <- names(which(apply(kras.mat,1,sum)>3 & apply(kras.mat,1,sum)<107))
# kras.mat <- kras.mat[tmp,]
# KRAS_G12V <- ifelse(KRAS_LUAD[KRAS_LUAD!="WT"] %in% c("G12V","G12C"),1,0)
# table(KRAS_G12V)
# res <- apply(kras.mat,1,function(x){fisher.test(x,KRAS_G12V,alternative="greater")$p.value})
# names(res) <- rownames(kras.mat)
# gold <- res[which(res<.05)]
# print(gold)
# paste(names(gold),collapse=" ")

# let's focus on the cell cycle genes: ATM MSH6 RAD21
ATM_LUAD <- luad_mut[which(rownames(luad_mut)=="ATM"),]
MSH6_LUAD <- luad_mut[which(rownames(luad_mut)=="MSH6"),]
RAD21_LUAD <- luad_mut[which(rownames(luad_mut)=="RAD21"),]

# all of them are mutually exclusive with STK11 and TP53 in KRAS mutant samples and are also mutually exclusive with each other !
table("ATM"=ATM_LUAD[KRAS_LUAD!="WT"],"MSH6"=MSH6_LUAD[KRAS_LUAD!="WT"])
table("ATM"=ATM_LUAD[KRAS_LUAD!="WT"],"RAD21"=RAD21_LUAD[KRAS_LUAD!="WT"])
table("MSH6"=MSH6_LUAD[KRAS_LUAD!="WT"],"RAD21"=RAD21_LUAD[KRAS_LUAD!="WT"])

# see how they are distributed among kras mutations
table("KRAS"=KRAS_LUAD[KRAS_LUAD!="WT"],"ATM"=ATM_LUAD[KRAS_LUAD!="WT"])
table("KRAS"=KRAS_LUAD[KRAS_LUAD!="WT"],"RAD21"=RAD21_LUAD[KRAS_LUAD!="WT"])
table("KRAS"=KRAS_LUAD[KRAS_LUAD!="WT"],"MSH6"=MSH6_LUAD[KRAS_LUAD!="WT"])

# visualize it
super.mat <- rbind(ALL_KRAS_LUAD[KRAS_LUAD!="WT"],ATM_LUAD[KRAS_LUAD!="WT"],STK11_LUAD[KRAS_LUAD!="WT"],TP53_LUAD[KRAS_LUAD!="WT"],MSH6_LUAD[KRAS_LUAD!="WT"],RAD21_LUAD[KRAS_LUAD!="WT"])
rownames(super.mat) <- c("KRAS","ATM","STK11","TP53","MSH6","RAD21")
super.mat <- super.mat[c("KRAS","TP53","STK11","ATM","MSH6","RAD21"),]
library(gplots)
heatmap.2(x=super.mat,col=greenred(3),trace="none")

# make the datasets coherent between mut and exp
colnames(luad_exp) <- substr(colnames(luad_exp),1,12)
colnames(luad_mut) <- substr(colnames(luad_mut),1,12)
tmp <- intersect(colnames(luad_exp),colnames(luad_mut))
luad_exp <- luad_exp[,tmp]
luad_mut <- luad_mut[,tmp]
KRAS_LUAD <- KRAS_LUAD[tmp]
STK11_LUAD <- luad_mut[which(rownames(luad_mut)=="STK11"),]
TP53_LUAD <- luad_mut[which(rownames(luad_mut)=="TP53"),]
ATM_LUAD <- luad_mut[which(rownames(luad_mut)=="ATM"),]
MSH6_LUAD <- luad_mut[which(rownames(luad_mut)=="MSH6"),]
RAD21_LUAD <- luad_mut[which(rownames(luad_mut)=="RAD21"),]
rm(tmp)

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

# let's use the ActArea as ccle_drug: "_IC50Norm" or "_ActAreaNorm"
ccle_drug <- ccle_drug_ActAreaNorm

#modify the ccle_mut into a binary matrix
ccle_mut[ccle_mut!="0"] <- 1
tmp <- rownames(ccle_mut)
ccle_mut <- apply(ccle_mut,2,as.numeric)
rownames(ccle_mut) <- tmp
rm(tmp)

# make ccle gene expression and mutation coherent
tmp <- intersect(colnames(ccle_mut),colnames(ccle_exp))
tmp <- tmp[grep(pattern="LUNG",tmp)]
tmp <- intersect(tmp,ccle_info$CCLE.name[ccle_info$Hist.Subtype1!="small_cell_carcinoma"])
ccle_exp <- ccle_exp[,tmp]
ccle_mut <- ccle_mut[,tmp]

# load the sanger data
source("/home/cferte/FELLOW/cferte/KRAS_Analysis/data_input/sanger_load_data.R")

# make coherent the tcga and the ccle genes and sanger
tmp <- intersect(rownames(luad_exp), rownames(ccle_exp))
tmp <- intersect(tmp,rownames(sanger_exp))
luad_exp <- luad_exp[tmp,]
ccle_exp <- ccle_exp[tmp,]
sanger_exp <- sanger_exp[tmp,]

# normalize the gene expression data from ccle and sanger to have the same mean and variance than luad
# Justin Guinney's function to rescale the validation data to get the same mean/var than the training set
normalize_to_X <- function(mean.x, sd.x, Y){
  m.y <- rowMeans(Y)
  sd.y <- apply(Y, 1, sd)
  Y.adj <- (Y - m.y) * sd.x / sd.y  + mean.x 
  Y.adj[sd.y == 0] <- mean.x[sd.y==0]
  Y.adj
}

ccle_exp <- normalize_to_X(rowMeans(luad_exp), apply(luad_exp, 1, sd), ccle_exp)
sanger_exp <- normalize_to_X(rowMeans(luad_exp), apply(luad_exp, 1, sd), sanger_exp)




# define the key variables for ccle_mut
kras_ccle <- ccle_mut[which(rownames(ccle_mut)=="KRAS"),]
stk11_ccle <- ccle_mut[which(rownames(ccle_mut)=="STK11"),]
tp53_ccle <- ccle_mut[which(rownames(ccle_mut)=="TP53"),]
atm_ccle <- ccle_mut[which(rownames(ccle_mut)=="ATM"),]
msh6_ccle <- ccle_mut[which(rownames(ccle_mut)=="MSH6"),]
rad21_ccle <- ccle_mut[which(rownames(ccle_mut)=="RAD21"),]

kras_stk11_ccle <- ifelse(kras_ccle==1 & stk11_ccle==1,1,0)
kras_tp53_ccle <- ifelse(kras_ccle==1 & tp53_ccle==1,1,0)
kras_atm_ccle <- ifelse(kras_ccle==1 & atm_ccle==1,1,0)
kras_msh6_ccle <- ifelse(kras_ccle==1 & msh6_ccle==1,1,0)
kras_rad21_ccle <- ifelse(kras_ccle==1 & rad21_ccle==1,1,0)

table(kras_stk11_ccle,kras_ccle)
table(kras_tp53_ccle,kras_ccle)
table(kras_atm_ccle,kras_ccle)
table(kras_msh6_ccle,kras_ccle)
table(kras_rad21_ccle,kras_ccle)



##############################################################################################################################
# build gene expression models on the two variables KRAS_STK11 and KRAS_TP53 among all KRAS
##############################################################################################################################

N <- 20
fit.tp53 <- c()
fit.stk11 <- c()
fit.neg <- c()
for(i in 1:N){
  train <- sample(names(KRAS_LUAD[KRAS_LUAD!="WT"]),replace=TRUE)
  x <- luad_exp[,train]

  y <- factor(ifelse(STK11_LUAD[train]==1,1,0))
cv.fit <- cv.glmnet(x=t(x), y=y, nfolds=3, alpha=.1, family="binomial")
  fit.stk11 <- c(fit.stk11,list(glmnet(x=t(x),y=y,family="binomial",alpha=.1,lambda=cv.fit$lambda.1se)))

  y <- factor(ifelse(TP53_LUAD[train]==1,1,0))
cv.fit <- cv.glmnet(x=t(x), y=y, nfolds=3, alpha=.1, family="binomial")
  fit.tp53 <- c(fit.tp53,list(glmnet(x=t(x),y=y,family="binomial",alpha=.1,lambda=cv.fit$lambda.1se)))
  
  y <- factor(ifelse(TP53_LUAD[train]==0 & STK11_LUAD[train]==0 ,1,0))
  cv.fit <- cv.glmnet(x=t(x), y=y, nfolds=3, alpha=.1, family="binomial")
  fit.neg <- c(fit.neg,list(glmnet(x=t(x),y=y,family="binomial",alpha=.1,lambda=cv.fit$lambda.1se)))
}

##############################################################################################################################
# apply this model in the TCGA RPPA data
# first run a sanity check and plot the AUCs
##############################################################################################################################
luad_rppa <- loadEntity('syn464306')
myfile <- paste(luad_rppa$cacheDir,list.files(luad_rppa$cacheDir),sep="/")
foo <- read.delim(file=myfile)
rownames(foo) <- foo$X.probe
colnames(foo) <- gsub(pattern=".",replacement="-",fixed=TRUE,x=colnames(foo))
foo <- foo[,grep(pattern="TCGA",x=colnames(foo))]
colnames(foo) <- substr(x=colnames(foo),start=1,stop=12)
luad_rppa <- foo

rm(foo,myfile)

tmp <- intersect(colnames(luad_rppa),names(which(KRAS_LUAD!="WT")))

yhat.stk11.rppa <- sapply(fit.stk11,function(x){ predict(x, t(luad_exp[,tmp]),type="response")})
pred <- prediction(apply(yhat.stk11.rppa,1,median), STK11_LUAD[tmp])
perf <- performance(pred, measure = "tpr", x.measure = "fpr") 
plot(perf, col=rainbow(10), main="KRAS_STK11 model")
auc <- performance(pred,"auc")
auc <- format(unlist(slot(auc, "y.values")),digits=2)
text(x=.6,y=.4,labels=paste("auc=",auc))

# correlate the yhats.stk11.rppa and the rppa data
foo <- apply(yhat.stk11.rppa,1,median)
names(foo) <- tmp
pval <- apply(luad_rppa[,tmp],1,function(x){cor.test(x,foo[tmp],method="spearman",use="pairwise.complete.obs")$p.value})
coeff <- apply(luad_rppa[,tmp],1,function(x){cor.test(x,foo[tmp],method="spearman",use="pairwise.complete.obs")$estimate})
coeff[names(sort(pval))[1:20]]

##############################################################################################################################
# apply this model in the ccle data
# first run a sanity check and plot the AUCs
##############################################################################################################################

library(ROCR)
par(mfrow=c(1,3))
xj <- ccle_exp[,kras_ccle==1]

yhat.stk11 <- sapply(fit.stk11,function(x){ predict(x, t(xj),type="response")})
pred <- prediction(apply(yhat.stk11,1,median), kras_stk11_ccle[kras_ccle==1])
perf <- performance(pred, measure = "tpr", x.measure = "fpr") 
plot(perf, col=rainbow(10), main="KRAS_STK11 model")
auc <- performance(pred,"auc")
auc <- format(unlist(slot(auc, "y.values")),digits=2)
text(x=.6,y=.4,labels=paste("auc=",auc))
rownames(yhat.stk11) <- colnames(xj)

yhat.tp53 <- sapply(fit.tp53,function(x){ predict(x, t(xj),type="response")})
pred <- prediction(apply(yhat.tp53,1,median), kras_tp53_ccle[kras_ccle==1])
perf <- performance(pred, measure = "tpr", x.measure = "fpr") 
plot(perf, col=rainbow(10), main="KRAS_TP53 model")
auc <- performance(pred,"auc")
auc <- format(unlist(slot(auc, "y.values")),digits=2)
text(x=.6,y=.4,labels=paste("auc=",auc))
rownames(yhat.tp53) <- colnames(xj)

kras_neg_ccle  <- ifelse(kras_ccle==1 & stk11_ccle==0 & tp53_ccle==0,1,0)
yhat.neg <- sapply(fit.neg,function(x){ predict(x, t(xj),type="response")})
pred <- prediction(apply(yhat.neg,1,median), kras_neg_ccle[kras_ccle==1])
perf <- performance(pred, measure = "tpr", x.measure = "fpr") 
plot(perf, col=rainbow(10), main="KRAS_neg model")
auc <- performance(pred,"auc")
auc <- format(unlist(slot(auc, "y.values")),digits=2)
text(x=.6,y=.4,labels=paste("auc=",auc))
rownames(yhat.neg) <- colnames(xj)


##############################################################################################################################
# apply this model in the sanger data
# first run a sanity check and plot the AUCs
##############################################################################################################################
kras_stk11_sanger <- ifelse(sanger_mut["KRAS",]!="wt" & sanger_mut["STK11",]!="wt",1,0)
xj <- sanger_exp[,kras_sanger!="WT"]
yhat.stk11.sanger <- sapply(fit.stk11,function(x){ predict(x, t(xj),type="response")})
pred <- prediction(apply(yhat.stk11.sanger,1,median), kras_stk11_sanger[kras_sanger!="WT"])
perf <- performance(pred, measure = "tpr", x.measure = "fpr") 
plot(perf, col=rainbow(10), main="KRAS_STK11 model")
auc <- performance(pred,"auc")
auc <- format(unlist(slot(auc, "y.values")),digits=2)
text(x=.6,y=.4,labels=paste("auc=",auc))
rownames(yhat.stk11.sanger) <- colnames(xj)


###########################################################################
#correlate the yhat.tp53 with drug sensitivity
###########################################################################

# not all the ccle cell lines have been tested for the drugs...
tmp <- intersect(rownames(yhat.tp53),rownames(ccle_drug))

a <- c()
b <- c()
for(i in 1:ncol(ccle_drug)) {
  mat1 <- yhat.tp53[tmp,]
  vec1 <- ccle_drug[tmp,i]
  a <- cbind(a,apply(mat1,2,function(X){cor.test(X,vec1,method="spearman")$p.value}))
  b <- cbind(b,apply(mat1,2,function(X){cor.test(X,vec1,method="spearman")$estimate}))
}

b <- apply(b,2,median,na.rm=TRUE)
a <- 10^(apply(log10(a),2,median,na.rm=TRUE))

names(a) <- names(b) <- colnames(ccle_drug)
cor.tp53 <- rbind(b,a)
rownames(cor.tp53) <- c("correlation","p.value")
cor.tp53 <- cor.tp53[,names(sort(b))]





###########################################################################
#correlate the yhat.stk11 with drug sensitivity
###########################################################################

# not all the ccle cell lines have been tested for the drugs...
tmp <- intersect(rownames(yhat.stk11),rownames(ccle_drug))

a <- c()
b <- c()
for(i in 1:ncol(ccle_drug)) {
  mat1 <- yhat.stk11[tmp,]
  vec1 <- ccle_drug[tmp,i]
  a <- cbind(a,apply(mat1,2,function(X){cor.test(X,vec1,method="spearman")$p.value}))
  b <- cbind(b,apply(mat1,2,function(X){cor.test(X,vec1,method="spearman")$estimate}))
}

b <- apply(b,2,median,na.rm=TRUE)
a <- 10^(apply(log10(a),2,median,na.rm=TRUE))

names(a) <- names(b) <- colnames(ccle_drug)
cor.stk11 <- rbind(b,a)
rownames(cor.stk11) <- c("correlation","p.value")
cor.stk11 <- cor.stk11[,names(sort(b))]

# tmp2 <- ifelse(ccle_mut["STK11",tmp]==1 & ccle_mut["KRAS",tmp]==1,1,0)
# boxplot(ccle_drug_ActAreaNorm[tmp,16]~tmp2, main=paste(colnames(ccle_drug_ActAreaNorm)[16]),xlab="KRAS mut & STK11 mut",ylab="actarea")
# boxplot(ccle_drug_ActAreaNorm[tmp,17]~tmp2, main=paste(colnames(ccle_drug_ActAreaNorm)[17]),xlab="KRAS mut & STK11 mut",ylab="actarea")
# 
# 
# 
# boxplot(ccle_drug_ActAreaNorm[intersect(colnames(ccle_mut),rownames(ccle_drug_ActAreaNorm)),16]~ccle_mut["KRAS",intersect(colnames(ccle_mut),rownames(ccle_drug_ActAreaNorm))])
# ccle_drug_ActAreaNorm[intersect(colnames(ccle_mut),rownames(ccle_drug_ActAreaNorm)),16])

###########################################################################
#correlate the yhat.neg with drug sensitivity
###########################################################################

# not all the ccle cell lines have been tested for the drugs...
tmp <- intersect(rownames(yhat.neg),rownames(ccle_drug))

a <- c()
b <- c()
for(i in 1:ncol(ccle_drug)) {
  mat1 <- yhat.neg[tmp,]
  vec1 <- ccle_drug[tmp,i]
  a <- cbind(a,apply(mat1,2,function(X){cor.test(X,vec1,method="spearman")$p.value}))
  b <- cbind(b,apply(mat1,2,function(X){cor.test(X,vec1,method="spearman")$estimate}))
}

b <- apply(b,2,median,na.rm=TRUE)
a <- 10^(apply(log10(a),2,median,na.rm=TRUE))

names(a) <- names(b) <- colnames(ccle_drug)
cor.neg <- rbind(b,a)
rownames(cor.neg) <- c("correlation","p.value")
cor.neg <- cor.neg[,names(sort(b))]

###########################################################################
# plot the drug sensitivity correlations as volcano plots !
###########################################################################

par(mfrow=c(1,1))

col.drug <- topo.colors(ncol(ccle_drug))
names(col.drug) <- colnames(cor.tp53)

plot(cor.stk11["correlation",],-log10(cor.stk11["p.value",]),xlab="correlation with drug sensitivity (Spearman rho)",
     ylab="significance:  -log10(p.value)",pch=19,
     col=col.drug[colnames(cor.stk11)], xlim=c(-.7,.7),ylim=c(0,2),
     main="kras.stk11 model ~ drug sensitivity",cex=-log10(cor.stk11["p.value",]))
abline(h=-log10(.05),col="red",lty=2)
text(x=cor.stk11["correlation",]+.2,y=-log10(cor.stk11["p.value",]),labels=colnames(cor.stk11),cex=.6)

plot(cor.tp53["correlation",],-log10(cor.tp53["p.value",]),xlab="correlation with drug sensitivity (Spearman rho)",
     ylab="significance:  -log10(p.value)",pch=19,
     col=col.drug[colnames(cor.tp53)], xlim=c(-.7,.7),ylim=c(0,2),
     main="kras.tp53 model ~ drug sensitivity",cex=-log10(cor.tp53["p.value",]))
abline(h=-log10(.05),col="red",lty=2)
text(x=cor.tp53["correlation",]+.2,y=-log10(cor.tp53["p.value",]),labels=colnames(cor.tp53),cex=.6)

plot(cor.neg["correlation",],-log10(cor.neg["p.value",]),xlab="correlation with drug sensitivity (Spearman rho)",
     ylab="significance:  -log10(p.value)",pch=19,
     col=col.drug[colnames(cor.neg)], xlim=c(-.7,.7),ylim=c(0,2),
     main="kras.neg model ~ drug sensitivity",cex=-log10(cor.neg["p.value",]))
abline(h=-log10(.05),col="red",lty=2)
text(x=cor.neg["correlation",]+.2,y=-log10(cor.neg["p.value",]),labels=colnames(cor.neg),cex=.6)


#######################################################

###########################################################################
#correlate the yhat.stk11.sanger with drug sensitivity
###########################################################################

# not all the sanger cell lines have been tested for the drugs...
tmp <- intersect(rownames(yhat.stk11.sanger),rownames(sanger_drug))

a <- c()
b <- c()
for(i in 1:ncol(sanger_drug)) {
  mat1 <- yhat.stk11.sanger[tmp,]
  vec1 <- sanger_drug[tmp,i]
  a <- cbind(a,apply(mat1,2,function(X){cor.test(X,vec1,method="spearman")$p.value}))
  b <- cbind(b,apply(mat1,2,function(X){cor.test(X,vec1,method="spearman")$estimate}))
}

b <- apply(b,2,median,na.rm=TRUE)
a <- 10^(apply(log10(a),2,median,na.rm=TRUE))

names(a) <- names(b) <- colnames(sanger_drug)
cor.stk11.sanger <- rbind(b,a)
rownames(cor.stk11.sanger) <- c("correlation","p.value")
cor.stk11.sanger <- cor.stk11.sanger[,names(sort(b))]

par(mfrow=c(1,1))
col.drug <- topo.colors(ncol(sanger_drug))
names(col.drug) <- colnames(cor.stk11.sanger)
plot(-cor.stk11.sanger["correlation",],-log10(cor.stk11.sanger["p.value",]),xlab="correlation with drug sensitivity (Spearman rho)",
     ylab="significance:  -log10(p.value)",pch=19,
     col=col.drug[colnames(cor.stk11.sanger)], xlim=c(-.6,.9),ylim=c(0,2),
     main="kras.stk11 model ~ drug sensitivity",cex=-log10(cor.stk11.sanger["p.value",]))
abline(h=-log10(.05),col="red",lty=2)
text(x=-cor.stk11.sanger["correlation",]+.2,y=-log10(cor.stk11.sanger["p.value",]),labels=colnames(cor.stk11.sanger),cex=.6)

