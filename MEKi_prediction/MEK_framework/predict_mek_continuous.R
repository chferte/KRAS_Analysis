# Charles Ferté
# Sage Bionetworks
# 14 Feb 2012

#############################
# load the data
#############################

# load the mek data from ccle
source("/home/cferte/FELLOW/cferte/KRAS_Analysis/MEKi_prediction/MEK_framework/load_mek_ccle.R")


# load the  cell status and probabilities
ccle_probs_status <- loadEntity("syn1709732")
ccle_probs_status <- ccle_probs_status$objects$ccle_probs_status
all.prob <- ccle_probs_status[[1]]
cell.status <- ccle_probs_status[[2]]

#######################################################
# load the validation set and make it coherent with the training data 
#######################################################


# foo  <-  loadEntity("syn418003")
# val_exp <- read.table(list.files(foo$cacheDir,full.names=TRUE),
#                       row.names=1,comment="",quote="",sep="\t",header=TRUE)
# tmp <- sapply(strsplit(x=rownames(val_exp),split="|",fixed=TRUE),function(x){x[[1]]})
# tmp1 <- unique(c(which(tmp=="?"), which(duplicated(tmp)==TRUE)))
# val_exp <- val_exp[-tmp1,]
# rownames(val_exp) <- tmp[-tmp1]
# m <- apply(val_exp, 1, mean)
# val_exp <- val_exp[m > 1,]
# val_exp <- log(val_exp) + 1
# rm(tmp,tmp1,m)

# load the crc gene expression tcga data
source("/home/cferte/FELLOW/cferte/KRAS_Analysis/MEKi_prediction/MEK_framework/load_crc_tcga_data.R")
#val.set <- crc.g

# load the breast gene expression tcga data
#source("/home/cferte/FELLOW/cferte/KRAS_Analysis/MEKi_prediction/MEK_framework/load_breast_tcga_data.R")

# load the luad gene expression tcga data
#source("/home/cferte/FELLOW/cferte/KRAS_Analysis/MEKi_prediction/MEK_framework/load_luad_tcga_data.R")


#load the AML tcga data
#source("/home/cferte/FELLOW/cferte/KRAS_Analysis/MEKi_prediction/MEK_framework/load_laml_tcga_data.R")

#load the Melanoma tcga data
#source("/home/cferte/FELLOW/cferte/KRAS_Analysis/MEKi_prediction/MEK_framework/load_laml_tcga_data.R")

# load the crc pdx nki gene expression data
#source("/home/cferte/FELLOW/cferte/KRAS_Analysis/MEKi_prediction/MEK_framework/validation_PDX_NKI.R")  


# # load the tcga lung mutation data as val set
# foo  <-  loadEntity("syn1676707")
# foo <- foo$objects$luad_data
# val_mut <- foo[[2]]
# colnames(val_mut) <- substr(colnames(val_mut),1,12)
# colnames(val_mut) <- gsub(pattern="-",replacement=".",x=colnames(val_mut),fixed=TRUE)
# rm(foo)
# val.set <- val_mut

#######################################################
# predictive modeling
# we predict the ActArea 
# training in the gene expression data of all cells. 
#######################################################

# define globalmatrix (gene expression only)
global.matrix <- ccle_exp

# define globalmatrix (mutations only) for the lung model
#global.matrix <- ccle_mut

#######################################################
# Make the val set coherent with the training data 
#######################################################
tmp <- intersect(rownames(val.set),rownames(global.matrix))
val.set <- val.set[tmp,]
global.matrix <- global.matrix[tmp,]
rm(tmp)

#######################################################
# create a penalty vector
#######################################################
# pen.vec <- rep(x=1,times=nrow(global.matrix))
# pen.vec[which(rownames(global.matrix) %in% c("STK11","TP53","KRAS"))] <- 0
# pen.vec[which(rownames(global.matrix) %in% c("BRAF_mut","PTEN_mut","CTNNB1_mut","PIK3CA_mut","TP53_mut","KRAS_mut"))] <- 0
# names(pen.vec) <- rownames(global.matrix)


#######################################################
# create a penalty vector for 
#######################################################
# pen.vec <- rep(x=1,times=nrow(global.matrix))
# pen.vec[which(rownames(global.matrix) %in% c("STK11","TP53","KRAS"))] <- 0
# pen.vec[which(rownames(global.matrix) %in% c("BRAF_mut","PTEN_mut","CTNNB1_mut","PIK3CA_mut","TP53_mut","KRAS_mut"))] <- 0
# names(pen.vec) <- rownames(global.matrix)



#######################################################
# start the statistical learning process
#######################################################
require(multicore)

N=100
i <- 0

PARAL <- mclapply(X=1:N,FUN=function(x){
  print(i)
  i <- 1+1
  
  #standard sampling
  train <- sample(mek.cells,replace=TRUE)
  
  # training for the lung model idx
  #train <- sample(cell.names.idx,replace=TRUE)
  
  # carcinoma model (exclude hemato melanoma glioma)
  #train <- sample(mek.cells[-which(mek.cells %in% c(hemal.mek.cells,glioma.mek.cells,melanoma.mek.cells))],replace=TRUE)
  
  # tissue specific models
  #train <- sample(mek.cells,replace=TRUE)
    
  #################################
  # set up the vector of response
  #################################
  vec.train <-apply(ccle_drug[train,mek.inhib],1,mean)
  
 
  #################################
  # fitting the model
  #################################
  # standard fitting
  cv.fit <- cv.glmnet(t(global.matrix[,train]), y=vec.train,nfolds=3, alpha=.1)
  fit <- glmnet(x=t(global.matrix[,train]),y=vec.train,alpha=.1,lambda=cv.fit$lambda.1se)
  
  # fitting with penalty factor
  #cv.fit <- cv.glmnet(t(global.matrix[,train]), y=vec.train,nfolds=3, alpha=.1,penalty.factor=pen.vec)
  #fit <- glmnet(x=t(global.matrix[,train]),y=vec.train,alpha=.1,lambda=cv.fit$lambda.1se,penalty.factor=pen.vec)
  
  return(list(fit,train)) },mc.set.seed=TRUE,mc.cores=6)
  

#######################################################
# get the yhats for INTERNAL Validation
#######################################################

yhat.all <- c()
yhat.breast <- c()
yhat.nsclc <- c()
yhat.crc <- c()
yhat.pancreas <- c()
yhat.ovary <- c()
yhat.melanoma <- c()
yhat.hemal  <- c()
selected <- c()

for(i in c(1:N)){
  train <- PARAL[[i]][[2]]
  fit <- PARAL[[i]][[1]]
  val <-mek.cells[-which(mek.cells %in% train)]
  selected <- c(selected,list(fit$beta))
  yhat.all <- c(yhat.all,list(predict(fit, t(global.matrix[,val]))))
  yhat.breast <- c(yhat.breast,list(predict(fit,t(global.matrix[,breast.mek.cells[-which(breast.mek.cells %in% train)]]))))
  yhat.nsclc <- c(yhat.nsclc,list(predict(fit,t(global.matrix[,nsclc.mek.cells[-which(nsclc.mek.cells %in% train)]]))))
  yhat.crc <- c(yhat.crc,list(predict(fit,t(global.matrix[,crc.mek.cells[-which(crc.mek.cells %in% train)]]))))
  #   yhat.melanoma <- c(yhat.melanoma,list(predict(fit,t(global.matrix[,melanoma.mek.cells]))))
  #   yhat.hemal <- c(yhat.hemal,list(predict(fit,t(global.matrix[,hemal.mek.cells]))))
  
  yhat.pancreas <- c(yhat.pancreas,list(predict(fit,t(global.matrix[,pancreas.mek.cells[-which(pancreas.mek.cells %in% train)]]))))
  yhat.ovary <- c(yhat.ovary,list(predict(fit,t(global.matrix[,ovary.mek.cells[-which(ovary.mek.cells %in% train)]]))))
  yhat.melanoma <- c(yhat.melanoma,list(predict(fit,t(global.matrix[,melanoma.mek.cells[-which(melanoma.mek.cells %in% train)]]))))
  yhat.hemal <- c(yhat.hemal,list(predict(fit,t(global.matrix[,hemal.mek.cells[-which(hemal.mek.cells %in% train)]]))))
  print(i)}  

# asess the perf of each model 

obs.sens <- apply(ccle_drug[,mek.inhib],1,mean)

# crc
abc <- matrix(NA,nrow=length(crc.mek.cells),ncol=N)
rownames(abc) <- crc.mek.cells
for(i in 1:N){abc[rownames(yhat.crc[[i]]),i] <- yhat.crc[[i]]}
boxplot(t(abc),col="red")
COR.crc <- c()
for(i in 1:N){
predictor <- as.numeric(abc[which(!is.na(abc[,i])),i])
response <- as.numeric(obs.sens[rownames(abc)[which(!is.na(abc[,i]))]])
roc <- roc(predictor=predictor, response=response)
COR.crc <- c(COR.crc,cor.test(predictor, response,method="spearman")$p.value)
}  
a <- COR.crc<.1
print(sum(a))
retained <- names(!is.na(apply(abc[,a],1,function(x) median(x,na.rm=TRUE))))
print("normal")
cor.test(apply(abc[,a],1,function(x) median(x,na.rm=TRUE))[retained],obs.sens[retained],method="spearman")
print("turbo")
cor.test(apply(abc,1,function(x) median(x,na.rm=TRUE)),obs.sens[crc.mek.cells],method="spearman")

# breast
abc <- matrix(NA,nrow=length(breast.mek.cells),ncol=N)
rownames(abc) <- breast.mek.cells
for(i in 1:N){abc[rownames(yhat.breast[[i]]),i] <- yhat.breast[[i]]}
boxplot(t(abc),col="red")
COR.breast <- c()
for(i in 1:N){
  predictor <- as.numeric(abc[which(!is.na(abc[,i])),i])
  response <- as.numeric(obs.sens[rownames(abc)[which(!is.na(abc[,i]))]])
  roc <- roc(predictor=predictor, response=response)
  COR.breast <- c(COR.breast,cor.test(predictor, response,method="spearman")$p.value)
}  
a <- COR.breast<.1
print(sum(a))
retained <- names(!is.na(apply(abc[,a],1,function(x) median(x,na.rm=TRUE))))
print("turbo")
cor.test(apply(abc[,a],1,function(x) median(x,na.rm=TRUE))[retained],obs.sens[retained],method="spearman")
print("normal")
cor.test(apply(abc,1,function(x) median(x,na.rm=TRUE)),obs.sens[breast.mek.cells],method="spearman")


# nsclc
abc <- matrix(NA,nrow=length(nsclc.mek.cells),ncol=N)
rownames(abc) <- nsclc.mek.cells
for(i in 1:N){abc[rownames(yhat.nsclc[[i]]),i] <- yhat.nsclc[[i]]}
boxplot(t(abc),col="red")
COR.nsclc <- c()
for(i in 1:N){
  predictor <- as.numeric(abc[which(!is.na(abc[,i])),i])
  response <- as.numeric(obs.sens[rownames(abc)[which(!is.na(abc[,i]))]])
  roc <- roc(predictor=predictor, response=response)
  COR.nsclc <- c(COR.nsclc,cor.test(predictor, response,method="spearman")$p.value)
}  

a <- COR.nsclc <.1
print(sum(a))
retained <- names(!is.na(apply(abc[,a],1,function(x) median(x,na.rm=TRUE))))
print("turbo")
cor.test(apply(abc[,a],1,function(x) median(x,na.rm=TRUE))[retained],obs.sens[retained],method="spearman")
print("normal")
cor.test(apply(abc,1,function(x) median(x,na.rm=TRUE)),obs.sens[nsclc.mek.cells],method="spearman")

# ovary
abc <- matrix(NA,nrow=length(ovary.mek.cells),ncol=N)
rownames(abc) <- ovary.mek.cells
for(i in 1:N){abc[rownames(yhat.ovary[[i]]),i] <- yhat.ovary[[i]]}
boxplot(t(abc),col="red")
COR.ovary <- c()
for(i in 1:N){
  predictor <- as.numeric(abc[which(!is.na(abc[,i])),i])
  response <- as.numeric(obs.sens[rownames(abc)[which(!is.na(abc[,i]))]])
  roc <- roc(predictor=predictor, response=response)
  COR.ovary <- c(COR.ovary,cor.test(predictor, response,method="spearman")$p.value)
}  

a <- COR.ovary <.1
print(sum(a))
retained <- names(!is.na(apply(abc[,a],1,function(x) median(x,na.rm=TRUE))))
print("turbo")
cor.test(apply(abc[,a],1,function(x) median(x,na.rm=TRUE))[retained],obs.sens[retained],method="spearman")
print("normal")
cor.test(apply(abc,1,function(x) median(x,na.rm=TRUE)),obs.sens[ovary.mek.cells],method="spearman")



