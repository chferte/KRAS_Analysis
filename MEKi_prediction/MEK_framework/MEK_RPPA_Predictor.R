# charles ferté
# march 21
# predict MEKi based on pERK

# load rppa
e <- loadEntity("syn464306")
LUAD_RPPA <- read.table(paste(e$cacheDir,e$files,sep="/"),header=TRUE,row.names=1,comment="",quote="",sep="\t")
colnames(LUAD_RPPA) <- gsub(pattern=".",replacement="-",x=colnames(LUAD_RPPA),fixed=TRUE)

# load tcga luad 
load("/home/jguinney/projects/AZRasPaper/data~/luad/luad_rnaseq_v3.1.6.rda")
load("/home/cferte/FELLOW/cferte/KRAS_Analysis/mutations_LUAD.RData")
LUAD_EXP <- RPKM.cqn
LUAD_MUT <- MATMUT_LUAD

# clean up - everybody does his share !
rm(MATMUT_LUAD, Counts,RPKM.cqn,e)

# make the names of the samples coherent
colnames(LUAD_EXP) <- substr(x=colnames(LUAD_EXP),start=1,stop=12)
colnames(LUAD_MUT) <- substr(x=colnames(LUAD_MUT),start=1,stop=12)
colnames(LUAD_RPPA) <- substr(x=colnames(LUAD_RPPA),start=1,stop=12)

# make the datasets coherents
tmp <- intersect(colnames(LUAD_RPPA), colnames(LUAD_MUT))
tmp <- intersect(tmp,colnames(LUAD_EXP))
LUAD_EXP <- LUAD_EXP[,tmp]
LUAD_MUT <- LUAD_MUT[,tmp]
LUAD_RPPA <- LUAD_RPPA[,tmp]


# se what looks like the RPPA data
s <- svd(LUAD_RPPA)
plot(s$v[,1],s$v[,2])
rownames(LUAD_RPPA)

# target
rownames(LUAD_RPPA)
target <- "K-Ras-M-C"

# plot density target
plot(density(as.numeric(LUAD_RPPA[target,])),main=paste(target, "RPPA in tcga luad"))

# build luad.matrix
luad.matrix <- rbind(LUAD_EXP,LUAD_MUT)
rownames(luad.matrix) <- c(paste(rownames(LUAD_EXP),"_exp",sep=""),paste(rownames(LUAD_MUT),"_mut",sep=""))

####################################################################################################################
# load the mek data from ccle
####################################################################################################################
source("/home/cferte/FELLOW/cferte/KRAS_Analysis/MEKi_prediction/MEK_framework/load_mek_ccle.R")

# load the  cell status and probabilities
ccle_probs_status <- loadEntity("syn1709732")
ccle_probs_status <- ccle_probs_status$objects$ccle_probs_status
all.prob <- ccle_probs_status[[1]]
cell.status <- ccle_probs_status[[2]]

####################################################################################################################
# set up the nsclc.ccle.matrix
####################################################################################################################

nsclc.ccle.matrix <- rbind(ccle_exp,ccle_mut)
rownames(nsclc.ccle.matrix) <- c(paste(rownames(ccle_exp),"_exp",sep=""),paste(rownames(ccle_mut),"_mut",sep=""))

####################################################################################################################
# make coherent the TCGA and ccle matrix
####################################################################################################################
tmp <- intersect(rownames(luad.matrix),rownames(nsclc.ccle.matrix))
luad.matrix <- luad.matrix[tmp,]
nsclc.ccle.matrix <- nsclc.ccle.matrix[tmp,]

####################################################################################################################
# Penalized regression
# train a model to predict ERK RPPA in TCGA
####################################################################################################################

require(multicore)
yhat.luad <- c()
selected <- c()
N <- 50
PARAL <- mclapply(X=1:N,FUN=function(x){
  print(i)
  train <- sample(colnames(luad.matrix),replace=TRUE)
  vec.train <-as.numeric(LUAD_RPPA[target,train])
  names(vec.train) <- train
cv.fit <- cv.glmnet(t(luad.matrix[,train]), y=vec.train,nfolds=3, alpha=.1)
fit <- glmnet(x=t(luad.matrix[,train]),y=vec.train,alpha=.1,lambda=cv.fit$lambda.1se)
  return(list(fit,train)) },mc.set.seed=TRUE,mc.cores=6)


for(i in c(1:N)){
  train <- PARAL[[i]][[2]]
  fit <- PARAL[[i]][[1]]
  val <-colnames(LUAD_EXP)[-which(colnames(LUAD_EXP) %in% train)]
  selected <- c(selected,list(fit$beta))
  yhat.luad <- c(yhat.luad,list(predict(fit, t(luad.matrix[,val]))))
  print(i)}  


# first transform the list of yhats into a matrix 
abc <- c()
abc <- as.data.frame(matrix(NA,nrow=ncol(LUAD_EXP),ncol=N))
rownames(abc) <- colnames(LUAD_EXP)
colnames(abc) <- c(1:N)
for(l in c(1:N)){abc[rownames(yhat.luad[[l]]),l] <- yhat.luad[[l]]}
abc <- apply(abc,1,median,na.rm=TRUE)

# internal validation
plot(abc,as.numeric(LUAD_RPPA[target,names(abc)]),pch=20,xlab="pred",ylab="obs")
cor.test(abc,as.numeric(LUAD_RPPA[target,names(abc)]),method="spearman")

# let's build a consensus predictor by aggregating the selected matrix and training a ridge like model
abc <- c()
abc <- as.data.frame(matrix(NA,nrow=nrow(fit$beta),ncol=N))
rownames(abc) <- rownames(fit$beta)
colnames(abc) <- c(1:N)
for(l in c(1:N)){abc[rownames(selected[[l]]),l] <- as.numeric(selected[[l]])}
abc[abs(abc)>0] <- 1
rownames(abc) <- rownames(fit$beta)
hist(log10(apply(abc,1,sum)),col="red")
abline(v=log10(quantile(apply(abc,1,sum),probs=.95)),col="blue")
tmp <- names(which(apply(abc,1,sum)>quantile(apply(abc,1,sum),probs=.95)))
cv.fit2 <- cv.glmnet(t(luad.matrix[tmp,]), y=as.numeric(LUAD_RPPA[target,]),nfolds=3, alpha=0)
fit2 <- glmnet(x=t(luad.matrix[,train]),y=as.numeric(LUAD_RPPA[target,]),alpha=0,lambda=cv.fit2$lambda.1se)

# let's see how this work in lung cell lines to predict MEK inhibition
pred <- predict(fit2, t(nsclc.ccle.matrix))
obs <- apply(mek.ActArea,1,mean)

plot(pred[nsclc.mek.cells,],obs[nsclc.mek.cells])
plot(pred,obs)
cor.test(pred[nsclc.mek.cells,],obs[nsclc.mek.cells])

