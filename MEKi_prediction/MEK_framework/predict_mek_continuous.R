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
# load the validation data and make it coherent with the training data 
#######################################################
foo  <-  loadEntity("syn418003")
val_exp <- read.table(list.files(foo$cacheDir,full.names=TRUE),row.names=1,comment="",quote="",sep="\t",header=TRUE)
tmp <- sapply(strsplit(x=rownames(val_exp),split="|",fixed=TRUE),function(x){x[[1]]})
tmp1 <- unique(c(which(tmp=="?"), which(duplicated(tmp)==TRUE)))
val_exp <- val_exp[-tmp1,]
rownames(val_exp) <- tmp[-tmp1]
rm(tmp,tmp1)

tmp <- intersect(rownames(val_exp),rownames(ccle_exp))
val_exp <- val_exp[tmp,]
ccle_exp <- ccle_exp[tmp,]
rm(tmp)

#######################################################
# predictive modeling
# we predict the ActArea 
# training in the gene expression data of all cells. 
#######################################################

# define globalmatrix (without eigengenes)
#global.matrix <- rbind(ccle_exp,ccle_cnv,ccle_mut)
#rownames(global.matrix) <- c(paste(rownames(ccle_exp),"_exp",sep=""),paste(rownames(ccle_cnv),"_cnv",sep=""),paste(rownames(ccle_mut),"_mut",sep=""))

# define globalmatrix (gene expression only)
global.matrix <- ccle_exp
rownames(global.matrix) <- c(paste(rownames(ccle_exp),"_exp",sep=""))

# define globalmatrix (gene expression + mutations)
#global.matrix <- rbind(ccle_exp,ccle_mut[c("STK11","TP53","KRAS"),])
#rownames(global.matrix) <- c(paste(rownames(ccle_exp),"_exp",sep=""),paste(rownames(ccle_mut[c("STK11","TP53","KRAS"),]),"_mut",sep=""))

# define globalmatrix (mutations only)
#global.matrix <- ccle_mut
#rownames(global.matrix) <- paste(rownames(ccle_mut),"_mut",sep="")

# define globalmatrix (the "3" lung mutations)
#global.matrix <- ccle_mut[c("STK11","TP53","KRAS"),]
#rownames(global.matrix) <- paste(rownames(ccle_mut[c("STK11","TP53","KRAS"),]),"_mut",sep="")

# create a vector for penalty
#pen.vec <- rep(x=1,times=nrow(global.matrix))
#pen.vec[which(rownames(global.matrix) %in% c("STK11_mut","TP53_mut","KRAS_mut"))] <- 0
#pen.vec[which(rownames(global.matrix) %in% c("BRAF_mut","PTEN_mut","CTNNB1_mut","PIK3CA_mut","TP53_mut","KRAS_mut"))] <- 0
#names(pen.vec) <- rownames(global.matrix)



# # define globalmatrix (with tissue specificity)
# abc <- model.matrix(~tissue.origin)
# tmp <- colnames(abc)
# abc <- t(abc)
# rownames(abc) <- tmp
# colnames(abc) <- names(tissue.origin)
# rm(tmp)
# global.matrix <- rbind(ccle_exp,ccle_cnv,ccle_mut,abc)
# rownames(global.matrix) <- c(paste(rownames(ccle_exp),"_exp",sep=""),paste(rownames(ccle_cnv),"_cnv",sep=""),paste(rownames(ccle_mut),"_mut",sep=""),rownames(abc))

# not penalize the tissue spec matrix
#pen.vec <- c(rep(1,times=nrow(global.matrix)-nrow(abc)),rep(0,times=nrow(abc)))


# # assign weights base on the density of ActArea
# tmp <- all.prob[[1]]
# prob.weights <- as.numeric(apply(tmp,1,function(x){max(x)}))
# prob.weights <- prob.weights
# names(prob.weights) <- rownames(tmp)
# plot(density(prob.weights), main="posterior probabilities weights")

# define globalmatrix2 (with eigengenes)
# eigengenes <- t(eigengenes[colnames(global.matrix),c(1:30)])
# global.matrix2 <- rbind(global.matrix,eigengenes)
# rownames(global.matrix2) <- c(rownames(global.matrix),paste("PC",c(1:30),sep=""))

# assign weights base on the density of ActArea
# mean.mek.sens <- apply(ccle_drug[mek.cells,mek.inhib],1,mean)
# density.id <- sapply(c(1:length(mean.mek.sens)),function(x){min(which(density(mean.mek.sens)$x>mean.mek.sens[x]))})
# density.weights <- log(1/(density(mean.mek.sens)$y[density.id]))
# names(density.weights) <- names(mean.mek.sens)
# rm(density.id,mean.mek.sens)
# density.weights  <- rep(1,times=length(density.weights))                    
# 
# summary(density.weights)
# # boxplot(density.weights)
# density.weights <- rep(0.5,times=length(mek.cells))
# names(density.weights) <- mek.cells
# density.weights[names(new.weight)] <- new.weight

require(multicore)

N=150
#models <- 0
i <- 0

# # 
# # set up the Q1:Q4 for running balanced models
# q25 <- quantile(apply(ccle_drug[mek.cells,mek.inhib],1,mean),probs=.35) 
# q50 <- quantile(apply(ccle_drug[mek.cells,mek.inhib],1,mean),probs=.5)
# q75 <- quantile(apply(ccle_drug[mek.cells,mek.inhib],1,mean),probs=.65)
# Q1 <- names(which(apply(ccle_drug[mek.cells,mek.inhib],1,mean) < q25))
# Q2 <- names(which(apply(ccle_drug[mek.cells,mek.inhib],1,mean) > q25 & apply(ccle_drug[mek.cells,mek.inhib],1,mean) < q50 ))
# Q3 <- names(which(apply(ccle_drug[mek.cells,mek.inhib],1,mean) > q50 & apply(ccle_drug[mek.cells,mek.inhib],1,mean) < q75 ))
# Q4 <- names(which(apply(ccle_drug[mek.cells,mek.inhib],1,mean) > q75))
# rm(q25,q50,q75)



PARAL <- mclapply(X=1:N,FUN=function(x){
  print(i)
  i <- 1+1
#while(models<N)
#{
  
  #standard sampling
  #train <- sample(mek.cells,replace=TRUE)
  #train <- sample(mek.cells[-which(mek.cells %in% c(hemal.mek.cells,glioma.mek.cells,melanoma.mek.cells))],replace=TRUE)
  
  # balanced model
  #train <- c(sample(x=c(Q1),replace=TRUE,size=220),sample(Q4,replace=TRUE,size=220))
  
  # weighted model on the distribution
  #train <- sample(mek.cells,replace=TRUE)
  
  
  # tissue specific models
  train <- sample(mek.cells,replace=TRUE)
  
  # weighted lung models
  #train <- c(sample(mek.cells,replace=TRUE,size=144),sample(nsclc.mek.cells,replace=TRUE))
  
  # mixed balanced & tissue specific model
  #train  <- c(sample(nsclc.mek.cells,replace=TRUE),sample(x=c(Q1),replace=TRUE,size=36),sample(Q4,replace=TRUE,size=36))
  
  vec.train <-apply(ccle_drug[train,mek.inhib],1,mean)
  
  # standard training
  cv.fit <- cv.glmnet(t(global.matrix[,train]), y=vec.train,nfolds=3, alpha=.1)
  fit <- glmnet(x=t(global.matrix[,train]),y=vec.train,alpha=.1,lambda=cv.fit$lambda.1se)
  
  # weighted models
  #cv.fit <- cv.glmnet(t(global.matrix[,train]), y=vec.train,nfolds=3, alpha=.1,weights=prob.weights[train])
  #fit <- glmnet(x=t(global.matrix[,train]),y=vec.train,alpha=.1,lambda=cv.fit$lambda.1se,weights=prob.weights[train])
  
  # penalty factor
  #cv.fit <- cv.glmnet(t(global.matrix[,train]), y=vec.train,nfolds=3, alpha=.1,penalty.factor=pen.vec)
  #fit <- glmnet(x=t(global.matrix[,train]),y=vec.train,alpha=.1,lambda=cv.fit$lambda.1se,penalty.factor=pen.vec)
  
  return(list(fit,train)) },mc.set.seed=TRUE,mc.cores=6)
  


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
  
#   yhat.glioma <- c(yhat.glioma,list(predict(fit,t(global.matrix[,glioma.mek.cells]))))
#   yhat.melanoma <- c(yhat.melanoma,list(predict(fit,t(global.matrix[,melanoma.mek.cells]))))
#   yhat.hemal <- c(yhat.hemal,list(predict(fit,t(global.matrix[,hemal.mek.cells]))))
  yhat.pancreas <- c(yhat.pancreas,list(predict(fit,t(global.matrix[,pancreas.mek.cells[-which(pancreas.mek.cells %in% train)]]))))
  yhat.ovary <- c(yhat.ovary,list(predict(fit,t(global.matrix[,ovary.mek.cells[-which(ovary.mek.cells %in% train)]]))))
  yhat.melanoma <- c(yhat.melanoma,list(predict(fit,t(global.matrix[,melanoma.mek.cells[-which(melanoma.mek.cells %in% train)]]))))
  yhat.hemal <- c(yhat.hemal,list(predict(fit,t(global.matrix[,hemal.mek.cells[-which(hemal.mek.cells %in% train)]]))))
  print(i)}  



#####################################################################################################################
# save the fit objects
#####################################################################################################################

#global_model_yhats <- list(yhat.all,yhat.nsclc,yhat.breast,yhat.crc,yhat.hemal,yhat.melanoma,yhat.pancreas,yhat.ovary)
#save(global_model_yhats,file="/home/cferte/RESULTS/MEKi/GLOBAL_MODEL/ROBJECTS/global_model_yhats.Rda")

# names(pure_tissue_models_yhats) <- c("yhat.all","yhat.nsclc","yhat.breast","yhat.crc","yhat.hemal","yhat.glioma","yhat.melanoma")
# save(pure_tissue_models_yhats,file="/home/cferte/RESULTS/MEKi/GLOBAL_MODEL/ROBJECTS/pure_tissue.models.Rda")

# stk11 kras tp53 mutation models
#stk11_kras_tp53_yhats <- list(yhat.all,yhat.nsclc,yhat.breast,yhat.crc,yhat.hemal,yhat.melanoma,yhat.pancreas,yhat.ovary)
#save(stk11_kras_tp53_yhats,file="/home/cferte/RESULTS/MEKi/GLOBAL_MODEL/ROBJECTS/stk11_kras_tp53_yhats.Rda")

# gene expression + stk11 kras tp53 mutation models
#stk11_kras_tp53_genexp_yhats <- list(yhat.all,yhat.nsclc,yhat.breast,yhat.crc,yhat.hemal,yhat.melanoma,yhat.pancreas,yhat.ovary)
#save(stk11_kras_tp53_genexp_yhats,file="/home/cferte/RESULTS/MEKi/GLOBAL_MODEL/ROBJECTS/stk11_kras_tp53_genexp_yhats.Rda")

# all mutations + stk11 kras tp53 mutation models
#stk11_kras_tp53_mutations_yhats <- list(yhat.all,yhat.nsclc,yhat.breast,yhat.crc,yhat.hemal,yhat.melanoma,yhat.pancreas,yhat.ovary)
#save(stk11_kras_tp53_mutations_yhats,file="/home/cferte/RESULTS/MEKi/GLOBAL_MODEL/ROBJECTS/stk11_kras_tp53_mutations_yhats.Rda")


#lkb1_mutations_model_yhats <- list(yhat.all,yhat.nsclc,yhat.breast,yhat.crc,yhat.hemal,yhat.melanoma,yhat.pancreas,yhat.ovary)
#save(lkb1_mutations_model_yhats,file="/home/cferte/RESULTS/MEKi/GLOBAL_MODEL/ROBJECTS/lkb1_mutations_model_yhats.Rda")

#ovarian_mutations_model_yhats <- list(yhat.all,yhat.nsclc,yhat.breast,yhat.crc,yhat.hemal,yhat.melanoma,yhat.pancreas,yhat.ovary)
#save(ovarian_mutations_model_yhats,file="/home/cferte/RESULTS/MEKi/GLOBAL_MODEL/ROBJECTS/ovarian_mutations_model_yhats.Rda")

# global_model_tissue_adjusted_yhats <- list(yhat.all,yhat.nsclc,yhat.breast,yhat.crc,yhat.hemal,yhat.glioma,yhat.melanoma)
# names(global_model_tissue_adjusted_yhats) <-  c("yhat.all","yhat.nsclc","yhat.breast","yhat.crc","yhat.hemal","yhat.glioma","yhat.melanoma")
# save(global_model_tissue_adjusted_yhats,file="/home/cferte/RESULTS/MEKi/GLOBAL_MODEL_TISSUE_ADJUSTED/ROBJECTS/global_model_tissue_adjusted_yhats.Rda")


# # which are the most frequently selected mutations ? (linear aggregation)
# abc <- matrix(data=NA,nrow=nrow(selected[[1]]),ncol=length(selected))
# colnames(abc) <- 1:N
# rownames(abc) <- rownames(global.matrix)
# for(i in 1:N){abc[,i] <- as.numeric(selected[[i]])}
# abc[abc!=0] <- 1
# abc <- apply(abc,1,function(x){sum(abs(x))})
# names(abc) <-rownames(global.matrix) 
# hist(abc,col="red",breaks=50)
# foo <- names(sort(abc,decreasing=TRUE)[1:50])
# paste(gsub(pattern="_mut",replacement="",x=foo),collapse=" ")


