# Charles Ferté
# Sage Bionetworks
# 14 Feb 2012

#############################
# load the data
#############################

# load the mek data from ccle
source("/home/cferte/FELLOW/cferte/KRAS_Analysis/MEKi_prediction/MEK_framework/load_mek_ccle.R")

# load the  cell statsu and probabilities
ccle_probs_status <- loadEntity("syn1709732")
ccle_probs_status <- ccle_probs_status$objects$ccle_probs_status
all.prob <- ccle_probs_status[[1]]
cell.status <- ccle_probs_status[[2]]



#######################################################
# predictive modeling
# we predict the ic50 
# training all cells global matrix without eigengenes and with eigengenes in parallel
#######################################################

# define globalmatrix (without eigengenes)
global.matrix <- rbind(ccle_exp,ccle_cnv,ccle_mut)
rownames(global.matrix) <- c(paste(rownames(ccle_exp),"_exp",sep=""),paste(rownames(ccle_cnv),"_cnv",sep=""),paste(rownames(ccle_mut),"_mut",sep=""))

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

N=100
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
  train <- sample(mek.cells[-which(mek.cells %in% c(hemal.mek.cells,glioma.mek.cells,melanoma.mek.cells))],replace=TRUE)
  
  # balanced model
  #train <- c(sample(x=c(Q1),replace=TRUE,size=220),sample(Q4,replace=TRUE,size=220))
  
  # weighted model on the distribution
  #train <- sample(mek.cells,replace=TRUE)
  
  
  # tissue specific models
  #train <- sample(melanoma.mek.cells,replace=TRUE)
  
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
yhat.glioma <- c()
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
#   #yhat.glioma <- c(yhat.glioma,list(predict(fit,t(global.matrix[,glioma.mek.cells[-which(glioma.mek.cells %in% train)]]))))
#   #yhat.melanoma <- c(yhat.melanoma,list(predict(fit,t(global.matrix[,melanoma.mek.cells[-which(melanoma.mek.cells %in% train)]]))))
#   #yhat.hemal <- c(yhat.hemal,list(predict(fit,t(global.matrix[,hemal.mek.cells[-which(hemal.mek.cells %in% train)]]))))
#   print(i)}  

#####################################################################################################################
# save the objects
#####################################################################################################################

#global_model_yhats <- list(yhat.all,yhat.nsclc,yhat.breast,yhat.crc,yhat.hemal,yhat.glioma,yhat.melanoma)
#save(global_model_yhats,file="/home/cferte/RESULTS/MEKi/GLOBAL_MODEL/ROBJECTS/global_model_yhats.Rda")


# # save it in pure_tissue_models_yhats
# pure_tissue_models_yhats <- c(pure_tissue_models_yhats,list(yhat.melanoma))
# # names as yhat.all,yhat.nsclc,yhat.breast,yhat.crc,yhat.hemal,yhat.glioma,yhat.melanoma)
# names(pure_tissue_models_yhats) <- c("yhat.all","yhat.nsclc","yhat.breast","yhat.crc","yhat.hemal","yhat.glioma","yhat.melanoma")
# #save it on belltown
# save(pure_tissue_models_yhats,file="/home/cferte/RESULTS/MEKi/GLOBAL_MODEL/ROBJECTS/pure_tissue.models.Rda")

# global_model_tissue_adjusted_yhats <- list(yhat.all,yhat.nsclc,yhat.breast,yhat.crc,yhat.hemal,yhat.glioma,yhat.melanoma)
# names(global_model_tissue_adjusted_yhats) <-  c("yhat.all","yhat.nsclc","yhat.breast","yhat.crc","yhat.hemal","yhat.glioma","yhat.melanoma")
# save(global_model_tissue_adjusted_yhats,file="/home/cferte/RESULTS/MEKi/GLOBAL_MODEL_TISSUE_ADJUSTED/ROBJECTS/global_model_tissue_adjusted_yhats.Rda")

# balanced_model_yhats <- list(yhat.all,yhat.nsclc,yhat.breast,yhat.crc,yhat.hemal,yhat.glioma,yhat.melanoma)
# names(balanced_model_yhats) <-  c("yhat.all","yhat.nsclc","yhat.breast","yhat.crc","yhat.hemal","yhat.glioma","yhat.melanoma")
# save(balanced_model_yhats,file="/home/cferte/RESULTS/MEKi/balanced_model_yhats.Rda")

#models without hemal and glioma
#global_carcinoma_model <- list(yhat.all,yhat.nsclc,yhat.breast,yhat.crc,yhat.hemal,yhat.glioma,yhat.melanoma)
#save(global_carcinoma_model,file="/home/cferte/RESULTS/MEKi/GLOBAL_MODEL/ROBJECTS/global_carcinoma_model_yhats.Rda")


