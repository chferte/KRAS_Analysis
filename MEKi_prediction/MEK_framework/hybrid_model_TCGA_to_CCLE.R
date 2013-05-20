# charles ferte
# 9 mai 2013
# Sage bionetworks

################################################################################################
# apply the hybrid model into ccle for validation
################################################################################################

# we want to focus on ntop genes
ntop <- 100

# pull out the ccle cnv data
library(cgdsr)
getCBIO_CCLECalls <- function(genes){
  mycgds = CGDS("http://www.cbioportal.org/public-portal/")
  ccle_muts <- getProfileData(mycgds, genes, "ccle_broad_mutations","ccle_broad_complete")
  tmp <- as.matrix(data.frame(lapply(ccle_muts, as.character), stringsAsFactors=FALSE))
  tmp[tmp=="NaN"] = NA
  mutM <- matrix(!is.na(tmp),nrow=nrow(tmp),dimnames=dimnames(ccle_muts))
  ccle_cna <- getProfileData(mycgds, genes, "ccle_broad_CNA","ccle_broad_complete")
  return (list(muts=mutM, cna=ccle_cna))
}

super.genes <- sapply(strsplit(top.selected.names[1:ntop],split="_"),function(x){x[[3]]})

ccle_cna <- t(getCBIO_CCLECalls(super.genes)[[2]])

# make samples coherent btw ccle mut and ccle cna
tmp <- intersect(colnames(ccle_cna),colnames(ccle_mut))
ccle_cna <- ccle_cna[,tmp]

# make coherent the  genes with ccle mut and ccle cna
tmp <- intersect(rownames(ccle_mut),rownames(ccle_cna))
ccle_cna <- ccle_cna[tmp,]

# transform into a mut_Amp, mut_del, mut, del, amp 
datmut <- ccle_mut[rownames(ccle_cna),colnames(ccle_cna)]
datcnv <- ccle_cna
gene.idx <- rownames(datmut)

mut_amp <- c()
for(i in gene.idx){
  mut_amp <- rbind(mut_amp,ifelse(datmut[i,]==1 & datcnv[i,]==1,1,0))
}
rownames(mut_amp) <- paste("mut_amp_",gene.idx,sep="")

mut_del <- c()
for(i in gene.idx){
  mut_del <- rbind(mut_del,ifelse(datmut[i,]==1 & datcnv[i,]==-1,1,0))
}
rownames(mut_del) <- paste("mut_del_",gene.idx,sep="")

mut_only <- c()
for(i in gene.idx){
  mut_only <- rbind(mut_only,ifelse(datmut[i,]==1 & datcnv[i,]==0,1,0))
}
rownames(mut_only) <- paste("mut_only_",gene.idx,sep="")

del_only <- c()
for(i in gene.idx){
  del_only <- rbind(del_only,ifelse(datmut[i,]==0 & datcnv[i,]==-1,1,0))
}
rownames(del_only) <- paste("del_only_",gene.idx,sep="")

amp_only <- c()
for(i in gene.idx){
  amp_only <- rbind(amp_only,ifelse(datmut[i,]==0 & datcnv[i,]==1,1,0))
}
rownames(amp_only) <- paste("amp_only_",gene.idx,sep="")


a <- rbind(mut_amp,mut_del)
b <- rbind(del_only,mut_only)
global.matrix <- rbind(rbind(a,b),amp_only)

rm(a,b,datmut,datcnv,amp_only,del_only,mut_only,mut_amp,mut_del)

# combine the top aberrations recovered to the gene expression matrix (ccle data)
tmp <- intersect(top.selected.names[1:ntop], rownames(global.matrix))
global.matrix <- global.matrix[tmp,]

# get rid of the NA's
naz <- unique(as.numeric(unlist(apply(global.matrix,2,function(x){which(is.na(x))}))))
tmp <- tmp[-naz]
global.matrix <- global.matrix[tmp,]


#concatenate the moelcular aberrations found in TCGA with the gene expression
ccle_exp <- ccle_exp[,colnames(global.matrix)]

# reduce the genbe expression matrix to the top variant probes
tmp2 <- apply(ccle_exp,1,var)
tmp2 <- names(which(tmp2>quantile(tmp2,probs=.8)))
global.matrix <- rbind(global.matrix,ccle_exp[tmp2,])
rownames(global.matrix)[1:100]

# set up the penalty factor (we don't want to penalize the new "tmp" features)
pen.vec <- c(rep(0,times=length(tmp)),rep(1,times=nrow(global.matrix)-length(tmp)))
rm(tmp,naz)  

################################################################################################
# generate a new model of drug sensitivity in cell lines based on gene expression 
# & the features recovered from the tcga dsv model
################################################################################################

  
N <- 100
PARAL2 <- mclapply(X=1:N,FUN=function(x){
    print(i)
    i <- 1+1
    train <- sample(colnames(global.matrix),replace=TRUE)
    vec.train <-apply(ccle_drug[train,mek.inhib],1,mean)
    cv.fit <- cv.glmnet(t(global.matrix[,train]), y=vec.train,nfolds=3, alpha=.1,penalty.factor=pen.vec)
    fit <- glmnet(x=t(global.matrix[,train]),y=vec.train,alpha=.1,lambda=cv.fit$lambda.1se,penalty.factor=pen.vec)
    return(list(fit,train)) },mc.set.seed=TRUE,mc.cores=6)

mek.cells <- intersect(colnames(global.matrix),mek.cells)
melanoma.mek.cells <- intersect(colnames(global.matrix),melanoma.mek.cells)
hemal.mek.cells <- intersect(colnames(global.matrix),hemal.mek.cells)

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
  train <- PARAL2[[i]][[2]]
  fit <- PARAL2[[i]][[1]]
  val <-mek.cells[-which(mek.cells %in% train)]
  selected <- c(selected,list(fit$beta))
  yhat.all <- c(yhat.all,list(predict(fit, t(global.matrix[,val]),)))
  yhat.breast <- c(yhat.breast,list(predict(fit,t(global.matrix[,breast.mek.cells[-which(breast.mek.cells %in% train)]]))))
  yhat.nsclc <- c(yhat.nsclc,list(predict(fit,t(global.matrix[,nsclc.mek.cells[-which(nsclc.mek.cells %in% train)]]))))
  yhat.crc <- c(yhat.crc,list(predict(fit,t(global.matrix[,crc.mek.cells[-which(crc.mek.cells %in% train)]]))))
  yhat.pancreas <- c(yhat.pancreas,list(predict(fit,t(global.matrix[,pancreas.mek.cells[-which(pancreas.mek.cells %in% train)]]))))
  yhat.ovary <- c(yhat.ovary,list(predict(fit,t(global.matrix[,ovary.mek.cells[-which(ovary.mek.cells %in% train)]]))))
  yhat.melanoma <- c(yhat.melanoma,list(predict(fit,t(global.matrix[,melanoma.mek.cells[-which(melanoma.mek.cells %in% train)]]))))
  yhat.hemal <- c(yhat.hemal,list(predict(fit,t(global.matrix[,hemal.mek.cells[-which(hemal.mek.cells %in% train)]]))))
  print(i)}  

# save it !
washing_machine_lung_model_yhats <- list(yhat.all,yhat.nsclc,yhat.breast,yhat.crc,yhat.hemal,yhat.melanoma,yhat.pancreas,yhat.ovary)
save(washing_machine_lung_model_yhats,file="/home/cferte/RESULTS/MEKi/GLOBAL_MODEL/ROBJECTS/washing_machine_lung_model_yhats.Rda")



