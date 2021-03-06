# charles ferté
# Sage Bionetworks
# 07 mai 2013

#####################################################################################################################
# predict the virtual drug sensitivity (vds) into the val_exp and correlate it with response (best.res)
#####################################################################################################################


normalize_to_X <- function(mean.x, sd.x, Y){
  m.y <- rowMeans(Y)
  sd.y <- apply(Y, 1, sd)
  Y.adj <- (Y - m.y) * sd.x / sd.y  + mean.x 
  Y.adj[sd.y == 0] <- mean.x[sd.y==0]
  Y.adj
}

val.set2 <- normalize_to_X(rowMeans(global.matrix), apply(global.matrix, 1, sd), val.set)

# # see what is the structure of this
 s <- svd(val.set2)
 plot(s$d^2/(sum(s$d^2)))
plot(s$v[,1],s$v[,2],pch=19,col=c(rep(1,times=19),rep(2,times=19)))


# let's remove the three outliers

val.set2 <- val.set2[,which(s$v[,1]<.1)]

vds <- c()
for(i in c(1:N)){
  fit <- PARAL[[i]][[1]]
  vds<- cbind(vds,as.numeric(predict(fit, t(val.set2))))
}

rownames(vds) <- substr(colnames(val.set2),1,12)

vds <- apply(vds,1,median)
names(vds) <- substr(colnames(val.set2),1,12)

plot(vds,pch=19,ylab="tumor response (%)",100*best.res[names(vds)],xlab="mek sensitivity score",xlim=c(1,3),main="tumor response ~ mek sensitivity score")
text(2.5,50,paste("spearman rho=" ,format(cor.test(vds,best.res[names(vds)],method="spearman")$estimate,digits=2)),adj=0)
text(2.5,45,paste("p value=" ,format(cor.test(vds,best.res[names(vds)],method="spearman")$p.value,digits=2)),adj=0)

fit <- smooth.spline(vds,100*best.res[names(vds)])
lines(fit,col="royalblue",lwd=3,lty=2)

#####################################################
# treated vs untreated comparison (on the 38 samples)
#####################################################

# do some clean up becausese we now want to use the n=38 PDX samples
rm(val.set,val.set2,val.set3,response,annot,best.res,
   design,nki_pdx_cqnNormalized,s,tmp,vds,i)

# load the MEK PDX data (n=38)
load("/home/jguinney/projects/AZRasPaper/data~/NKIpdx/NKI_PDX_CQNnormalized.rda")
response <- read.table("/home/jguinney/projects/AZRasPaper/data~/NKIpdx/PDXTumorRegression.txt",header=TRUE)
rownames(response) <- response$Days
annot <- read.table("/home/cferte/NKI_MEK_data/NKI_PDX_annot..txt",header=TRUE)
annot$PDX_ID <- gsub(pattern="-",replacement="_",x=annot$PDX_ID)
rownames(annot) <- paste(annot$ID,annot$tumor_sample_ID,sep="_")
colnames(nki_pdx_cqnNormalized) <- substr(colnames(nki_pdx_cqnNormalized),1,6)
tmp <- intersect(colnames(nki_pdx_cqnNormalized),rownames(annot))
nki_pdx_cqnNormalized <- nki_pdx_cqnNormalized[,tmp]
annot <- annot[tmp,]
colnames(nki_pdx_cqnNormalized) <- annot$PDX_ID

val.set <- normalize_to_X(rowMeans(global.matrix), apply(global.matrix, 1, sd), nki_pdx_cqnNormalized[rownames(global.matrix),])
vds <- c()
for(i in c(1:N)){
  fit <- PARAL[[i]][[1]]
  vds<- cbind(vds,as.numeric(predict(fit, t(val.set))))
}
vds <- apply(vds,1,median)
names(vds) <- substr(colnames(val.set),1,12)


boxplot(vds[paste("PDX_",1:38,sep="")]~ c(rep("naive tumors",times=19),rep("mek resisting tumors",times=19)),
        outline=FALSE,ylab="mek sensitivity score", main="mek sensitivity score ~ treatment status")
pval <- format(wilcox.test(vds[paste("PDX_",1:38,sep="")]~ c(rep("untreated",times=19),rep("treated",times=19)),paired=TRUE)$p.value,digits=2)

stripchart(list("treated"=vds[paste("PDX_",20:38,sep="")],"untreated"=vds[paste("PDX_",1:19,sep="")]),method="jitter",
           pch=19,col=c("red","black"),vertical=TRUE,add=TRUE)
text(1.5,1.1,paste("wilcoxon pairwise test \np value=",pval))

#summary(vds[paste("PDX_",20:38,sep="")])
#summary(vds[paste("PDX_",1:19,sep="")])
# 
# vds.null <- c()
# for(i in c(1:N)){
#   fit <- PARAL[[i]][[1]]
#   vds.null <- cbind(vds.null,as.numeric(predict(fit, t(val.set2[sample(rownames(val.set2),replace=FALSE),]))))
# }
# vds.null <- apply(vds.null,1,mean)
# names(vds.null) <- names(vds)
# boxplot(vds.null[paste("PDX_",1:38,sep="")]~ c(rep("untreated",times=19),rep("treated",times=19)))
# wilcox.test(vds.null[paste("PDX_",1:38,sep="")]~ c(rep("untreated",times=19),rep("treated",times=19)),paired=TRUE)
# 
# # load the E-MEXP-3557 xenografts data 
# load("/home/jguinney/projects/AZRasPaper/data~/MEXP3557/MEXP3557_eset.rda")
# phenoData(eset)
# eset@phenoData@data
# 
# # plot 1:19 and 20 38 in differet colors
# nTT <- paste("PDX_",1:19,sep="")
# TT <- names(vds)[names(vds)!=nTT]
# nTT
# vec.res <- response$X28
# names(vec.res) <- response$Days
# tmp <- names(sort(vec.res[names(vds)]))
# 
# plot(100*vec.res[tmp],vds[tmp],pch=19,type="p",
#      xlab=" Day 21 response (%)", ylab="predicted MEK sensitivity score",cex=1.2)
# q <- lowess(100*vec.res[tmp],vds[tmp])
# lines(q$x,q$y,col="royalblue",lwd=4,lty=4,type="l")
# cor <- cor.test(vec.res[names(vds)],vds,method="spearman",pch=19)
# text(x=30,y=1.5,labels=paste("p.value:",format(cor$p.value,digit=3),sep=" "),adj=0)
# text(x=30,y=1.6,labels=paste("Spearman rho:",format(cor$estimate,digit=3),sep=" "),adj=0)
# 
# 
# # try to aggregate all the correlations from day 1 to day 
# par(mfrow=c(2,3))
# rho <- c()
# pval <- c()
# for(i in 2:ncol(response)){
#   rho <- c(rho,cor.test(response[names(vds),i],vds,method="spearman",pch=19)$estimate)
#   pval <- c(pval,cor.test(response[names(vds),i],vds,method="spearman",pch=19)$p.value)
#   plot(vds,response[names(vds),i], pch=19,
#        xlab="predited MEK sensitivity", ylab="tumor response")
# }
# # save the vds (virtual drug sensitivity)
# #save(vds,file="/home/cferte/RESULTS/vds_lung_mut.Rda")
# #save(vds,file="/home/cferte/RESULTS/vds_crc_mut.Rda")
# #save(vds,file="/home/cferte/RESULTS/vds_breast_mut.Rda")
# #save(vds,file="/home/cferte/RESULTS/vds_laml_mut.Rda")
# 
# #####################################################################################################################
# # generate bootstrapped sparse (lasso) models using the gistic and the mutations calls of the same data that predict the vds
# #####################################################################################################################
# 
# # load the gistic calls
# foo  <-  loadEntity("syn1687610")
# val_gistic <- read.table(list.files(foo$cacheDir,full.names=TRUE),row.names=1,comment="",quote="",sep="\t",header=TRUE)
# colnames(val_gistic) <- substr(colnames(val_gistic),1,12)
# rm(foo)
# sample.idx <- intersect(colnames(val_gistic),names(vds))
# 
# # load the mutation data
# foo  <-  loadEntity("syn1676707")
# foo <- foo$objects$luad_data
# val_mut <- foo[[2]]
# colnames(val_mut) <- substr(colnames(val_mut),1,12)
# colnames(val_mut) <- gsub(pattern="-",replacement=".",x=colnames(val_mut),fixed=TRUE)
# rm(foo)
# sample.idx <- intersect(colnames(val_mut),sample.idx)
# 
# # feature selection on the genes using KEGG:
# gsets <- loadEntity("syn1679661")
# gsets_all <- gsets$objects$gsets
# assign(x="gsets",gsets_all[[1]])
# rm(gsets_all)
# tmp <- unique(c(grep(pattern="CANCER",names(gsets)), grep(pattern="LEUKEM",names(gsets)), 
#                 grep(pattern="MELANOM",names(gsets)), grep(pattern="GLIOM",
#                                                            names(gsets)),grep(pattern="SCLC",names(gsets)), 
#                 grep(pattern="CARCIN",names(gsets))))
# 
# new.gsets <- gsets[tmp]
# genes <- unique(unlist(sapply(1:length(new.gsets),function(x){new.gsets[[x]]})))
# 
# genes2 <- read.delim2(file="/home/cferte/cancer_gene_census.txt",header=TRUE)
# genes2 <- genes2$Symbol
# gene.idx <- unique(c(genes,genes2))
# gene.idx <- intersect(gene.idx,rownames(val_gistic))
# gene.idx <- intersect(gene.idx,rownames(val_mut))
# rm(genes,genes2,new.gsets,gsets)
# 
# # transform into a mut_Amp, mut_del, mut, del, amp 
# datmut <- val_mut[gene.idx,sample.idx]
# datcnv <- val_gistic[gene.idx,sample.idx]
# 
# mut_amp <- c()
# for(i in gene.idx){
#   mut_amp <- rbind(mut_amp,ifelse(datmut[i,]==1 & datcnv[i,]==1,1,0))
# }
# rownames(mut_amp) <- paste("mut_amp_",gene.idx,sep="")
# 
# mut_del <- c()
# for(i in gene.idx){
#   mut_del <- rbind(mut_del,ifelse(datmut[i,]==1 & datcnv[i,]==-1,1,0))
# }
# rownames(mut_del) <- paste("mut_del_",gene.idx,sep="")
# 
# mut_only <- c()
# for(i in gene.idx){
#   mut_only <- rbind(mut_only,ifelse(datmut[i,]==1 & datcnv[i,]==0,1,0))
# }
# rownames(mut_only) <- paste("mut_only_",gene.idx,sep="")
# 
# del_only <- c()
# for(i in gene.idx){
#   del_only <- rbind(del_only,ifelse(datmut[i,]==0 & datcnv[i,]==-1,1,0))
# }
# rownames(del_only) <- paste("del_only_",gene.idx,sep="")
# 
# amp_only <- c()
# for(i in gene.idx){
#   amp_only <- rbind(amp_only,ifelse(datmut[i,]==0 & datcnv[i,]==1,1,0))
# }
# rownames(amp_only) <- paste("amp_only_",gene.idx,sep="")
# 
# 
# a <- rbind(mut_amp,mut_del)
# b <- rbind(del_only,mut_only)
# global.matrix <- rbind(rbind(a,b),amp_only)
# 
# rm(a,b,datmut,datcnv,amp_only,del_only,mut_only,mut_amp,mut_del)
# 
# # train sparse models
# require(multicore)
# 
# N=500
# #models <- 0
# 
# i <- 0
# 
# 
# NEW.PARAL <- mclapply(X=1:N,FUN=function(x){
#   print(i)
#   i <- 1+1
#   
#   train <- sample(sample.idx,replace=TRUE)
#   cv.fit <- cv.glmnet(x=t(global.matrix[,train]),y=vds[train],alpha=1)
#   fit <-glmnet(x=t(global.matrix[,train]),y=vds[train],alpha=1,lambda=cv.fit$lambda.1se)
#   
#   return(list(fit))},mc.set.seed=TRUE,mc.cores=6)
# 
# 
# # return the most selected features
# abc <- matrix(data=NA, nrow=nrow(global.matrix),ncol=N)
# abc <- sapply(1:N, function(x){
#   foo <- NEW.PARAL[[x]][[1]]
#   abc[,x] <- as.numeric(foo$beta)
# })
# 
# rownames(abc) <- rownames(global.matrix)
# 
# top.features <- apply(abc,1,function(x){mean(x,na.rm=TRUE)})
# names(top.features) <- rownames(global.matrix)
# order <- names(sort(abs(top.features),decreasing=TRUE))[1:20]
# top.features[order]
# 
# 
# top.incorporated <- apply(abc,1,function(x){length(which(x!=0))})
# names(top.incorporated) <- rownames(global.matrix)
# sort(top.incorporated,decreasing=TRUE)[1:30]
# hist(log(top.incorporated),col="red")
# 
# 
# 
# # play with the the selected features
# 
# 
# abcd <- abc
# abcd[abcd!=0] <- 1
# abcd[abcd==0] <- 0
# 
# i <- 5
# plot(density(abcd[top.selected.names[i],]),main=paste(top.selected.names[i]))
# 
# top.selected.names <- names(sort(top.incorporated,decreasing=TRUE))
# par(oma=c(2,2,4,8))
# M <- abc[top.selected.names[1:150],]
# heatmap.2(M,trace="none",col=greenred(10),tracecol=NULL,main="fit$beta values")
# 
# top.selected.names <- names(sort(top.incorporated,decreasing=TRUE)[1:50])
# par(oma=c(2,2,4,8))
# M1 <- abcd[top.selected.names,]
# heatmap.2(M1,trace="none",col=greenred(2),tracecol=NULL,main="binary matrix")
# 
# 
# 
# 
# 
# 
