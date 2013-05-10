# charles ferté
# Sage Bionetworks
# 07 mai 2013

#####################################################################################################################
# predict the virtual drug sensitivity (vds) into the val_exp
#####################################################################################################################

vds <- c()
for(i in c(1:N)){
  fit <- PARAL[[i]][[1]]
  vds<- cbind(vds,as.numeric(predict(fit, t(val_exp))))
  print(i)
}

vds <- apply(vds,1,mean)
names(vds) <- substr(colnames(val_exp),1,12)
vds

#####################################################################################################################
# generate bootstrapped sparse (lasso) models using the gistic and the mutations calls of the same data that predict the vds
#####################################################################################################################

# load the gistic calls
foo  <-  loadEntity("syn1687610")
val_gistic <- read.table(list.files(foo$cacheDir,full.names=TRUE),row.names=1,comment="",quote="",sep="\t",header=TRUE)
colnames(val_gistic) <- substr(colnames(val_gistic),1,12)
rm(foo)
sample.idx <- intersect(colnames(val_gistic),names(vds))

# load the mutation data
foo  <-  loadEntity("syn1676707")
foo <- foo$objects$luad_data
val_mut <- foo[[2]]
colnames(val_mut) <- substr(colnames(val_mut),1,12)
colnames(val_mut) <- gsub(pattern="-",replacement=".",x=colnames(val_mut),fixed=TRUE)
rm(foo)
sample.idx <- intersect(colnames(val_mut),sample.idx)

# feature selection on the genes using KEGG:
gsets <- loadEntity("syn1679661")
gsets_all <- gsets$objects$gsets
assign(x="gsets",gsets_all[[1]])
rm(gsets_all)
tmp <- unique(c(grep(pattern="CANCER",names(gsets)), grep(pattern="LEUKEM",names(gsets)), 
         grep(pattern="MELANOM",names(gsets)), grep(pattern="GLIOM",
        names(gsets)),grep(pattern="SCLC",names(gsets)), 
         grep(pattern="CARCIN",names(gsets))))

new.gsets <- gsets[tmp]
genes <- unique(unlist(sapply(1:length(new.gsets),function(x){new.gsets[[x]]})))

genes2 <- read.delim2(file="/home/cferte/cancer_gene_census.txt",header=TRUE)
genes2 <- genes2$Symbol
gene.idx <- unique(c(genes,genes2))
gene.idx <- intersect(gene.idx,rownames(val_gistic))
gene.idx <- intersect(gene.idx,rownames(val_mut))
rm(genes,genes2,new.gsets,gsets)

# transform into a mut_Amp, mut_del, mut, del, amp 
datmut <- val_mut[gene.idx,sample.idx]
datcnv <- val_gistic[gene.idx,sample.idx]

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

# train sparse models
require(multicore)

N=500
#models <- 0

i <- 0


NEW.PARAL <- mclapply(X=1:N,FUN=function(x){
  print(i)
  i <- 1+1
  
train <- sample(sample.idx,replace=TRUE)
cv.fit <- cv.glmnet(x=t(global.matrix[,train]),y=vds[train],alpha=1)
fit <-glmnet(x=t(global.matrix[,train]),y=vds[train],alpha=1,lambda=cv.fit$lambda.1se)
  
return(list(fit))},mc.set.seed=TRUE,mc.cores=6)


# return the most selected features
abc <- matrix(data=NA, nrow=nrow(global.matrix),ncol=N)
abc <- sapply(1:N, function(x){
  foo <- NEW.PARAL[[x]][[1]]
  abc[,x] <- as.numeric(foo$beta)
})

rownames(abc) <- rownames(global.matrix)

top.features <- apply(abc,1,function(x){mean(x,na.rm=TRUE)})
names(top.features) <- rownames(global.matrix)
order <- names(sort(abs(top.features),decreasing=TRUE))[1:20]
top.features[order]


top.incorporated <- apply(abc,1,function(x){length(which(x!=0))})
names(top.incorporated) <- rownames(global.matrix)
sort(top.incorporated,decreasing=TRUE)[1:30]
hist(log(top.incorporated),col="red")



# play with the the selected features


abcd <- abc
abcd[abcd!=0] <- 1
abcd[abcd==0] <- 0

i <- 5
plot(density(abcd[top.selected.names[i],]),main=paste(top.selected.names[i]))

top.selected.names <- names(sort(top.incorporated,decreasing=TRUE))
par(oma=c(2,2,4,8))
M <- abc[top.selected.names[1:150],]
heatmap.2(M,trace="none",col=greenred(10),tracecol=NULL,main="fit$beta values")

top.selected.names <- names(sort(top.incorporated,decreasing=TRUE)[1:50])
par(oma=c(2,2,4,8))
M1 <- abcd[top.selected.names,]
heatmap.2(M1,trace="none",col=greenred(2),tracecol=NULL,main="binary matrix")


