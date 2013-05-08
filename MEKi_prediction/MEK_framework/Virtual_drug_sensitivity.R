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
fit.original <- fit

#####################################################################################################################
# generate bootstrapped sparse (lasso) models using the gistic and the mutations calls of the same data that predict the vds
#####################################################################################################################

# load the gistic calls
foo  <-  loadEntity("syn1687610")
val_gistic <- read.table(list.files(foo$cacheDir,full.names=TRUE),row.names=1,comment="",quote="",sep="\t",header=TRUE)
colnames(val_gistic) <- substr(colnames(val_gistic),1,12)
rm(foo)
sample.idx <- intersect(colnames(val_gistic),names(vds))

table(is.na(val_gistic))

# feature selection on the genes using KEGG:
gsets <- loadEntity("syn1679661")
gsets_all <- gsets$objects$gsets
assign(x="gsets",gsets_all[[2]])
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
rm(genes,genes2,new.gsets,gsets)


# train sparse models
require(multicore)

N=150
#models <- 0

i <- 0


PARAL <- mclapply(X=1:N,FUN=function(x){
  print(i)
  i <- 1+1
  
train <- sample(sample.idx,replace=TRUE)
cv.fit <- cv.glmnet(x=t(val_gistic[gene.idx,train]),y=vds[train],alpha=1)
fit <-glmnet(x=t(val_gistic[gene.idx,train]),y=vds[train],alpha=1,lambda=cv.fit$lambda.1se)
  
return(list(fit))},mc.set.seed=TRUE,mc.cores=6)


# return the most selected features
abc <- matrix(data=NA, nrow=length(gene.idx),ncol=N)
rownames(abc) <- gene.idx
abc <- sapply(1:N, function(x){
  foo <- PARAL[[x]][[1]]
  abc[,x] <- as.numeric(foo$beta)
})

top.features <- apply(abc,1,function(x){mean(x,na.rm=TRUE)})
names(top.features) <- gene.idx
order <- names(sort(abs(top.features),decreasing=TRUE))[1:200]
top.features[order]
