# charles ferté
# sage bionetworks
# mid july 2013

# load the MEK PDX data
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

# # see the latent structure in it
# val.set <- nki_pdx_cqnNormalized
# s <- svd(val.set)
# plot(s$d^2/sum(s$d^2),pch=19)
# plot(s$v[,1],s$v[,2],pch=19)

# restrict to the treated samples
mask <- paste("PDX_",c(1:19),sep="")
val.set <- nki_pdx_cqnNormalized
val.set <- val.set[,mask]

# get a vector of the best response
response$X0 <- NULL
best.res <- apply(response[,2:7],1,function(x){min(x,na.rm=TRUE)})
rm(mask,tmp)

# differential expression
require("limma")
design <- model.matrix(~best.res)
fit <- eBayes(lmFit(val.set,design))
hist(fit$p.value[,2],breaks=50)