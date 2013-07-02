
library(synapseClient)
library(cgdsr)
library(glmnet)
library(parallel)
library(ROCR)
library(hgu133plus2.db)

# function combine probe 2 gene
combine_probes_2_gene <- function(expr, genes, method="svd"){
  
  if(is.list(genes)) genes <- unlist(genes)
  
  stopifnot(dim(expr)[1] ==  length(genes))
  ugenes <- unique(genes)
  ugenes <- sort(ugenes[!is.na(ugenes)])
  M <- matrix(NaN, ncol=dim(expr)[2],nrow=length(ugenes),
              dimnames=list(ugenes, colnames(expr)))
  
  for(gene in ugenes){
    sub.expr <- as.matrix(expr[which(genes == gene),])
    if(dim(sub.expr)[2] == 1){
      M[gene,] <- sub.expr
    }else{
      tmp <- svd(sub.expr - rowMeans(sub.expr))$v[,1]
      tmp.c <- mean(cor(tmp, t(sub.expr)))
      #cat(gene," ", tmp.c, "\n")
      multiplier <- ifelse(tmp.c < 0, -1, 1)
      M[gene,] <- tmp * multiplier
    }
  }
  M
}

# load the luad cancer tcga
e <- loadEntity("syn418003")
luad.rnaseq <- read.table(paste(e$cacheDir,e$files,sep="/"),header=TRUE,row.names=1,comment="",quote="",sep="\t")

is.tumor <- as.numeric(gsub("TCGA\\.\\w{2}\\.\\w{4}\\.(\\d{2}).*","\\1",colnames(luad.rnaseq))) < 10
tmp <- luad.rnaseq[, is.tumor]
m <- apply(tmp, 1, mean)
tmp <- tmp[m > 1,]
tmp <- log(tmp + 1)
genes <- gsub("(.*?)\\|.*","\\1",rownames(tmp))
mask <- genes != "?" & !duplicated(genes)
luad.rnaseq.g <- tmp[mask,]
rownames(luad.rnaseq.g) <- genes[mask]
colnames(luad.rnaseq.g) <- gsub("(TCGA\\.\\w{2}\\.\\w{4}).*","\\1", colnames(tmp))

val.set <- luad.rnaseq.g


