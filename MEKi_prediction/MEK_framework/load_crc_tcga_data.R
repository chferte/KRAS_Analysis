
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

# load the colorectal data
e <- loadEntity("syn1446274")
read <- read.table(paste(e$cacheDir,e$files,sep="/"),header=TRUE,row.names=1,comment="")
e <- loadEntity("syn1446195")
coad <- read.table(paste(e$cacheDir,e$files,sep="/"),header=TRUE,row.names=1,comment="")

all(rownames(read) == rownames(coad))
crc <- as.matrix(cbind(read,coad))
m <- apply(crc, 1, mean)
crc <- crc[m > 1,]
crc <- log(crc + 1)
crc <- crc[,grepl("TCGA",colnames(crc))]
pat.ids <- gsub("(TCGA\\.\\w{2}\\.\\w{4}).*","\\1", colnames(crc))
genes <- gsub("(.*?)\\|.*","\\1",rownames(crc))
mask <- genes != "?"
crc.g <- combine_probes_2_gene(crc[mask,],genes[mask])


