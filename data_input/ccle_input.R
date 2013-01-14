# charles ferté
# sage bionetworks
# November 12th


#load the different packages
options(stringsAsFactors=FALSE)

library(affy)
library(limma)
library(caret)
library(glmnet)
library(snm)
synapseLogin("charles.ferte@sagebase.org","charles")

##################################################################################
# Input ccle moelcular data
##################################################################################
ccle_kras <- loadEntity('syn1443160')
ccle_kras <- ccle_kras$objects$KRAS_CCLE
ccle_kras <- substr(ccle_kras,3,nchar(ccle_kras))
ccle_kras[ccle_kras==""] <- "WT"
ccle_kras <- ifelse(ccle_kras %in% c("G12A","G12S","G13D","G12D","G13C","Q61H","Q61K","Q61L"),"rare",ccle_kras)
table(ccle_kras)

liste <- list.celfiles(path="/home/bmecham/ccle/celFiles",full.names=T)
map <- read.delim("/home/cferte/FELLOW/cferte/KRAS_Project/CCLE/CCLE_map.txt",header=T)
map <- map[ map$CCLE.name %in% names(ccle_kras),]
liste <- liste[match(paste(map$ID,".CEL",sep=""),substr(liste,29,nchar(liste)))]

ccle_EXP <- ReadAffy(filenames=liste)
sampleNames(ccle_EXP) <- map$CCLE.name

# ##################################################################################
# # rma
# ##################################################################################
# ccle_rma <- rma(ccle_EXP)
# # retrieve the feature names 
# 
# library(hgu133plus2.db)
# tmp <- unlist(mget(x=featureNames(ccle_rma),hgu133plus2SYMBOL,ifnotfound=NA))
# 
# combine_probes_2_gene <- function(expr, genes, method="svd"){
#   
#   if(is.list(genes)) genes <- unlist(genes)
#   
#   stopifnot(dim(expr)[1] ==  length(genes))
#   ugenes <- unique(genes)
#   ugenes <- sort(ugenes[!is.na(ugenes)])
#   M <- matrix(NaN, ncol=dim(expr)[2],nrow=length(ugenes),
#               dimnames=list(ugenes, colnames(expr)))
#   
#   for(gene in ugenes){
#     sub.expr <- as.matrix(expr[which(genes == gene),])
#     if(dim(sub.expr)[2] == 1){
#       M[gene,] <- sub.expr
#     }else{
#       tmp <- svd(sub.expr - rowMeans(sub.expr))$v[,1]
#       tmp.c <- mean(cor(tmp, t(sub.expr)))
#       #cat(gene," ", tmp.c, "\n")
#       multiplier <- ifelse(tmp.c < 0, -1, 1)
#       M[gene,] <- tmp * multiplier
#     }
#   }
#   M
# }
# 
# ccle_rma <- combine_probes_2_gene(expr=ccle_rma,genes=tmp)
# colnames(ccle_rma) <- sampleNames(ccle_EXP)
# CCLE_RMA <- ccle_rma
# rm(ccle_rma)
# 
# 
####################################################################################
# perform supervised normalization using snm (and rma summarization)
###################################################################################
rawdata <- ccle_EXP
map$Batch <- as.factor(map$Batch)
n <- length(table(map$Batch))

# the batch variable is a confounder of the gene expression and provides unwanted latent structure to the data
s <- svd(exprs(rawdata))
plot(s$v[,1],s$v[,2],col=rainbow(n)[map$Batch],pch=20,cex=2)
plot(s$v[,3],s$v[,4],col=rainbow(n)[map$Batch],pch=20,cex=2)


#SCANBATCH
a <- rawdata@protocolData@data$ScanDate
a <- substr(a,1,nchar(a)-9)
a <- gsub(pattern="-",replacement="/",x=a)
a <- gsub(pattern="T",replacement="",x=a)
a <- gsub(pattern="2009",replacement="09",x=a)
a <- gsub(pattern="2010",replacement="09",x=a)
n <- length(table(a))
a <- as.factor(a)
# the batch variable is a confounder of the gene expression and provides unwanted latent structure to the data
s <- svd(exprs(rawdata))
plot(s$v[,1],s$v[,2],col=rainbow(n)[a],pch=20,cex=2)
plot(s$v[,3],s$v[,4],col=rainbow(n)[a],pch=20,cex=2)

# we know that KRAS are biological & study variables of interest 
bio.var <- model.matrix(~ ccle_kras)

# map$Batch is a batch variable provided by the ccle group
adj.var <- model.matrix(~ a)

myobject <- log2(pm(rawdata))
snm.fit <- snm(myobject, 
               bio.var=bio.var, 
               adj.var=adj.var, 
               rm.adj=TRUE)
new.expr <- snm.fit$norm.dat
pm(rawdata) <- 2^new.expr
myNormSummarized <- rma(rawdata, background=F, normalize=F)
dim(myNormSummarized)
tmp <- "ccle_snm"
assign(tmp,myNormSummarized)

# retrieve the feature names 
library(hgu133plus2.db)
tmp <- unlist(mget(x=featureNames(ccle_snm),hgu133plus2SYMBOL,ifnotfound=NA))

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

ccle_snm <- combine_probes_2_gene(expr=ccle_snm,genes=tmp)
colnames(ccle_snm) <- sampleNames(rawdata)

# check is the batch variable still induce latent structure in the data
s <- svd(ccle_snm)
plot(s$v[,1],s$v[,2],col=rainbow(n)[map$Batch],pch=20,cex=2)
plot(s$v[,3],s$v[,4],col=rainbow(n)[map$Batch],pch=20,cex=2)

# check the influence of KRAS to the latent structure:
n <-length(table(as.factor(ccle_kras)))
plot(s$v[,1],s$v[,2],col=rainbow(n)[as.factor(ccle_kras)],pch=20,cex=2)
plot(s$v[,3],s$v[,4],col=rainbow(n)[as.factor(ccle_kras)],pch=20,cex=2)

CCLE_SNM <- ccle_snm
rm(ccle_snm)

############################################################################################
# save this in synapse
############################################################################################
# save ccle exp
ccle_snm <- Data(list(name = "CCLE_SNM", parentId = 'syn1337457'))
ccle_snm <- createEntity(ccle_snm)

# add object into the data entity
ccle_snm <- addObject(ccle_snm,CCLE_SNM)

# push the raw data into this entity
ccle_snm <- storeEntity(entity=ccle_snm)


# save ccle rma
ccle_rma <- Data(list(name = "CCLE_RMA", parentId = 'syn1337457'))
ccle_rma <- createEntity(ccle_rma)

# add object into the data entity
ccle_rma <- addObject(ccle_rma,CCLE_RMA)

# push the raw data into this entity
ccle_rma <- storeEntity(entity=ccle_rma)

