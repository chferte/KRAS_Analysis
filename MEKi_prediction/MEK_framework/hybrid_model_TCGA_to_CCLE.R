# charles ferte
# 9 mai 2013
# Sage bionetworks


# apply the hybrid model into ccle for validation

#load the data

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

super.genes <- sapply(strsplit(top.selected.names[1:150],split="_"),function(x){x[[3]]})

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
