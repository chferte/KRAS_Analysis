## charles ferté,MD
## Sage Bionetworks 
## "2012-dec-06"

## KRAS project LUAD
## input the RNASeq data

require(synapseClient)
require(Biobase)
## synapse Login
synapseLogin("charles.ferte@sagebase.org","charles")

### magic option
options(stringsAsFactors=FALSE)

#######################################################################################################
# 1. read the RNASeq information file of TCGA LUAD normalized by Justin Guinney  and load the and see the intersect
#######################################################################################################

source("/home/cferte/FELLOW/cferte/KRAS_Analysis/data_input/TCGA_LUAD_input.R")

load("/home/cferte/FELLOW/cferte/KRAS_Analysis/MUT_TYPE_LUAD.RData")
load("/home/cferte/FELLOW/cferte/KRAS_Analysis/mutations_LUAD.RData")

# create the KRAS_LUAD vector
tmp <- MATMUT_LUAD[grep(pattern="KRAS", rownames(MATMUT_LUAD)),]
rownames(tmp) <- substr(x=rownames(tmp),8,nchar(rownames(tmp)))
tmp1 <- c()
for(i in colnames(tmp)){
  tmp1 <- c(tmp1,ifelse(sum(tmp[,i])==0,"WT",rownames(tmp)[which(tmp[,i]==1)]))
}
names(tmp1) <- colnames(tmp)
KRAS_LUAD <- tmp1
rm(tmp,tmp1)

# #######################################################################################################
# # 2. create a KRAS_LUAD variable 
# #######################################################################################################
# 
# G12C <- as.numeric(MATMUT_LUAD["KRAS_rs121913530_C_A_chr12:25398285-25398285",])
# G12V <- as.numeric(MATMUT_LUAD["KRAS_rs121913529_C_A_chr12:25398284-25398284",])
# G12D <- as.numeric(MATMUT_LUAD["KRAS_rs121913529_C_T_chr12:25398284-25398284",])
# G13D <- as.numeric(MATMUT_LUAD["KRAS_rs112445441_C_T_chr12:25398281-25398281",])
# tmp <- setdiff(c(rownames(MATMUT_LUAD)[grep("KRAS",rownames(MATMUT_LUAD))]),
#                c("KRAS_rs121913530_C_A_chr12:25398285-25398285",     #G12C
#                  "KRAS_rs121913529_C_A_chr12:25398284-25398284",     #G12V
#                  "KRAS_rs121913529_C_T_chr12:25398284-25398284",     #G12D
#                  "KRAS_rs112445441_C_T_chr12:25398281-25398281"))    #G13D
# RARE <- as.numeric(apply(MATMUT_LUAD[tmp,],2,sum))
# names(RARE) <- colnames(MATMUT_LUAD)
# 
# foo <- rbind(G12C,G12V)
# foo <- rbind(foo,G12D)
# foo <- rbind(foo,G13D)
# foo <- rbind(foo,RARE)
# rownames(foo) <- c("G12C","G12V","G12D","G13D","rare")
# foo2 <- as.data.frame(t(foo))
# WT <- ifelse(apply(foo2,1,sum)==0,1,0)
# foo2$WT <- WT
# KRAS_LUAD <- unlist(apply(foo2,1,function(x){colnames(foo2)[which(x==1)]}))
# rm(G12C,G12D,G12V,G13D,foo,foo2,RARE,WT,tmp)
# table(KRAS_LUAD)
# save(file="/home/cferte/FELLOW/cferte/KRAS_Analysis/KRAS_LUAD.RData",KRAS_LUAD)
# 
# ###################################################################################################
# # 3. plot the overall mutation rate in LUAD (y axis = mutations per megabase, x axis= samples )
# ###################################################################################################
# 
# MUTRATE <- apply(MATMUT_LUAD,2,sum)
# par(mfrow=c(1,1))
# x <- sort(MUTRATE,decreasing=TRUE)*95859/1000000
# summary(x)
# #plot(x,pch=20,cex=.01,xlab="249 LUAD TCGA samples",
# #ylab="Overall number of mutations per Mb")
# #abline(h=quantile(x=x,probs=.15),lty=2,col="red")
# #abline(h=quantile(x=x,probs=.75),lty=2,col="red")
# 
# plot(log10(sort(MUTRATE,decreasing=TRUE)*94771/1000000),pch=20,cex=.5,xlab="249 LUAD TCGA samples",
#      ylab="Overall number of mutations per Mb (log10 scale)",yaxt="n")
# axis(side=2,at=c(.5,1,1.5,2,2.5,3),labels=c(3,10,30,100,300,1000))
# abline(h=log10(quantile(x=x,probs=.75)),lty=2,col="red")
# abline(h=log10(quantile(x=x,probs=.25)),lty=2,col="red")
# rm(x)


# ###############################################################################################################################
# # 3. compute the mutations that are associated (overlapping or exclusive) with KRAS mutations
# ###############################################################################################################################
# kras <- apply(MATMUT_LUAD[grep("KRAS",rownames(MATMUT_LUAD)),],2,sum)
# kras <- ifelse(kras==0,0,1)
# kras.overlap <- apply(MATMUT_LUAD,1,function(x){ fisher.test(as.numeric(kras),as.numeric(x),alternative="greater")$p.value})
# kras.exclusive <- apply(MATMUT_LUAD,1,function(x){ fisher.test(as.numeric(kras),as.numeric(x),alternative="less")$p.value})
# 
# table(kras.overlap)
# 
# hist(kras.overlap, breaks=30)
# hist(kras.exclusive, breaks=30)
# 
# table(kras.overlap<.01)
# which(kras.overlap<.05)
# table(kras.exclusive<.01)
# which(kras.exclusive<.05)

################################################################################################################################
# 4. create a new matrix per gene instead of per specific mutation
###############################################################################################################################

# gene <- unique(sapply(strsplit(x=rownames(MATMUT_LUAD),split="_"), function(x){x[[1]]}))
# j <- c()
# new.matmut <- matrix(0,nrow=length(gene),ncol=length(colnames(MATMUT_LUAD)))
# rownames(new.matmut) <- gene
# colnames(new.matmut) <- colnames(MATMUT_LUAD)
# 
# for(i in rownames(new.matmut)){
# new.matmut[i,c(names(which(apply(MATMUT_LUAD[grep(pattern=i,x=rownames(MATMUT_LUAD)),],2,sum)!=0)))] <- 1
# }
# 
# MATMUT_GENE_LUAD <- new.matmut
# save(file="/home/cferte/FELLOW/cferte/KRAS_Analysis/MATMUT_GENE_LUAD.RData",MATMUT_GENE_LUAD)


# ###############################################################################################################################
# # 5. compute the overlap between KRAS specific mutations and specific mutations
# ###############################################################################################################################
# j <- which(KRAS_LUAD %in% c("WT","G12C"))
# G12C<- ifelse(KRAS_LUAD[j]=="G12C",1,0)
# k <- which(apply(MATMUT_LUAD[,j],1,sum)>0)
# G12C.overlap <- apply(MATMUT_LUAD[k,j],1,function(x){fisher.test(as.numeric(G12C),as.numeric(x),alternative="greater")$p.value})
# G12C.exclusive <- apply(MATMUT_LUAD[k,j],1,function(x){ fisher.test(as.numeric(G12C),as.numeric(x),alternative="less")$p.value})
# 
# 
# j <- which(KRAS_LUAD %in% c("WT","G12V"))
# G12V <- ifelse(KRAS_LUAD[j]=="G12V",1,0)
# k <- which(apply(MATMUT_LUAD[,j],1,sum)>0)
# G12V.overlap <- apply(MATMUT_LUAD[k,j],1,function(x){ fisher.test(as.numeric(G12V),as.numeric(x),alternative="greater")$p.value})
# G12V.exclusive <- apply(MATMUT_LUAD[k,j],1,function(x){ fisher.test(as.numeric(G12V),as.numeric(x),alternative="less")$p.value})
# 
# 
# # display the results
# 
# hist(G12C.overlap, breaks=30)
# hist(G12C.exclusive, breaks=30)
# table(G12C.overlap<.05)
# table(G12C.exclusive<.05)
# paste(names(which(G12C.overlap<.05)),collapse=" ")
# paste(names(which(G12C.exclusive<.05)),collapse=" ")
# table(G12V.overlap<.05)
# table(G12V.exclusive<.05)
# hist(G12V.overlap, breaks=30)
# hist(G12V.exclusive, breaks=30)
# paste(names(which(G12V.overlap<.05)),collapse=" ")
# paste(names(which(G12V.exclusive<.05)),collapse=" ")
# 
# # create G12Coverlap.rnk object to be analyzed in the GSEA java (pre ranked test)
# tmp <- G12C.overlap
# tmp <- tmp[-c(grep(pattern="KRAS",x=names(tmp)))]
# rnames <- sapply(strsplit(names(tmp),split="_"),function(x){x[[1]]})
# j <- c()
# vec <- as.numeric(rep(1,times=length(unique(rnames))))
# names(vec) <- unique(rnames)
# for(i in names(vec))
# { 
#   j <- which(rnames==i)
#   vec[i] <- min(tmp[j])
# }
# table(vec)
# plot(sort(vec))
# plot(sort(-log10(vec)))
# foo <- as.data.frame(cbind(names(vec),as.numeric(-log10(as.numeric(vec)))))
# write.table(foo,file="G12Coverlap.rnk",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
# 
# 
# 
# # create G12Voverlap.rnk object to be analyzed in the GSEA java (pre ranked test)
# tmp <- G12V.overlap
# tmp <- tmp[-c(grep(pattern="KRAS",x=names(tmp)))]
# rnames <- sapply(strsplit(names(tmp),split="_"),function(x){x[[1]]})
# j <- c()
# vec <- as.numeric(rep(1,times=length(unique(rnames))))
# names(vec) <- unique(rnames)
# for(i in names(vec))
# { 
#   j <- which(rnames==i)
#   vec[i] <- min(tmp[j])
# }
# table(vec)
# plot(sort(vec))
# plot(sort(-log10(vec)))
# foo <- as.data.frame(cbind(names(vec),as.numeric(-log10(as.numeric(vec)))))
# write.table(foo,file="G12Voverlap.rnk",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
# 
# 
# ###############################################################################################################################
# # 7. compute the overlap between KRAS G12C and G12V mutations
# ###############################################################################################################################
# j <- which(KRAS_LUAD %in% c("G12V","G12C"))
# G12C<- ifelse(KRAS_LUAD[j]=="G12C",1,0)
# k <- which(apply(MATMUT_LUAD[,j],1,sum)>0)
# G12C.G12V.overlap <- apply(MATMUT_LUAD[k,j],1,function(x){fisher.test(as.numeric(G12C),as.numeric(x),alternative="greater")$p.value})
# G12C.G12V.exclusive <- apply(MATMUT_LUAD[k,j],1,function(x){ fisher.test(as.numeric(G12C),as.numeric(x),alternative="less")$p.value})
# 
# # create G12C.G12V.overlap.rnk object to be analyzed in the GSEA java (pre ranked test)
# tmp <- G12C.G12V.overlap
# tmp <- tmp[-c(grep(pattern="KRAS",x=names(tmp)))]
# rnames <- sapply(strsplit(names(tmp),split="_"),function(x){x[[1]]})
# j <- c()
# vec <- as.numeric(rep(1,times=length(unique(rnames))))
# names(vec) <- unique(rnames)
# for(i in names(vec))
# { 
#   j <- which(rnames==i)
#   vec[i] <- min(tmp[j])
# }
# table(vec)
# plot(sort(vec))
# plot(sort(-log10(vec)))
# foo <- as.data.frame(cbind(names(vec),as.numeric(-log10(as.numeric(vec)))))
# write.table(foo,file="G12C.G12V.genes.overlap.rnk",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
# 
# # create G12C.G12V.exclusive.rnk object to be analyzed in the GSEA java (pre ranked test)
# tmp <- G12C.G12V.exclusive
# tmp <- tmp[-c(grep(pattern="KRAS",x=names(tmp)))]
# rnames <- sapply(strsplit(names(tmp),split="_"),function(x){x[[1]]})
# j <- c()
# vec <- as.numeric(rep(1,times=length(unique(rnames))))
# names(vec) <- unique(rnames)
# for(i in names(vec))
# { 
#   j <- which(rnames==i)
#   vec[i] <- min(tmp[j])
# }
# table(vec)
# plot(sort(vec))
# plot(sort(-log10(vec)))
# foo <- as.data.frame(cbind(names(vec),as.numeric(-log10(as.numeric(vec)))))
# write.table(foo,file="G12C.G12V.genes.exclusive.rnk",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")


###############################################################################################################################
# 8. compute the overlap between KRAS G12C and G12V mutations AND genes
###############################################################################################################################

# load the MATMUT_GENE_LUAD (mutation matrix per gene)
load(file="/home/cferte/FELLOW/cferte/KRAS_Analysis/MATMUT_GENE_LUAD.RData")
new.matmut <- MATMUT_GENE_LUAD

j <- which(KRAS_LUAD %in% c("WT","G12C"))
G12C<- ifelse(KRAS_LUAD[j]=="G12C",1,0)
k <- names(which(apply(new.matmut[,j],1,sum)>1 & apply(new.matmut[,j],1,sum)<length(j)))
G12C.overlap.genes <- apply(new.matmut[k,j],1,function(x){fisher.test(as.numeric(G12C),as.numeric(x),alternative="greater")$p.value})
G12C.exclusive.genes <- apply(new.matmut[k,j],1,function(x){ fisher.test(as.numeric(G12C),as.numeric(x),alternative="less")$p.value})


j <- which(KRAS_LUAD %in% c("WT","G12V"))
G12V <- ifelse(KRAS_LUAD[j]=="G12V",1,0)
k <- names(which(apply(new.matmut[,j],1,sum)>1 & apply(new.matmut[,j],1,sum)<length(j)))
G12V.overlap.genes <- apply(new.matmut[k,j],1,function(x){fisher.test(as.numeric(G12V),as.numeric(x),alternative="greater")$p.value})
G12V.exclusive.genes <- apply(new.matmut[k,j],1,function(x){ fisher.test(as.numeric(G12V),as.numeric(x),alternative="less")$p.value})

#j <- which(KRAS_LUAD %in% c("G12V"))
G12C <- ifelse(KRAS_LUAD=="G12C",1,0)
k <- names(which(apply(new.matmut,1,sum)>1 & apply(new.matmut,1,sum)<401))
G12C.REST.overlap.genes <- apply(new.matmut[k,],1,function(x){fisher.test(as.numeric(G12C),as.numeric(x),alternative="greater")$p.value})
G12C.REST.exclusive.genes <- apply(new.matmut[k,],1,function(x){ fisher.test(as.numeric(G12C),as.numeric(x),alternative="less")$p.value})




foo <- as.data.frame(cbind(names(G12C.overlap.genes),as.numeric(-log10(as.numeric(G12C.overlap.genes)))))
foo <- foo[-which(foo$V1=="KRAS"),]
write.table(foo,file="G12Cgenesoverlap.rnk",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
foo <- as.data.frame(cbind(names(G12C.exclusive.genes),as.numeric(-log10(as.numeric(G12C.exclusive.genes)))))
foo <- foo[-which(foo$V1=="KRAS"),]
write.table(foo,file="G12Cgenesexclusive.rnk",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
foo <- as.data.frame(cbind(names(G12V.overlap.genes),as.numeric(-log10(as.numeric(G12V.overlap.genes)))))
foo <- foo[-which(foo$V1=="KRAS"),]
write.table(foo,file="G12Vgenesoverlap.rnk",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
foo <- as.data.frame(cbind(names(G12V.exclusive.genes),as.numeric(-log10(as.numeric(G12V.exclusive.genes)))))
foo <- foo[-which(foo$V1=="KRAS"),]
write.table(foo,file="G12Vgenesexclusive.rnk",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")


