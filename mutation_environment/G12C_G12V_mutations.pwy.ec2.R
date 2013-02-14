# charles ferté & justin guinney
# feb 8th 2013
# explore the dependencies that exist with specific RAS mutations with other mutations
# (make it on ec2 !)

library(multicore)
library(synapseClient)
synapseLogin(username="charles.ferte@sagebase.org",password="charles")

# load the data

gsets <- loadEntity("syn1679661")
gsets_all <- gsets$objects$gsets
assign(x="gsets",gsets_all[[1]])
pid_gsets <- gsets[grepl("^PID_",names(gsets))]

luad_all <- loadEntity("syn1676707")
luad_all <- luad_all$objects$luad_data
luad_mut <- assign(names(luad_all)[2],luad_all[[2]])
luad_kras <- assign(names(luad_all)[3],luad_all[[3]])

# compute the distribution of the gene mutations per gene in order to weight the importance of each gene
distri.mut <-apply(luad_mut,1,sum)/ncol(luad_mut)
hist(distri.mut,breaks=200, main="distribution of the number of mutations per gene \n in the luad population (n=401)",xlab="number of times a gene is mutated", col="orange")
freq.weight <- log(1/distri.mut)
freq.weight <- (freq.weight-min(freq.weight))/max(freq.weight-min(freq.weight))
hist(freq.weight,breaks=80, main="distribution of the weights across all genes \n (weights are ponderated on the frequency  \nof the mutations per gene in the entire luad population)",col="orange")

# assign a "good weight" (freq.weight =1) on the cancer related genes
cancerIdxs <- c()
cancer.pwy <- names(gsets)[grep(pattern="NON_SMALL_CELL_LUNG_CANCER",x=names(gsets))]
for(i in cancer.pwy){
cancerIdxs <- unique(c(cancerIdxs,unlist(gsets[i])))
}
freq.weight[names(freq.weight) %in% cancerIdxs] <- 1

# adjust for the numbe of mutations per Mb / per sample
sample.mut <- colSums(luad_mut)*1e+06/nrow(luad_mut)
plot(sort(sample.mut))
plot(sort(log10(sample.mut),decreasing=TRUE),ylab="overall rate of mutations (log10 number per Mb)",main="overall mutation rate per sample",xlab="n = 401 samples",pch=20)

boxplot(log10(sample.mut)~luad_kras,las=2,ylab="overall rate of mutations (log10 number per Mb)",main="overall mutation rate per KRAS mutation")
stripchart(log10(sample.mut)~luad_kras,col="aquamarine4",add=TRUE,method="jitter",vertical=TRUE,pch=20)

hist(log10(sample.mut), breaks=100,col="orange",xlab="% mutation per sample",main="distribution of the overall mutation rate per sample")

sample.weight <- 1/log10(sample.mut)
plot(sort(sample.weight))

mask <- luad_kras=="G12C" | luad_kras=="G12V"
tmp <- luad_mut[, mask]
tmp <- tmp[-which(rownames(tmp)=="KRAS"),]

factor <- rep("G12C", ncol(tmp))
factor[colnames(tmp) %in% names(luad_kras)[luad_kras %in% "G12V"]] <- "G12V"
f <- factor(factor)


# load the function to convert any gene set into the gene ids of the genes comprised in this gene set
convertGsetToGIdxs <- function(geneVec, gsets){
  lapply(gsets, function(x){
    which(geneVec %in% x)
  })
}

gsetIdxs <- convertGsetToGIdxs(rownames(tmp), gsets)

# compute the test itself 
test.mut.pathways <- function(MUTtbl, gsets, classFactor, countThreshold=1){
  stopifnot(ncol(MUTtbl) == length(classFactor))
  
  sapply(gsets, function(gsetIdxs){
    if(length(gsetIdxs) < 10){ return (NA)}
    cs <- colSums(t(t(MUTtbl[gsetIdxs,])*sample.weight[colnames(MUTtbl)])*freq.weight[gsetIdxs])
    #stat <- factor(cs >= countThreshold )
    abs(mean(cs[classFactor=="G12C"]) - mean(cs[classFactor=="G12V"]))
    #if(length(levels(stat)) < 2){ return(NA) }
    #fisher.test(stat,classFactor,alternative="less")$p.value
    #wilcox.test(cs ~ classFactor)$p.value
  })
}

R<- test.mut.pathways(MUTtbl=tmp, gsets=gsetIdxs, classFactor=f)
names(R) <- names(gsetIdxs)

hist(R,breaks=100)
sort(R,decreasing=TRUE)[1:10]

# run the permutation test and parallelize it
n.permut <- 1000
#R.null <- replicate(n.permut,test.mut.pathways(MUTtbl=tmp, gsets=gsetIdxs, classFactor=factor(factor)[sample(ncol(tmp))]))
R.null <- mclapply(X= c(1:n.permut), FUN= function(X){ test.mut.pathways(MUTtbl=tmp, gsets=gsetIdxs, classFactor=factor(factor)[sample(ncol(tmp))])} , mc.cores = 10,mc.set.seed=TRUE)
tmp.obj <- c()
R.null <- sapply(c(1:length(R.null)),function(x){ tmpo <- cbind(tmp.obj,R.null[[x]])})
epval <- sapply(c(1:length(R)),function(x){sum(R.null[x,]<=R[x])/n.permut})
names(epval) <- names(R)
sort(epval)[1:20]

PWY <- names(sort(epval)[1:20])

####################################################################################################################
# explore the combined pathways that are selected from grat
####################################################################################################################

#PWY <- names(gsets)[grep(pattern="MTOR",x=names(gsets))]

for (i in PWY){

#grat <- rat[c(20,6,30)]
par(oma=c(8,0,2,1))
pwy <- i
pwy1 <- pwy[1]
pwy2 <- pwy[1]
pwy3 <- pwy[1]

intersect(rownames(tmp)[gsetIdxs[[pwy1]]],rownames(tmp)[gsetIdxs[[pwy2]]])
intersect(rownames(tmp)[gsetIdxs[[pwy3]]],rownames(tmp)[gsetIdxs[[pwy1]]])
intersect(rownames(tmp)[gsetIdxs[[pwy3]]],rownames(tmp)[gsetIdxs[[pwy2]]])

g12c.erk.muts <- rbind(rbind(tmp[gsetIdxs[[pwy1]],f=="G12C"],tmp[gsetIdxs[[pwy2]],f=="G12C"]),tmp[gsetIdxs[[pwy3]],f=="G12C"])
g12c.erk.muts <- g12c.erk.muts[-which(duplicated(rownames(g12c.erk.muts))),]
g12c.erk.muts <- g12c.erk.muts[order(rowSums(g12c.erk.muts), decreasing=T),]

g12c.erk.muts <- g12c.erk.muts[,order(g12c.erk.muts[1, ], g12c.erk.muts[2, ], g12c.erk.muts[3, ], g12c.erk.muts[4, ], g12c.erk.muts[5, ], g12c.erk.muts[6, ], g12c.erk.muts[7, ],g12c.erk.muts[8, ], decreasing=T)]  

g12v.erk.muts <- rbind(rbind(tmp[gsetIdxs[[pwy1]],f=="G12V"],tmp[gsetIdxs[[pwy2]],f=="G12V"]),tmp[gsetIdxs[[pwy3]],f=="G12V"])
g12v.erk.muts <- g12v.erk.muts[-which(duplicated(rownames(g12v.erk.muts))),]
g12v.erk.muts <- g12v.erk.muts[,order(g12v.erk.muts[1, ], g12v.erk.muts[2, ], g12v.erk.muts[3, ], g12v.erk.muts[4, ], g12v.erk.muts[5, ], g12v.erk.muts[6, ], g12v.erk.muts[7, ],g12v.erk.muts[8, ], decreasing=T)]

  require(gplots)
  M <- cbind(g12c.erk.muts,g12v.erk.muts)

    vec <- c(rep(x=1,times=length(colnames(g12c.erk.muts))),rep(x=2,times=length(colnames(g12v.erk.muts))))
  heatmap.2(M,Rowv=FALSE,Colv="Rowv",dendrogram="none",scale="none",trace="none",
            colsep = 1:ncol(M), rowsep = 1:nrow(M), sepcolor = 'white', sepwidth=c(0.05, 0.05),
            keysize = 1.25, density.info = 'none',key=FALSE,
            ColSideColors=c("royalblue","orange")[vec],col = colorpanel(2, 'gray60', 'brown3'),main=paste(pwy1,"\n according to the G12C vs G12V status"))
cat("mutations in",pwy1,"that co-occur with G12C mutations are:","\n",names(which(apply(g12c.erk.muts,1,sum)>0)),"\n")
}
#
PWY

