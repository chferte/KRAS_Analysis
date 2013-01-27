# Justin Guinney and Charles Ferté
# Sage Bionetworks
# Jan 23 2013

# yield information about the pathways associated with the mutations that cooccur 
# or ar mutually exclusive with a particular mutation

source("/home/jguinney/projects/h3/analysis/JGLibrary.R")
gsets <- load.gmt.data("/home/jguinney/projects/h3/resources/c2.cp.v3.1.symbols.gmt")
pid_gsets <- gsets[grepl("^PID_",names(gsets))]

convertGsetToGIdxs <- function(geneVec, gsets){
  lapply(gsets, function(x){
    which(geneVec %in% x)
  })
}


test.mut.pathways <- function(MUTtbl, gsets, classFactor, countThreshold=1){
  stopifnot(ncol(MUTtbl) == length(classFactor))
  
  sapply(gsets, function(gsetIdxs){
    if(length(gsetIdxs) < 10){ return (NA)}
    cs <- colSums(MUTtbl[gsetIdxs,])
    stat <- as.numeric(cs >= countThreshold )
    if(length(levels(stat)) < 2){ return(NA) }
    wilcox.test(stat ~ classFactor,alternative="less")$p.value
  })
}

load("/home/cferte/FELLOW/cferte/KRAS_Analysis/KRAS_LUAD.RData")
load(file="/home/cferte/FELLOW/cferte/KRAS_Analysis/MATMUT_LUAD.RData")

mask <- KRAS_LUAD=="G12C" | KRAS_LUAD=="G12V"
tmp <- MATMUT_LUAD[, mask]
tmp <- tmp[-which(rownames(tmp)=="KRAS"),]

factor <- rep("G12C", ncol(tmp))
factor[colnames(tmp) %in% names(KRAS_LUAD)[KRAS_LUAD %in% "G12V"]] <- "G12V"

gsetIdxs <- convertGsetToGIdxs(rownames(tmp), gsets)
#system.time(replicate(10,test.mut.pathways(MUTtbl=tmp, gsets=gsetIdxs, classFactor=factor(factor))))
#system.time(replicate(10,test.mut.pathways.original(MUTtbl=tmp, gsets=gsetIdxs, classFactor=factor(factor))))

R <- test.mut.pathways(MUTtbl=tmp, gsets=gsetIdxs, classFactor=factor(factor))
sort(R)[1:10]

# compute the permutation test
R.null <- replicate(10000,test.mut.pathways(MUTtbl=tmp, gsets=gsets, factor(factor)[sample(ncol(tmp))]))


emp.pval <- c()
res <- sapply(names(R),function(x){
  emp.val <- c(emp.pval,length(which(R.null[x,]<R[x]))/10000)
  })

sort(res)[1:10]
