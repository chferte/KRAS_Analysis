# Justin Guinney and Charles Ferté
# Sage Bionetworks
# Jan 23 2013

# extract information about the pathways of the mutations that cooccur 
# or ar mutually exclusive with a particular mutation


source("/home/jguinney/projects/h3/analysis/JGLibrary.R")
gsets <- load.gmt.data("/home/jguinney/projects/h3/resources/c2.cp.v3.1.symbols.gmt")
pid_gsets <- gsets[grepl("^PID_",names(gsets))]


test.mut.pathways <- function(MUTtbl, gsets, classFactor, countThreshold=1){
  stopifnot(ncol(MUTtbl) == length(classFactor))
  
  sapply(gsets, function(gset){
    idxs <- rownames(MUTtbl) %in% gset
    if(sum(idxs) < 10){ return (NA)}
    cs <- colSums(MUTtbl[idxs,])
    stat <- as.numeric(cs >= countThreshold )
    if(sum(stat) < 3 | sum(stat)==length(stat)){ return(NA) }
    fisher.test(stat,classFactor)$p.value
  })
}

load("/home/cferte/FELLOW/cferte/KRAS_Analysis/KRAS_LUAD.RData")
load(file="/home/cferte/FELLOW/cferte/KRAS_Analysis/MATMUT_GENE_LUAD.RData")

mask <- KRAS_LUAD=="G12C" | KRAS_LUAD=="G12V"
tmp <- MATMUT_GENE_LUAD[, mask]
tmp <- tmp[-which(rownames(tmp)=="KRAS"),]

factor <- rep("G12C", ncol(tmp))
factor[colnames(tmp) %in% names(KRAS_LUAD)[KRAS_LUAD %in% "G12V"]] <- "G12V"

R <- test.mut.pathways(MUTtbl=tmp, gsets=gsets, classFactor=factor(factor))
R.null <- replicate(1000,test.mut.pathways(MUTtbl=tmp, gsets=gsets, factor(factor)[sample(ncol(tmp))]))

emp.pval <- c()
res <- sapply(names(R),function(x){
  emp.val <- c(emp.pval,length(which(R.null[x,]<R[x]))/1000)
  })

sort(res)[1:50]
