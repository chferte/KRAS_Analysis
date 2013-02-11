# charles ferté & justin guinney
# explore the dependencies that exist with specific RAS mutations with other mutations

library(synapseClient)
synapseLogin(username="charles.ferte@sagebase.org",password="charles")

# # gsets <- load.gmt.data("/Volumes/users/jguinney/projects/h3/resources/c2.cp.v3.1.symbols.gmt")
# source("/home/cferte/FELLOW/cferte/KRAS_Analysis/mutation_environment/JGLibrary.R")
# gsets_c2.cp.v3.1 <- load.gmt.data("/home/jguinney/projects/h3/resources/c2.cp.v3.1.symbols.gmt")
# gsets_c2.cgp.v3.1 <- load.gmt.data("/home/jguinney/projects/h3/resources/c2.cgp.v3.1.symbols.gmt")
# gsets_c6.all.v3.1.symbols.gmt <- load.gmt.data("/home/jguinney/projects/h3/resources/c6.all.v3.1.symbols.gmt")
# gsets <- list(gsets_c2.cp.v3.1=gsets_c2.cp.v3.1,gsets_c2.cgp.v3.1 =gsets_c2.cgp.v3.1,gsets_c6.all.v3.1.symbols.gmt=gsets_c6.all.v3.1.symbols.gmt)
# 
# # save this list of gene sets into synapse
# # ...and save it into synapse
# gsets_all <- Data(list(name = "gene_sets", parentId = 'syn1670945'))
# gsets_all <- createEntity(gsets_all)
# 
# # # add object into the data entity
# gsets_all <- addObject(gsets_all,gsets)
# # 
# # # push the raw data into this entity
# gsets <- storeEntity(entity=gsets_all) 

gsets <- loadEntity("syn1679661")
gsets_all <- gsets$objects$gsets
assign(x="gsets",gsets_all[[1]])

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
    stat <- factor(cs >= countThreshold )
    if(length(levels(stat)) < 2){ return(NA) }
    fisher.test(stat,classFactor,alternative="less")$p.value
    #wilcox.test(cs ~ classFactor)$p.value
  })
}

#load("/Volumes/cferte/FELLOW/cferte/KRAS_Analysis/luad_kras.RData")
# load("/Volumes/cferte/FELLOW/cferte/KRAS_Analysis/mutations_LUAD.RData")
# load("/home/cferte/FELLOW/cferte/KRAS_Analysis/luad_kras.RData")
# load("/home/cferte/FELLOW/cferte/KRAS_Analysis/mutations_LUAD.RData")

luad_all <- loadEntity("syn1676707")
luad_all <- luad_all$objects$luad_data
luad_mut <- assign(names(luad_all)[2],luad_all[[2]])
luad_kras <- assign(names(luad_all)[3],luad_all[[3]])

# top pathway in G12C vs. G12V  
# rat<-read.delim("/home/cferte/restricted_pwy4.txt",header=F)
# rat <- as.character(rat$V1)
# rat <- gsub(pattern=" ",replacement="",x=rat)
# paste(rat,collapse=" ")

rat <- "KEGG_TERPENOID_BACKBONE_BIOSYNTHESIS KEGG_BIOSYNTHESIS_OF_UNSATURATED_FATTY_ACIDS KEGG_HEMATOPOIETIC_CELL_LINEAGE KEGG_OLFACTORY_TRANSDUCTION KEGG_ASTHMA BIOCARTA_RANKL_PATHWAY BIOCARTA_IL17_PATHWAY MIPS_60S_RIBOSOMAL_SUBUNIT_CYTOPLASMIC MIPS_ANTI_BHC110_COMPLEX MIPS_DDB2_COMPLEX MIPS_12S_U11_SNRNP MIPS_SMN_COMPLEX MIPS_CENP_A_NAC_CAD_COMPLEX MIPS_POLYCYSTIN_1_MULTIPROTEIN_COMPLEX REACTOME_GENERIC_TRANSCRIPTION_PATHWAY REACTOME_P75NTR_SIGNALS_VIA_NFKB REACTOME_SIGNALING_BY_GPCR REACTOME_METABOLISM_OF_POLYAMINES REACTOME_ADP_SIGNALLING_THROUGH_P2RY1 REACTOME_GPCR_DOWNSTREAM_SIGNALING REACTOME_PD1_SIGNALING REACTOME_SIGNAL_AMPLIFICATION REACTOME_ZINC_TRANSPORTERS REACTOME_THROMBOXANE_SIGNALLING_THROUGH_TP_RECEPTOR REACTOME_ADP_SIGNALLING_THROUGH_P2RY12 REACTOME_ACTIVATION_OF_THE_AP1_FAMILY_OF_TRANSCRIPTION_FACTORS REACTOME_CALNEXIN_CALRETICULIN_CYCLE REACTOME_N_GLYCAN_TRIMMING_IN_THE_ER_AND_CALNEXIN_CALRETICULIN_CYCLE REACTOME_ADVANCED_GLYCOSYLATION_ENDPRODUCT_RECEPTOR_SIGNALING REACTOME_HEMOSTASIS REACTOME_EARLY_PHASE_OF_HIV_LIFE_CYCLE"
rat <- unlist(strsplit(x=rat,split=" "))

mask <- luad_kras=="G12C" | luad_kras=="G12V"
tmp <- luad_mut[, mask]
tmp <- tmp[-which(rownames(tmp)=="KRAS"),]

factor <- rep("G12C", ncol(tmp))
factor[colnames(tmp) %in% names(luad_kras)[luad_kras %in% "G12V"]] <- "G12V"
f <- factor(factor)

mutex <- apply(tmp, 1, function(x){
  x <- factor(x)
  if(length(levels(x)) < 2) { return (NA)}
  fisher.test(x,f,alternative="less")$p.value
})

gsetIdxs <- convertGsetToGIdxs(rownames(tmp), gsets)
R <- test.mut.pathways(MUTtbl=tmp, gsets=gsetIdxs, classFactor=f)

# run permutation test
R.null <- replicate(100,test.mut.pathways(MUTtbl=tmp, gsets=gsetIdxs, classFactor=factor(factor)[sample(ncol(tmp))]))
permuted.ptw <- sapply(names(R)[!is.na(R)],function(x){length(which(R.null[x,]<R[x]))/100})

# restrict on the patways of interest
restricted.pwy <- names(which(permuted.ptw<=.01))

R1 <- test.mut.pathways(MUTtbl=tmp, gsets=gsetIdxs[restricted.pwy], classFactor=f)
identical(R[names(R1)],R1)

# run new permutation test on the restricted pathways of interest
R.null <- replicate(1000,test.mut.pathways(MUTtbl=tmp, gsets=gsetIdxs[restricted.pwy], classFactor=factor(factor)[sample(ncol(tmp))]))
new.permuted.pwy <- sapply(names(R1),function(x){length(which(R.null[x,]<R[x]))/1000})
new.restricted.pwy <- names(which(new.permuted.pwy<=.001))


R2 <- test.mut.pathways(MUTtbl=tmp, gsets=gsetIdxs[new.restricted.pwy], classFactor=f)
identical(R[names(R2)],R2)

# run new permutation test on the restricted pathways of interest
R.null <- replicate(5000,test.mut.pathways(MUTtbl=tmp, gsets=gsetIdxs[new.restricted.pwy], classFactor=factor(factor)[sample(ncol(tmp))]))
new.permuted.pwy1 <- sapply(names(R2),function(x){length(which(R.null[x,]<R[x]))/5000})
new.restricted.pwy1 <- names(which(new.permuted.pwy1<.0002))


R3 <- test.mut.pathways(MUTtbl=tmp, gsets=gsetIdxs[new.restricted.pwy1], classFactor=f)
identical(R[names(R3)],R3)
# run new permutation test on the restricted pathways of interest
R.null <- replicate(10000,test.mut.pathways(MUTtbl=tmp, gsets=gsetIdxs[new.restricted.pwy1], classFactor=factor(factor)[sample(ncol(tmp))]))
new.permuted.pwy2 <- sapply(names(R3),function(x){length(which(R.null[x,]<R[x]))/10000})
new.restricted.pwy2 <- names(which(new.permuted.pwy2<.0001))

R4 <- test.mut.pathways(MUTtbl=tmp, gsets=gsetIdxs[new.restricted.pwy2], classFactor=f)
identical(R[names(R4)],R4)
R.null <- replicate(50000,test.mut.pathways(MUTtbl=tmp, gsets=gsetIdxs[new.restricted.pwy2], classFactor=factor(factor)[sample(ncol(tmp))]))
new.permuted.pwy3 <- sapply(names(R4),function(x){length(which(R.null[x,]<R[x]))/50000})
new.restricted.pwy3 <- names(which(new.permuted.pwy3<.00001))

# run new permutation test on the restricted pathways of interest (n=100 000)
R5 <- test.mut.pathways(MUTtbl=tmp, gsets=gsetIdxs[new.restricted.pwy3], classFactor=f)
identical(R[names(R5)],R5)
R.null <- replicate(100000,test.mut.pathways(MUTtbl=tmp, gsets=gsetIdxs[new.restricted.pwy3], classFactor=factor(factor)[sample(ncol(tmp))]))
new.permuted.pwy4 <- sapply(names(R5),function(x){length(which(R.null[x,]<R[x]))/100000})
new.restricted.pwy4 <- names(which(new.permuted.pwy4<1e-05))
new.permuted.pwy4[new.restricted.pwy4]

# run new permutation test on the restricted pathways of interest (n=1000 000)
new.restricted.pwy4 <- rat
R6<- test.mut.pathways(MUTtbl=tmp, gsets=gsetIdxs[new.restricted.pwy4], classFactor=f)
identical(R[names(R6)],R6)
R.null <- replicate(1000000,test.mut.pathways(MUTtbl=tmp, gsets=gsetIdxs[new.restricted.pwy4], classFactor=factor(factor)[sample(ncol(tmp))]))
new.permuted.pwy5 <- sapply(names(R6[which(!is.na(R6))]),function(x){length(which(R.null[x,]<R[x]))/1000000})
new.restricted.pwy5 <- names(which(new.permuted.pwy5<1e-06))
new.permuted.pwy5[new.restricted.pwy5]


system.time(replicate(10,test.mut.pathways.original(MUTtbl=tmp, gsets=gsets, classFactor=factor(factor))))

# explore the ERK pahway
par(oma=c(8,2,2,2))
pwy <- "BIOCARTA_ERK_PATHWAY"
g12c.erk.muts <- tmp[gsetIdxs[[pwy]],f=="G12C"]
g12v.erk.muts <- tmp[gsetIdxs[[pwy]],f=="G12V"]

wilcox.test(apply(g12c.erk.muts,2,sum),apply(g12v.erk.muts,2,sum))
#boxplot(list(G12C=apply(g12c.erk.muts,2,sum),G12V=apply(g12v.erk.muts,2,sum)), outline=FALSE)
#stripchart(list(G12C=apply(g12c.erk.muts,2,sum),G12V=apply(g12v.erk.muts,2,sum)),method="jitter",pch=20,vertical=TRUE,add=TRUE,col="royalblue")

library(gplots)
M <- cbind(g12c.erk.muts,g12v.erk.muts)
vec <- c(rep(x=1,times=length(colnames(g12c.erk.muts))),rep(x=2,times=length(colnames(g12v.erk.muts))))
heatmap.2(M,Rowv=FALSE,Colv="Rowv",scale="none",trace="none",ColSideColors=c("royalblue","red")[vec],col=greenred(2))
title(main=paste(pwy,"~ G12C vs G12V status"),outer=TRUE)

R[pwy]

# explore the SMN pahway
par(oma=c(8,2,2,2))
pwy <- rat[11]
pwy1 <- rat[12]
pwy2 <- rat[8]


g12c.erk.muts <- rbind(rbind(tmp[gsetIdxs[[pwy]],f=="G12C"],tmp[gsetIdxs[[pwy1]],f=="G12C"]),tmp[gsetIdxs[[pwy2]],f=="G12C"])
g12v.erk.muts <- rbind(rbind(tmp[gsetIdxs[[pwy]],f=="G12V"],tmp[gsetIdxs[[pwy1]],f=="G12V"]),tmp[gsetIdxs[[pwy2]],f=="G12V"])
g12c.erk.muts <- g12c.erk.muts[-which(duplicated(rownames(g12c.erk.muts))),]
g12v.erk.muts <- g12v.erk.muts[-which(duplicated(rownames(g12v.erk.muts))),]

wilcox.test(apply(g12c.erk.muts,2,sum),apply(g12v.erk.muts,2,sum))

library(gplots)
M <- cbind(g12c.erk.muts,g12v.erk.muts)
vec <- c(rep(x=1,times=length(colnames(g12c.erk.muts))),rep(x=2,times=length(colnames(g12v.erk.muts))))
heatmap.2(M,Rowv=FALSE,Colv="Rowv",scale="none",trace="none",ColSideColors=c("royalblue","red")[vec],col=greenred(2),main=paste(pwy,"~ G12C vs G12V status"))

R[pwy]

# explore the GPCR pahway
par(oma=c(8,2,2,2))
pwy <- rat[22]
pwy1 <- rat[24]
pwy2 <- rat[25]
pwy3 <- rat[19]


g12c.erk.muts <- rbind(rbind(rbind(tmp[gsetIdxs[[pwy]],f=="G12C"],tmp[gsetIdxs[[pwy1]],f=="G12C"]),tmp[gsetIdxs[[pwy2]],f=="G12C"]),tmp[gsetIdxs[[pwy3]],f=="G12C"])
g12v.erk.muts <- rbind(rbind(rbind(tmp[gsetIdxs[[pwy]],f=="G12V"],tmp[gsetIdxs[[pwy1]],f=="G12V"]),tmp[gsetIdxs[[pwy2]],f=="G12V"]),tmp[gsetIdxs[[pwy3]],f=="G12V"])
g12c.erk.muts <- g12c.erk.muts[-which(duplicated(rownames(g12c.erk.muts))),]
g12v.erk.muts <- g12v.erk.muts[-which(duplicated(rownames(g12v.erk.muts))),]

wilcox.test(apply(g12c.erk.muts,2,sum),apply(g12v.erk.muts,2,sum))

library(gplots)
M <- cbind(g12c.erk.muts,g12v.erk.muts)
vec <- c(rep(x=1,times=length(colnames(g12c.erk.muts))),rep(x=2,times=length(colnames(g12v.erk.muts))))
heatmap.2(M,Rowv=FALSE,dendrogram="none",Colv="Rowv",scale="none",trace="none",ColSideColors=c("royalblue","red")[vec],col=greenred(2),main=paste(pwy,"~ G12C vs G12V status"))

R[pwy]

# explore the PIK3CA pahway
par(oma=c(8,2,2,2))
pwy <- "REACTOME_PI3K_AKT_ACTIVATION"
pwy1 <- "BIOCARTA_MTOR_PATHWAY"
pwy2 <- "BIOCARTA_MTOR_PATHWAY"


g12c.erk.muts <- rbind(rbind(tmp[gsetIdxs[[pwy]],f=="G12C"],tmp[gsetIdxs[[pwy1]],f=="G12C"]),tmp[gsetIdxs[[pwy2]],f=="G12C"])
g12v.erk.muts <- rbind(rbind(tmp[gsetIdxs[[pwy]],f=="G12V"],tmp[gsetIdxs[[pwy1]],f=="G12V"]),tmp[gsetIdxs[[pwy2]],f=="G12V"])
g12c.erk.muts <- g12c.erk.muts[-which(duplicated(rownames(g12c.erk.muts))),]
g12v.erk.muts <- g12v.erk.muts[-which(duplicated(rownames(g12v.erk.muts))),]

wilcox.test(apply(g12c.erk.muts,2,sum),apply(g12v.erk.muts,2,sum))

library(gplots)
M <- cbind(g12c.erk.muts,g12v.erk.muts)
vec <- c(rep(x=1,times=length(colnames(g12c.erk.muts))),rep(x=2,times=length(colnames(g12v.erk.muts))))
heatmap.2(M,Rowv=FALSE,Colv="Rowv",dendrogram="none",scale="none",trace="none",ColSideColors=c("royalblue","red")[vec],col=greenred(2),main=paste(pwy,"~ G12C vs G12V status"))


R[pwy]
R[pwy1]

####################################################################################################################
####################################################################################################################

# explore the PD1 pahway
par(oma=c(8,2,2,2))
pwy <- rat[21]
pwy1 <- rat[21]
pwy2 <- rat[21]


g12c.erk.muts <- rbind(rbind(tmp[gsetIdxs[[pwy]],f=="G12C"],tmp[gsetIdxs[[pwy1]],f=="G12C"]),tmp[gsetIdxs[[pwy2]],f=="G12C"])
g12v.erk.muts <- rbind(rbind(tmp[gsetIdxs[[pwy]],f=="G12V"],tmp[gsetIdxs[[pwy1]],f=="G12V"]),tmp[gsetIdxs[[pwy2]],f=="G12V"])
g12c.erk.muts <- g12c.erk.muts[-which(duplicated(rownames(g12c.erk.muts))),]
g12v.erk.muts <- g12v.erk.muts[-which(duplicated(rownames(g12v.erk.muts))),]

wilcox.test(apply(g12c.erk.muts,2,sum),apply(g12v.erk.muts,2,sum))

library(gplots)
M <- cbind(g12c.erk.muts,g12v.erk.muts)
vec <- c(rep(x=1,times=length(colnames(g12c.erk.muts))),rep(x=2,times=length(colnames(g12v.erk.muts))))
heatmap.2(M,Rowv=FALSE,Colv="Rowv",dendrogram="none",scale="none",trace="none",ColSideColors=c("royalblue","red")[vec],col=greenred(2),main=paste(pwy,"~ G12C vs G12V status"))


R[pwy]
R[pwy1]

####################################################################################################################
# explore the RANKL TF via AP1
par(oma=c(8,2,2,2))
pwy <- rat[6]
pwy1 <- rat[26]
pwy2 <- rat[21]


g12c.erk.muts <- rbind(rbind(tmp[gsetIdxs[[pwy]],f=="G12C"],tmp[gsetIdxs[[pwy1]],f=="G12C"]),tmp[gsetIdxs[[pwy2]],f=="G12C"])
g12v.erk.muts <- rbind(rbind(tmp[gsetIdxs[[pwy]],f=="G12V"],tmp[gsetIdxs[[pwy1]],f=="G12V"]),tmp[gsetIdxs[[pwy2]],f=="G12V"])
g12c.erk.muts <- g12c.erk.muts[-which(duplicated(rownames(g12c.erk.muts))),]
g12v.erk.muts <- g12v.erk.muts[-which(duplicated(rownames(g12v.erk.muts))),]

wilcox.test(apply(g12c.erk.muts,2,sum),apply(g12v.erk.muts,2,sum))

library(gplots)
M <- cbind(g12c.erk.muts,g12v.erk.muts)
vec <- c(rep(x=1,times=length(colnames(g12c.erk.muts))),rep(x=2,times=length(colnames(g12v.erk.muts))))
heatmap.2(M,Rowv=FALSE,Colv="Rowv",dendrogram="none",scale="none",trace="none",ColSideColors=c("royalblue","red")[vec],col=greenred(2),main=paste(pwy,"~ G12C vs G12V status"))


R[pwy]
R[pwy1]

####################################################################################################################
####################################################################################################################

# explore all the significant pathways from rat
for(i in 1:length(rat))

{
par(oma=c(8,2,2,2))
pwy <- rat[i]
pwy1 <- pwy
pwy2 <- pwy


g12c.erk.muts <- rbind(rbind(tmp[gsetIdxs[[pwy]],f=="G12C"],tmp[gsetIdxs[[pwy1]],f=="G12C"]),tmp[gsetIdxs[[pwy2]],f=="G12C"])
g12v.erk.muts <- rbind(rbind(tmp[gsetIdxs[[pwy]],f=="G12V"],tmp[gsetIdxs[[pwy1]],f=="G12V"]),tmp[gsetIdxs[[pwy2]],f=="G12V"])
g12c.erk.muts <- g12c.erk.muts[-which(duplicated(rownames(g12c.erk.muts))),]
g12v.erk.muts <- g12v.erk.muts[-which(duplicated(rownames(g12v.erk.muts))),]

print(i)
wilcox.test(apply(g12c.erk.muts,2,sum),apply(g12v.erk.muts,2,sum))

library(gplots)
M <- cbind(g12c.erk.muts,g12v.erk.muts)
vec <- c(rep(x=1,times=length(colnames(g12c.erk.muts))),rep(x=2,times=length(colnames(g12v.erk.muts))))
heatmap.2(M,Rowv=FALSE,Colv="Rowv",dendrogram="none",scale="none",trace="none",ColSideColors=c("royalblue","red")[vec],col=greenred(2),main=paste(pwy,"~ G12C vs G12V status"))


R[pwy]
}

####################################################################################################################
####################################################################################################################
N <- nrow(g12v.erk.muts)
M <- matrix(NA,nrow=N,ncol=N)
for(i in 1:(N-1)){
  for(j in (i+1):N){
    a <- factor(g12v.erk.muts[i,])
    b <- factor(g12v.erk.muts[j,])
    if(length(levels(a)) < 2 | length(levels(b)) < 2) { next; }
    M[i,j] <- fisher.test(a,b,alternative="greater")$p.value
  }
}

g12c.erk.muts[c(26,16),]

# compute the permutation test
R.null <- replicate(10000,test.mut.pathways(MUTtbl=tmp, gsets=gsets, factor(factor)[sample(ncol(tmp))]))