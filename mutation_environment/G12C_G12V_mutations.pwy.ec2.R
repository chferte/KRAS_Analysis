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

# rat is the character list of the top pathways we want to specifically test
rat <- "KEGG_TERPENOID_BACKBONE_BIOSYNTHESIS KEGG_BIOSYNTHESIS_OF_UNSATURATED_FATTY_ACIDS KEGG_OLFACTORY_TRANSDUCTION KEGG_ASTHMA BIOCARTA_RANKL_PATHWAY BIOCARTA_IL17_PATHWAY MIPS_60S_RIBOSOMAL_SUBUNIT_CYTOPLASMIC MIPS_ANTI_BHC110_COMPLEX MIPS_DDB2_COMPLEX MIPS_12S_U11_SNRNP MIPS_SMN_COMPLEX MIPS_CENP_A_NAC_CAD_COMPLEX MIPS_POLYCYSTIN_1_MULTIPROTEIN_COMPLEX REACTOME_GENERIC_TRANSCRIPTION_PATHWAY REACTOME_P75NTR_SIGNALS_VIA_NFKB REACTOME_SIGNALING_BY_GPCR REACTOME_METABOLISM_OF_POLYAMINES REACTOME_ADP_SIGNALLING_THROUGH_P2RY1 REACTOME_GPCR_DOWNSTREAM_SIGNALING REACTOME_PD1_SIGNALING REACTOME_SIGNAL_AMPLIFICATION REACTOME_ZINC_TRANSPORTERS REACTOME_THROMBOXANE_SIGNALLING_THROUGH_TP_RECEPTOR REACTOME_ADP_SIGNALLING_THROUGH_P2RY12 REACTOME_ACTIVATION_OF_THE_AP1_FAMILY_OF_TRANSCRIPTION_FACTORS REACTOME_CALNEXIN_CALRETICULIN_CYCLE REACTOME_N_GLYCAN_TRIMMING_IN_THE_ER_AND_CALNEXIN_CALRETICULIN_CYCLE REACTOME_ADVANCED_GLYCOSYLATION_ENDPRODUCT_RECEPTOR_SIGNALING REACTOME_HEMOSTASIS REACTOME_EARLY_PHASE_OF_HIV_LIFE_CYCLE"
rat <- unlist(strsplit(x=rat,split=" "))

# or do it in that specific pathway
grat <- names(gsets)[grep(pattern="ERK",x=names(gsets))]

# load the functions

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

# start the process !

mask <- luad_kras=="G12C" | luad_kras=="G12V"
tmp <- luad_mut[, mask]
tmp <- tmp[-which(rownames(tmp)=="KRAS"),]

factor <- rep("G12C", ncol(tmp))
factor[colnames(tmp) %in% names(luad_kras)[luad_kras %in% "G12V"]] <- "G12V"
f <- factor(factor)

gsetIdxs <- convertGsetToGIdxs(rownames(tmp), gsets)

# run new permutation test on the restricted pathways of interest
new.restricted.pwy <- rat
n.permut <- 1000

R<- test.mut.pathways(MUTtbl=tmp, gsets=gsetIdxs[new.restricted.pwy], classFactor=f)

# this is what we want to parallelize
#R.null <- replicate(n.permut,test.mut.pathways(MUTtbl=tmp, gsets=gsetIdxs[new.restricted.pwy], classFactor=factor(factor)[sample(ncol(tmp))]))
R.null <- mclapply(X= c(1:n.permut), FUN= function(X){ test.mut.pathways(MUTtbl=tmp, gsets=gsetIdxs[new.restricted.pwy], classFactor=factor(factor)[sample(ncol(tmp))])} , mc.cores = 5,mc.set.seed=TRUE)
tmp.obj <- c()
tmp.name <- names(R.null[[1]])
R.null <- sapply(c(1:length(R.null)),function(x){ tmpo <- cbind(tmp.obj,R.null[[x]])})
rownames(R.null) <- tmp.name
rm(tmp.obj,tmp.name)
new.permuted.pwy <- sapply(names(R[which(!is.na(R))]),function(x){length(which(R.null[x,]<R[x]))/n.permut})
new.permuted.pwy
new.restricted.pwy <- names(which(new.permuted.pwy<(1/n.permut)))

# display the results
cat("the significant pathways for ", n.permut,"permutations are:","\n")
new.permuted.pwy[new.restricted.pwy]


####################################################################################################################
####################################################################################################################

#grat <- names(gsets)[grep(pattern="SPLIC",x=names(gsets))]
grat <- rat
# explore all the significant pathways from grat
for(i in 1:length(grat))
  
{
  par(oma=c(8,0,2,1))
  pwy <- grat[i]
  pwy1 <- pwy
  pwy2 <- pwy
  
  
  g12c.erk.muts <- rbind(rbind(tmp[gsetIdxs[[pwy]],f=="G12C"],tmp[gsetIdxs[[pwy1]],f=="G12C"]),tmp[gsetIdxs[[pwy2]],f=="G12C"])
  g12v.erk.muts <- rbind(rbind(tmp[gsetIdxs[[pwy]],f=="G12V"],tmp[gsetIdxs[[pwy1]],f=="G12V"]),tmp[gsetIdxs[[pwy2]],f=="G12V"])
  g12c.erk.muts <- g12c.erk.muts[-which(duplicated(rownames(g12c.erk.muts))),]
  g12v.erk.muts <- g12v.erk.muts[-which(duplicated(rownames(g12v.erk.muts))),]

  
  require(gplots)
  M <- cbind(g12c.erk.muts,g12v.erk.muts)
  vec <- c(rep(x=1,times=length(colnames(g12c.erk.muts))),rep(x=2,times=length(colnames(g12v.erk.muts))))
  heatmap.2(M,Rowv=FALSE,Colv="Rowv",dendrogram="none",scale="none",trace="none",
            colsep = 1:ncol(M), rowsep = 1:nrow(M), sepcolor = 'white', sepwidth=c(0.05, 0.05),
            keysize = 1.25, density.info = 'none',key=FALSE,
            ColSideColors=c("royalblue","orange")[vec],col = colorpanel(2, 'gray70', 'brown3'),main=paste(pwy,"\n~ G12C vs G12V status"))
  
  R[pwy]
}




