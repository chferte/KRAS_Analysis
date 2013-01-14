# Charles Ferté
# Sage Bionetworks
# 14 Sept 2012

###############################################################
# load the LUAD data  for the KRAS analysis
###############################################################

load("/home/jguinney/projects/AZRasPaper/data~/luad/luad_rnaseq_v3.1.6.rda")
load("/home/cferte/FELLOW/cferte/KRAS_Analysis/mutations_LUAD.RData")
load("/home/cferte/FELLOW/cferte/KRAS_Analysis/KRAS_LUAD.RData")
LUAD_EXP <- RPKM.cqn

###############################################################
# make the sample names coherent between gene expr and mutation
###############################################################

colnames(LUAD_EXP) <- substr(x=colnames(RPKM.cqn),1,16)
names(KRAS_LUAD) <- substr(x=names(KRAS_LUAD),1,16)
tmp <- intersect(colnames(LUAD_EXP),names(KRAS_LUAD))
LUAD_EXP <- LUAD_EXP[,tmp]
KRAS_LUAD <- KRAS_LUAD[tmp]
rm(tmp,RPKM.cqn,Counts,MATMUT_LUAD)

