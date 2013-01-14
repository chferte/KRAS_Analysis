# Charles Ferté
# Sage Bionetworks
# Dec 8th 2012


###############################################################
# load the CHEMORES data  
###############################################################

load("/home/cferte/FELLOW/cferte/KRAS_Project/OBJECTS/CHEMORES_EXP.RData")
load("/home/cferte/FELLOW/cferte/KRAS_Project/OBJECTS/CHEMORES_CLIN.RData")
load("/home/cferte/FELLOW/cferte/KRAS_Project/OBJECTS/KRAS_CHEMORES_STATUS.RData")
tmp <- ifelse(KRAS_CHEMORES %in% c("G12A","Q61H"),"rare",KRAS_CHEMORES)
names(tmp) <- names(KRAS_CHEMORES)
KRAS_CHEMORES <- tmp
rm(tmp)

###############################################################
# get rid of SCC and SCLC samples
###############################################################

CHEMORES_CLIN <- CHEMORES_CLIN[ CHEMORES_CLIN$Histology!="SCC",]
CHEMORES_CLIN <- CHEMORES_CLIN[ CHEMORES_CLIN$Histology!="SCLC",]
CHEMORES_EXP <- CHEMORES_EXP[,rownames(CHEMORES_CLIN)]
KRAS_CHEMORES <- KRAS_CHEMORES[rownames(CHEMORES_CLIN)]

