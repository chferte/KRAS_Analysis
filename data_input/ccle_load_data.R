# charles ferté

###############################################################
# load the CCLE data (snm normalized) 
###############################################################
ccle_eset <- loadEntity('syn1467966')
ccle_eset <- ccle_eset$objects$CCLE_RMA

# ccle_eset <- loadEntity('syn1354643') # previous eset
ccle_drug <- loadEntity('syn1354656')
ccle_kras <- loadEntity('syn1443160')
KRAS_CCLE <- ccle_kras$objects$KRAS_CCLE
KRAS_CCLE <- substr((KRAS_CCLE),3,nchar(KRAS_CCLE))
KRAS_CCLE[KRAS_CCLE==""] <- "WT"

# make  ccle_drug and ccle_EXP coherent objects
ccle_kras <- ccle_kras$objects$KRAS_CCLE
ccle_drug <- ccle_drug$objects$ccle_drug
ccle_drug <- ccle_drug[names(ccle_kras),]
ccle_drug <- ccle_drug@data
ccle_EXP <- ccle_eset
CCLE_EXP <- ccle_EXP[,names(ccle_kras)]
tmp <- ifelse(KRAS_CCLE %in% c("G12A","G12S","G13D","G13C","Q61H","Q61K","Q61L"),"rare",KRAS_CCLE)
names(tmp) <- names(KRAS_CCLE)
KRAS_CCLE <- tmp
rm(tmp,ccle_eset,ccle_EXP,ccle_kras)

