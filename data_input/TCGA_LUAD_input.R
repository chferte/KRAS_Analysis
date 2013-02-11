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


# luad_exp <- LUAD_EXP
# luad_mut <- MATMUT_LUAD
# luad_kras <- KRAS_LUAD
# luad_data <- list(luad_exp=luad_exp,luad_mut=luad_mut, luad_kras=luad_kras)
# 
# # save the object luad_data in synapse
# luad_all <- Data(list(name = "luad_all", parentId = 'syn1670945'))
# luad_all <- createEntity(luad_all)
# 
# # add object into the data entity
# luad_all <- addObject(luad_all,luad_data)
# 
# # push the raw data into this entity
# luad_all <- storeEntity(entity=luad_all)


###############################################################
# make the sample names coherent between gene expr and mutation
###############################################################

colnames(LUAD_EXP) <- substr(x=colnames(RPKM.cqn),1,16)
names(KRAS_LUAD) <- substr(x=names(KRAS_LUAD),1,16)
tmp <- intersect(colnames(LUAD_EXP),names(KRAS_LUAD))
LUAD_EXP <- LUAD_EXP[,tmp]
KRAS_LUAD <- KRAS_LUAD[tmp]
rm(tmp,RPKM.cqn,Counts,MATMUT_LUAD)

