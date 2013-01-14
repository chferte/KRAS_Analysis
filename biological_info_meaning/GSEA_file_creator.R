#charles ferté
# jan 7 2013
# sage bionetworks

# create data files for GSEA java

###############################################################
# load the LUAD data  
###############################################################

source("/home/cferte/FELLOW/cferte/KRAS_Analysis/data_input/TCGA_LUAD_input.R")
tmp <- names(KRAS_LUAD)
KRAS_LUAD[KRAS_LUAD %in% c("A146P","D33E","G12A","G12D","G12F","G12R","G12S","G12Y","G13C","G13D","L19F","Q61H","Q61L")] <- "rare"
table(KRAS_LUAD)
KRAS_LUAD[KRAS_LUAD=="G12C"] <- 1
KRAS_LUAD[KRAS_LUAD=="G12V"] <- 2
KRAS_LUAD[KRAS_LUAD=="G12A"] <- 3
KRAS_LUAD[KRAS_LUAD=="G12D"] <- 3
KRAS_LUAD[KRAS_LUAD=="rare"] <- 3
KRAS_LUAD[KRAS_LUAD=="WT"] <- 0
KRAS_LUAD <- as.numeric(KRAS_LUAD)
names(KRAS_LUAD) <- tmp
table(KRAS_LUAD)

# make the GSEA files
library(ArrayTools)
blah <- as.data.frame(KRAS_LUAD)
A <- ExpressionSet(LUAD_EXP,phenoData=AnnotatedDataFrame(blah))

# make GSEA files
createGSEAFiles(mydir = "/gluster/home/cferte/FELLOW/cferte/KRAS_Analysis/biological_info_meaning/",eSet=A, catVar="KRAS_LUAD")



