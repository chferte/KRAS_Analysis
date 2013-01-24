#charles ferté
# jan 7 2013
# sage bionetworks

# create data files for GSEA java

###############################################################
# load the LUAD data  
###############################################################

source("/home/cferte/FELLOW/cferte/KRAS_Analysis/data_input/TCGA_LUAD_input.R")
tmp <- names(KRAS_LUAD[KRAS_LUAD %in% c("G12C","G12V")])

KRAS_LUAD <- KRAS_LUAD[tmp]
KRAS_LUAD[KRAS_LUAD=="G12C"] <- 1
KRAS_LUAD[KRAS_LUAD=="G12V"] <- 0

KRAS_LUAD <- as.numeric(KRAS_LUAD)
names(KRAS_LUAD) <- tmp
table(KRAS_LUAD)
LUAD_EXP <- LUAD_EXP[,names(KRAS_LUAD)]

# make the GSEA files
library(ArrayTools)
blah <- as.data.frame(KRAS_LUAD)
A <- ExpressionSet(LUAD_EXP,phenoData=AnnotatedDataFrame(blah))

# make GSEA files
createGSEAFiles(mydir = "/gluster/home/cferte/FELLOW/cferte/KRAS_Analysis/biological_info_meaning/GSEA_gene_expression/",eSet=A, catVar="KRAS_LUAD")



