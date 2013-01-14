# charles ferté
# sage bionetworks
# dec 8th 2012

##########################################
# load the Sanger data 
##########################################
source("~/FELLOW/cferte/KRAS_Analysis/data_input/CellLine_input.R")

SangerExpr <- getSangerExpr()

# retrieve the names of the enes of the Sanger dataset
library(org.Hs.eg.db)
rownames(SangerExpr) <- sapply(strsplit(rownames(SangerExpr),split="_"),function(x){x[[1]]})
tmp <- unlist(mget(x=rownames(SangerExpr),org.Hs.egSYMBOL,ifnotfound=NA))
rownames(SangerExpr) <- tmp

# restrict the analysis on the Lung datsets
SangerExpr <- SangerExpr[,grep(pattern="LUNG",colnames(SangerExpr))]
colnames(SangerExpr) <- sapply(strsplit(colnames(SangerExpr),split="_"),function(x){x[[1]]})

# get Sanger Drug sensitivity information and make it coherent with he gene expression dataset
SangerDrug <- getSangerDrug()
SangerDrug <- SangerDrug[intersect(rownames(SangerDrug),colnames(SangerExpr)),]
SangerExpr <- SangerExpr[,intersect(rownames(SangerDrug),colnames(SangerExpr))]
SangerInfo$CellName <- gsub(pattern="-",replacement="",x=SangerInfo$SampleName)
tmp <- intersect(rownames(SangerDrug),SangerInfo$CellName)
SangerLung <- SangerInfo[SangerInfo$CellName %in% tmp,]
SangerLung <- SangerLung[which(duplicated(SangerLung$CellName)!=TRUE),]
SangerDrug <- SangerDrug[tmp,]
SangerExpr <- SangerExpr[,tmp]

# retrieve the KRAS mutation status information
SangerMut <- as.data.frame(getSangerMut())
SangerMut <- SangerMut[ SangerMut$Cell.Line %in% SangerLung$SampleName,]
SangerMut$CellName <- gsub(pattern="-",replacement="",x=SangerMut$Cell.Line)
rownames(SangerMut) <- SangerMut$CellName
SangerMut <- SangerMut[tmp,]
KRAS_Sanger <- gsub(pattern="p.",replacement="",x=SangerMut$KRAS)
names(KRAS_Sanger) <- rownames(SangerMut)
KRAS_Sanger[KRAS_Sanger=="wt"] <- "WT"
KRAS_Sanger[KRAS_Sanger %in% c("G12A","G12F","G12S","G13C","G13D","Q61H")] <- "rare"
table(KRAS_Sanger)

# clean up !
SANGER_EXP <- SangerExpr
KRAS_SANGER <- KRAS_Sanger
rm(SangerInfo,SangerExpr,SangerLung,SangerMut,tmp,KRAS_Sanger)