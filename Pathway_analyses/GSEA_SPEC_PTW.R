## Charles Ferté
### july 2012

# adaptation from the code of Inn Sock Jang
options(stringsAsFactors=FALSE)

### load the data
library(affy)
library(corpcor)
library(lattice)
library(limma)
require("foreach")
library(synapseClient)
synapseLogin("charles.ferte@sagebase.org","charles")

# load the G12C G12V LUAD pval data from synapse
rawdata <- loadEntity('syn1419785')
rawdata <- rawdata$objects$LESS

#define PVAL and CVAL
PVAL <- rawdata
CVAL <- rep(x=1,times=length(PVAL))

### load GSEA by ISJ
# load the miSig DB (these code for MsigDB for geneSet from pathways)
mSigDB_annotations <- loadEntity("syn105363")
mSigDB_symbolID <- loadEntity("syn105350")
DB<-mSigDB_symbolID$objects$MsigDB_symbolID

# load the x pathways that were significant for G12C LUAD during the first pass
# PTW <- read.delim("/gluster/home/cferte/FELLOW/cferte/KRAS_Project/OBJECTS/MUT_OVERLAP/LUAD_MUT_OVERLAP/G12V_MUT_LUAD_PROFILE.txt")
# PTW$Enriched.Pathway
# FOCP <- paste(PTW$Enriched.Pathway,collapse=" ")
# #FOpaste(sapply(strsplit(x=paste(PTW$p.value)[1:28],split=" "),function(x){x[[1]]}),collapse=" ")

# enter the names of the pathways as a unique character string sepoarted by spaces
#FOCP  <- "REACTOME_ACTIVATION_OF_ATR_IN_RESPONSE_TO_REPLICATION_STRESS REACTOME_CELL_CYCLE_CHECKPOINTS REACTOME_DOWNSTREAM_SIGNALING_OF_ACTIVATED_FGFR REACTOME_FORMATION_AND_MATURATION_OF_MRNA_TRANSCRIPT REACTOME_FORMATION_OF_A_POOL_OF_FREE_40S_SUBUNITS REACTOME_GENE_EXPRESSION KEGG_RIBOSOME  REACTOME_OLFACTORY_SIGNALING_PATHWAY REACTOME_PEPTIDE_CHAIN_ELONGATION REACTOME_TRANSCRIPTION REACTOME_VIRAL_MRNA_TRANSLATION REACTOME_INFLUENZA_VIRAL_RNA_TRANSCRIPTION_AND_REPLICATION KEGG_OXIDATIVE_PHOSPHORYLATION KEGG_CHRONIC_MYELOID_LEUKEMIA REACTOME_ELONGATION_AND_PROCESSING_OF_CAPPED_TRANSCRIPTS"
#FOCP <- unlist(strsplit(x=FOCP,split=" "))
#FOCP <- FOCP[-8]

# select C2 and C5 (you can select other object for your pathway)
#allPathways <- c(mSigDB_annotations$objects$C2$KEGG,mSigDB_annotations$objects$C2$CGP,mSigDB_annotations$objects$C2$BIOCARTA,mSigDB_annotations$objects$C2$REACTOME)
allPathways <- c(mSigDB_annotations$objects$C2$KEGG,mSigDB_annotations$objects$C2$BIOCARTA,mSigDB_annotations$objects$C2$REACTOME)
#allPathways <- FOCP

# geneset elements
geneAllSetList <-DB$genesets[is.element(DB$geneset.names,allPathways)]

source("~/KRAS_Analysis/myPathwayAnalysis1.R")
source("~/KRAS_Analysis/preRankedTest.R")

referenceSet <- -log10(PVAL)*CVAL


analyticResult <-foreach (curPathway = allPathways) %do%{
  mSigDB_index <- which(DB$geneset.names == curPathway)
  curPathwayGenes <- DB$genesets[[mSigDB_index]]
  curPathwayGenes <- intersect(names(PVAL),curPathwayGenes)
  
  # R5 class call
  ## enter the dataset = referenceSet as a txt file with fiurst comunm= gene 2 nd column= statistic value, 
  
  pathwayAnalysis<-myPathwayAnalysis1$new()
  pathwayAnalysis$gsea(referenceSet,curPathwayGenes,np=1000000,w =1)
  #pathwayAnalysis$fet(AllGenes,curPathwayGenes,testSet)
  return(pathwayAnalysis)
  
}



## get the pthw that are associated with p< .05 (SIGPTW)
k <- c()
for(i in (1:length(allPathways))){
  k <- c(k, analyticResult[[i]]$gseaResult$p.value)
}
names(k) <- allPathways
result <- list(k=k,analyticResult=analyticResult)

# store the k and analyticResult in Synapse (here the project id = syn1194693)
#load the data in the workspace ...here it is G12C_OVERLAP_PVAL your project
# then, store the data under the Synapse id project
library(synapseClient)
synapseLogin("charles.ferte@sagebase.org","charles")
newLayer <- Data(list(name = "G12V_G12C_LUAD_MUT_GSEA", parentId = "syn1194693"))
newLayer <- addObject(newLayer, result )
newLayer <- storeEntity(newLayer)
