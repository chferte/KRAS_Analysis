## Charles Ferté
### july 2012
#### adaptation of the gsea r code

options(stringsAsFactors=FALSE)

### load the data
library(affy)
library(corpcor)
library(lattice)
library(limma)
require("foreach")
library(synapseClient)
synapseLogin("charles.ferte@sagebase.org","charles")

## load the different lung gene expression files
# CAUTION: make sure to assign the rank.features and coeff.features values
# rank.features should be the ranking (decreasing =TRUE) of the list of features
# default: coeff.features = 1 for all the rank.features 
# coeff.features <- rep(1,times=length(rank.features))


### load GSEA by ISJ
# load the miSig DB (these code for MsigDB for geneSet from pathways)
mSigDB_annotations <- loadEntity("syn105363")
mSigDB_symbolID <- loadEntity("syn105350")
DB<-mSigDB_symbolID$objects$MsigDB_symbolID

# select C2 and C5 (you can select other object for your pathway)
#allPathways <- c(mSigDB_annotations$objects$C2$KEGG,mSigDB_annotations$objects$C2$CGP,mSigDB_annotations$objects$C2$BIOCARTA,mSigDB_annotations$objects$C2$REACTOME)
allPathways <- c(mSigDB_annotations$objects$C2$KEGG,mSigDB_annotations$objects$C2$BIOCARTA,mSigDB_annotations$objects$C2$REACTOME)
# geneset elements
geneAllSetList <-DB$genesets[is.element(DB$geneset.names,allPathways)]

source("/home/cferte/FELLOW/cferte/TOOLS/myPathwayAnalysis1.R")
source("/home/cferte/FELLOW/cferte/TOOLS/preRankedTest.R")

referenceSet <- rank.features*coeff.features


analyticResult <-foreach (curPathway = allPathways) %do%{
  mSigDB_index <- which(DB$geneset.names == curPathway)
  curPathwayGenes <- DB$genesets[[mSigDB_index]]
  curPathwayGenes <- intersect(names(rank.features),curPathwayGenes)
  
  # R5 class call
  ## enter the dataset = referenceSet as a txt file with fiurst comunm= gene 2 nd column= statistic value, 
  
  pathwayAnalysis<-myPathwayAnalysis1$new()
  pathwayAnalysis$gsea(referenceSet,curPathwayGenes,np=1000,w =1)
  #pathwayAnalysis$fet(AllGenes,curPathwayGenes,testSet)
  return(pathwayAnalysis)
  
}



## get the pthw that are associated with p< .05 (SIGPTW)
k <- c()
for(i in (1:length(allPathways))){
  k <- c(k, analyticResult[[i]]$gseaResult$p.value)
}
names(k) <- allPathways
SIGPTW <- k[which(k<.05)]



save(k,analyticResult,allPathways,file=final_file)

## you can plot the ES:
# analyticResult[[i]]$gseaPlot()

