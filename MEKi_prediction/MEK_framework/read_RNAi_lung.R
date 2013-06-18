# charles ferte
# read the RNAi screens from NKI

# select the top 30 features using the bioinformatic pipeline
vdx  <- gsub(pattern=c("mut_"),replacement="",fixed=TRUE,rownames(drug_fmat[[1]])[1:30])
vdx  <- gsub(pattern=c("amp_"),replacement="",fixed=TRUE,vdx)
vdx  <- gsub(pattern=c("del_"),replacement="",fixed=TRUE,vdx)
vdx


# select the top features from nki rnai screens
cancer.census <- read.delim(file="/home/cferte/cancer_gene_census.txt")
cancer.census <- cancer.census$Symbol
foo <- read.csv(file="/external-data/DAT_108__MEK_NKI/2013_05/MEK_screens_SAGE/results_cell lines/sw480/ut_azd/A=ut B=azd.csv",header=TRUE)
foo <- foo[which(foo$GeneSymbol %in% cancer.census),]
foo$padj2 <- p.adjust(foo$pval,method="BH")
bar <- read.csv(file="/external-data/DAT_108__MEK_NKI/2013_05/MEK_screens_SAGE/results_cell lines/sw620/ut_azd/A=ut B=azd.csv",header=TRUE)
bar <- bar[which(bar$GeneSymbol %in% cancer.census),]
bar$padj2 <- p.adjust(bar$pval,method="BH")

o1 <- as.character(unique(foo$GeneSymbol[which( foo$padj2<.05 & abs(foo$log2FoldChange)>1.5)]))
o2 <- as.character(unique(bar$GeneSymbol[which( bar$padj2<.05 & abs(bar$log2FoldChange)>1.5)]))
nki.features <- unique(c(o1,o2))
intersect(nki.features,vdx)


# translate it into the 
