# charles ferte
# read the RNAi screens from NKI

# select the top 30 features using the bioinformatic pipeline
vdx.lung  <- gsub(pattern=c("mut_"),replacement="",fixed=TRUE,rownames(drug_fmat[[1]])[1:30])
vdx.lung  <- gsub(pattern=c("amp_"),replacement="",fixed=TRUE,vdx.lung)
vdx.lung  <- gsub(pattern=c("del_"),replacement="",fixed=TRUE,vdx.lung)
vdx.lung


# select the top features from nki rnai screens
cancer.census <- read.delim(file="/home/cferte/cancer_gene_census.txt")
cancer.census <- cancer.census$Symbol
foo <- read.csv(file="/external-data/DAT_108__MEK_NKI/2013_05/MEK_screens_SAGE/results_cell lines/h358/ut_azd/A=ut B=azd.csv",header=TRUE)
foo <- foo[which(foo$GeneSymbol %in% cancer.census),]
foo$padj2 <- p.adjust(foo$pval,method="BH")
bar <- read.csv(file="/external-data/DAT_108__MEK_NKI/2013_05/MEK_screens_SAGE/results_cell lines/h1792/ut_azd/A=ut B=azd.csv",header=TRUE)
bar <- bar[which(bar$GeneSymbol %in% cancer.census),]
bar$padj2 <- p.adjust(bar$pval,method="BH")

o1 <- unique(foo$GeneSymbol[which( foo$padj2<.05 & abs(foo$log2FoldChange)>1.5)])
o2 <- unique(bar$GeneSymbol[which( bar$padj2<.05 & abs(bar$log2FoldChange)>1.5)])
nki.features <- unique(c(o1,o2))

intersect(nki.features,vdx.lung)

foo[ foo$GeneSymbol =="LCK" ,]