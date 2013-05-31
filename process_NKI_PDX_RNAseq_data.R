library(Rsamtools)
library(cqn)
library(org.Hs.eg.db)

library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

bamDir <- "/external-data/DAT_108__MEK_NKI/2013_04/"

txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
gnModel <- exonsBy(txdb, "gene")
bamFls <- list.files(bamDir, "bam$", full=TRUE)
names(bamFls) <- sub("\\..*", "", basename(bamFls))

# ras mutations
which <- RangesList("12"=IRanges(25398280L, 25398285L))
which <- RangesList("12"=IRanges(25358180L,25403854L))
what <- c("rname", "strand", "pos", "qwidth", "seq")
param <- ScanBamParam(which=which, what=what)
bam <- scanBam(bamFls[2], param=param)

counter <- function(fl, gnModel){
  cat("processing ",fl,"\n")
  aln <- readGappedAlignments(fl)
  seqnames(aln) <- factor(paste("chr", as.character(seqnames(aln)),sep=""))
  strand(aln) <- "*" # for strand-blind sample prep protocol
  hits <- countOverlaps(aln, gnModel)
  counts <- countOverlaps(gnModel, aln[hits==1])
  names(counts) <- names(gnModel)
  counts
}
counts <- sapply(bamFls, counter, gnModel)
save(counts, file="/home/cferte/NKIcounts.rda")

load("/home/cferte/NKIcounts.rda")
genes <- mget(rownames(counts), org.Hs.egSYMBOL, ifnotfound=NA)
counts <- counts[!is.na(genes),]
genes <- genes[!is.na(genes)]
rownames(counts) <- genes



# load the annotation file
annot <- read.table("/home/cferte/NKI_PDX_annot..txt",header=TRUE)
annot$treatment[annot$treatment!="Untreated"] <- "AZD6244"
annot$True_ID <- paste(annot$ID,annot$tumor_sample_ID,sep="_")
colnames(counts) <- substr(colnames(counts),1,6)
tmp <- intersect(colnames(counts),annot$True_ID)
counts <- counts[,tmp]
rownames(annot) <- annot$True_ID
annot <- annot[tmp,]


# Principal componnet analysis
s <- svd(counts)
plot(s$d^2/sum(s$d^2),pch=19)
color <- as.numeric(factor(annot$treatment))
plot(s$v[,1],s$v[,2],pch=19,col=color)
boxplot(log.counts)
boxplot(nki_pdx_cqnNormalized)

# compute the RPKM from the counts matrix
load("/home/cferte/SPV_ERCC1_RNAseq/gene.annot.rda")
gene.annot <- gene.annot[!duplicated(gene.annot$hgnc_symbol),]
rownames(gene.annot) <- gene.annot$hgnc_symbol
tmp <- intersect(gene.annot$hgnc_symbol,rownames(counts))
counts <- counts[tmp,]
gene.annot <- gene.annot[tmp,]


final.result <- c()
for(sample.idx in colnames(counts)){
  result <- c()
  for(gene.idx in rownames(counts))
  {
    RPK <- counts[gene.idx,sample.idx]/(gene.annot$length[gene.annot$hgnc_symbol==gene.idx]/1000)
    RPKM <- RPK/(colSums(counts)[sample.idx]/1e6)
    #names(RPKM) <- gene.idx
    result <- c(result,RPKM)
  }
  final.result <- cbind(final.result,result)
}
rownames(final.result) <- rownames(counts) 
colnames(final.result) <- colnames(counts)
boxplot(log(final.result+1))
RPKM <- final.result

### counts to CQN normalized
gene.annot <- read.table("data~/gene.annot.txt",sep="\t",header=TRUE)
idxs <- match(genes.nona, gene.annot$hgnc_symbol)
counts.m <- counts.nona[!is.na(idxs),]
gene.annot.m <- gene.annot[na.omit(idxs),]

cqn <- cqn(counts.m, lengths=gene.annot.m$length, x=gene.annot.m$gcpercent,sizeFactors=colSums(counts.m),verbose=TRUE)
nki_pdx_cqnNormalized <- cqn$y + cqn$offset
rownames(nki_pdx_cqnNormalized) <- gene.annot.m$hgnc_symbol
save(nki_pdx_cqnNormalized, file="~/data/NKI_PDX_CQNnormalized.rda")

