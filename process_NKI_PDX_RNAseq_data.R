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
save(counts, file="NKIcounts.rda")

genes <- mget(rownames(counts), org.Hs.egSYMBOL, ifnotfound=NA)

counts.nona <- counts[!is.na(genes),]
genes.nona <- genes[!is.na(genes)]


### counts to CQN normalized
gene.annot <- read.table("data~/gene.annot.txt",sep="\t",header=TRUE)
idxs <- match(genes.nona, gene.annot$hgnc_symbol)
counts.m <- counts.nona[!is.na(idxs),]
gene.annot.m <- gene.annot[na.omit(idxs),]

cqn <- cqn(counts.m, lengths=gene.annot.m$length, x=gene.annot.m$gcpercent,sizeFactors=colSums(counts.m),verbose=TRUE)
nki_pdx_cqnNormalized <- cqn$y + cqn$offset
rownames(nki_pdx_cqnNormalized) <- gene.annot.m$hgnc_symbol
save(nki_pdx_cqnNormalized, file="~/data/NKI_PDX_CQNnormalized.rda")

