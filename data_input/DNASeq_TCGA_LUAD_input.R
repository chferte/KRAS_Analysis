## charles ferté,MD
## Sage Bionetworks 
## "2012-dec-06"

## KRAS project LUAD
## input the TCGA DNASeq data and build a matrix of specific mutations (rownames) for each sample (colnames)

require(synapseClient)
require(Biobase)
## synapse Login
synapseLogin("charles.ferte@sagebase.org","charles")

### magic option
options(stringsAsFactors=FALSE)

#######################################################################################################
# 1. read the maf file and name it MAF 
#######################################################################################################
#MAF <- read.delim(file="/home/jguinney/projects/AZRasPaper/data~/luad/LUAD.exome.cleaned.somatic.maf",header=T,skip=1)
MAF <- read.delim(file="/home/jguinney/projects/tcgaras/data~/PR_TCGA_LUAD_PAIR_Capture_All_Pairs_QCPASS.aggregated.capture.tcga.uuid.somatic.maf",header=T,skip=2)
MAF$Hugo_Symbol <- as.character(MAF$Hugo_Symbol)
MAF$Tumor_Sample_Barcode <- as.character(MAF$Tumor_Sample_Barcode)
MAF$Variant_Classification <- as.character(MAF$Variant_Classification)
MAF$Reference_Allele <- as.character(MAF$Reference_Allele)
MAF$Tumor_Seq_Allele1 <- as.character(MAF$Tumor_Seq_Allele1)
MAF$Tumor_Seq_Allele2 <- as.character(MAF$Tumor_Seq_Allele2)
MAF$Protein_Change <- as.character(MAF$Protein_Change)

#######################################################################################################
# 2. work with MAF1 from now to build the MATMUT matrix
#######################################################################################################

MAF1 <- MAF
MAF1 <- MAF1[which(MAF1$Hugo_Symbol!="Unknown"),]

MAF1$CHR_POS <- paste("chr",MAF1$Chromosome,":",MAF1$Start_Position,"-",MAF1$End_Position,sep="")
MAF1$CODE <- as.character(MAF1$dbSNP_RS)
MAF1$CODE[which(nchar(MAF1$CODE)==0)] <- "NA"

MAF1$MUTID <- paste(MAF1$Hugo_Symbol,MAF1$Protein_Change,sep="_")

MAF1 <- MAF1[,c("MUTID","Tumor_Sample_Barcode","Hugo_Symbol","Variant_Classification","CHR_POS","Reference_Allele","Tumor_Seq_Allele1","cDNA_Change", "Protein_Change")]

# create a bianry matrix full of zeros with colnames equal to the sampel names and rownames equal the MUTIDs

MATMUT<-matrix(0,nrow=length(unique(MAF1$MUTID)),ncol=length(unique(MAF1$Tumor_Sample_Barcode)))
colnames(MATMUT) <- unique(MAF1$Tumor_Sample_Barcode)
rownames(MATMUT) <- unique(MAF1$MUTID)

# assign 1 to any sample mutated for any MUTID
for(i in rownames(MATMUT)){
  MATMUT[i,c(MAF1$Tumor_Sample_Barcode[which(MAF1$MUTID==i)])] <- 1  
}

foo <- MAF1[-which(duplicated(MAF1$MUTID)),c("MUTID","Variant_Classification")]
rownames(foo) <- foo$MUTID
foo$frequency <- apply(MATMUT,1,sum)
MUT_TYPE <- foo[rownames(MATMUT),]

MATMUT_LUAD <- as.data.frame(MATMUT)
MUT_TYPE_LUAD <- as.data.frame(MUT_TYPE)

rm(foo,MAF,MAF1,i,MATMUT)

#######################################################################################################
# 2. save the objects MATMUT and MUT_TYPE in Belltown
#######################################################################################################

save(MUT_TYPE_LUAD,file="/home/cferte/FELLOW/cferte/KRAS_Analysis/MUT_TYPE_LUAD.RData")
save(MATMUT_LUAD,file="/home/cferte/FELLOW/cferte/KRAS_Analysis/mutations_LUAD.RData")


#######################################################################################################
# 3. create KRAS_LUAD vector
#######################################################################################################
tmp <- MATMUT_LUAD[grep(pattern="KRAS", rownames(MATMUT_LUAD)),]
rownames(tmp) <- substr(x=rownames(tmp),8,nchar(rownames(tmp)))
tmp1 <- c()
for(i in colnames(tmp)){
tmp1 <- c(tmp1,ifelse(sum(tmp[,i])==0,"WT",rownames(tmp)[which(tmp[,i]==1)]))
}
names(tmp1) <- colnames(tmp)
KRAS_LUAD <- tmp1
save(KRAS_LUAD,file="/home/cferte/FELLOW/cferte/KRAS_Analysis/KRAS_LUAD.RData")

