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

# get rid of the unknwon gene s and the Sielnt mutations
MAF1 <- MAF1[which(MAF1$Hugo_Symbol!="Unknown"),]
MAF1 <- MAF1[which(MAF1$Variant_Classification!="Silent"),]
gene <- unique(MAF1$Hugo_Symbol)
samples <- unique(MAF1$Tumor_Sample_Barcode)

j <- c()

MAF1 <- MAF1[,c( "Tumor_Sample_Barcode","Hugo_Symbol","Protein_Change","Genome_Change")]
MAF1$Change <- MAF1$Protein_Change
MAF1$Change[which(MAF1$Protein_Change=="")] <- MAF1$Genome_Change[which(MAF1$Protein_Change=="")] 
MAF1 <- MAF1[,-c(3,4)]
which(colnames(MAF1) %in% c("Protein_Change","Genome_Change"))
MAF1$MUTID <- paste(MAF1$Hugo_Symbol,MAF1$Change,sep="_")

# create a bianry matrix full of zeros with colnames equal to the sampel names and rownames equal the MUTIDs
MATMUT<-matrix(0,nrow=length(unique(MAF1$MUTID)),ncol=length(unique(MAF1$Tumor_Sample_Barcode)))
colnames(MATMUT) <- unique(MAF1$Tumor_Sample_Barcode)
rownames(MATMUT) <- unique(MAF1$MUTID)

# assign 1 to any sample mutated for any MUTID
for(i in rownames(MATMUT)){MATMUT[i,c(MAF1$Tumor_Sample_Barcode[which(MAF1$MUTID==i)])] <- 1}

idx <- sapply(strsplit(x=rownames(MATMUT),split="_"), function(x){x[[1]]})
gene <- unique(idx)
j <- c()
new.matmut <- matrix(0,nrow=length(gene),ncol=length(colnames(MATMUT)))
rownames(new.matmut) <- gene
colnames(new.matmut) <- colnames(MATMUT)
 
for(i in rownames(new.matmut))
  {
if(length(which(idx==i))>1) 
  {
  new.matmut[i,c(names(which(apply(MATMUT[which(idx==i),],2,sum)!=0)))] <- 1 } 
else 
  {new.matmut[i,names(which(MATMUT[which(idx==i),]==1))] <- 1
}
}
 
MATMUT_LUAD <- new.matmut

# sanity check
for(i in c(1:length(colnames(new.matmut))))
{
a <- rownames(new.matmut)[which(new.matmut[,i]==1)]
b <- unique(MAF1$Hugo_Symbol[ MAF1$Tumor_Sample_Barcode ==colnames(new.matmut)[i]])
print(table(is.na(match(b,a))))
}
    

save(file="/home/cferte/FELLOW/cferte/KRAS_Analysis/MATMUT_LUAD.RData",MATMUT_LUAD)
rm(MAF,MAF1,i,MATMUT)

#######################################################################################################
# 2. save the objects MATMUT and MUT_TYPE in Belltown
#######################################################################################################

save(MUT_TYPE_LUAD,file="/home/cferte/FELLOW/cferte/KRAS_Analysis/MUT_TYPE_LUAD.RData")
save(MATMUT_LUAD,file="/home/cferte/FELLOW/cferte/KRAS_Analysis/mutations_LUAD.RData")
save(file="/home/cferte/FELLOW/cferte/KRAS_Analysis/MATMUT_LUAD.RData",MATMUT_LUAD)

# #######################################################################################################
# # 3. create KRAS_LUAD vector
# #######################################################################################################
# tmp <- MATMUT_LUAD[grep(pattern="KRAS", rownames(MATMUT_LUAD)),]
# rownames(tmp) <- substr(x=rownames(tmp),8,nchar(rownames(tmp)))
# tmp1 <- c()
# for(i in colnames(tmp)){
# tmp1 <- c(tmp1,ifelse(sum(tmp[,i])==0,"WT",rownames(tmp)[which(tmp[,i]==1)]))
# }
# names(tmp1) <- colnames(tmp)
# KRAS_LUAD <- tmp1
# save(KRAS_LUAD,file="/home/cferte/FELLOW/cferte/KRAS_Analysis/KRAS_LUAD.RData")
# 
