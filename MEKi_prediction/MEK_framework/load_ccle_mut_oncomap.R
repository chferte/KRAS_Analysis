oncomap <- read.delim("/home/cferte/cell_line_data/CCLE_Oncomap3_2012-04-09.maf",header=TRUE)
oncomap <- oncomap[which(oncomap$Hugo_Symbol!="Unknown"),]
oncomap <- oncomap[,c("Hugo_Symbol","Tumor_Sample_Barcode","Protein_Change", "Genome_Change")]
oncomap$Protein_Change[oncomap$Protein_Change==""] <- oncomap$Genome_Change[ oncomap$Protein_Change==""]
MATMUT<-matrix(0,nrow=length(unique(oncomap$Hugo_Symbol)),ncol=length(unique(oncomap$Tumor_Sample_Barcode)))
colnames(MATMUT) <- unique(oncomap$Tumor_Sample_Barcode)
rownames(MATMUT) <- unique(oncomap$Hugo_Symbol)
for(i in rownames(MATMUT)){
  MATMUT[i,c(oncomap$Tumor_Sample_Barcode[which(oncomap$Hugo_Symbol==i)])] <- c(oncomap$Protein_Change[which(oncomap$Hugo_Symbol==i)])  
}
unique(oncomap$Hugo_Symbol)
ccle_oncomap <- MATMUT
rm(MATMUT,oncomap,i)
