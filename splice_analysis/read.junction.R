# charles ferté
# feb 8th 2013
# read junction data


library(synapseClient)
synapseLogin(username="charles.ferte@sagebase.org",password="charles")

# load the data
setwd("/home/cferte/FELLOW/cferte/KRAS_Analysis/splice_luad_data/unc.edu_LUAD.IlluminaHiSeq_RNASeqV2.Level_3.1.7.0/")
myfiles <- list.files(pattern="junction")


# create a matrix of myfiles columns by the junction features as rows
tmp <- read.delim(file=myfiles[1],header=TRUE)
tmp <- unique(tmp$junction)
junction.matrix <- matrix(NA,nrow=length(tmp),ncol=length(myfiles))
rownames(junction.matrix) <- tmp
colnames(junction.matrix) <- myfiles

# read the junction file and make a list

tmp <- c()
tmp1 <- c()
for(k in c(1:length(myfiles))){
  tmp <- read.delim(file=myfiles[k],header=TRUE)
  duplicates <- tmp$junction[which(duplicated(tmp$junction)==TRUE)]
  
  # get rid of the duplicated features by averadging them
  n <- length(duplicates)
  for(i in c(1:n))
  { tmp$raw_counts[ tmp$junction==duplicates[i]] <- mean(tmp$raw_counts[ tmp$junction==duplicates[i]])
  }
  tmp <- tmp[-which(duplicated(tmp$junction)==TRUE),]
  rownames(tmp) <- tmp$junction
  tmp$junction <- NULL
  tmp1 <- rownames(tmp)
  tmp <- as.numeric(tmp$raw_counts)
  names(tmp) <- tmp1
  junction.matrix[names(tmp),k] <- tmp 
  rm(tmp1,tmp,duplicates,n)
}

save(junction.matrix,file="/home/cferte/FELLOW/cferte/KRAS_Analysis/splice_luad_data/junction.matrix.RDa")
load("/home/cferte/FELLOW/cferte/KRAS_Analysis/splice_luad_data/junction.matrix.RDa")

# read the mage file
mage <- read.delim("/home/cferte/FELLOW/cferte/KRAS_Analysis/splice_luad_data/unc.edu_LUAD.IlluminaHiSeq_RNASeqV2.mage-tab.1.8.0/unc.edu_LUAD.IlluminaHiSeq_RNASeqV2.1.8.0.sdrf.txt",header=TRUE)

myfiles1 <- sub(pattern="unc.edu.",replacement="",x=myfiles)
myfiles1 <- sub(pattern=".junction_quantification.txt",replacement="",x=myfiles1)
myfiles1 <- unlist(sapply(strsplit(x=myfiles1,split=".",fixed=TRUE),function(x){x[[1]]}))
mage$assay_edit  <-substr(x=mage$Assay.Name,start=1,stop=36) 
mage <- mage[-which(duplicated(mage$Assay.Name)),]
mage <- mage[ mage$assay_edit==myfiles1,]

colnames(junction.matrix) <- mage$Comment..TCGA.Barcode.

# save this in synapse

junction_matrix <- Data(list(name = "junction_matrix_luad", parentId = 'syn1670945'))
junction_matrix <- createEntity(junction_matrix)

# add object into the data entity
junction_matrix <- addObject(junction_matrix,junction.matrix)

# push the raw data into this entity
junction_matrix <- storeEntity(entity=junction_matrix)


