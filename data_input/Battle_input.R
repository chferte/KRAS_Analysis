# charles ferté
# sage bionetworks
# November 9th


##################################################################################
# Input BATTLE moelcular data
##################################################################################
library(affy)
liste <- list.celfiles(path="/home/cferte/FELLOW/cferte/KRAS_Project/BATTLE2/BATTLE2_CEL/",full.names=T)
battle_EXP <- ReadAffy(filenames=liste)
sampleNames(battle_EXP) <- substr(sampleNames(battle_EXP),1,9)

rawdata <- battle_EXP
battle_EXP <- rma(battle_EXP)
library(hugene10sttranscriptcluster.db)
tmp <- unlist(mget(x=featureNames(battle_EXP),hugene10sttranscriptclusterSYMBOL,ifnotfound=NA))

combine_probes_2_gene <- function(expr, genes, method="svd"){
  
  if(is.list(genes)) genes <- unlist(genes)
  
  stopifnot(dim(expr)[1] ==  length(genes))
  ugenes <- unique(genes)
  ugenes <- sort(ugenes[!is.na(ugenes)])
  M <- matrix(NaN, ncol=dim(expr)[2],nrow=length(ugenes),
              dimnames=list(ugenes, colnames(expr)))
  
  for(gene in ugenes){
    sub.expr <- as.matrix(expr[which(genes == gene),])
    if(dim(sub.expr)[2] == 1){
      M[gene,] <- sub.expr
    }else{
      tmp <- svd(sub.expr - rowMeans(sub.expr))$v[,1]
      tmp.c <- mean(cor(tmp, t(sub.expr)))
      #cat(gene," ", tmp.c, "\n")
      multiplier <- ifelse(tmp.c < 0, -1, 1)
      M[gene,] <- tmp * multiplier
    }
  }
  M
}

BATTLE_EXP <- combine_probes_2_gene(expr=battle_EXP,genes=tmp)


colnames(BATTLE_EXP) <- sampleNames(battle_EXP)
rm(tmp,liste,mapCdfName)



##################################################################################
# save the data in Synapse

##################################################################################
# Input BATTLE clinical data
##################################################################################

# read the matrix series file
battle_CLIN <- as.data.frame(read.delim("/home/cferte/FELLOW/cferte/KRAS_Project/BATTLE2/GSE33072_series_matrix.txt",header=TRUE,skip=36,nrow=30))

# get rid of the KRAS status NA
KRAS_battle <- unlist(apply(battle_CLIN,2,function(x) {x[grep(pattern="kras mutation:",x=x)]}))
table(KRAS_battle)
names(table(KRAS_battle))
battle_CLIN <- battle_CLIN[,names(which(KRAS_battle!="kras mutation: NA"))]
KRAS_battle <- unlist(apply(battle_CLIN,2,function(x) {x[grep(pattern="kras mutation:",x=x)]}))
table(KRAS_battle)

blah <- battle_CLIN[,names(which(KRAS_battle=="kras mutation: Mutant"))]
codon <- unlist(apply(blah,2,function(x) {x[grep(pattern="codon",x=x)]}))
table(codon)

aa <-  unlist(apply(blah,2,function(x) {x[grep(pattern="replaced",x=x)]}))
table(aa)

tmp <- rbind(codon,aa)
tmp[1,] <- substr(tmp[1,],nchar(tmp[1,])-1,nchar(tmp[1,]))
tmp[2,] <- substr(tmp[2,],nchar(tmp[2,]),nchar(tmp[2,]))

KRAS_battle[colnames(tmp)] <- apply(tmp,2,function(x){paste(x[1],x[2],sep="")})
rm(tmp)
KRAS_battle[KRAS_battle=="kras mutation: WT"] <- "WT"
KRAS_battle <- KRAS_battle[KRAS_battle!="kras mutation: Mutant"]
tmp <- ifelse(KRAS_battle %in% c("10r","12A","12r","13D","61r"),"rare",KRAS_battle)
names(tmp) <- names(KRAS_battle)
KRAS_battle <- tmp
rm(tmp)
KRAS_battle[KRAS_battle=="12C"] <- "G12C"
KRAS_battle[KRAS_battle=="12V"] <- "G12V"
KRAS_battle[KRAS_battle=="12D"] <- "G12D"
table(KRAS_battle)
rm(blah,codon,aa,battle_CLIN)

# make coherent the battle_EXP and the KRAS_battle
tmp <- intersect(colnames(BATTLE_EXP),names(KRAS_battle))
BATTLE_EXP <- BATTLE_EXP[,tmp]
KRAS_BATTLE <- KRAS_battle[tmp]
rm(tmp,battle_EXP,KRAS_battle)


###################################################################################
# perform supervised normalization using snm (and rma summarization)
###################################################################################
rawdata <- rawdata[,names(KRAS_BATTLE)]

#tmp <- match(names(KRAS_BATTLE),sampleNames(rawdata))
SCANBATCH <- substr(rawdata@protocolData@data$ScanDate,1,8)
SCANBATCH <- as.factor(as.Date(SCANBATCH,format="%m/%d/%y"))
table(SCANBATCH)

# we have a  problem because the date of processing of the cel files is big confounder of the expression
s <- svd(exprs(rawdata))
plot(s$v[,1],s$v[,2],col=rainbow(15)[SCANBATCH],pch=20,cex=3)
plot(s$v[,3],s$v[,4],col=rainbow(15)[SCANBATCH],pch=20,cex=3)

# we know that KRAS are biological & study variables of interest 
bio.var <- model.matrix(~ KRAS_BATTLE)

# SCANBATCH is a variable concatenating the study name and the probable batch (grouped according to the cel files date of production)
adj.var <- model.matrix(~ SCANBATCH )

myobject <- log2(pm(rawdata))
snm.fit <- snm(myobject, 
               bio.var=bio.var, 
               adj.var=adj.var, 
               rm.adj=TRUE)
new.expr <- snm.fit$norm.dat
pm(rawdata) <- 2^new.expr
myNormSummarized <- rma(rawdata, background=F, normalize=F)
dim(myNormSummarized)
tmp <- "Battle_snm"
assign(tmp,myNormSummarized)

# retrieve the feature names 
library(hugene10sttranscriptcluster.db)
tmp <- unlist(mget(x=featureNames(Battle_snm),hugene10sttranscriptclusterSYMBOL,ifnotfound=NA))

combine_probes_2_gene <- function(expr, genes, method="svd"){
  
  if(is.list(genes)) genes <- unlist(genes)
  
  stopifnot(dim(expr)[1] ==  length(genes))
  ugenes <- unique(genes)
  ugenes <- sort(ugenes[!is.na(ugenes)])
  M <- matrix(NaN, ncol=dim(expr)[2],nrow=length(ugenes),
              dimnames=list(ugenes, colnames(expr)))
  
  for(gene in ugenes){
    sub.expr <- as.matrix(expr[which(genes == gene),])
    if(dim(sub.expr)[2] == 1){
      M[gene,] <- sub.expr
    }else{
      tmp <- svd(sub.expr - rowMeans(sub.expr))$v[,1]
      tmp.c <- mean(cor(tmp, t(sub.expr)))
      #cat(gene," ", tmp.c, "\n")
      multiplier <- ifelse(tmp.c < 0, -1, 1)
      M[gene,] <- tmp * multiplier
    }
  }
  M
}

Battle_snm <- combine_probes_2_gene(expr=Battle_snm,genes=tmp)
colnames(Battle_snm) <- sampleNames(rawdata)

############################################################################################
# save this in synapse
############################################################################################
# save battle exp
battle_snm <- Data(list(name = "BATTLE_SNM", parentId = 'syn1337457'))
battle_snm <- createEntity(battle_snm)

# add object into the data entity
battle_snm <- addObject(battle_snm,Battle_snm)

# push the raw data into this entity
battle_snm <- storeEntity(entity=battle_snm)


# save kras battle
kras_BATTLE <- Data(list(name = "KRAS_BATTLE", parentId = 'syn1337457'))
kras_BATTLE <- createEntity(kras_BATTLE)

# add object into the data entity
kras_BATTLE <- addObject(kras_BATTLE,KRAS_BATTLE)

# push the raw data into this entity
kras_BATTLE <- storeEntity(entity=kras_BATTLE)
