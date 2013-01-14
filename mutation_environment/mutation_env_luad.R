## charles ferté,MD
## Sage Bionetworks 
## "2012-07-01"

## KRAS project LUAD

## mutational profile in LUAD & CRC
require(synapseClient)
require(HDclassif)
require(mclust)
require(Biobase)
require(ggplot2)
require(corpcor)
require(survival)
library(affy)
library(corpcor)
library(lattice)
library(snm)
library(WGCNA)
## synapse Login
synapseLogin("charles.ferte@sagebase.org","charles")

### magic option
options(stringsAsFactors=FALSE)


#########################################################################################################################
#########compute significant mutual exclusivity and  co-occurring mutations with KRAS  in LUAD (fischer exact test)
#########################################################################################################################
## load the mutations LUAD
load("/home/cferte/FELLOW/cferte/KRAS_Project/LUAD/mutations_LUAD/mutations_LUAD.RData")
KRAS_LUAD <- as.numeric(ifelse(apply(MATMUT_LUAD[grep("KRAS",rownames(MATMUT_LUAD)),],2,sum)==0,0,1))
# 
# GREATER <- apply(MATMUT_LUAD,1,function(x){ fisher.test(KRAS_LUAD,as.numeric(x),alternative="greater")$p.value})
# LESS <- apply(MATMUT_LUAD,1,function(x){ fisher.test(KRAS_LUAD,as.numeric(x),alternative="less")$p.value})
# 
# G_LUAD <- names(GREATER)[which(GREATER<.05)]
# L_LUAD <- names(LESS)[which(LESS<.05)]
# 
# L_LUAD <- substr(L_LUAD,1,nchar(L_LUAD)-7)
# paste(L_LUAD,collapse=" ")
# 
# G_LUAD <- substr(G_LUAD,1,nchar(G_LUAD)-4)
# G_LUAD <- sapply(strsplit(G_LUAD,"_"),function(x){x[[1]]})
# paste(setdiff(G_LUAD,"KRAS"),collapse=" ")

###############################################################################################
# Work specifically the KRAS_G12C in LUAD
###############################################################################################

rownames(MATMUT_LUAD)[grep("KRAS",rownames(MATMUT_LUAD))]
apply(MATMUT_LUAD[grep("KRAS",rownames(MATMUT_LUAD)),],1,sum)

G12C <- MATMUT_LUAD["KRAS_rs121913530_C_A",]
G12V <- MATMUT_LUAD["KRAS_rs121913529_C_A",]
G12D <- MATMUT_LUAD["KRAS_rs121913529_C_T",]
G13D <- MATMUT_LUAD["KRAS_rs112445441_C_T",]
RARE <- apply(MATMUT_LUAD[c("KRAS_rs121913240_T_A","KRAS_rs121913535_C_A","KRAS_NA_T_C","KRAS_rs121913529_C_G","KRAS_rs121913529_CC_AA","KRAS_rs121913527_C_G"),],2,sum)

foo <- rbind(G12C,G12V)
foo <- rbind(foo,G12D)
foo <- rbind(foo,G13D)
foo <- rbind(foo,RARE)
rownames(foo) <- c("G12C","G12V","G12D","G13D","rare")

foo2 <- as.data.frame(t(foo))
WT <- ifelse(apply(foo2,1,sum)==0,1,0)
foo2$WT <- WT

foo2[which(apply(foo2,1,sum)==2),"rare"] <- 1

KRAS_LUAD <- unlist(apply(foo2,1,function(x){colnames(foo2)[which(x==1)]}))
rm(G12C,G12D,G12V,G13D,foo,foo2,RARE,WT)
save(KRAS_LUAD,file="/home/cferte/FELLOW/cferte/KRAS_Project/OBJECTS/KRAS_LUAD_STATUS.RData")

MUTRATE <- apply(MATMUT_LUAD,2,sum)

###################################################################################################
# plot the overall mutation rate in LUAD (y axis = mutations per megabase, x axis= )
###################################################################################################

par(mfrow=c(1,2))
x <- sort(MUTRATE,decreasing=TRUE)*94771/1000000
summary(x)
#plot(x,pch=20,cex=.01,xlab="249 LUAD TCGA samples",
     #ylab="Overall number of mutations per Mb")
#abline(h=quantile(x=x,probs=.15),lty=2,col="red")
#abline(h=quantile(x=x,probs=.75),lty=2,col="red")

plot(log10(sort(MUTRATE,decreasing=TRUE)*94771/1000000),pch=20,cex=.5,xlab="249 LUAD TCGA samples",
     ylab="Overall number of mutations per Mb (log10 scale)",yaxt="n")
axis(side=2,at=c(.5,1,1.5,2,2.5,3),labels=c(3,10,30,100,300,1000))
abline(h=log10(quantile(x=x,probs=.75)),lty=2,col="red")
abline(h=log10(quantile(x=x,probs=.25)),lty=2,col="red")

rm(x)

###################################################################################################
# KRAS mutation is associated with higher overall mutation rate
###################################################################################################

#par(mfrow=c(1,2),oma=c(0,0,3,0))
par(mfrow=c(1,1))
tmp <- ifelse(KRAS_LUAD=="WT","Wild-type 
(n=189)","Mutant
  (n=60)")
boxplot(MUTRATE~tmp,outline=FALSE,ylab="Overall mutation rate (N)")
title(main= "KRAS mutant samples are associated with 
higher overall mutation rate in lung adenocarcinoma",sub=paste("P =",format(wilcox.test(MUTRATE[which(KRAS_LUAD!="WT")],MUTRATE[which(KRAS_LUAD=="WT")])$p.value,digits=2), "(wilcoxon test)"))

# ### KRAS G12C mutation is associated with the highest overall mutation rate
# boxplot(MUTRATE~as.factor(KRAS),outline=FALSE,cex.axis=.9,ylab="mutation load (N)")
# title("mutation load in lung adenocarcinoma
# according to KRAS mutation status",outer=TRUE,cex.main=1.3)
# blah <- table(KRAS)
# blah <- paste("n=",blah)
# blah <- paste(blah,collapse="        ")
# title(sub=blah)

### KRAS G12C mutation is associated with the highest overall mutation rate
### as compared to G12V and to the rare variants
plot.new()
tmp <- ifelse(KRAS_LUAD=="G12C","G12C",ifelse(KRAS_LUAD=="WT","WT",ifelse(KRAS_LUAD=="G12V","G12V","rare variants")))

boxplot(MUTRATE~tmp,outline=FALSE)
a <- wilcox.test(MUTRATE[which(tmp=="G12C")],MUTRATE[which(tmp=="WT")])
b <- wilcox.test(MUTRATE[which(tmp=="G12C")],MUTRATE[which(tmp=="G12V")])
c <- wilcox.test(MUTRATE[which(tmp=="G12C")],MUTRATE[which(tmp=="rare variants")])
d <- wilcox.test(MUTRATE[which(tmp=="G12V")],MUTRATE[which(tmp=="WT")])
e <- wilcox.test(MUTRATE[which(tmp=="rare variants")],MUTRATE[which(tmp=="WT")])
title(main="KRAS G12C mutation is associated 
with the highest overall mutation rate")
title(sub=paste(paste("G12C vs. WT","p =",format(a$p.value,digits=2)),
                paste("G12V vs. WT","p =",format(d$p.value,digits=2)),
                #paste("rare vs. WT","p =",format(e$p.value,digits=2)),
                sep= "     "),
                cex.sub=.75)

# plot.new()    
# par(mfrow=c(1,1), oma= c(0,0,4,0))
# tmp <- ifelse(KRAS=="G12C","G12C",ifelse(KRAS=="WT","WT","Other KRAS mutation"))
# boxplot(MUTRATE~tmp,outline=FALSE)
# title("mutation load in lung adenocarcinoma
# according to KRAS mutation status 
#   (G12C vs. Other KRAS mutation vs. WT",outer=TRUE,cex.main=1.3)
# ab <- wilcox.test(MUTRATE[which(tmp=="G12C")],MUTRATE[which(tmp=="Other KRAS mutation")])
# bc <- wilcox.test(MUTRATE[which(tmp=="WT")],MUTRATE[which(tmp=="Other KRAS mutation")])
# ac <- wilcox.test(MUTRATE[which(tmp=="G12C")],MUTRATE[which(tmp=="WT")])

summary(MUTRATE[which(tmp=="G12C")])
summary(MUTRATE[which(tmp=="G12V")])
summary(MUTRATE[which(tmp=="rare variants")])
summary(MUTRATE[which(tmp=="WT")])
par(mfrow=c(2,2))
plot(density(MUTRATE[which(tmp=="G12C")]),main="G12C")
plot(density(MUTRATE[which(tmp=="G12V")]),main="G12V")
plot(density(MUTRATE[which(tmp=="WT")]),main="WT")
rm(a,b,c,d,e,tmp)

###################################################################################################
# compute the mutations that are associated with G12C specifically
###################################################################################################

load("/home/cferte/FELLOW/cferte/KRAS_Project/LUAD/mutations_LUAD/mutations_LUAD.RData")


G12C <- MATMUT_LUAD["KRAS_rs121913530_C_A",]
#BOTH <- apply(MATMUT_LUAD,1,function(x){ fisher.test(as.numeric(G12C),as.numeric(x))$p.value})
GREATER <- apply(MATMUT_LUAD,1,function(x){ fisher.test(as.numeric(G12C),as.numeric(x),alternative="greater")$p.value})
LESS <- apply(MATMUT_LUAD,1,function(x){ fisher.test(as.numeric(G12C),as.numeric(x),alternative="less")$p.value})

G12C_OVERLAP_PVAL <- GREATER
names(G12C_OVERLAP_PVAL) <-  sapply(strsplit(x=names(GREATER),split="_"),function(x){x[[1]]})
setwd("/gluster/home/cferte/FELLOW/cferte/KRAS_Project/OBJECTS/MUT_OVERLAP/")
save(G12C_OVERLAP_PVAL,file="G12C_OVERLAP_PVAL.RData")

## enter the list of G12C mutations in the GSEA crawler (in sock)

## which of the overlapping mutations have a prevalence of more than 5 occurence out of the 249 pts
G_G12C <- GREATER[which(GREATER<.01)]
sort(G_G12C)
names(G_G12C) <-  sapply(strsplit(x=names(G_G12C),split="_"),function(x){x[[1]]})
paste(names(G_G12C),collapse=" ")

L_G12C <- names(LESS)[which(LESS<.05)]
Try <- G_G12C[which(G>5)]
Try <- Try[-grep("KRAS",Try)]

###################################################################################################
# compute the mutations that are associated with G12V specifically in LUAD
###################################################################################################

G12V <- MATMUT_LUAD["KRAS_rs121913529_C_A",]
GREATER <- apply(MATMUT_LUAD,1,function(x){ fisher.test(as.numeric(G12V),as.numeric(x),alternative="greater")$p.value})
LESS <- apply(MATMUT_LUAD,1,function(x){ fisher.test(as.numeric(G12V),as.numeric(x),alternative="less")$p.value})

G12V_OVERLAP_PVAL <- GREATER
names(G12V_OVERLAP_PVAL) <-  sapply(strsplit(x=names(GREATER),split="_"),function(x){x[[1]]})
setwd("/gluster/home/cferte/FELLOW/cferte/KRAS_Project/OBJECTS/MUT_OVERLAP/LUAD_MUT_OVERLAP/")
save(G12V_OVERLAP_PVAL,file="G12V_OVERLAP_PVAL.RData")

G_G12V <- names(GREATER)[which(GREATER<.05)]
L_G12V <- names(LESS)[which(LESS<.05)]

G_G12V <- substr(G_G12V,1,nchar(G_G12V)-4)
G_G12V <- sapply(strsplit(G_G12V,"_"),function(x){x[[1]]})
paste(G_G12V,collapse=" ")

###################################################################################################
# compute the diffrences between G12C and G12V in the LUAD mutational environment
###################################################################################################
G12C <- as.numeric(MATMUT_LUAD["KRAS_rs121913530_C_A",])
G12V <- as.numeric(MATMUT_LUAD["KRAS_rs121913529_C_A",])
tmp <- as.matrix(MATMUT_LUAD[,c(which(G12C==1),which(G12V==1))])
tmp <- tmp[which(apply(tmp,1,sum)>1),]
#tmp <- tmp[-grep("KRAS",rownames(tmp)),]
vec <- c(rep(1,times=30),rep(0,times=15))

#GREATER <- apply(tmp,1,function(x){ fisher.test(vec,as.numeric(x),alternative="greater")$p.value})
LESS <- apply(tmp,1,function(x){ fisher.test(vec,as.numeric(x))$p.value})
#names(LESS) <- sapply(strsplit(x=names(LESS),split="_"),function(x){x[[1]]})
par(mfrow=c(1,1))
plot(density(LESS))

# take the top 5% point mutations (first 300) and compute the overlap (fisher exact test) on the gsea website
paste(sapply(strsplit(x=names(sort(LESS)[1:300]),split="_"),function(x){x[[1]]}),collapse=" ")

# # store LESS in Synapse
# library(synapseClient)
# synapseLogin("charles.ferte@sagebase.org","charles")
# newLayer <- Data(list(name = "G12V_G12C_LUAD_MUT_PVAL", parentId = "syn1203262"))
# newLayer <- addObject(newLayer, LESS)
# newLayer <- storeEntity(newLayer)

# draw a heatmap with the top 300 features
tmp <- as.matrix(MATMUT_LUAD[,c(which(G12C==1),which(G12V==1))])
tmp <- tmp[which(apply(tmp,1,sum)>1),]
#tmp <- tmp[-grep("KRAS",rownames(tmp)),]
par(oma=c(5,2,1,2))
tmp2 <- tmp
tmp2 <- tmp2[sort(LESS,index.return=T)$ix[1:400],]
rownames(tmp2) <- sapply(strsplit(x=rownames(tmp2),split="_"),function(x){x[[1]]})
require(gplots)
heatmap.2(tmp2,trace="none",col=greenred(100),dendrogram="row",ColSideColors=c("blue","red")[vec+1],Colv=NULL)

# draw the heatmap with the first pathways associated with G12V
par(oma=c(5,2,1,10))
tmp3 <- tmp
tmp3 <- tmp3[sort(LESS,index.return=T)$ix[1:400],]
tmp3 <- tmp3[c(grep("TP53",rownames(tmp3))[1:2],grep("ATM",rownames(tmp3)),grep("MYT1",rownames(tmp3)),grep("DAXX",rownames(tmp3)),grep("PML",rownames(tmp3))),]
#rownames(tmp3) <- sapply(strsplit(x=rownames(tmp3),split="_"),function(x){x[[1]]})
require(gplots)
heatmap.2(tmp3,trace="none",col=greenred(10),dendrogram="none",ColSideColors=c("blue","red")[vec+1],Colv=NULL)

# draw the heatmap with the first pathways associated with G12C
par(oma=c(5,2,1,10))
tmp <- as.matrix(MATMUT_LUAD[,c(which(G12C==1),which(G12V==1))])
tmp <- tmp[which(apply(tmp,1,sum)>1),]
tmp3 <- tmp
#tmp3 <- tmp3[sort(LESS,index.return=T)$ix[1:1000],]
tmp3 <- tmp3[c(grep("SHC2",rownames(tmp3)),grep("SOS",rownames(tmp3)),grep("GRB2",rownames(tmp3)),grep("MAP3K1",rownames(tmp3)),grep("TP53",rownames(tmp3))[1:2],grep("ATM",rownames(tmp3)),grep("MYT1",rownames(tmp3)),grep("DAXX",rownames(tmp3)),grep("PML",rownames(tmp3))),]
tmp3 <- tmp3[which(MUT_TYPE_LUAD$Variant_Classification[match(rownames(tmp3),rownames(tmp))] %in% c("Missense_Mutation", "Silent","Nonsense_Mutation")),]
require(gplots)
heatmap.2(tmp3,trace="none",col=greenred(10),ColSideColors=c("blue","red")[vec+1],Colv=NULL,Rowv=NULL)

### retrievez the data on functionnal or not mutations
load("/home/cferte/FELLOW/cferte/KRAS_Project/LUAD/mutations_LUAD/MUT_TYPE_LUAD.RData")




###################################################################################################
# sort out the first pathways associated with G12C  specifically in LUAD
###################################################################################################

load("/home/cferte/FELLOW/cferte/KRAS_Project/OBJECTS/MUT_OVERLAP/LUAD_MUT_OVERLAP/G12C_OVERLAP_MUT_GSEA.RData")
G12C_LUAD_GSEA <- analyticResult
G12C_PTW <- k
rm(analyticResult,k)

# load the additional subsequent data (10^6 permut) from synapse and replace them into the list
G12C_MILLION_PERMUT <- loadEntity('syn1419605')
ADD <- G12C_MILLION_PERMUT$objects$result$k
G12C_PTW[match(names(ADD),names(G12C_PTW))]  <- ADD
hist(G12C_PTW)

#  p adjust (BH) on the pathways
S <- p.adjust(G12C_PTW,method="BH")
hist(S, main="adjusted BH")
GOLD <- sort(S[which(S<.01)])
names(GOLD)

G12C_MUT_LUAD_PROFILE <- as.data.frame(sort(S[which(S<.01)]))
colnames(G12C_MUT_LUAD_PROFILE) <- "p value"

setwd("/gluster/home/cferte/FELLOW/cferte/KRAS_Project/OBJECTS/MUT_OVERLAP/LUAD_MUT_OVERLAP/")
write.table(G12C_MUT_LUAD_PROFILE,file="G12C_MUT_LUAD_PROFILE.txt")

##################################################################################
# what are the mutations overalapping with G12C within these pathways
##################################################################################
load("/home/cferte/FELLOW/cferte/KRAS_Project/LUAD/mutations_LUAD/mutations_LUAD.RData")
load("/home/cferte/FELLOW/cferte/KRAS_Project/OBJECTS/MUT_OVERLAP/LUAD_MUT_OVERLAP/G12C_OVERLAP_PVAL.RData")

# #restrict to the potentially fonctionnal mutations
# load("/home/cferte/FELLOW/cferte/KRAS_Project/LUAD/mutations_LUAD/mutations_LUAD.RData")
# load("/home/cferte/FELLOW/cferte/KRAS_Project/LUAD/mutations_LUAD/MUT_TYPE_LUAD.RData")
# names(table(MUT_TYPE_LUAD$Variant_Classification))
# NONFON <- which(MUT_TYPE_LUAD$Variant_Classification %in% c("3'UTR","5'UTR","Silent","Intron"))
# MATMUT_LUAD <- MATMUT_LUAD[-NONFON,]

G12C <- as.numeric(MATMUT_LUAD["KRAS_rs121913530_C_A",])
G12C <- as.numeric(MATMUT_LUAD["KRAS_rs121913530_C_A",])
G12V <- as.numeric(MATMUT_LUAD["KRAS_rs121913529_C_A",])
G12D <- as.numeric(MATMUT_LUAD["KRAS_rs121913529_C_T",])
G13D <- as.numeric(MATMUT_LUAD["KRAS_rs112445441_C_T",])

# # small code to explore the pathways (usefull)
# blah <- strsplit(split="\n",x="paste the genes of the pathway in txt")
# PTT <- unlist(blah)
# 
# blah <- c()
# for(i in PTT){
#   blah <- apply(MATMUT_LUAD[which(names(G12C_OVERLAP_PVAL)==i),],2,sum)
#   print(i)
#   print(table(G12C,as.numeric(ifelse(blah==0,0,1))))
#   #fisher.test(G12C,as.numeric(ifelse(blah==0,0,1)))
# }

# G12C mutations are concomitant with SOS1 mutations (interesting)
blah <- apply(MATMUT_LUAD[which(names(G12C_OVERLAP_PVAL)=="SOS1"),],2,sum)
SOS1 <- as.numeric(ifelse(blah==0,0,1))
table(G12C,SOS1)
fisher.test(G12C,SOS1,alternative="greater")

# G12C mutations are concomitant with RAPGEF1 mutations (interesting)
blah <- apply(MATMUT_LUAD[which(names(G12C_OVERLAP_PVAL)=="RAPGEF1"),],2,sum)
RAPGEF1 <- as.numeric(ifelse(blah==0,0,1))
table(G12C,RAPGEF1)
fisher.test(G12C,RAPGEF1,alternative="greater")

# G12C mutations are concomitant with HRAS mutations (interesting)
blah <- apply(MATMUT_LUAD[which(names(G12C_OVERLAP_PVAL)=="HRAS"),],2,sum)
HRAS <- as.numeric(ifelse(blah==0,0,1))
table(G12C,HRAS)
fisher.test(G12C,HRAS,alternative="greater")

# G12C mutations are concomitant with NRAS mutations (interesting)
blah <- apply(MATMUT_LUAD[which(names(G12C_OVERLAP_PVAL)=="NRAS"),],2,sum)
NRAS <- as.numeric(ifelse(blah==0,0,1))
table(G12C,NRAS)
fisher.test(G12C,NRAS,alternative="greater")

# G12C mutations are concomitant with RAP1A mutations (interesting)
blah <- apply(MATMUT_LUAD[which(names(G12C_OVERLAP_PVAL)=="RAP1A"),],2,sum)
RAP1A <- as.numeric(ifelse(blah==0,0,1))
table(G12C,RAP1A)
fisher.test(G12C,RAP1A,alternative="greater")

# G12C mutations are concomitant with DUSP4 mutations (interesting)
blah <- apply(MATMUT_LUAD[which(names(G12C_OVERLAP_PVAL)=="DUSP4"),],2,sum)
DUSP4 <- as.numeric(ifelse(blah==0,0,1))
table(G12C,DUSP4)
fisher.test(G12C,DUSP4,alternative="greater")

# G12C mutations are concomitant with DUSP6 mutations (interesting)
blah <- apply(MATMUT_LUAD[which(names(G12C_OVERLAP_PVAL)=="DUSP6"),],2,sum)
DUSP6 <- as.numeric(ifelse(blah==0,0,1))
table(G12C,DUSP6)
fisher.test(G12C,DUSP6,alternative="greater")

# G12C mutations are concomitant with NGF mutations (interesting)
blah <- apply(MATMUT_LUAD[which(names(G12C_OVERLAP_PVAL)=="NGF"),],2,sum)
NGF <- as.numeric(ifelse(blah==0,0,1))
table(G12C,NGF)
fisher.test(G12C,NGF,alternative="greater")

# G12C mutations are concomitant with RIT2 mutations (interesting)
blah <- apply(MATMUT_LUAD[which(names(G12C_OVERLAP_PVAL)=="RIT2"),],2,sum)
RIT2 <- as.numeric(ifelse(blah==0,0,1))
table(G12C,RIT2)
fisher.test(G12C,RIT2,alternative="greater")

# G12C mutations are concomitant with RIT1 mutations (interesting)
blah <- apply(MATMUT_LUAD[which(names(G12C_OVERLAP_PVAL)=="RIT1"),],2,sum)
RIT1 <- as.numeric(ifelse(blah==0,0,1))
table(G12C,RIT1)
fisher.test(G12C,RIT1,alternative="less")

# G12C mutations are concomitant with BRAF mutations (interesting and not expected)
blah <- apply(MATMUT_LUAD[which(names(G12C_OVERLAP_PVAL)=="BRAF"),],2,sum)
BRAF <- as.numeric(ifelse(blah==0,0,1))
table(G12C,BRAF)
fisher.test(G12C,BRAF,alternative="greater")

# G12C mutations are concomitant with CCNB2 mutations (interesting)
blah <- apply(MATMUT_LUAD[which(names(G12C_OVERLAP_PVAL)=="CCNB2"),],2,sum)
CCNB2 <- as.numeric(ifelse(blah==0,0,1))
table(G12C,CCNB2)
fisher.test(G12C,CCNB2,alternative="greater")

# G12C mutations are concomitant with MAP2K1 mutations (interesting)
blah <- apply(MATMUT_LUAD[which(names(G12C_OVERLAP_PVAL)=="MAP2K1"),],2,sum)
MAP2K1 <- as.numeric(ifelse(blah==0,0,1))
table(G12C,MAP2K1)
fisher.test(G12C,MAP2K1,alternative="greater")

# G12C mutations are concomitant with MAP2K2 mutations (interesting)
blah <- apply(MATMUT_LUAD[which(names(G12C_OVERLAP_PVAL)=="MAP2K2"),],2,sum)
MAP2K2 <- as.numeric(ifelse(blah==0,0,1))
table(G12C,MAP2K2)
fisher.test(G12C,MAP2K2,alternative="greater")

# G12C mutations are concomitant with MAP3K1 mutations (interesting)
blah <- apply(MATMUT_LUAD[which(names(G12C_OVERLAP_PVAL)=="MAP3K1"),],2,sum)
MAP3K1 <- as.numeric(ifelse(blah==0,0,1))
table(G12C,MAP3K1)
fisher.test(G12C,MAP3K1,alternative="greater")

# G12C mutations are concomitant with MAPK1 mutations (interesting)
blah <- apply(MATMUT_LUAD[which(names(G12C_OVERLAP_PVAL)=="MAPK1"),],2,sum)
MAPK1 <- as.numeric(ifelse(blah==0,0,1))
table(G12C,MAPK1)
fisher.test(G12C,MAPK1,alternative="greater")

# G12C mutations are concomitant with MAPK3 mutations (interesting)
blah <- apply(MATMUT_LUAD[which(names(G12C_OVERLAP_PVAL)=="MAPK3"),],2,sum)
MAPK3 <- as.numeric(ifelse(blah==0,0,1))
table(G12C,MAPK3)
fisher.test(G12C,MAPK3,alternative="greater")

# G12C mutations are concomitant with SHC2 mutations (interesting)
blah <- apply(MATMUT_LUAD[which(names(G12C_OVERLAP_PVAL)=="SHC2"),],2,sum)
SHC2 <- as.numeric(ifelse(blah==0,0,1))
table(G12C,SHC2)
fisher.test(G12C,SHC2,alternative="greater")

# G12C mutations are concomitant with SHC1 mutations (interesting)
blah <- apply(MATMUT_LUAD[which(names(G12C_OVERLAP_PVAL)=="SHC1"),],2,sum)
SHC1 <- as.numeric(ifelse(blah==0,0,1))
table(G12C,SHC1)
fisher.test(G12C,SHC1,alternative="greater")

# G12C mutations appear mutually exclusive with GRB2 mutations (interesting)
blah <- apply(MATMUT_LUAD[which(names(G12C_OVERLAP_PVAL)=="GRB2"),],2,sum)
GRB2 <- as.numeric(ifelse(blah==0,0,1))
table(G12C,GRB2)
fisher.test(G12C,GRB2)

# G12C mutations appear mutually exclusive with IRS2 mutations (interesting)
blah <- apply(MATMUT_LUAD[which(names(G12C_OVERLAP_PVAL)=="IRS2"),],2,sum)
IRS2 <- as.numeric(ifelse(blah==0,0,1))
table(G12C,IRS2)
fisher.test(G12C,IRS2,alternative="less")

# G12C mutations appear mutually exclusive with IRS1 mutations (interesting)
blah <- apply(MATMUT_LUAD[which(names(G12C_OVERLAP_PVAL)=="IRS1"),],2,sum)
IRS1 <- as.numeric(ifelse(blah==0,0,1))
table(G12C,IRS1)
fisher.test(G12C,IRS1,alternative="less")

# G12C mutations appear mutually exclusive with FRS2 mutations (interesting)
blah <- apply(MATMUT_LUAD[which(names(G12C_OVERLAP_PVAL)=="FRS2"),],2,sum)
FRS2 <- as.numeric(ifelse(blah==0,0,1))
table(G12C,FRS2)
fisher.test(G12C,FRS2,alternative="less")

# G12C mutations appear mutually exclusive with RAF1 mutations (interesting)
blah <- apply(MATMUT_LUAD[which(names(G12C_OVERLAP_PVAL)=="RAF1"),],2,sum)
RAF1 <- as.numeric(ifelse(blah==0,0,1))
table(G12C,RAF1)
fisher.test(G12C,RAF1,alternative="less")

# G12C mutations appear mutually exclusive with AKT21 mutations (interesting)
blah <- apply(MATMUT_LUAD[which(names(G12C_OVERLAP_PVAL)=="AKT2"),],2,sum)
AKT2 <- as.numeric(ifelse(blah==0,0,1))
table(G12C,AKT2)
fisher.test(G12C,AKT2,alternative="less")

# G12C mutations appear mutually exclusive with MTOR mutations (interesting)
blah <- apply(MATMUT_LUAD[which(names(G12C_OVERLAP_PVAL)=="MTOR"),],2,sum)
MTOR <- as.numeric(ifelse(blah==0,0,1))
table(G12C,MTOR)
fisher.test(G12C,MTOR,alternative="less")

# G12C mutations appear mutually exclusive with RPTOR mutations (interesting)
blah <- apply(MATMUT_LUAD[which(names(G12C_OVERLAP_PVAL)=="RPTOR"),],2,sum)
RPTOR <- as.numeric(ifelse(blah==0,0,1))
table(G12C,RPTOR)
fisher.test(G12C,RPTOR,alternative="less")

# G12C mutations appear mutually exclusive with STK11 mutations (interesting)
blah <- apply(MATMUT_LUAD[which(names(G12C_OVERLAP_PVAL)=="STK11"),],2,sum)
STK11 <- as.numeric(ifelse(blah==0,0,1))
table(G12C,STK11)
fisher.test(G12C,STK11,alternative="less")

# G12C mutations appear mutually exclusive with PIK3CA mutations (interesting)
blah <- apply(MATMUT_LUAD[which(names(G12C_OVERLAP_PVAL)=="PIK3CA"),],2,sum)
PIK3CA <- as.numeric(ifelse(blah==0,0,1))
table(G12C,PIK3CA)
fisher.test(G12C,PIK3CA,alternative="less")

# G12C mutations appear mutually exclusive with PIK3CB mutations (interesting)
blah <- apply(MATMUT_LUAD[which(names(G12C_OVERLAP_PVAL)=="PIK3CB"),],2,sum)
PIK3CB <- as.numeric(ifelse(blah==0,0,1))
table(G12C,PIK3CB)
fisher.test(G12C,PIK3CB,alternative="less")

# G12C mutations appear mutually exclusive with FGFR1 mutations (interesting)
blah <- apply(MATMUT_LUAD[which(names(G12C_OVERLAP_PVAL)=="FGFR1"),],2,sum)
FGFR1 <- as.numeric(ifelse(blah==0,0,1))
table(G12C,FGFR1)
fisher.test(G12C,FGFR1,alternative="less")


names(G12C_OVERLAP_PVAL)[grep(pattern="MTOR",names(G12C_OVERLAP_PVAL))]
table(IRS1,G12C) ## mutations of SOS1 and SHC2 are exclusive

table(GRB2,G12C) ## mutations of SOS1 and SHC2 are exclusive
table(SOS1,SHC2) ## mutations of SOS1 and SHC2 are exclusive
table(SOS1,MAP3K1) ## mutations of SOS1 and MAP3K1 are nearly exclusive
table(SHC2,MAP3K1) ## mutations of MAP3K1 and SHC2 are nearly exclusive

# co mutations in the RAS pathway
TMP <- data.frame(G12C,SHC2,GRB2,SOS1)
TMP <- TMP[TMP$G12C==1,]
image(as.matrix(t(TMP)))


########################

TMP <- cbind(G12C,G12V,G12D,G13D,SOS1,SHC2,RIT2,RIT1,GRB2,IRS1,FRS2,IRS2,MAP3K1,MAPK1,MAPK3,MAP2K1,MAP2K2,DUSP4,DUSP6,RAPGEF1,RAP1A,RAF1,NRAS,HRAS,BRAF,STK11,RPTOR,MTOR,AKT2,PIK3CA,PIK3CB,FGFR1)


setwd("/gluster/home/cferte/FELLOW/cferte/KRAS_Project/OBJECTS/MUT_OVERLAP/LUAD_MUT_OVERLAP/")
write.csv(TMP,file="ALL_KRAS_SOS_SHC_GRB2.csv")

library(e1071)
heatmap(hamming.distance(t(TMP)))
hamming.distance(t(TMP))


###################################################################################################
# sort out the first pathways associated with G12V  specifically in LUAD
###################################################################################################

load("/home/cferte/FELLOW/cferte/KRAS_Project/OBJECTS/MUT_OVERLAP/LUAD_MUT_OVERLAP/G12V_OVERLAP_MUT_GSEA.RData")
G12V_LUAD_GSEA <- analyticResult
G12V_PTW <- k
rm(analyticResult,k)

# load the additional subsequent data (10^6 permut) from synapse and replace them into the list
G12V_MILLION_PERMUT <- loadEntity('syn1419601')
ADD <- G12V_MILLION_PERMUT$objects$result$k
G12V_PTW[match(names(ADD),names(G12V_PTW))]  <- ADD
hist(G12V_PTW)

#  p adjust (BH) on the pathways
S <- p.adjust(G12V_PTW,method="BH")
hist(S, main="adjusted BH")
sort(S[which(S<.05)])

G12V_MUT_LUAD_PROFILE <- as.data.frame(sort(S[which(S<.05)]))
colnames(G12V_MUT_LUAD_PROFILE) <- "p value"

setwd("/gluster/home/cferte/FELLOW/cferte/KRAS_Project/OBJECTS/MUT_OVERLAP/LUAD_MUT_OVERLAP/")
write.table(G12V_MUT_LUAD_PROFILE,file="G12V_MUT_LUAD_PROFILE.txt")

##################################################################################
# what are the mutations overalapping with G12V within these pathways
##################################################################################
load("/home/cferte/FELLOW/cferte/KRAS_Project/LUAD/mutations_LUAD/mutations_LUAD.RData")
load("/home/cferte/FELLOW/cferte/KRAS_Project/OBJECTS/MUT_OVERLAP/LUAD_MUT_OVERLAP/G12V_OVERLAP_PVAL.RData")
length(which(G12V_OVERLAP_PVAL<.05))
test <- read.delim("/home/cferte/FELLOW/cferte/KRAS_Project/OBJECTS/MUT_OVERLAP/LUAD_MUT_OVERLAP/geneset.txt",skip=1)


# G12V mutations are concomitant with mutations in the ATR pathway (interesting)
for(i in test[1:dim(test)[1],]){
  blah <-  apply(MATMUT_LUAD[which(names(G12V_OVERLAP_PVAL)==i),],2,sum)
  blah <- ifelse(blah==0,0,1)
  print(i)
  print(table(as.numeric(MATMUT_LUAD["KRAS_rs121913530_C_A",]),as.numeric(blah)))
  print(fisher.test(as.numeric(MATMUT_LUAD["KRAS_rs121913530_C_A",]),as.numeric(blah),alternative="greater")$p.value)
  }

# G12V mutations are concomitant with mutations in the FGFR pathway (interesting)
for(i in test[1:dim(test)[1],]){
  blah <-  apply(MATMUT_LUAD[which(names(G12V_OVERLAP_PVAL)==i),],2,sum)
  blah <- ifelse(blah==0,0,1)
  print("###################")
  print(i)
  d <- table(as.numeric(MATMUT_LUAD["KRAS_rs121913530_C_A",]),as.numeric(blah))
  print(d)
  ifelse(dim(d)[2]==2,print(fisher.test(as.numeric(MATMUT_LUAD["KRAS_rs121913530_C_A",]),as.numeric(blah))$p.value),"NA")
}


## co muttaions with SOS1, PIK3CA, PIK3C3,FGFR1
## mutual exclusive with SRC, INSR,INS,IRS1,IRS2

blah <- apply(MATMUT_LUAD[which(names(G12V_OVERLAP_PVAL)=="SOS1"),],2,sum)
SOS1 <- ifelse(blah==0,0,1)
table(as.numeric(MATMUT_LUAD["KRAS_rs121913530_C_A",]),as.numeric(SOS1))
fisher.test(as.numeric(MATMUT_LUAD["KRAS_rs121913530_C_A",]),as.numeric(SOS1),alternative="greater")





# G12C mutations are concomitant with SOS1 mutations (interesting)
blah <- apply(MATMUT_LUAD[which(names(G12V_OVERLAP_PVAL)=="SOS1"),],2,sum)
SOS1 <- ifelse(blah==0,0,1)
table(as.numeric(MATMUT_LUAD["KRAS_rs121913530_C_A",]),as.numeric(SOS1))
fisher.test(as.numeric(MATMUT_LUAD["KRAS_rs121913530_C_A",]),as.numeric(SOS1),alternative="greater")








#########################################################################################################################
# load the mutations CRC
#########################################################################################################################
load("/home/cferte/FELLOW/cferte/KRAS_Project/CRC/mutations_CRC.RData")
apply(MATMUT_CRC[grep("KRAS",rownames(MATMUT_CRC)),],1,sum)
KRAS_CRC <- as.numeric(ifelse(apply(MATMUT_CRC[grep("KRAS",rownames(MATMUT_CRC)),],2,sum)==0,0,1))




########################################################################################################################
# correlation between the G12D G13D G12V G12C mutations with overall mutation rate in CRC
########################################################################################################################

apply(MATMUT_CRC[grep("KRAS",rownames(MATMUT_CRC)),],1,sum)
apply(MATMUT_CRC[grep("KRAS",rownames(MATMUT_CRC)),],2,sum)

G12D <- MATMUT_CRC["KRAS_p.G12D",]
G13D <- MATMUT_CRC["KRAS_p.G13D",]
G12V <- MATMUT_CRC["KRAS_p.G12V",]
G12C <- MATMUT_CRC["KRAS_p.G12C",]
A146T <- MATMUT_CRC["KRAS_p.A146T",]
RARE <- apply(MATMUT_CRC[c("KRAS_.","KRAS_p.R68S","KRAS_p.K117N","KRAS_p.G12S","KRAS_p.G12A","KRAS_p.Q22K","KRAS_p.Q61L","KRAS_p.E98X"),],2,sum)

foo <- rbind(G12D,G13D)
foo <- rbind(foo,G12V)
foo <- rbind(foo,G12C)
foo <- rbind(foo,A146T)
foo <- rbind(foo,RARE)
rownames(foo) <- c("G12D","G13D","G12V","G12C","A146T","RARE")

foo2 <- as.data.frame(t(foo))
WT <- ifelse(apply(foo2,1,sum)==0,1,0)
foo2$WT <- WT

KRAS_CRC <- unlist(apply(foo2,1,function(x){colnames(foo2)[which(x==1)]}))
rm(G12V,G12C,G13D,G12D,A146T,RARE,foo,foo2,WT)
setwd("/home/cferte/FELLOW/cferte/KRAS_Project/OBJECTS/")
names(KRAS_CRC) <- substr(names(KRAS_CRC),1,12)
save(KRAS_CRC,file="KRAS_CRC_STATUS.RData")
MUTRATE <- apply(MATMUT_CRC,2,sum)


###################################################################################################
# plot the overall mutation rate in CRC (y axis = mutations per megabase, x axis= 140 samples)
###################################################################################################


plot(log10(sort(MUTRATE,decreasing=TRUE)*69450/1000000),pch=20,cex=.5,xlab=" 140 CRC TCGA samples",
     ylab="Overall number of mutations per Mb (log10 scale)",yaxt="n")
axis(side=2,at=c(.5,1,1.5,2,2.5,3),labels=c(3,10,30,100,300,1000))
x <- sort(MUTRATE,decreasing=TRUE)*69450/1000000
abline(h=log10(quantile(x=x,probs=.75)),lty=2,col="red")
abline(h=log10(quantile(x=x,probs=.25)),lty=2,col="red")


summary(sort(MUTRATE,decreasing=TRUE)*69450/1000000)

hyperMutCRC <- ifelse(MUTRATE*69450/1000000>12,1,0)

###################################################################################################
# KRAS mutation is not associated with higher overall mutation rate in CRC
###################################################################################################

#par(mfrow=c(1,2),oma=c(0,0,3,0))
par(mfrow=c(1,1))
tmp <- ifelse(KRAS_CRC=="WT","Wild-type 
(n=79)","Mutant
  (n=61)")
boxplot(MUTRATE~tmp,outline=FALSE,ylab="Overall mutation rate (N)")
title(main= "KRAS mutant samples are not associated with 
higher overall mutation rate in colorectal cancer",sub=paste("P =",format(wilcox.test(MUTRATE[which(KRAS_CRC!="WT")],MUTRATE[which(KRAS_CRC=="WT")])$p.value,digits=2), "(wilcoxon test)"))

### KRAS mutation is not associated with the highest overall mutation rate in CRC
plot.new()
boxplot(MUTRATE~as.factor(KRAS_CRC),outline=FALSE,cex.axis=.9,ylab="mutation load (N)")
title("mutation load in colorectal
according to KRAS mutation status",outer=TRUE,cex.main=1.3)
blah <- table(KRAS_CRC)
blah <- paste("n=",blah)
blah <- paste(blah,collapse="        ")
title(sub=blah)

### none of the KRAS mutation is significantly associated with higher overall mutation rate

plot.new()
tmp <- ifelse(KRAS_CRC=="G12D","G12D",ifelse(KRAS_CRC=="WT","WT",ifelse(KRAS_CRC=="G12V","G12V",ifelse(KRAS_CRC=="G13D","G13D","rare variants"))))
#tmp <- ifelse(KRAS_CRC=="G12D","G12D",ifelse(KRAS_CRC=="WT","WT",ifelse(KRAS_CRC=="G12V","G12V","rare variants")))
boxplot(MUTRATE~tmp,outline=TRUE)
#axis(1,at=c(1,2,3,4,5),labels=c(paste(names(table(tmp)),"(n=",table(tmp),")")))
a <- wilcox.test(MUTRATE[which(tmp=="G12D")],MUTRATE[which(tmp=="WT")])
b <- wilcox.test(MUTRATE[which(tmp=="rare variants")],MUTRATE[which(tmp=="WT")])
c <- wilcox.test(MUTRATE[which(tmp=="G13D")],MUTRATE[which(tmp=="WT")])
d <- wilcox.test(MUTRATE[which(tmp=="G12V")],MUTRATE[which(tmp=="WT")])
e <- wilcox.test(MUTRATE[which(tmp=="rare variants")],MUTRATE[which(tmp=="WT")])
title(main="None of the common KRAS specific mutation 
is associated with higher overall mutation rate")
title(sub=paste(paste("G12D vs. WT","p =",format(a$p.value,digits=2)),
                paste("G12V vs. WT","p =",format(d$p.value,digits=2)), 
                paste("G13D vs. WT","p =",format(c$p.value,digits=2)),
                #paste("rare vs. WT","p =",format(b$p.value,digits=2)),
                sep= "     "),
      cex.sub=.75)

rm(a,b,c,d,e,blah,tmp)

#############################################################################################################################
### correlate the KRAS mutation with the hypermutated status in CRC
#############################################################################################################################

KRAS <- ifelse(KRAS_CRC=="WT",0,1)
table(KRAS,hyperMutCRC)
fisher.test(KRAS,hyperMutCRC)

table(KRAS_CRC,hyperMutCRC)



#############################################################################################################################
### compute the mutations that are associated with G12D specifically in CRC
#############################################################################################################################

G12D <- MATMUT_CRC["KRAS_p.G12D",]
GREATER <- apply(MATMUT_CRC,1,function(x){ fisher.test(as.numeric(G12D),as.numeric(x),alternative="greater")$p.value})
LESS <- apply(MATMUT_CRC,1,function(x){ fisher.test(as.numeric(G12D),as.numeric(x),alternative="less")$p.value})

G12D_CRC_OVERLAP_PVAL <- GREATER
names(G12D_CRC_OVERLAP_PVAL) <-  sapply(strsplit(x=names(GREATER),split="_"),function(x){x[[1]]})
setwd("/gluster/home/cferte/FELLOW/cferte/KRAS_Project/OBJECTS/MUT_OVERLAP/CRC_MUT_OVERLAP/")
save(G12D_CRC_OVERLAP_PVAL,file="G12D_CRC_OVERLAP_PVAL.RData")


G_G12D <- names(GREATER)[which(GREATER<.05)]
L_G12D <- names(LESS)[which(LESS<.05)]

G_G12D <- sapply(strsplit(G_G12D,"_"),function(x){x[[1]]})
paste(G_G12D,collapse=" ")


#############################################################################################################################
### compute the mutations that are associated with G13D specifically in CRC
#############################################################################################################################
G13D <- MATMUT_CRC["KRAS_p.G13D",]

GREATER <- apply(MATMUT_CRC,1,function(x){ fisher.test(as.numeric(G13D),as.numeric(x),alternative="greater")$p.value})
LESS <- apply(MATMUT_CRC,1,function(x){ fisher.test(as.numeric(G13D),as.numeric(x),alternative="less")$p.value})

G13D_CRC_OVERLAP_PVAL <- GREATER
names(G13D_CRC_OVERLAP_PVAL) <-  sapply(strsplit(x=names(GREATER),split="_"),function(x){x[[1]]})
setwd("/gluster/home/cferte/FELLOW/cferte/KRAS_Project/OBJECTS/MUT_OVERLAP/")
save(G13D_CRC_OVERLAP_PVAL,file="G13D_CRC_OVERLAP_PVAL.RData")

G_G13D <- names(GREATER)[which(GREATER<.05)]
L_G13D <- names(LESS)[which(LESS<.05)]

G_G13D <- sapply(strsplit(G_G13D,"_"),function(x){x[[1]]})
paste(G_G13D,collapse=" ")

#############################################################################################################################
### compute the mutations that are associated with G12V specifically in CRC
#############################################################################################################################
G12V <- MATMUT_CRC["KRAS_p.G12V",]
GREATER <- apply(MATMUT_CRC,1,function(x){ fisher.test(as.numeric(G12V),as.numeric(x),alternative="greater")$p.value})
LESS <- apply(MATMUT_CRC,1,function(x){ fisher.test(as.numeric(G12V),as.numeric(x),alternative="less")$p.value})

G12V_CRC_OVERLAP_PVAL <- GREATER
names(G12V_CRC_OVERLAP_PVAL) <-  sapply(strsplit(x=names(GREATER),split="_"),function(x){x[[1]]})
setwd("/gluster/home/cferte/FELLOW/cferte/KRAS_Project/OBJECTS/MUT_OVERLAP/")
save(G12V_CRC_OVERLAP_PVAL,file="G12V_CRC_OVERLAP_PVAL.RData")

G_G12V <- names(GREATER)[which(GREATER<.05)]
L_G12V <- names(LESS)[which(LESS<.05)]

G_G12V <- sapply(strsplit(G_G12V,"_"),function(x){x[[1]]})
paste(G_G12V,collapse=" ")


###################################################################################################
# sort out the first pathways associated with G12D  specifically in CRC
###################################################################################################

#load the data
load("/home/cferte/FELLOW/cferte/KRAS_Project/OBJECTS/MUT_OVERLAP/CRC_MUT_OVERLAP/G12D_CRC_OVERLAP_GSEA.RData")
G12D_CRC_GSEA <- analyticResult
G12D_PTW <- k
rm(analyticResult,k)


# load the additional subsequent data (10^6 permut) from synapse and replace them into the list
G12D_MILLION_PERMUT <- loadEntity('syn1359990')
ADD <- G12D_MILLION_PERMUT$objects$result$k
G12D_PTW[match(names(ADD),names(G12D_PTW))]  <- ADD
hist(G12D_PTW)

#  p adjust (BH) on the pathways
S <- p.adjust(G12D_PTW,method="BH")
hist(S, main="adjusted BH")
sort(S)[1:10]

G12D_MUT_CRC_PROFILE <- as.data.frame(sort(S)[1:10])
colnames(G12D_MUT_CRC_PROFILE) <- "p value"

setwd("/gluster/home/cferte/FELLOW/cferte/KRAS_Project/OBJECTS/MUT_OVERLAP/CRC_MUT_OVERLAP/")
write.table(G12D_MUT_CRC_PROFILE,file="G12D_MUT_CRC_PROFILE.txt")

###################################################################################################
# sort out the first pathways associated with G12V  specifically in CRC
###################################################################################################

# load the data from synapse
load("/home/cferte/FELLOW/cferte/KRAS_Project/OBJECTS/MUT_OVERLAP/CRC_MUT_OVERLAP/G12V_CRC_OVERLAP_GSEA.RData")
G12V_CRC_GSEA <- analyticResult
G12V_PTW <- k
rm(analyticResult,k)


# load the subsequent data from synapse (10^6 permutations)
G12V_MILLION_PERMUT <- loadEntity('syn1391161')
ADD <- G12V_MILLION_PERMUT$objects$result$k

# input inside G12V_PTW the data from the subsequent analysis and assign e-10 to the 6 very significant probes
G12V_PTW[match(names(ADD),names(G12V_PTW))]  <- ADD
G12V_PTW[G12V_PTW==0.000000000] <- 0.0000000001
hist(G12V_PTW)

# adjust for multiple testing with BH
S <- p.adjust(p=G12V_PTW,method="BH")
hist(S,main="G12V adjusted BH",breaks=20)
sort(S)[1:10]

G12V_MUT_CRC_PROFILE <- as.data.frame(sort(S)[1:10])
colnames(G12V_MUT_CRC_PROFILE) <- "p value"

setwd("/gluster/home/cferte/FELLOW/cferte/KRAS_Project/OBJECTS/MUT_OVERLAP/CRC_MUT_OVERLAP/")
write.table(G12V_MUT_CRC_PROFILE,file="G12V_MUT_CRC_PROFILE.txt")

###################################################################################################
# sort out the first pathways associated with G13D  specifically in CRC
###################################################################################################

# load the data
load("/home/cferte/FELLOW/cferte/KRAS_Project/OBJECTS/MUT_OVERLAP/CRC_MUT_OVERLAP/G13D_CRC_OVERLAP_GSEA.RData")
G13D_CRC_GSEA <- analyticResult
G13D_PTW <- k
rm(analyticResult,k)
S <- G13D_PTW[which(G13D_PTW<.05)]
sort(S)[1:10]

# adjust for multiple testing with BH
S <- p.adjust(p=G13D_PTW,method="BH")
hist(S,main="G12V adjusted BH",breaks=20)
sort(S)[1:10]

G13D_MUT_CRC_PROFILE  <- as.data.frame(sort(G13D_PTW[S])[1:10])
colnames(G13D_MUT_CRC_PROFILE) <- "p value"
setwd("/gluster/home/cferte/FELLOW/cferte/KRAS_Project/OBJECTS/MUT_OVERLAP/CRC_MUT_OVERLAP/")
write.table(G13D_MUT_CRC_PROFILE,file="G13D_MUT_CRC_PROFILE.txt")



