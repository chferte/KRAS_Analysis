# charles ferté 
# feb 11th 2013
# explore the spliceosome junction data

#
library(synapseClient)
library(limma)
synapseLogin(username="charles.ferte@sagebase.org",password="charles")

# load data
junction.matrix <- loadEntity("syn1680085")
junction.matrix <- junction.matrix$objects$junction.matrix

luad_data <- loadEntity("syn1676707")
luad_data <- luad_data$objects$luad_data
names(luad_data)
luad_exp <- luad_data[[1]]
luad_mut <- luad_data[[2]]
luad_kras <- luad_data[[3]]

# make coherent colnames(junction.matrix) with luad_kras
names(luad_kras) <- substr(names(luad_kras),1,12)
colnames(junction.matrix) <- substr(colnames(junction.matrix),1,12)
tmp <- intersect(colnames(junction.matrix),names(luad_kras))
luad_kras <- luad_kras[tmp]
junction.matrix <- junction.matrix[,tmp]


# see if the G12C and G12V phenotypes are associated with different overall junctions
boxplot(apply(junction.matrix,2,sum)~luad_kras)
stripchart(apply(junction.matrix,2,sum)~luad_kras,method="jitter",vertical=TRUE,pch=20,col="aquamarine4",add=TRUE)

#
s.G12C <- svd(t(junction.matrix[,luad_kras=="G12C"]))
s.G12V <- svd(t(junction.matrix[,luad_kras=="G12V"]))
par(mfrow=c(1,2))
plot((s.G12C$d)^2/sum((s.G12C$d)^2))
plot((s.G12V$d)^2/sum((s.G12V$d)^2))


# restrict to KRAS G12C and G12V
tmp <- luad_kras == "G12C" | luad_kras=="G12V"

luad_kras <- luad_kras[tmp==TRUE]
setwd("/home/cferte/FELLOW/cferte/KRAS_Analysis/splice_analysis/")
blah <- junction.matrix
blah[blah!=0] <- 1
boxplot(apply(blah[,tmp],2,sum)~luad_kras)

blah <- junction.matrix[,tmp==TRUE]

blah <- blah[,sort(luad_kras,index.return=TRUE)$ix]


blah <- blah[-which(apply(blah,1,function(x){length(which(x!=0))})<5),]
n <- sample(x=rownames(blah),size=20000,replace=FALSE)
boxplot(blah[n,],outline=FALSE)
it(junction.matrix,model.matrix(~luad_kras)))
hist(fit$p.value[,2],breaks=100, main="junction matrix ~ kras G12C/G12V")



