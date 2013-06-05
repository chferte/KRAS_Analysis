# charles ferté
# sage bionetworks
# June 3rd 2013

# load the mek data from ccle
source("/home/cferte/FELLOW/cferte/KRAS_Analysis/MEKi_prediction/MEK_framework/load_mek_ccle.R")
source("/home/cferte/FELLOW/cferte/KRAS_Analysis/MEKi_prediction/MEK_framework/load_ccle_mut_oncomap.R")
ccle_oncomap <- ifelse(ccle_oncomap=="0",0,1)
tmp <- intersect(colnames(ccle_oncomap),nsclc.mek.cells)

# boxplots MEK sensitivity ~ p53 stk11 & kras
foo <- ccle_drug_ActAreaNorm[tmp,mek.inhib]
foo <- apply(foo,1,mean)

par(mfrow=c(2,3), oma=c(0,0,0,0), mar=c(3,3,3,3))
for (i in c("KRAS","STK11","TP53")){
boxplot(foo~ccle_oncomap[i,tmp],main=paste(i,"mutation"),ylab="MEK sensitivity",outline=FALSE,ylim=c(0,5))
stripchart(list("wt"=foo[ccle_oncomap[i,tmp]==0],"mut"=foo[ccle_oncomap[i,tmp]==1]),method="jitter",vertical=TRUE,add=TRUE,pch=19,col="red")
print(wilcox.test(foo~ccle_oncomap[i,tmp])$p.value)
}

foo <- ccle_drug_ActAreaNorm[nsclc.mek.cells,mek.inhib]
foo <- apply(foo,1,mean)

for (i in c("KRAS","STK11","TP53")){
  boxplot(foo~ccle_mut[i,nsclc.mek.cells],main=paste(i,"mutation"),ylab="MEK sensitivity",outline=FALSE,ylim=c(0,5))
  stripchart(list("wt"=foo[ccle_mut[i,nsclc.mek.cells]==0],"mut"=foo[ccle_mut[i,nsclc.mek.cells]==1]),method="jitter",vertical=TRUE,add=TRUE,pch=19,col="red")
  print(wilcox.test(foo~ccle_mut[i,nsclc.mek.cells])$p.value)
}

# plot the concommitant mutantions
par(mfrow=c(2,1), oma=c(0,0,0,0))

foo <- ccle_drug_ActAreaNorm[tmp,mek.inhib]
foo <- apply(foo,1,mean)
gold <- ccle_oncomap["KRAS",tmp]
gold[gold==1] <- "KRAS only"
gold[gold==0] <- "WT"
table(gold)
gold[names(which(ccle_oncomap["KRAS",tmp]==1 & ccle_oncomap["TP53",tmp]==1))] <- "KRAS TP53"
gold[names(which(ccle_oncomap["KRAS",tmp]==1 & ccle_oncomap["STK11",tmp]==1))] <- "KRAS STK11"
table(gold)
boxplot(foo~gold,cex.lab=.8)

foo <- ccle_drug_ActAreaNorm[nsclc.mek.cells,mek.inhib]
foo <- apply(foo,1,mean)
gold <- ccle_mut["KRAS",nsclc.mek.cells]
gold[gold==1] <- "KRAS only"
gold[gold==0] <- "WT"
table(gold)
gold[names(which(ccle_mut["KRAS",nsclc.mek.cells]==1 & ccle_mut["TP53",nsclc.mek.cells]==1))] <- "KRAS TP53"
gold[names(which(ccle_mut["KRAS",nsclc.mek.cells]==1 & ccle_mut["STK11",nsclc.mek.cells]==1))] <- "KRAS STK11"
table(gold)
boxplot(foo~gold,cex.lab=.8)


par(mfrow=c(1,1))
boxplot(foo~KRAS_STK11[tmp],main=" KRAS mut & STK11 mut",ylab="MEK sensitivity (Act. Area)")
boxplot(foo~KRAS_TP53[tmp],main="KRAS mut & TP53 mut",ylab="MEK sensitivity (Act. Area)")

