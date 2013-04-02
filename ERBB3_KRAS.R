# assess wether ERBB3 is univariately associated with MEK sensitivity in the entire ccle dataset

source("/home/cferte/FELLOW/cferte/KRAS_Analysis/MEKi_prediction/MEK_framework/load_mek_ccle.R")
selected <- c(nsclc.mek.cells,pancreas.mek.cells,crc.mek.cells)
mek.response <- apply(ccle_drug_ActAreaNorm[selected,mek.inhib],1,mean)

tmp1 <- grep(pattern="ERBB3",x=rownames(ccle_exp))
tmp2 <- grep(pattern="ERBB3",x=rownames(ccle_cnv))

erbb3_median <- ifelse(ccle_exp[tmp1,selected]>median(ccle_exp[tmp1,selected]),"over-expressed","under-expressed")
boxplot(mek.response ~erbb3_median)
stripchart(list("high erbb3 expression"=mek.response[erbb3_median="over-expressed"],"low erbb3 expression"=mek.response[erbb3_median="under-expressed"]),vertical=TRUE,add=TRUE,method="jitter",col="red",pch=20)
wilcox.test(mek.response ~erbb3_median)

# plot gene expression erbb3 and MEK response
erbb3_exp <- ccle_exp[tmp1,selected]
kras.col <- ifelse(ccle_mut["KRAS",selected]==1,"green","blue")
plot(erbb3_exp,mek.response, pch=20, col=kras.col,
     main="sensitivity to MEKi ~ erbb3 gene expression \nin NSCLC, CRC & Pancreatic carcinoma",ylab="Activity Area")
fit <- loess(mek.response~erbb3_exp)
lines(fit$x,fit$fitted,col="red",lwd=3,type="p",cex=.3)
cor.test(erbb3_exp,mek.response, method="spearman")

#
erbb3_highexp <- ifelse(ccle_exp[tmp1,selected]>0,"high erbb3 expression","low erbb3 expression")
kras.status <- as.factor(ifelse(ccle_mut["KRAS",selected]==1,"kras mutant","kras wt"))

par(mfrow=c(1,2))
boxplot(mek.response[erbb3_highexp=="low erbb3 expression"]~kras.status[erbb3_highexp=="low erbb3 expression"],ylim=c(0,5),main="erbb3 low expression", ylab="Activity Area")
boxplot(mek.response[erbb3_highexp=="high erbb3 expression"]~kras.status[erbb3_highexp=="high erbb3 expression"],ylim=c(0,5), main="erbb3 high expression", ylab="Activity Area")
wilcox.test(mek.response[erbb3_highexp=="low erbb3 expression"]~kras.status[erbb3_highexp=="low erbb3 expression"])
wilcox.test(mek.response[erbb3_highexp=="high erbb3 expression"]~kras.status[erbb3_highexp=="high erbb3 expression"])
# plot cnv erbb3 and MEK response
erbb3_cnv <- ccle_cnv[tmp2,]
plot(erbb3_cnv,mek.response, pch=20, main="sensitivity to MEKi ~ erbb3 cnv",ylab="Activity Area")
fit <- loess(mek.response~erbb3_cnv)
lines(fit$x,fit$fitted,col="red",lwd=3,type="p")
cor.test(erbb3_cnv,mek.response, method="spearman")

# boxplot the MEK response according to erbb3 mutation status
ERBB3_mut <- ifelse(ccle_mut["ERBB3",mek.cells]==1,"erbb3.mutant","erbb3.wt")
boxplot(mek.response~ERBB3_mut, main="MEKi response according to erbb3 mutation status",outline=FALSE)
stripchart(list("erbb3.mutant"=mek.response[ERBB3_mut=="erbb3.mutant"],"erbb3.wt"=mek.response[ERBB3_mut=="erbb3.wt"]),method="jitter",add=TRUE,vertical=TRUE,col="red",pch=20)
wilcox.test(mek.response~ERBB3_mut)

