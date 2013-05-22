# charles ferté
# 20 mai 2013

# see if the logical  association of  mutations induce changes in drug sensitivity
mek.nsclc.actarea <-apply(ccle_drug[nsclc.mek.cells,mek.inhib],1,mean)
TP53 <- ccle_mut["TP53",names(mek.nsclc.actarea)]
KRAS <- ccle_mut["KRAS",names(mek.nsclc.actarea)]
STK11 <- ccle_mut["STK11",names(mek.nsclc.actarea)]
KRAS_TP53_CCLE <- ifelse(KRAS==1 & TP53 ==1,1,0)
KRAS_STK11_CCLE <- ifelse(KRAS==1 & STK11 ==1,1,0)
table(KRAS_TP53_CCLE)
table(KRAS,TP53)
cor.test(TP53,mek.nsclc.actarea,method="spearman")
cor.test(KRAS,mek.nsclc.actarea,method="spearman")
cor.test(KRAS_TP53_CCLE,mek.nsclc.actarea,method="spearman")
par(mfrow=c(2,2))
boxplot(mek.nsclc.actarea~KRAS,main="KRAS mut")
boxplot(mek.nsclc.actarea~TP53,main="TP53 mut")
boxplot(mek.nsclc.actarea~KRAS_TP53_CCLE,main="KRAS mut and TP53 mut")
boxplot(mek.nsclc.actarea~KRAS_STK11_CCLE,main="KRAS mut and STK11 mut")
wilcox.test(mek.nsclc.actarea~KRAS_TP53_CCLE)