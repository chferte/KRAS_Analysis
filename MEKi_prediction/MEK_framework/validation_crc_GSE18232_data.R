# charles ferté
# Sage Bionetworks
# 19 juillet 2013

#####################################################################################################################
# predict the virtual drug sensitivity (vds) into the val.set
#####################################################################################################################

normalize_to_X <- function(mean.x, sd.x, Y){
  m.y <- rowMeans(Y)
  sd.y <- apply(Y, 1, sd)
  Y.adj <- (Y - m.y) * sd.x / sd.y  + mean.x 
  Y.adj[sd.y == 0] <- mean.x[sd.y==0]
  Y.adj
}

val.set2 <- normalize_to_X(rowMeans(global.matrix), apply(global.matrix, 1, sd), val.set)

# # # see what is the structure of this
# s <- svd(val.set2)
# plot(s$d^2/(sum(s$d^2)))
# plot(s$v[,1],s$v[,2],pch=19,col=c(rep(1,times=19),rep(2,times=19)))
# 
# # let's remove the outliers
# val.set2 <- val.set2[,which(s$v[,1]> -.8)]

# determine the predicted values in this val.set2
vds <- c()
for(i in c(1:N)){
  fit <- PARAL[[i]][[1]]
  vds<- cbind(vds,as.numeric(predict(fit, t(val.set2))))
}

vds <- apply(vds,1,median)
names(vds) <- colnames(val.set2)

pheno <- pheno[colnames(val.set2),]
par(mfrow=c(2,2))
boxplot(vds~pheno$treatment, main="MEK sens scrore ~ treatment")
boxplot(vds~pheno$mek.inhib, main="MEK sens scrore ~ mek treatment")
boxplot(vds~pheno$cell.lines, main="MEK sens scrore ~ cell lines")
boxplot(vds~pheno$kras.mutation, main="MEK sensitivity score ~ KRAS mutation")


