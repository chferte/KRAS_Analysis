# charles ferté
# Sage Bionetworks
# 07 mai 2013

#####################################################################################################################
# predict the virtual drug sensitivity (vds) into the val_exp
#####################################################################################################################

load("/home/cferte/Bootstrp_300_gene_exp_mek_models.rda")
N <- length(PARAL)
val.set <- crc[[1]][rownames(global.matrix),]

normalize_to_X <- function(mean.x, sd.x, Y){
  m.y <- rowMeans(Y)
  sd.y <- apply(Y, 1, sd)
  Y.adj <- (Y - m.y) * sd.x / sd.y  + mean.x 
  Y.adj[sd.y == 0] <- mean.x[sd.y==0]
  Y.adj
}

val.set2 <- normalize_to_X(rowMeans(global.matrix), apply(global.matrix, 1, sd), val.set)

seq <- c(1:N)
#seq <- as.numeric(gsub(pattern="m",replacement="",hemal.models))

vds <- c()
for(i in seq){
  fit <- PARAL[[i]][[1]]
  vds<- cbind(vds,as.numeric(predict(fit, t(val.set2))))
}

rownames(vds) <- substr(colnames(val.set2),1,12)
vds.mean <- apply(vds,1,median)
order <- names(sort(vds.mean))
boxplot(t(vds[order,]), main="VDS",col="orangered")
abline(h=quantile(vds,probs=c(.25,.75)),col="red",lty=2,lwd=2)
save(vds,file="/home/cferte/RESULTS/vds_crc.Rda")
