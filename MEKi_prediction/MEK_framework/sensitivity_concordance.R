# charles ferté
# feb 28
# sage bionetworks

# script to plot the concordance between the IC50 and ActArea from Sanger and CCLE for the 2 MEK inhibs.

###########################################################
# load the data
###########################################################
source("/home/cferte/FELLOW/cferte/KRAS_Analysis/MEKi_prediction/MEK_framework/load_mek_sanger.R")
source("/home/cferte/FELLOW/cferte/KRAS_Analysis/MEKi_prediction/MEK_framework/load_mek_ccle.R")

###################################################################################################
# plot the concordant cell lines between ccle and sanger
###################################################################################################
tmp <- intersect(rownames(mek.sanger.ic50),rownames(mek.ActArea))
x <- rank(apply(mek.sanger.ic50[tmp,],1,mean))
y <- rank(-1*apply(mek.ActArea[tmp,],1,mean),ties.method=)

par(mfrow=c(1,1))
plot(x,y,pch=20,xlab="sanger ranks",ylab= "ccle ranks")
plot(density(abs(x-y)[tmp]))
hist(abs(x-y)[tmp],breaks=70,col="blue")
abline(v=quantile(abs(x-y),probs=c(.25,.5,.75)),col="red")
plot(x[tmp],y[tmp],pch=20)
fit <- rlm(y~x)
summary(fit)
abline(a=fit$coefficients[1],b=fit$coefficients[2],col="orangered",lty=2,lwd=4)
cor.test(x,y)
quantile(fit$residuals,probs=c(.25,.75))
plot(fit$residuals)
hist(fit$residuals,breaks=30,xlim=c(-2,2))
plot(x,y,pch=20,col=rainbow(7)[exp((abs(fit$residuals)))])

boxplot(abs(fit$residuals))
#plot(fit$residuals)
plot(x,y,pch=20,col=rainbow(10)[exp(abs(fit$residuals))])
abline(a=fit$coefficients[1],b=fit$coefficients[2],col="gray50",lty=2,lwd=4)

tissue.color <- as.numeric(factor(gsub("^.*?_(.*)","\\1",tmp)))
plot(x,y,pch=20,col=rainbow(23)[as.numeric(tissue.color)])
     
cor.test(x,-log(y))

new.weight <- 1/exp(abs(fit$residuals))

new.weight  <- (new.weight-min(new.weight))/max(new.weight)
min(new.weight)
max(new.weight)
hist(new.weight,breaks=30)
final.weight  <- rep(0.5,times=length(mek.cells))
names(final.weight) <- names(mek.cells)
final.weight[names(new.weight)] <- new.weight
plot(density(final.weight))

# plot the concordance of these IC50
plot(density(x))
plot(density(-log(y)))
range <- seq(from=0,to=1,by=.01)

res1 <- c()  
a <- sort(x)
  b <- sort(y)  
for(i in c(1: length(a))){
  

a1 <- a
  a1 <- rep(0,times=length(a1))
  a1[(length(a1)-i):length(a1)] <- 1
  names(a1) <- names(a)
  res <- c()
  for(j in 1:(length(b)-2)){
  b1=b
  b1 <- rep(0,times=length(b1))
  b1[(length(b1)-j):length(b1)] <- 1
  names(b1) <- names(b)
  res <- c(res,-log10(fisher.test(a1,b1,alternative="greater")$p.value))
  }
res1 <- c(res1,res)
} 

  plot(density(res1))
  
  for(i in range){
    ab <-c(ab,length(intersect(a[i:72],b[j:72])))
  }
  concordant <- c(concordant,ab)
}
plot(density(concordant))
z <- rep(range,times=length(range))
w <- rep(range,each=length(range))


library(scatterplot3d)
scatterplot3d(z,w,concordant,pch=20)


which(a)