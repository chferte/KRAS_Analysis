# plot the concordant cell lines


range <- seq(from=0,to=1,by=.01)
range
concordant <- c()
for(j in range){
  


  # plot the correlation between mek inhibitors in ccle
  x <- ccle_drug[mek.cells,mek.inhib[1]]
  y <- ccle_drug[mek.cells,mek.inhib[2]]
  cor(x,y)
  
  
  abline(a=0,b=1)
  a <- sort(ccle_drug[mek.cells,mek.inhib[1]])
  b <- sort(ccle_drug[mek.cells,mek.inhib[2]])
  i=100
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
  plot(res,type="l")
  
  for(i in range){
    ab <-c(ab,length(intersect(a[i:72],b[j:72])))
  }
  concordant <- c(concordant,ab)
}
plot(density(concordant))
x <- rep(range,times=length(range))
y <- rep(range,each=length(range))


library(scatterplot3d)
scatterplot3d(x,y,concordant,pch=20)


which(a)