# plot the concordant cell lines


range <- seq(from=0,to=1,by=.01)
range
concordant <- c()
for(j in range){
  


  ab <- c()
  a <- names(sort(ccle_drug[nsclc.mek.cells,mek.inhib[1]]))
  b <- names(sort(ccle_drug[nsclc.mek.cells,mek.inhib[2]]))
  
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