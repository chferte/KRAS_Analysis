# charles ferté
# sage bionetworks
# Feb 28th, 2013


#############################
# load the data
#############################
source("/home/cferte/FELLOW/cferte/KRAS_Analysis/MEKi_prediction/MEK_framework/load_mek_ccle.R")


################################################################################################################
# compute mixture models to assign sensitive and resistant cells
################################################################################################################

require(flexmix)

Nclust <- 2

cells <- list(mek.cells,nsclc.mek.cells,breast.mek.cells,crc.mek.cells,hemal.mek.cells,melanoma.mek.cells,pancreas.mek.cells,ovary.mek.cells)
cell.names <- list("ALL CELLS","NSCLC","BREAST","CRC","Hematologic\nMalignancies","MELANOMA","PANCREAS","OVARY")


#################
all.prob <- c()
all.prob.sens <- c()
all.prob.res <- c()
total.cell.status <-c()
for(i in c(1:length(cells))){
  
  blah <- apply(ccle_drug[cells[[i]],mek.inhib],1,mean)
  blah <- blah[names(sort(blah))]
  rank(blah)
  Nclust <- 2
  
  tissue.post.prob <- c()
  K <- 150
  while(length(tissue.post.prob)<K) {
    sample.out <- as.numeric(format(length(names(blah))*.07,digits=1))
    iter.cells <- sample(names(blah),replace=FALSE,size=length(names(blah))-sample.out)
    #iter.cells <- sample(rownames(blah),replace=FALSE,size=length(rownames(blah))-2)
    #tryclust <- flexmix(blah[iter.cells,]~1,k=Nclust, model=FLXMCmvnorm(diag=FALSE))
    tryclust <- flexmix(blah[iter.cells]~1,k=Nclust)
    clust <- tryclust@cluster
    table(clust)
    if(length(table(clust))<Nclust){} else{
      names(clust) <- iter.cells
      
      post.prob <- tryclust@posterior$scaled
      if(mean(blah[names(which(clust==1))]) > mean(blah[names(which(clust==2))])) {
        colnames(post.prob) <- c("sens","res")
      } else { colnames(post.prob) <- c("res","sens")}
      post.prob <- post.prob[,c("sens","res")]
      rownames(post.prob) <- iter.cells
      tissue.post.prob  <- c(tissue.post.prob,list(post.prob))
    }
  }
  
  
  final.prob.sens <- matrix(NA,ncol=length(tissue.post.prob),nrow=length(cells[[i]])) 
  rownames(final.prob.sens) <- names(blah)
  for(x in  1:length(tissue.post.prob)){final.prob.sens[rownames(tissue.post.prob[[x]]),x]<-(tissue.post.prob[[x]][,"sens"])}
  
  all.prob.sens <- c(all.prob.sens,list(final.prob.sens))
  
  final.prob.res <- matrix(NA,ncol=length(tissue.post.prob),nrow=length(cells[[i]])) 
  rownames(final.prob.res) <- names(blah)
  for(x in  1:length(tissue.post.prob)){final.prob.res[rownames(tissue.post.prob[[x]]),x]<-(tissue.post.prob[[x]][,"res"])}
  all.prob.res <- c(all.prob.res,list(final.prob.res))
  
  final.prob <- cbind(apply(final.prob.sens,1,function(x){exp(mean(log(x[!is.na(x)])))}),apply(final.prob.res,1,function(x){exp(mean(log(x[!is.na(x)])))}))
  colnames(final.prob) <- c("sens","res")
  
  all.prob <- c(all.prob,list(final.prob))
  
}


###################
# plot the ActArea density plot and the distribution of status (sens, res, inter)
##################
par(mfrow=c(2,4),oma=c(0,0,0,0))
for(i in 1:length(all.prob))
{
  #col.vec <- (apply(all.prob[[i]],1,function(x){names(which(x==max(x)))}))
  #col.vec <- ifelse(col.vec=="sens",1,2)
  col.vec <- factor(-1*(all.prob[[i]][,1]+1))
  density.object <- density(apply(ccle_drug[cells[[i]],mek.inhib],1,mean))
  plot(density.object,lwd=3,
       main=paste(cell.names[[i]]),
       xlim=c(-1.5,8),ylim=c(0,1),
       xlab=paste("Activity Area (Normalized data)","\n(N= ",density.object$n,", Bandwidth= ",format(density.object$bw,digits=2),")",sep=""))
  ypos <- abs(rnorm(sd=.04,mean=.8,n=length(rownames(all.prob[[i]]))))
  points(apply(ccle_drug[rownames(all.prob[[i]]),mek.inhib],1,mean),ypos,col=greenred(length(col.vec))[col.vec],pch=19,cex=1)
}


# cell status: sensitive are the ones >.66, res are <.33, inter are the rest !
cell.status <- sapply(1:length(cell.names),function(x){ifelse(all.prob[[x]][,1]>.66,"sens",ifelse(all.prob[[x]][,1]<.33,"res","inter"))})
names(cell.status) <- cell.names
names(all.prob) <- cell.names
ccle_probs_status <- list(Posterior.probs=all.prob,cell.status=cell.status)

# tabulate the distribution
sapply(names(cell.status),function(x){table(cell.status[[x]])})

#save the cell status in synapse

#ccle_sens_status <- Data(list(name = "ccle_sens_status", parentId = 'syn1670945'))
#ccle_sens_status <- createEntity(ccle_sens_status)

ccle_sens_status <- loadEntity("syn1709732")
ccle_sens_status <- deleteObject("ccle_probs_status",owner=ccle_sens_status)
ccle_sens_status <- addObject(ccle_sens_status,ccle_probs_status)
ccle_sens_status <- storeEntity(entity=ccle_sens_status)




