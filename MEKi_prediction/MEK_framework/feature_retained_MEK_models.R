# charles ferté
# Feb 22th 2013
# sage bionetworks

################################################################################################################
# explore what is retained in the models
################################################################################################################

abc <- matrix(NA,ncol=length(selected),nrow=nrow(global.matrix))
rownames(abc) <- rownames(global.matrix)
for(i in c(1:length(selected)))
{  abc[,i]<- as.numeric(selected[[i]])  
}

fit.betas.mek.models <- selected

# # store it into synapse
# fits <- Data(list(name = "fitbetas_mek_models", parentId = 'syn1670945'))
# fits<- createEntity(fits)
# 
# # add object into the data entity
# fits <- addObject(fits,fit.betas.mek.models)
# 
# # push the raw data into this entity
# fits <- storeEntity(entity=fits)

par(mfrow=c(1,1))
abc[abc!=0] <-1 
hist(log10(rowSums(abs(abc))),breaks=50,col="red")
h <- rowSums(abs(abc))
sort(h,decreasing=TRUE)[1:100]
which(h>quantile(h,probs=.99))
