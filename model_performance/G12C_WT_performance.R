# Charles Ferté
# Sage Bionetworks
# Dec 14th 2012


##############################################################
# assess the performance of the G12C_WT model in terms of AUC
###############################################################

#######################################################
#evaluate the performance of the model in BATTLE
#######################################################
tmp <- ifelse(KRAS_BATTLE=="G12C",1,0)
names(tmp) <- names(KRAS_BATTLE)
G12C_BATTLE <- tmp
rm(tmp)


########################################################
# plot a ROC curve to asses the performance of our model in BATTLE
########################################################
AUC_BATTLE <- c()
#par(mfrow=c(2,2))
require(ROCR)
for (i in c(1:N)){
  Pred <- prediction(as.numeric(yhat_BATTLE[,i]),as.numeric(G12C_BATTLE))
  Perf <- performance(prediction.obj=Pred,"tpr","fpr")
  AUC <- performance(prediction.obj=Pred,"auc")
  AUC_BATTLE <- c(AUC_BATTLE,as.numeric(AUC@y.values))
}

#plot(Perf, col="royalblue",main="predicting KRAS G12C in BATTLE 
#(modele trained in TCGA)")
#text(x=.7,y=.4,labels=paste("AUC=",format(x=AUC@y.values,digits=2)),col="royalblue")

################################################################
# apply the model in CHEMORES
################################################################



#######################################################
#evaluate the performance of the model in CHEMORES
#######################################################
tmp <- ifelse(KRAS_CHEMORES=="G12C",1,0)
names(tmp) <- names(KRAS_CHEMORES)
G12C_CHEMORES <- tmp
rm(tmp)

#boxplot(yhat_CHEMORES~G12C_CHEMORES,xlab=c("KRAS G12C mutational status"),ylab="model of G12C",main="predicting KRAS G12C in CHEMORES
#(modele trained in TCGA)")
#stripchart(yhat_CHEMORES~G12C_CHEMORES,vertical=TRUE,pch=20,add=T,col="royalblue",cex=.7 )

########################################################
# plot a ROC curve to asses the performance of our model in CHEMORES
########################################################

AUC_CHEMORES <- c()
for (i in c(1:N)){
  Pred <- prediction(as.numeric(yhat_CHEMORES[,i]),as.numeric(G12C_CHEMORES))
  Perf <- performance(prediction.obj=Pred,"tpr","fpr")
  AUC <- performance(prediction.obj=Pred,"auc")
  AUC_CHEMORES <- c(AUC_CHEMORES,as.numeric(AUC@y.values))
}
boxplot(AUC_CHEMORES)


#######################################################
#evaluate the performance of the model in ccle
#######################################################
tmp <- ifelse(KRAS_CCLE=="G12C",1,0)
names(tmp) <- names(KRAS_CCLE)
G12C_CCLE <- tmp
rm(tmp)


#boxplot(yhat_CCLE~G12C_CCLE,xlab=c("KRAS G12C mutational status"),ylab="model of G12C",main="predicting KRAS G12C in ccle
#(modele trained in TCGA)")
#stripchart(yhat_CCLE~G12C_CCLE,vertical=TRUE,pch=20,add=T,col="royalblue",cex=.7,method="jitter")

##################################################################
# plot a ROC curve to asses the performance of our model in ccle
##################################################################

AUC_CCLE <- c()
for (i in c(1:N)){
  Pred <- prediction(as.numeric(yhat_CCLE[,i]),as.numeric(G12C_CCLE))
  Perf <- performance(prediction.obj=Pred,"tpr","fpr")
  AUC <- performance(prediction.obj=Pred,"auc")
  AUC_CCLE <- c(AUC_CCLE,as.numeric(AUC@y.values))
}
boxplot(AUC_CCLE)
# 
# plot(Perf, col="royalblue",main="predicting KRAS G12C in ccle
# (modele trained in TCGA)")
# text(x=.7,y=.4,labels=paste("AUC=",format(x=AUC@y.values,digits=2)),col="royalblue")


#######################################################
#evaluate the performance of the model in Sanger
#######################################################
tmp <- ifelse(KRAS_SANGER=="G12C",1,0)
names(tmp) <- names(KRAS_SANGER)
G12C_SANGER <- tmp
rm(tmp)


#boxplot(yhat_SANGER~G12C_SANGER,xlab=c("KRAS G12C mutational status"),ylab="model of G12C",main="predicting KRAS G12C in Sanger
#(modele trained in TCGA)")
#stripchart(yhat_SANGER~G12C_SANGER,vertical=TRUE,pch=20,add=T,col="royalblue",cex=.7, method="jitter")

##################################################################
# plot a ROC curve to asses the performance of our model in SANGER
##################################################################

AUC_SANGER <- c()
for (i in c(1:N)){
  Pred <- prediction(as.numeric(yhat_SANGER[,i]),as.numeric(G12C_SANGER))
  Perf <- performance(prediction.obj=Pred,"tpr","fpr")
  AUC <- performance(prediction.obj=Pred,"auc")
  AUC_SANGER <- c(AUC_SANGER,as.numeric(AUC@y.values))
}
boxplot(AUC_SANGER)
#
plot(Perf, col="royalblue",main="predicting KRAS G12C in Sanger
(modele trained in TCGA)")
text(x=.7,y=.4,labels=paste("AUC=",format(x=AUC@y.values,digits=2)),col="royalblue")

#######################
par(mfrow=c(1,1))
tmp <- cbind(AUC_BATTLE,AUC_CHEMORES,AUC_CCLE,AUC_SANGER)
rownames(tmp) <- c(1:dim(tmp)[1])
boxplot(tmp,ylab="AUC",outline=FALSE,ylim=c(.5,1))
stripchart(list(BATTLE=tmp[,1], CHEMORES=tmp[,2],CCLE=tmp[,3],SANGER=tmp[,4]), add=T,vertical=TRUE,method="jitter",col="red",pch=20)
abline(h=c(.5,.6,.7,.8,.9),lty=2,lwd=.7)

