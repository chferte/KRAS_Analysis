# Charles Ferté
# Sage Bionetworks
# Dec 14th 2012


##############################################################
# assess the performance of the G12CV_WT model in terms of AUC
###############################################################


#######################################################
#evaluate the performance of the model in BATTLE
#######################################################
tmp <- ifelse(KRAS_BATTLE %in% c("G12C","G12V"),1,0)
names(tmp) <- names(KRAS_BATTLE)
G12CV_BATTLE <- tmp
rm(tmp)


#par(mfrow=c(1,1))
#boxplot(yhat_BATTLE~G12CV_BATTLE,xlab=c("KRAS G12CV mutational status"),ylab="model of G12CV",main="predicting KRAS G12CV in BATTLE
#(modele trained in TCGA+CHEMORES)")
#stripchart(yhat_BATTLE~G12CV_BATTLE,vertical=TRUE,pch=20,add=T,col="royalblue",cex=.7,method="jitter")

########################################################
# plot a ROC curve to asses the performance of our model in BATTLE
########################################################
AUC_BATTLE <- c()
#par(mfrow=c(2,2))
require(ROCR)
for (i in c(1:N)){
  Pred <- prediction(as.numeric(yhat_BATTLE[,i]),as.numeric(G12CV_BATTLE))
  Perf <- performance(prediction.obj=Pred,"tpr","fpr")
  AUC <- performance(prediction.obj=Pred,"auc")
  AUC_BATTLE <- c(AUC_BATTLE,as.numeric(AUC@y.values))
}

#plot(Perf, col="royalblue",main="predicting KRAS G12CV in BATTLE 
#(modele trained in TCGA)")
#text(x=.7,y=.4,labels=paste("AUC=",format(x=AUC@y.values,digits=2)),col="royalblue")

################################################################
# apply the model in CHEMORES
################################################################



#######################################################
#evaluate the performance of the model in CHEMORES
#######################################################
tmp <- ifelse(KRAS_CHEMORES %in% c("G12C","G12V"),1,0)
names(tmp) <- names(KRAS_CHEMORES)
G12CV_CHEMORES <- tmp
rm(tmp)

#boxplot(yhat_CHEMORES~G12CV_CHEMORES,xlab=c("KRAS G12CV mutational status"),ylab="model of G12CV",main="predicting KRAS G12CV in CHEMORES
#(modele trained in TCGA)")
#stripchart(yhat_CHEMORES~G12CV_CHEMORES,vertical=TRUE,pch=20,add=T,col="royalblue",cex=.7 )

########################################################
# plot a ROC curve to asses the performance of our model in CHEMORES
########################################################

AUC_CHEMORES <- c()
for (i in c(1:N)){
  Pred <- prediction(as.numeric(yhat_CHEMORES[,i]),as.numeric(G12CV_CHEMORES))
  Perf <- performance(prediction.obj=Pred,"tpr","fpr")
  AUC <- performance(prediction.obj=Pred,"auc")
  AUC_CHEMORES <- c(AUC_CHEMORES,as.numeric(AUC@y.values))
}
boxplot(AUC_CHEMORES)


#######################################################
#evaluate the performance of the model in ccle
#######################################################
tmp <- ifelse(KRAS_CCLE %in% c("G12C","G12V"),1,0)
names(tmp) <- names(KRAS_CCLE)
G12CV_CCLE <- tmp
rm(tmp)


#boxplot(yhat_CCLE~G12CV_CCLE,xlab=c("KRAS G12CV mutational status"),ylab="model of G12CV",main="predicting KRAS G12CV in ccle
#(modele trained in TCGA)")
#stripchart(yhat_CCLE~G12CV_CCLE,vertical=TRUE,pch=20,add=T,col="royalblue",cex=.7,method="jitter")

##################################################################
# plot a ROC curve to asses the performance of our model in ccle
##################################################################

AUC_CCLE <- c()
for (i in c(1:N)){
  Pred <- prediction(as.numeric(yhat_CCLE[,i]),as.numeric(G12CV_CCLE))
  Perf <- performance(prediction.obj=Pred,"tpr","fpr")
  AUC <- performance(prediction.obj=Pred,"auc")
  AUC_CCLE <- c(AUC_CCLE,as.numeric(AUC@y.values))
}
boxplot(AUC_CCLE)
# 
# plot(Perf, col="royalblue",main="predicting KRAS G12CV in ccle
# (modele trained in TCGA)")
# text(x=.7,y=.4,labels=paste("AUC=",format(x=AUC@y.values,digits=2)),col="royalblue")


#######################################################
#evaluate the performance of the model in Sanger
#######################################################
tmp <- ifelse(KRAS_SANGER %in% c("G12C","G12V"),1,0)
names(tmp) <- names(KRAS_SANGER)
G12CV_SANGER <- tmp
rm(tmp)


#boxplot(yhat_SANGER~G12CVV_SANGER,xlab=c("KRAS G12CV mutational status"),ylab="model of G12CVV",main="predicting KRAS G12CV in Sanger
#(modele trained in TCGA)")
#stripchart(yhat_SANGER~G12CV_SANGER,vertical=TRUE,pch=20,add=T,col="royalblue",cex=.7, method="jitter")

##################################################################
# plot a ROC curve to asses the performance of our model in SANGER
##################################################################

AUC_SANGER <- c()
for (i in c(1:N)){
  Pred <- prediction(as.numeric(yhat_SANGER[,i]),as.numeric(G12CV_SANGER))
  Perf <- performance(prediction.obj=Pred,"tpr","fpr")
  AUC <- performance(prediction.obj=Pred,"auc")
  AUC_SANGER <- c(AUC_SANGER,as.numeric(AUC@y.values))
}
boxplot(AUC_SANGER)
#
plot(Perf, col="royalblue",main="predicting KRAS G12CV in Sanger
(modele trained in TCGA)")
text(x=.7,y=.4,labels=paste("AUC=",format(x=AUC@y.values,digits=2)),col="royalblue")

#######################
tmp <- cbind(AUC_BATTLE,AUC_CHEMORES,AUC_CCLE,AUC_SANGER)
rownames(tmp) <- c(1:dim(tmp)[1])
boxplot(tmp,ylab="AUC",outline=FALSE)
stripchart(list(BATTLE=tmp[,1], CHEMORES=tmp[,2],CCLE=tmp[,3],SANGER=tmp[,4]), add=T,vertical=TRUE,method="jitter",col="red",pch=20)
abline(h=c(.4,.5,.6,.7,.8),lty=2,lwd=.7)
