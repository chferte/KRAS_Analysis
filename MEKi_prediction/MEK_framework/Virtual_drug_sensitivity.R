# charles ferté
# Sage Bionetworks
# 07 mai 2013

#####################################################################################################################
# predict the virtual drug sensitivity (vds) into the val_exp
#####################################################################################################################

vds <- c()
for(i in c(1:N)){
  fit <- PARAL[[i]][[1]]
  vds<- cbind(vds,as.numeric(predict(fit, t(val_exp))))
  print(i)
}

vds <- apply(vds,1,mean)
names(vds) <- colnames(val_exp)


#####################################################################################################################
# generate bootstrapped sparse (lasso) models using the gistic and the mutations calls of the same data that predict the vds
#####################################################################################################################

# load the gistic calls
syn1687610

