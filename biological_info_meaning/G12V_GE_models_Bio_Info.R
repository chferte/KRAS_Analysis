# Charles Ferté
# Sage Bionetworks
# 14 Sept 2012


# infer the biological information from the G12V_WT gene expression model

# first load the data and 
source("/home/cferte/FELLOW/cferte/KRAS_Analysis/codon_specific_models/G12V_model_GOLD.R")

# identify the top selected features (selected > 10 times)
tmp <- ifelse(abs(selected)>0,1,0)
tmp <- apply(tmp,1,sum)
names(tmp) <- rownames(selected)

table(tmp>10)
G12V_GE_features <- paste(names(tmp)[which(tmp>10)],collapse=" ")
hist(log10(tmp),breaks=100,col="red")
save(G12V_GE_features,file="/home/cferte/FELLOW/cferte/KRAS_Analysis/biological_info_meaning/G12V_GE_features.rda")
#rank.features <-as.numeric(sort(tmp,index.return=T)$ix)
#names(rank.features) <- rownames(selected)



# ResultBS <- c()
# for(k in c(1:dim(selected)[2])){
# ResultBS <- cbind(ResultBS,rank(as.numeric(selected[,k]),ties.method="min")/(dim(selected)[1]))
# }
# rank.features <- apply(ResultBS,1,sum)
# names(rank.features) <- rownames(selected)

# assign rank.features and coeff.features for GSEA_ISJ
#rank.features <- tmp
# coeff.features <- rep(x=1,times=length(rank.features))
# #rm(tmp)
# 
# # path and name of the file where the data will be stored
# final_file <- "/home/cferte/FELLOW/cferte/KRAS_Analysis/biological_info_meaning/G12V_GE_GSEA_results.rda"
# 
# # load GSEA_ISJ
# source("/home/cferte/FELLOW/cferte/KRAS_Analysis/biological_info_meaning/GSEA_ISJ.R")
# 
# # load the results
# load("/home/cferte/FELLOW/cferte/KRAS_Analysis/biological_info_meaning/G12V_GE_GSEA_results.rda")
# table(k<.01)
# k[which(k<.05)]
# 
# significant.ptw <- p.adjust(k,method="BH")
# sort(significant.ptw)[1:100]
# 
# 
