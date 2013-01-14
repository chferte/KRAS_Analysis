
# compute the difference between the features selected in G12C and G12V models
load("/home/cferte/FELLOW/cferte/KRAS_Analysis/biological_info_meaning/G12C_GE_features.rda")
load("/home/cferte/FELLOW/cferte/KRAS_Analysis/biological_info_meaning/G12V_GE_features.rda")

G12C <- unlist(strsplit(x=G12C_GE_features,split=" "))
G12V <- unlist(strsplit(x=G12V_GE_features,split=" "))

inter <- intersect(G12C,G12V)

fisher.test(x=c(48,325),y=c(48,217))

# specific G12C model fetaures
paste(setdiff(G12C,inter),collapse=" ")

# specific G12V model fetaures
paste(setdiff(G12V,inter),collapse=" ")

# common G12V and G12C model fetaures
paste(inter,collapse=" ")