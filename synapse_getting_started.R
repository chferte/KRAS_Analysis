## Charles Ferté
### july 2012
#### adaptation of the gsea r code

options(stringsAsFactors=FALSE)

# load the synapse client and login !
library(synapseClient)
synapseLogin("charles.ferte@sagebase.org","charles")

#load the data you want to store in synapse

# you have to create a project in Synapse first (example of project id = syn1337457)
#load the data in the workspace ...here it is G12C_OVERLAP_PVAL your project
# then, store the data under the Synapse id project
newLayer <- Data(list(name = "ccle_drug", parentId = "syn1337457"))
newLayer <- addObject(newLayer, ccle_drug)
newLayer <- storeEntity(newLayer)

#return the synapse id name of the new data in the folder !
propertyValue(newLayer, "id")




