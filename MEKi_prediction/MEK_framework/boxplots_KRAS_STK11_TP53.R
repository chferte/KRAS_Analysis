foo <- ccle_drug_ActAreaNorm[cell.names.idx,mek.inhib]
foo <- apply(foo,1,mean)
par(mfrow=c(2,3))
boxplot(foo~ccle_oncomap["STK11",cell.names.idx],main="STK11 mut",ylab="MEK sensitivity (Act. Area)")
boxplot(foo~ccle_oncomap["TP53",cell.names.idx], main="TP53 mut",ylab="MEK sensitivity (Act. Area)")
boxplot(foo~ccle_oncomap["KRAS",cell.names.idx], main="KRAS mut",ylab="MEK sensitivity (Act. Area)")
KRAS_STK11 <- ifelse(ccle_oncomap["KRAS",]==1 & ccle_oncomap["STK11",]==1,1,0)
KRAS_TP53 <- ifelse(ccle_oncomap["KRAS",]==1 & ccle_oncomap["TP53",]==1,1,0)
gold <- ifelse(ccle_oncomap["KRAS",]==1,"KRAS mut alone","wt")
gold[names(which(KRAS_STK11==1))] <- "KRAS mut & STK11 mut"
gold[names(which(KRAS_TP53==1))] <- "KRAS mut & TP53 mut"
table(gold)
boxplot(foo~gold[cell.names.idx])
par(mfrow=c(1,1))
boxplot(foo~KRAS_STK11[cell.names.idx],main=" KRAS mut & STK11 mut",ylab="MEK sensitivity (Act. Area)")
boxplot(foo~KRAS_TP53[cell.names.idx],main="KRAS mut & TP53 mut",ylab="MEK sensitivity (Act. Area)")

