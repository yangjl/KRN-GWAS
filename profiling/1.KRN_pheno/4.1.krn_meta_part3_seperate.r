# Phenotype for KRN meta-analysis
# Jinliang YANG
# 5/14/2012

# set the working directory on server:.7
setwd("/Users/yangjl/Documents/GWAS2_KRN/pheno")


training <- read.table("krn_GenSel_training_053112.txt", header=T);
namril <- subset(training, pop.=="NAMRIL");
bxril <- subset(training, pop.=="BxRIL");
mxril <- subset(training, pop.=="MxRIL");
diallel <- subset(training, pop.=="Diallel");

write.table(namril, "krn_GenSel_training_namril.txt", sep="\t", row.names=FALSE, quote=FALSE);
write.table(bxril, "krn_GenSel_training_bxril.txt", sep="\t", row.names=FALSE, quote=FALSE);
write.table(mxril, "krn_GenSel_training_mxril.txt", sep="\t", row.names=FALSE, quote=FALSE);
write.table(diallel, "krn_GenSel_training_diallel.txt", sep="\t", row.names=FALSE, quote=FALSE);


validation <- read.table("krn_GenSel_validation_053112.txt", header=T)

namril2 <- subset(validation, pop.=="NAMRIL");
bxril2 <- subset(validation, pop.=="BxRIL");
mxril2 <- subset(validation, pop.=="MxRIL");
diallel2 <- subset(validation, pop.=="Diallel");

write.table(namril2, "krn_GenSel_validation_namril.txt", sep="\t", row.names=FALSE, quote=FALSE);
write.table(bxril2, "krn_GenSel_validation_bxril.txt", sep="\t", row.names=FALSE, quote=FALSE);
write.table(mxril2, "krn_GenSel_validation_mxril.txt", sep="\t", row.names=FALSE, quote=FALSE);
write.table(diallel2, "krn_GenSel_validation_diallel.txt", sep="\t", row.names=FALSE, quote=FALSE);

### add $ to pheno file and then run the following:
nohup ./GenSel4.1 KRN_namril.inp
nohup ./GenSel4.1 KRN_bxril.inp
nohup ./GenSel4.1 KRN_mxril.inp
nohup ./GenSel4.1 KRN_diallel.inp





