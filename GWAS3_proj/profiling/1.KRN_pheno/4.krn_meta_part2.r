# Phenotype for KRN meta-analysis
# Jinliang YANG
# 5/14/2012

# set the working directory on server:.7
setwd("/Users/yangjl/Documents/GWAS2_KRN/pheno")

###---------------- DIALLEL data ----------------------###
diallel <- read.csv("/Users/yangjl/Documents/Heterosis_GWAS/pheno2011/pheno_diallel_master_BLUE.csv")
diallel <- diallel[, 1:2]
diallel$FID <- "Z000"
diallel$pop <- "Diallel"

###----------------BxRIL data 2011-------------------####
bxril <- read.table("/Users/yangjl/Documents/Heterosis_GWAS/pheno2011/pheno_BxNAM445_KRN_032212.txt", header=T)
bxril$Genotype <- paste("B73x", bxril$Genotype, sep="")
names(bxril)[4] <- "FID"

###---------------- NAM RIL--------------------------###
krnblup <- read.table("pheno_krn_meta_9reps_PLUP.txt", header=T)
idx0 <- grep("Z027", krnblup$Genotype)
krnblup <- krnblup[-idx0,]

####
krnblup$Note <- NA;
krnblup$Genotype <- as.character(krnblup$Genotype)
for(i in 1:nrow(krnblup)){
  tem <- unlist(strsplit(krnblup$Genotype[i], split="E"));
  krnblup$Note[i] <- tem[1];
  print(i)
}
krnblup$pop <- "NAMRIL"
names(krnblup)[3] <- "FID"

###---------------IBM BxRIL data 2008-----------------###
ibmkrn <- read.table("~/Documents/GWAS2_KRN/pheno/meankrn.txt", header=T)
idx1 <- grep("B73x", ibmkrn$Genotype)
bx <- ibmkrn[idx1,]
bx <- bx[bx$Genotype!="B73xMo17",]
names(bx)[2] <- "KRN"
bx$FID <- "Z017"
bx$Genotype <- gsub("B73x", "", bx$Genotype)

leaf <- read.csv("~/Documents/VariationDB/pheno/IBM_leaf_traits_20101025_KSU.csv", header=T)
leaf$MO.Number <- gsub("O", "0", leaf$MO.Number)


bx <- merge(bx, leaf[,1:2], by.x="Genotype", by.y="MO.Number", all.x=T)
tem1 <- bx[is.na(bx$Geno_Code),]
tem1$Geno_Code <- tem1$Genotype
tem1$Geno_Code <- gsub("M", "Z017EM", tem1$Geno_Code)
tem2 <- bx[!is.na(bx$Geno_Code),]
bx <- rbind(tem1, tem2)
bx$Genotype <- paste("B73x", bx$Geno_Code, sep="")

bxril <- rbind(bxril[, c(1,2,4)], bx[, c(1,2,4)])
bxril$pop <- "BxRIL"
dim(bxril)
#734


### IBM MxRIL data 2008
idx2 <- grep("Mo17x", ibmkrn$Genotype)
mx <- ibmkrn[idx2,]
mx <- mx[mx$Genotype!="Mo17xB73",]
names(mx)[2] <- "KRN"
mx$FID <- "Z017"
mx$Genotype <- gsub("Mo17x", "", mx$Genotype)


mx <- merge(mx, leaf[,1:2], by.x="Genotype", by.y="MO.Number", all.x=T)
tem3 <- mx[is.na(mx$Geno_Code),]
tem3$Geno_Code <- tem3$Genotype
tem3$Geno_Code <- gsub("M", "Z017EM", tem3$Geno_Code)
tem4 <- mx[!is.na(mx$Geno_Code),]
mx <- rbind(tem3, tem4)
mx$Genotype <- paste("Mo17x", mx$Geno_Code, sep="")

mxril <- mx[, c(1,2,4)]
mxril$pop <- "MxRIL"
dim(mxril)
#290 4

ibmril <- ibmkrn[c(-idx1, -idx2),]
ibmril$Genotype <- as.character(ibmril$Genotype)
ibmril <- ibmril[ibmril$Genotype!="B73" & ibmril$Genotype!="Mo17", ]

ibmril <- merge(ibmril, leaf[,1:2], by.x="Genotype", by.y="MO.Number", all.x=T)
tem5 <- ibmril[is.na(ibmril$Geno_Code),]
tem5$Geno_Code <- tem5$Genotype
tem5$Geno_Code <- gsub("M", "Z017EM", tem5$Geno_Code)
ibmril <- tem5
ibmril$Genotype <- tem5$Geno_Code
ibmril <- ibmril[,c(1,2)]
names(ibmril)[2] <- "KRN"
ibmril$FID <- "Z017"
ibmril$pop <- "NAMRIL"
dim(ibmril)
# 132 4


#######-------------------------------------------------------------#########


##############################################################################################
krnall <- rbind(krnblup,ibmril, bxril, mxril, diallel) #6745
krnall$Fix <- 1;
krnall <- krnall[, c(1,2,5,4,3)]

write.table(krnall, "krn_GenSel_6745_061412.txt", sep="\t", quote=FALSE, row.names=FALSE)

set.seed(1234)
sub1 <- subset(krnall, pop=="NAMRIL")
idx1 <- sample(1:nrow(sub1), size=as.integer(0.8*nrow(sub1)))
training1 <- sub1[idx1, ]
validate1 <- sub1[-idx1, ]

sub2 <- subset(krnall, pop=="BxRIL")
idx2 <- sample(1:nrow(sub2), size=as.integer(0.8*nrow(sub2)))
training2 <- sub2[idx2, ]
validate2 <- sub2[-idx2, ]

sub3 <- subset(krnall, pop=="Diallel")
idx3 <- sample(1:nrow(sub3), size=as.integer(0.8*nrow(sub3)))
training3 <- sub3[idx3, ]
validate3 <- sub3[-idx3, ]

sub4 <- subset(krnall, pop=="MxRIL")
idx4 <- sample(1:nrow(sub4), size=as.integer(0.8*nrow(sub4)))
training4 <- sub4[idx4, ]
validate4 <- sub4[-idx4, ]


training <- rbind(training1, training2, training3, training4)
validation <- rbind(validate1, validate2, validate3, validate4)
write.table(training, "krn_GenSel_training_053112.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(validation, "krn_GenSel_validation_053112.txt", sep="\t", quote=FALSE, row.names=FALSE)

#################################################
plot(density(subset(krnall, pop=="NAMRIL")$KRN ), col="red", lwd=2, ylim=c(0,.5), main="Density Plot of KRN", xlab="KRN")
lines(density(subset(krnall, pop=="Diallel")$KRN), col="bisque4", lwd=2)
lines(density(subset(krnall, pop=="BxRIL")$KRN), col="cadetblue4", lwd=2)
lines(density(subset(krnall, pop=="MxRIL")$KRN), col="chocolate4", lwd=2)
abline(v=17.1, col="blue", lty=2, lwd=2)
temp <- legend("topright", legend =c(" ", " ", " ", " "), title="Population", text.width=strwidth("NAMRIL"),
               lty=1, lwd=3, col=c("red", "bisque4", "cadetblue4", "chocolate4"), xjust=1, yjust=1)
text(temp$rect$left + temp$rect$w, temp$text$y,
       c("NAMRIL", "Diallel", "BxRIL", "MxRIL"), pos=2)



######---FOR PLINK use
pheno <- read.table("krn_GenSel_6745_061412.txt", header=T)
names(pheno)[5]<- "Note";
pheno$Genotype <- as.character(pheno$Genotype);
pheno$FID <- pheno$Genotype;
pheno$IID <- 1;

for(i in 1:nrow(pheno)){
	if(pheno$pop[i] != "Diallel"){
		tem <- unlist(strsplit(pheno$Genotype[i], split="E"));
		pheno$FID[i] <- tem[1];
		pheno$IID[i] <- tem[2];
	}
	print(i);
}

krn <- pheno[,c("FID", "IID", "KRN")];
covar <- pheno[,c("FID", "IID", "pop", "Note")]

write.table(krn, "pheno_KRN_PLINK_meta.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(covar, "covar_KRN_PLINK_meta.txt", sep="\t", row.names=FALSE, quote=FALSE);

covar <- read.table("covar_KRN_PLINK_meta.txt", header=F);

covar2 <- matrix(data=0, nrow=6745, ncol=34)
covar2 <- as.data.frame(covar2)
names(covar2) <- c("FID", "IID","Cluster", "NAMRIL", "B73x", "Mo17x", "Diallel",  "B73", as.character(unique(covar$V4))[1:26])

covar2$FID <- covar$V1
covar2$IID <- covar$V2
covar2$Cluster <- covar$V3

covar2$FID <- as.character(covar2$FID);
for(i in 1:nrow(covar2)){
	if(covar2$Cluster[i] == "NAMRIL"){
		tem <- covar2$FID[i];
		covar2[i, tem] <- 1;
		covar2[i, "NAMRIL"] <- 1;
		
	}else if(covar2$Cluster[i] == "BxRIL"){
		tem <- gsub("B73x", "", covar2$FID[i]);
		covar2[i, tem] <- 1;
		covar2[i, "B73x"] <- 1;
	}else if(covar2$Cluster[i] == "MxRIL"){
		tem <- gsub("Mo17x", "",covar2$FID[i]);
		covar2[i, tem] <- 1;
		covar2[i, "B73x"] <- 1; 
	}else if(covar2$Cluster[i] == "Diallel"){
		tem <- unlist(strsplit(covar2$FID[i], split="x"));
		covar2[i, tem[1]] <- 1;
		covar2[i, tem[2]] <- 1;
		covar2[i, "Diallel"] <- 1; 
	}
	print(i);
}

write.table(covar2[,c(-3,-4,-8,-9)], "covar2_KRN_PLINK_meta.txt", sep="\t", row.names=FALSE, quote=FALSE)

covar <- covar[, c(1,2,4)];
ch <- as.character(unique(covar$V4));
covar$V5 <- -9;
for(i in 1:26){
	covar[covar$V4==ch[i],]$V5 = i; 
}
covar3 <- covar[, c(1,2,4)]
write.table(covar3, "covar3_KRN_PLINK_meta.txt", sep="\t", row.names=FALSE, quote=FALSE, col.names=FALSE)
