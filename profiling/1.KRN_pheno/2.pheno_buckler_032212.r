# Processing Buckler's NAM data
# Jinliang Yang
# update: 03.22.2012

# set the working directory on server:.7
setwd("/Users/yangjl/Documents/GWAS2_KRN/pheno")


###############################################################################
# Step 1: read in the data
###############################################################################
pheno <- read.table("/Users/yangjl/Documents/GWAS2_KRN/pheno/traitMatrix_maize282NAM_v10-120314.txt", header=T)

### Keep only the following traits
### EAR: KRN, AKW, KC, TKW=earWeight -CobWeight
### COB: CW, CD and CL

mynm <- names(pheno)
mynm1 <- c("Trait","X20KernelWeight","CobDiameter","CobDiameter.1","CobWeight","CobWeight.1","EarDiameter","EarDiameter.1","EarLength","EarLength.1",
           "EarRowNumber","EarRowNumber.1","EarWeight","EarWeight.1","TotalKernelVolume", "TotalKernelVolume.1")
mypheno <- pheno[, mynm1]
nampheno <- mypheno[-1:-283,]
nampheno[nampheno==-999] <- NA;
nampheno$FID <- NA;
names(nampheno)[1] <- "Genotype";
nampheno$Genotype <- as.character(nampheno$Genotype);
for(i in 1:nrow(nampheno)){
  tem <- unlist(strsplit(nampheno$Genotype[i], split="E"));
  nampheno$FID[i] <- tem[1];
}
nampheno$Fix <- 1;

##### Read in the pheno2011 BxNAM and NAM data ##################
pheno2011 <- read.table("~/Documents/Heterosis_GWAS/pheno2011/pheno_BxNAM445_KRN_032212.txt", header=TRUE)

##### AKW #######################################################
akw <- nampheno[,c("Genotype", "X20KernelWeight","Fix","FID")]
names(akw)[2] <- "AKW";
akw$AKW <- as.numeric(as.character(akw$AKW))/20
akw <- akw[!is.na(akw$AKW),]
dim(akw)
#[1] 4363    4
akw2 <- subset(akw,!(Genotype %in% pheno2011$Genotype))
dim(akw2)
#[1] 4003    4
hist(akw, breaks=50, main="AKW")
write.table(akw, "nampheno_akw.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(akw2, "nampheno_akw_subset.txt", sep="\t", quote=FALSE, row.names=FALSE)

##### KRN #######################################################
mymean <- function(x) mean(x, na.rm=T)
krn <- nampheno[,c("Genotype", "EarRowNumber","EarRowNumber.1","Fix","FID")]
krn[,2] <- as.numeric(as.character(krn[,2]));
krn[,3] <- as.numeric(as.character(krn[,3]));
krn$KRN <- apply(krn[, 2:3], 1, mymean);
krn <- krn[!is.na(krn$KRN),]
krn <- krn[, c(1,6,4,5)]
dim(krn)
#[1] 4719    4
krn2 <- subset(krn,!(Genotype %in% pheno2011$Genotype) )
dim(krn2)
#[1] 4346    4
hist(krn$KRN, breaks=50, main="KRN")
write.table(krn, "nampheno_krn.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(krn2, "nampheno_krn_subset.txt", sep="\t", quote=FALSE, row.names=FALSE)

##### KC #######################################################
kc <- nampheno[,c("Genotype", "TotalKernelVolume", "TotalKernelVolume.1","Fix","FID")]
kc[,2] <- as.numeric(as.character(kc[,2]));
kc[,3] <- as.numeric(as.character(kc[,3]));
kc$KC <- apply(kc[, 2:3], 1, mymean);
kc <- kc[!is.na(kc$KC),]
kc <- kc[, c(1,6,4,5)]
dim(kc)
#[1] 4372    4
kc2 <- subset(kc,!(Genotype %in% pheno2011$Genotype) )
dim(kc2)
#[1] 4018    4
hist(kc$KC, breaks=50, main="KC")
write.table(kc, "nampheno_kc.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(kc2, "nampheno_kc_subset.txt", sep="\t", quote=FALSE, row.names=FALSE)

##### TKW #######################################################
tkw <- nampheno[,c("Genotype","EarWeight","EarWeight.1","CobWeight","CobWeight.1","Fix","FID")]
tkw[,2] <- as.numeric(as.character(tkw[,2]));
tkw[,3] <- as.numeric(as.character(tkw[,3]));
tkw[,4] <- as.numeric(as.character(tkw[,4]));
tkw[,5] <- as.numeric(as.character(tkw[,5]));

tkw[,2] <- tkw[,2]-tkw[,4]
tkw[,3] <- tkw[,3]-tkw[,5]
tkw <- tkw[, c(1:3,6:7)]

tkw$TKW <- apply(tkw[, 2:3], 1, mymean);
tkw <- tkw[!is.na(tkw$TKW),]
tkw <- tkw[order(tkw$TKW),]
tkw <- tkw[, c(1,6,4,5)]
tkw <- tkw[c(-1,-2,-nrow(tkw)),]
dim(tkw)
#[1] 4345    4
tkw2 <- subset(tkw,!(Genotype %in% pheno2011$Genotype) )
dim(tkw2)
#[1] 3997    4
hist(tkw$TKW, breaks=50, main="TKW")
write.table(tkw, "nampheno_tkw.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(tkw2, "nampheno_tkw_subset.txt", sep="\t", quote=FALSE, row.names=FALSE)


##### CW #######################################################
cw <- nampheno[,c("Genotype", "CobWeight","CobWeight.1","Fix","FID")]
cw[,2] <- as.numeric(as.character(cw[,2]));
cw[,3] <- as.numeric(as.character(cw[,3]));
cw$CW <- apply(cw[, 2:3], 1, mymean);
cw <- cw[!is.na(cw$CW),]
cw <- cw[, c(1,6,4,5)]
dim(cw)
#[1] 4719    4
cw <- cw[order(cw$CW, decreasing=T),]
cw <- cw[-1:-4,]
hist(cw$CW, breaks=50, main="CW")
cw2 <- subset(cw, !(Genotype %in% pheno2011$Genotype))

write.table(cw, "nampheno_cw.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(cw2, "nampheno_cw_subset.txt", sep="\t", quote=FALSE, row.names=FALSE)

##### CD=Ear Diameter #######################################################
cd <- nampheno[,c("Genotype", "EarDiameter","EarDiameter.1","Fix","FID")]
cd[,2] <- as.numeric(as.character(cd[,2]));
cd[,3] <- as.numeric(as.character(cd[,3]));
cd$CD <- apply(cd[, 2:3], 1, mymean);
cd <- cd[!is.na(cd$CD),]
cd <- cd[, c(1,6,4,5)]
dim(cd)
#[1] 4715    4
cd <- cd[order(cd$CD, decreasing=T),]
cd <- cd[-1:-2,]
hist(cd$CD, breaks=50, main="CD")
cd2 <- subset(cd, !(Genotype %in% pheno2011$Genotype))

write.table(cd, "nampheno_cd.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(cd2, "nampheno_cd_subset.txt", sep="\t", quote=FALSE, row.names=FALSE)

##### CD=Ear Diameter #######################################################
cl <- nampheno[,c("Genotype", "EarLength","EarLength.1","Fix","FID")]
cl[,2] <- as.numeric(as.character(cl[,2]));
cl[,3] <- as.numeric(as.character(cl[,3]));
cl$CL <- apply(cl[, 2:3], 1, mymean);
cl <- cl[!is.na(cl$CL),]
cl <- cl[, c(1,6,4,5)]
dim(cl)
#[1] 5019    4
hist(cl$CL, breaks=50, main="CL")
cl2 <- subset(cl, !(Genotype %in% pheno2011$Genotype))
dim(cl2)
#4622 4
write.table(cl, "nampheno_cl.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(cl2, "nampheno_cl_subset.txt", sep="\t", quote=FALSE, row.names=FALSE)

###########plotting
hist(krn$KRN, breaks=30, main="Buckler KRN", xlab="Kernel Row Number (row)")
hist(cw$CW, breaks=30, main="Buckler CW", xlab="Cob Weight (g)")

