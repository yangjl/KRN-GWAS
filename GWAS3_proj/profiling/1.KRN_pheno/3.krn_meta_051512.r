# Phenotype for KRN meta-analysis
# Jinliang YANG
# 5/14/2012

# set the working directory on server:.7
setwd("/Users/yangjl/Documents/GWAS2_KRN/pheno")
library(ggplot2)
library(nlme)

###############################################################################
# Step 1: read in the raw data of Ames
###############################################################################
### IBM KRN data
leaf <- read.csv("~/Documents/GWAS2_KRN/pheno/IBM_leaf_traits_20101025_KSU.csv")
ibmkrn <- read.table("~/Documents/GWAS2_KRN/pheno/meankrn.txt", header=T)
idx1 <- grep("B73x", ibmkrn$Genotype)
bx <- ibmkrn[idx1,]
bx$Genotype <- gsub("B73x", "", bx$Genotype)

idx2 <- grep("Mo17", ibmkrn$Genotype)
mx <- ibmkrn[idx2,]
mx$Genotype <- gsub("Mo17x", "", mx$Genotype)

leaf$MO.Number <- gsub("O", "0", leaf$MO.Number)
leaf <- leaf[, 1:2]
names(leaf) <- c("Genotype", "ID")

ibm <- merge(leaf, ibmkrn[, 1:2], by.x="ID", by.y="Genotype", all.x=T)
names(ibm)[3] <- "IBMRIL"
ibm <- merge(ibm, bx[,1:2], by.x="ID", by.y="Genotype", all.x=TRUE)
names(ibm)[4] <- "BxRIL"
ibm <- merge(ibm, mx[,1:2], by.x="ID", by.y="Genotype", all.x=TRUE)
names(ibm)[5] <- "MxRIL"

krnibm08 <- ibm[, 2:3]
krnibm08$Year <- 8;
krnibm08$Location <- "Johnson"
names(krnibm08)[2] <- "KRN"

########## 2008 and 2009
krn0809 <- read.csv("0809_fmt.nam_raw.csv")
dim(krn0809)
krn0809$KRN <- krn0809$KRN*2+2
head(krn0809)
krn0809$Location <- "Curtiss"

########## 2010
krn10 <- read.csv("10NAM_KRN_raw.csv")
krn10 <- krn10[!is.na(krn10$KRN),]
krn10 <- krn10[,2:7]
krn10$Year <- "Y10"
krn10$Location <- "Curtiss"
names(krn10) <- names(krn0809)
head(krn10)

############ 2010g
krn10g <- read.csv("10g_pheno_NAM_RIL_KRN_raw.csv")
krn10g <- krn10g[!is.na(krn10g$KRN),]
krn10g$Year <- "10g"
krn10g$Location <- "HawaiiRes"

########### 2011
krn11 <- read.csv("~/Documents/Heterosis_GWAS/pheno2011/krn11_NAMRIL_051612.csv")
idx <- grep("@", krn11$Genotype)
krn11 <- krn11[idx, ]
krn11$Genotype <- gsub("@", "", krn11$Genotype)
krn11$Year <- "Y11"
names(krn11)[8] <- "Location"

#########################################################
namlong <- rbind(krnibm08[, c("Genotype", "KRN", "Year", "Location") ], krn0809[, c("Genotype", "KRN", "Year", "Location") ], 
                 krn10[, c("Genotype", "KRN", "Year", "Location") ], 
                 krn10g[, c("Genotype", "KRN", "Year", "Location") ], krn11[, c("Genotype", "KRN", "Year", "Location") ]
                 )

#bct computes the Box-Cox family of transforms: y = (y^lambda - 1)/(lambda*gm^(lambda-1)), 
#where gm is the geometric mean of the y's. returns log(y)*gm when lambda equals 0.
bct <- function(y,lambda){
  gm <- exp( mean( log(y) ) )
  if(lambda==0) return( log(y)*gm )
  yt <- (y^lambda - 1)/( lambda * gm^(lambda-1) )
  return(yt)
}


namlong <- namlong[!is.na(namlong$KRN),]
#namlong$KRN <- round(bct(namlong$KRN, lambda=0.4), 0);
#namlong$KRN <- round(bcPower(namlong$KRN, lambda=0.4, jacobian.adjusted = TRUE), 0);


write.table(namlong, "pheno_NAMRIL_KRN_051612.csv", sep=",", row.names=FALSE, quote=FALSE)

lmeout <- lme(KRN~Genotype, data=namlong, random=~1|Year/Location)
ped.hat1 <- lmeout$coef$fixed;
ped.hat1[-1] <- ped.hat1[-1]+ped.hat1[1];
names(ped.hat1)[1]="Z001E0001";
names(ped.hat1) <- gsub("Genotype", "", names(ped.hat1));
tped <- data.frame(Genotype=names(ped.hat1), KRN=ped.hat1)

###############################################################################
# Step 2: My KRN trait
###############################################################################

#Note that we transformed BN and ERN.
#The raw data file is large. Here is a dropbox link:
#The abbreviations in this file are: db=TL, dc=SL, dg=BN,si=CL,sf=CD,so=ERN.
library(nlme)
inflo7 <- read.table("~/Documents/GWAS2_KRN/pheno/maize_inflo7_rawdata.txt", header=TRUE)
pbkrn <- subset(inflo7, trait=="so")

pbkrn$Genotype <- paste("ZZ", 1000+pbkrn$pop, "EE", 10000+pbkrn$entry, sep="")
pbkrn$Genotype <- gsub("Z1", "", pbkrn$Genotype)
pbkrn$Genotype <- gsub("E1", "", pbkrn$Genotype)
table(pbkrn$location)

pbBLUP <- list();
loc="08A";
for(loc in unique(pbkrn$location)){
  myloc <- subset(pbkrn, location == loc);
  dim(myloc)
  lmeout <- lme(value~Genotype, data=myloc, random=~1|pblock/block);
  ### stop here:
  
  ped.hat1 <- lmeout$coef$fixed;
  ped.hat1[-1] <- ped.hat1[-1]+ped.hat1[1];
  
  names(ped.hat1)[1]= myloc$Genotype[1];
  names(ped.hat1) <- gsub("Genotype", "", names(ped.hat1));
  tped <- data.frame(Genotype=names(ped.hat1), KRN=ped.hat1);
  
  pbBLUP[[loc]] <- tped;
  
}

write.table(tped, "tem9.csv", sep=",", row.names=FALSE, quote=FALSE)

##############################################################################################
#tem1 <- read.csv("tem1.csv") # Ames
tem2 <- read.csv("/Users/yangjl/Documents/GWAS2_KRN/pheno/tem2.csv") # 06A   06CL1 06FL1 06PR  06U   07A   07CL  08A 
tem3 <- read.csv("/Users/yangjl/Documents/GWAS2_KRN/pheno/tem3.csv") # 06A   06CL1 06FL1 06PR  06U   07A   07CL  08A 
tem4 <- read.csv("/Users/yangjl/Documents/GWAS2_KRN/pheno/tem4.csv") # 06A   06CL1 06FL1 06PR  06U   07A   07CL  08A 
tem5 <- read.csv("/Users/yangjl/Documents/GWAS2_KRN/pheno/tem5.csv") # 06A   06CL1 06FL1 06PR  06U   07A   07CL  08A 
tem6 <- read.csv("/Users/yangjl/Documents/GWAS2_KRN/pheno/tem6.csv") # 06A   06CL1 06FL1 06PR  06U   07A   07CL  08A 
tem7 <- read.csv("/Users/yangjl/Documents/GWAS2_KRN/pheno/tem7.csv") # 06A   06CL1 06FL1 06PR  06U   07A   07CL  08A 
tem8 <- read.csv("/Users/yangjl/Documents/GWAS2_KRN/pheno/tem8.csv") # 06A   06CL1 06FL1 06PR  06U   07A   07CL  08A 
tem9 <- read.csv("/Users/yangjl/Documents/GWAS2_KRN/pheno/tem9.csv") # 06A   06CL1 06FL1 06PR  06U   07A   07CL  08A 

names(tem2)[2] <- "L06A"
names(tem3)[2] <- "L06CL1"
names(tem4)[2] <- "L06FL1"
names(tem5)[2] <- "L06PR"
names(tem6)[2] <- "L06U"
names(tem7)[2] <- "L07A"
names(tem8)[2] <- "L07CL"
names(tem9)[2] <- "L08A"


temkrn <- merge(tped, tem2, all=TRUE, by="Genotype")
temkrn <- merge(temkrn, tem3, all=TRUE, by="Genotype");
temkrn <- merge(temkrn, tem4, all=TRUE, by="Genotype")
temkrn <- merge(temkrn, tem5, all=TRUE, by="Genotype")
temkrn <- merge(temkrn, tem6, all=TRUE, by="Genotype")
temkrn <- merge(temkrn, tem7, all=TRUE, by="Genotype")
temkrn <- merge(temkrn, tem8, all=TRUE, by="Genotype")
temkrn <- merge(temkrn, tem9, all=TRUE, by="Genotype")

write.table(temkrn, "krn_meta_9reps.csv", sep=",", row.names=FALSE, quote=FALSE)
pairs(temkrn[, 2:10], text.panel = diag, upper.panel=panel.smooth, 
 lower.panel=panel.cor, gap=0, main="KRN Correlation Plots of three obs")
  

###########################
#tem1 <- read.csv("tem1.csv") # Ames
tem2 <- read.csv("/Users/yangjl/Documents/GWAS2_KRN/pheno/tem2.csv") # 06A   06CL1 06FL1 06PR  06U   07A   07CL  08A 
tem3 <- read.csv("/Users/yangjl/Documents/GWAS2_KRN/pheno/tem3.csv") # 06A   06CL1 06FL1 06PR  06U   07A   07CL  08A 
tem4 <- read.csv("/Users/yangjl/Documents/GWAS2_KRN/pheno/tem4.csv") # 06A   06CL1 06FL1 06PR  06U   07A   07CL  08A 
tem5 <- read.csv("/Users/yangjl/Documents/GWAS2_KRN/pheno/tem5.csv") # 06A   06CL1 06FL1 06PR  06U   07A   07CL  08A 
tem6 <- read.csv("/Users/yangjl/Documents/GWAS2_KRN/pheno/tem6.csv") # 06A   06CL1 06FL1 06PR  06U   07A   07CL  08A 
tem7 <- read.csv("/Users/yangjl/Documents/GWAS2_KRN/pheno/tem7.csv") # 06A   06CL1 06FL1 06PR  06U   07A   07CL  08A 
tem8 <- read.csv("/Users/yangjl/Documents/GWAS2_KRN/pheno/tem8.csv") # 06A   06CL1 06FL1 06PR  06U   07A   07CL  08A 
tem9 <- read.csv("/Users/yangjl/Documents/GWAS2_KRN/pheno/tem9.csv") # 06A   06CL1 06FL1 06PR  06U   07A   07CL  08A 

tped$Location <- "Ames";
tem2$Location <- "L06A";
tem3$Location <- "L06CL1";
tem4$Location <- "L06FL1";
tem5$Location <- "L06PR";
tem6$Location <- "L06U";
tem7$Location <- "L07A";
tem8$Location <- "L07CL";
tem9$Location <- "L08A";

krn9rep <- rbind(tped, tem2, tem3, tem4, tem5, tem6, tem7, tem8, tem9)
krn9rep <- krn9rep[!is.na(krn9rep$KRN),]
hist(krn9rep$KRN, break=30)

lmeout2 <- lme(KRN~Genotype, data=krn9rep, random=~1|Location)
ped.hat2 <- lmeout2$coef$fixed;
ped.hat2[-1] <- ped.hat2[-1]+ped.hat2[1];
names(ped.hat2)[1]="Z001E0001";
names(ped.hat2) <- gsub("Genotype", "", names(ped.hat2));
tped2 <- data.frame(Genotype=names(ped.hat2), KRN=ped.hat2)


write.table(tped2, "pheno_krn_meta_9reps_PLUP.txt", sep="\t", row.names=FALSE, quote=FALSE)

hist(1:100, col="red")
hist(2:99, col="blue", add=TRUE)

##################################################################
krnblup <- read.table("pheno_krn_meta_9reps_PLUP.txt", header=T)
idx0 <- grep("Z027", krnblup$Genotype)
krnblup <- krnblup[-idx0,]

ninerep <- read.csv("krn_meta_9reps.csv")
idx <- grep("Z027", ninerep$Genotype)
ninerep <- ninerep[-idx, ]


names(krnblup)[2] <- "BLUP"
names(ninerep)[2] <- "Ames"
krn <- merge(krnblup, ninerep, by="Genotype")

source("~/Documents/Codes/reusableRcode.r")
pairs(krn[, 2:10], text.panel = diag, upper.panel=panel.smooth, 
      lower.panel=panel.cor, gap=0, main="KRN Correlation Plots")

############
krnblup$Note <- NA;
krnblup$Genotype <- as.character(krnblup$Genotype)
for(i in 1:nrow(krnblup)){
  tem <- unlist(strsplit(krnblup$Genotype[i], split="E"));
  krnblup$Note[i] <- tem[1];
  print(i)
}


###########################
#pdf("pheno_summary.pdf")
par (mfrow=c(2,1))
hist(krnblup$BLUP, xlab="KRN BLUP", breaks=50, col="darksalmon", 
     main="Histogram of KRN by genotype", axes=FALSE, xlim=c(8,24), ylim=c(0, 800))
axis(1, at=c(8, 12,16, 20, 24), labels= c(8, 12,16, 20, 24) )
axis(2, at=c(0,400,800), labels= c(0,400,800))

abline(v=17.1, lty=2, lwd=2, col="blue")
text(17.7, 600, "B73", col="blue")

boxplot(krnblup$BLUP ~ krnblup$Note, xaxt="n", col="darksalmon", ann=FALSE, axes = FALSE, 
        main="Boxplot of KRN by subpopulation", ylab="KRN BLUP")
abline(h=17.1, col="blue", lty=2, lwd=2)
axis(1, at=1:26, labels= as.character(unique(krnblup$Note)), las=2, cex.axis=0.8)
axis(2, at=c(10,15,20), labels= c(10, 15, 20))

#dev.off()

