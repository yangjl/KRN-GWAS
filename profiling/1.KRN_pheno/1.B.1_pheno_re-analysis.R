### Jinliang Yang
### Sep 1st, 2015
### re-analysis of the KRN phenotypic data

# set the working directory on server:.7
library(ggplot2)
library(nlme)

###############################################################################
# Step 1: read in the raw data of Ames
###############################################################################
### IBM KRN data
get_ibm_data <- function(){
    leaf <- read.csv("largedata/IBM_leaf_traits_20101025_KSU.csv")
    ibmkrn <- read.table("largedata/meankrn.txt", header=T)
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
    krnibm08$block <- 1
    krnibm08$plot <- 1:193
    return(krnibm08)
}

krnibm08 <- get_ibm_data()
########## 2008 and 2009
krn0809 <- read.csv("largedata/0809_fmt.nam_raw.csv")
dim(krn0809)
krn0809$KRN <- krn0809$KRN*2+2
head(krn0809)
krn0809$Location <- "Curtiss"
krn0809$block <- 2
krn0809$plot <- as.numeric(krn0809$Genotype)


########## 2010
krn10 <- read.csv("largedata/10NAM_KRN_raw.csv")
krn10 <- krn10[!is.na(krn10$KRN),]
krn10 <- krn10[,2:7]
krn10$Year <- "Y10"
krn10$Location <- "Curtiss"

krn10$block <- 3
krn10$plot <- as.numeric(krn10$Genotype)
names(krn10) <- names(krn0809)
head(krn10)

############ 2010g
krn10g <- read.csv("largedata/10g_pheno_NAM_RIL_KRN_raw.csv")
krn10g <- krn10g[!is.na(krn10g$KRN),]
krn10g$Year <- "10g"
krn10g$Location <- "HawaiiRes"
krn10g$block <- 4
krn10g$plot <- as.numeric(krn10g$Genotype)

########### 2011
krn11 <- read.csv("largedata/LL_KRN_test.csv")
idx <- grep("@", krn11$Genotype)
krn11 <- krn11[idx, ]
krn11$Genotype <- gsub("@", "", krn11$Genotype)
krn11$Year <- "Y11"
names(krn11)[8] <- "Location"
krn11$Row <- gsub("J", "", krn11$Row)
krn11$block <- round(as.numeric(as.character(krn11$Row))/50, 0)
krn11$plot <- as.numeric(krn11$Row)

#########################################################
namlong <- rbind(krnibm08[, c("Genotype", "KRN", "Year", "Location", "block", "plot") ], 
                 krn0809[, c("Genotype", "KRN", "Year", "Location", "block", "plot") ], 
                 krn10[, c("Genotype", "KRN", "Year", "Location", "block", "plot") ], 
                 krn10g[, c("Genotype", "KRN", "Year", "Location", "block", "plot") ], 
                 krn11[, c("Genotype", "KRN", "Year", "Location", "block", "plot") ]
)


namlong <- namlong[!is.na(namlong$KRN),]
#namlong$KRN <- round(bct(namlong$KRN, lambda=0.4), 0);
#namlong$KRN <- round(bcPower(namlong$KRN, lambda=0.4, jacobian.adjusted = TRUE), 0);


write.table(namlong, "largedata/pheno_NAMRIL_KRN_051612.csv", sep=",", row.names=FALSE, quote=FALSE)


namlong <- read.csv("largedata/pheno_NAMRIL_KRN_051612.csv")
namlong$KRN <- as.numeric(namlong$KRN)
namlong <- namlong[!is.na(namlong$KRN),]

library("lme4")
lmeout <- lmer(KRN ~ Genotype + (1 | Year) + + (1 | Location) + (1 | Genotype:Year)+ (1 | Genotype:Location), data=namlong)




