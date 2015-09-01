# Jinliang Yang
# Purpose: Stackplot
# start: 2.11.2012
# add the RNA-seq LM regression data

# location: 129.186.85.7
setwd("~/Documents/GWAS2_KRN/plot")
# a table of 10 chr length of the maize genome:


######################################################################################################
# This chromosome length table need to be updated:
plink_pval <- read.table("pval_AGPv2_2M_plink.txt", header=TRUE)

d1 <- plink_pval[,1:3] #CHR,SNP,BP
chrlength <- data.frame()
for (i in 1:10){
	idx <- which.max(subset(d1, CHR==i)$BP)
	tem <- subset(d1, CHR==i)[idx,]
	chrlength <- rbind(chrlength, tem)	
}

write.table(chrlength, "chr_length_B73v2.csv", sep=",", row.names=FALSE, quote=FALSE)


# RUN from here:
######################################################################################################
# function to change the chr position to a accumulative pos
# input: CHR, BP, P
# output: CHR, BP, P, pos
newpos <- function(dataframe, GAP=5000000){
	d <- dataframe
	
	if (!("CHR" %in% names(d) & "BP" %in% names(d) )) stop("Make sure your data frame contains columns CHR and BP")
	
	cl <- read.csv("chr_length_B73v2.csv")
	d$pos=d$BP
	cl$accumpos <- cl$BP
	cl<- cl[order(cl$CHR),]
	
	for (i in 2:10) {
		cl[cl$CHR==i,]$accumpos <- cl[cl$CHR==(i-1),]$accumpos + cl[cl$CHR==i,]$accumpos+GAP
		d[d$CHR==i, ]$pos=d[d$CHR==i, ]$BP+ cl[cl$CHR==(i-1),]$accumpos +GAP
	}
	return(d)
}



#######################################################################################################

gap=10000000
###########################
###--------Regular GWAS---
###########################
gwas <- read.table("pval_AGPv2_2M_plink.txt", header=TRUE)
gwas <- gwas[, c("CHR", "BP", "P", "SNP")]
gwas <- newpos(gwas, GAP=gap)
gwas$logp <- -log10(gwas$P)

gwas$snpid <- paste(gwas$CHR, gwas$BP, sep="_")
gwasbackup <- gwas
dim(gwas)
#2058739       7
gwas <- subset(gwas, logp > 2)

##### borrow the accumulated pos
map <- gwasbackup

###########################
###--------Bayesian GWAS---
###########################
bayes <- read.table("KRN_wmean_SNPv2.mrkRes1", header=TRUE)
#bayes$snpid <- paste(bayes$chr, bayes$pos, sep="_")
#bayes82 <- merge(kas82, bayes, by="snpid")
#dim(bayes82)
bayesbckup <- bayes
bayes <- bayes[, c("ID", "ModelFreq")]
bayes <- merge(bayes, map, by.x="ID", by.y="snpid")
bayes <- subset(bayes, ModelFreq!=0)

bayes <- bayes[, c(1,2,3, 7)]
#names(bayes)[1:2] <- c('CHR', 'BP')
#bayes <- newpos(bayes, GAP=gap)
dim(bayes)
#[1] 1265269       4

###########################################
###--------KAS info to version 2--------
###########################################

kas82 <- read.csv("subsnp82.csv")
v1 <- read.table("wmean_merged_GWAS.qassoc", header=T)
v1$snpid <- paste(v1$CHR,v1$BP, sep="_")
kas82 <- merge(kas82, v1, by="snpid");
kas82 <- kas82[,1:9]
dim(kas82)
#82  9

kas82 <- merge(kas82, gwasbackup, by="SNP");
kas82$newid <- paste(kas82$CHR.x, kas82$BP, sep="_")

kas82 <- merge(kas82, bayesbckup, by.x="newid", by.y="ID")


###########################
###------- QTL LOD---------
###########################

mrk <- read.table("~/Documents/GWAS2_KRN/QTL/NAMKRN_joint/mynam_krn.jz63", header=FALSE)
names(mrk) <- c("Chr", "Num", "Marker", "Pos", "LR", "Add")
mrk$UID <- paste(mrk$Chr, mrk$Num, sep="_")
mrk <- mrk[,c(1:4,7)]

#--read in the physical map (missing ones been imputed)
phymap <- read.delim("~/Documents/GWAS2_KRN/QTL/NAM_map_AGPv2.imputed.txt", header=TRUE)

myphymap <- merge(mrk, phymap, by.x="Marker", by.y="marker")
dim(myphymap)
#[1] 1864   10
myphymap <- myphymap[!duplicated(myphymap$UID),]

cim6 <- read.table("~/Documents/GWAS2_KRN/QTL/NAMKRN_joint/mynam_krn.z6", header=TRUE)
cim6 <- cim6[,1:8]
cim6 <- cim6[, c(1,2,3,4)]
cim6$ID <- paste(cim6$c, cim6$m, sep="_")
cim6$position <- 100* as.numeric(cim6$position)
names(cim6)[4] <- "lod"
cim6$c <- as.numeric(cim6$c)
cim6$lod <- 0.217*as.numeric(cim6$lod)

cim6 <- cim6[order(cim6$c, cim6$m, cim6$lod, decreasing=TRUE),]
cim6 <- cim6[!duplicated(cim6$ID),]

cim6 <- merge(cim6, myphymap, by.x="ID", by.y="UID")
cim6 <- cim6[,c("ch", "imputepos", "genetpos", "Marker", "lod")]
names(cim6)[1:2] <- c("CHR", "BP")

cim6 <- newpos(cim6, GAP=gap)
dim(cim6)
#[1] 886   5

###########################
###------- QTL Interval----
###########################
qtlthred <- read.table("~/Documents/GWAS2_KRN/QTL/NAMKRN_joint/mynam_krn.z6e", header=FALSE)
qtlcutoff <- quantile(qtlthred$V2, 0.95)*0.217

qtl <- read.table("~/Documents/GWAS2_KRN/QTL/Joint_CIM6_summary_v2.txt", header=T)
names(qtl)[c(1,6)] <- c("CHR", "BP")
qtl <- newpos(qtl, GAP=gap)
names(qtl)[c(6,16)] <- c("SigPhy","SigPos")
names(qtl)[9] <- "BP"
qtl <- newpos(qtl, GAP=gap)
names(qtl)[c(9, 17)] <- c("LeftPhy","LeftPos")
names(qtl)[12] <- "BP"
qtl <- newpos(qtl, GAP=gap)
names(qtl)[c(12, 18)] <- c("RightPhy", "RightPos")
dim(qtl)
#16	18

############################
# RNA-sea Regression
############################
rnaseq <- read.csv("RNAseq_lm.csv")
rnaseq <- subset(rnaseq, Chr!="chrMt"&Chr!="chrPt"&Chr!="chrUNKNOWN")
rnaseq$Chr <- gsub("chr", "", rnaseq$Chr)
names(rnaseq)[2:3] <- c("CHR", "BP")
rnaseq <- newpos(rnaseq, gap)

###################################################
save.image("data4stackplot_v2.RData")	











