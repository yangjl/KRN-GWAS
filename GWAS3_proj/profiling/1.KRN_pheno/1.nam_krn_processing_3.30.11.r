# NAM RIL KRN phenotypic Data Precessing
# Version 3: Merge all the 2008, 2009, 2010 and 2010g data together.
# 
# undate log:
# start: 03/26/2010
# update: 03/29/2010
# update: 09/05/2010
# update: 02/18/2011
# update: 04/01/2011
# update: 06/20/2011
# update: 09/08/2011


# set the working directory on server:.7
setwd("/Users/yangjl/Documents/KRN_GWAS/pheno")
# local directory:
#setwd("/Users/yangjl/Documents/workingSpace/NAMGWAS/Pheno/nam_krn")
# in the Linux Server .9
# setwd("/mnt/01/yangjl/Research/QTL/NAMkrn")
.libPaths("~/bin/Rlib")
library(ggplot2)



###############################################################################
# Step 1: read in the raw data
###############################################################################
krn0809 <- read.csv("0809_fmt.nam_raw.csv")
dim(krn0809)
krn0809$KRN <- krn0809$KRN*2+2
head(krn0809)

krn10 <- read.csv("10NAM_KRN_raw.csv")
krn10 <- krn10[!is.na(krn10$KRN),]
krn10 <- krn10[,2:7]
krn10$Year <- "10"
names(krn10) <- names(krn0809)
head(krn10)

krn10g <- read.csv("10g_pheno_NAM_RIL_KRN_raw.csv")
krn10g <- krn10g[!is.na(krn10g$KRN),]
names(krn10g) <- names(krn0809)
krn10g$Year <- "10g"
head(krn10g)

fmt.nam <- rbind(krn0809, krn10, krn10g)
fmt.nam[fmt.nam$Note=="MS71",]$Note <- "Ms71"

#------------------------------------------------------------------------------
# overview of the data
dim(fmt.nam)
#[1] 12050     7
head(fmt.nam)
#File  Genotype  Note KRN QC ST Year
#1 DSC04728.JPG Z011E0086 Il14H  16  3  S    8
#2 DSC04732.JPG Z024E0058   P39  18  4  T    8
#3 DSC04733.JPG Z024E0188   P39  18  5  S    8
#4 DSC04736.JPG Z024E0103   P39  16  4  S    8
#5 DSC04737.JPG Z024E0127   P39  18  5  S    8
#6 DSC04738.JPG Z024E0034   P39  16  3  S    8
tail(fmt.nam)
#File  Genotype   Note KRN QC ST Year
#17722 10g-2031-24@ Z008E0105  CML52  12  5  S  10g
#17742 10g-2032-21@ Z005E0071 CML277  16  5  S  10g
#17752 10g-2032-22@ Z005E0071 CML277  14  4  T  10g
#17762 10g-2033-21@ Z008E0050  CML52  10  5  S  10g
#17772 10g-2033-22@ Z008E0050  CML52  12  2  S  10g
#17782 10g-1819-23@ Z008E0047  CML52  14  4  S  10g

#[1] 4780
tem1 <- ddply(fmt.nam, .(Year, Note), "nrow") # summary by year and note(genotype)
tem2 <- ddply(fmt.nam, .(Year), "nrow") # summary by year
tem3 <- ddply(fmt.nam, .(Note), "nrow") # summary by note(genotype)

tem4 <- ddply(fmt.nam, .(Year), summarise, # summary by year, geno and note
	geno = length(unique(Genotype)),
	pop = length(unique(Note))
)

boxplot(KRN~as.factor(QC),data=fmt.nam, xlab="Quality Score",main="NAM KRN", ylab="KRN")

#Count the number in each QC scale.
QCcount <- ddply(fmt.nam, .(QC), "nrow")
QCcount

###############################################################################
# Step 2: ANOVA and Heritability
###############################################################################

fit <- aov(KRN ~ Genotype, data=fmt.nam)

# Broad Sense Heritability
48250/(48250+18491)
#[1] 0.7229439

fit2 <- aov(KRN ~ Genotype, weights=QC, data=fmt.nam)
183036/(183036+67424)
#0.731

#------------------------------------------------------------------------------
# fit by year
table_fityear <- data.frame()
for (year in unique(fmt.nam$Year)){
		subnam <- subset(fmt.nam, Year==year)
		fityear <- aov(KRN ~ Genotype, weights=QC, data=subnam)
		tem <- summary(fityear)[[1]]
		res <- tem[1,2]/(tem[1,2]+tem[2,2])
		table_res <- data.frame(Year=year, h2=res)
		table_fityear <- rbind(table_fityear, table_res)
		
}


#------------------------------------------------------------------------------
# fit by note
table_fitnote <- data.frame()
for (note in unique(fmt.nam$Note)){
	subnam <- subset(fmt.nam, Note==note)
	fitnote <- aov(KRN ~ Genotype, weights=QC, data=subnam)
	tem <- summary(fitnote)[[1]]
	res <- tem[1,2]/(tem[1,2]+tem[2,2])
	table_res <- data.frame(Note=note, h2=res)
	table_fitnote <- rbind(table_fitnote, table_res)
}

###############################################################################
# Step 3: Formatting and output for the following analysis
###############################################################################

# Get the weighted mean
wmean.nam <- ddply(fmt.nam, .(Genotype), summarise,
  Note = unique(Note)[1],
  Year = unique(Year)[1],
  mean = mean(KRN, na.rm=TRUE),
  wmean = weighted.mean(KRN, QC, na.rm=TRUE)
  ) 


#------------------------------------------------------------------------------
# Count the Twisted Ratio(TR) and Twisted Potential(TP)
ts.nam <- ddply(fmt.nam[!is.na(fmt.nam$ST),], .(Genotype), summarise,
  SC = length(ST[ST=="S"]),
  TC = length(ST[ST=="T"])
  )
ts.nam$TR <- with(ts.nam, TC/(TC+SC))
head(ts.nam)

ts.nam$TP <- NA
for(i in 1:nrow(ts.nam)){
if (ts.nam$TC[i]==0 & ts.nam$SC[i]!=0){
  ts.nam$TP[i] <- 0} 
  else if (ts.nam$TC[i] > 0){
  ts.nam$TP[i] <- 1}
  }
hist(ts.nam$TP, breaks=30, main="Distribution of Twisted Ratio", xlab="Twisted Ration", cex.lab=1.5)
table(ts.nam$TP)

#------------------------------------------------------------------------------
# Count the Cross Incompatibility Ratio(CIR->PSSR) and Potential(CIP ->PSSP).
# change to the poor seed set:
ci.nam <- ddply(fmt.nam, .(Genotype), summarise,
  CC = length(QC[QC>=3]),
  CI = length(QC[QC<=2])
  )
ci.nam$CIR <- with(ci.nam, CI/(CI+CC))

ci.nam$CIP <- NULL
for(i in 1:nrow(ci.nam)){
if(ci.nam$CI[i]==0 & ci.nam$CC[i]!=0){
  ci.nam$CIP[i] <- 0}
  else if(ci.nam$CI[i] > 0){
  ci.nam$CIP[i] <- 1}
  }
 table(ci.nam$CIP)

#------------------------------------------------------------------------------
# Merge All the data together
trait.nam <- merge(wmean.nam, ts.nam, by="Genotype", all=TRUE)
trait.nam <- merge(trait.nam, ci.nam, by="Genotype", all=TRUE)

d <- trait.nam
d$FID <- "Z"
d$IID <- "E"
for (i in 1:nrow(d)){
d[i,]$FID <- unlist(strsplit(as.character(d[i,]$Genotype), "E"))[1]
d[i,]$IID <- paste("E", unlist(strsplit(as.character(d[i,]$Genotype), "E"))[2], sep="") 
}
trait <- d[,c(1,14,15,2, 3:13)]
names(trait)[14:15] <- c("PSSR", "PSSP")
head(trait)
tail(trait)
dim(trait)

##### trait ###################################################################
write.table(trait, "nam_pheno_krn_040211.txt", sep="\t", row.names=FALSE, quote=FALSE)

###############################################################################
# Step 4: Plotting
###############################################################################
trait <- read.delim("nam_pheno_krn_040211.txt", header=TRUE)
dim(trait)

wmean.nam <- ddply(trait, .(Note, Year), "nrow")
wmean.nam <- ddply(trait, .(Note), "nrow")
wmean.nam <- wmean.nam[wmean.nam$nrow >50,]

###########################
pdf("pheno_summary.pdf")
par (mfrow=c(2,1))
hist(trait$wmean, ylab="Number of Ears", xlab="KRN Weighted Mean", breaks=50, col="darksalmon", 
main="Histogram of KRN by genotype", axes=FALSE)
axis(1, at=c(10, 15,20), labels= c(10, 15,20) )
axis(2, at=c(0,400,800), labels= c(0,400,800))

abline(v=17.1, lty=2, lwd=2, col="blue")
text(17.7, 600, "B73", col="blue")

boxplot(trait$wmean ~ trait$Note, xaxt="n", col="darksalmon", ann=FALSE, axes = FALSE, 
main="Boxplot of KRN by subpopulation", ylab="KRN Weighted Mean")
abline(h=17.1, col="blue", lty=2, lwd=2)
axis(1, at=1:25, labels= as.character(unique(trait$Note)), las=2, cex.axis=0.8)
axis(2, at=c(10,15,20), labels= c(10, 15, 20))

dev.off()
################################################################################


hist(trait$TR, ylab="", xlab="Twisted Ratio", col="red", main="Histogram of TR by Genotype")
hist(trait$PSSR, ylab="", xlab="Poor Seed Set Ratio", col="red", main="Histogram of PSSR by Genotype")

trait <- trait[with(trait, order(Year, Genotype)),]
q <- qplot(Note, wmean, data=trait, geom="boxplot", main="Boxplot of KRN by population", xlab="", fill="blue", ylab="Weighted Mean")
q + opts(axis.text.x=theme_text(angle=-90, hjust=0))+geom_hline(yintercept=17.1, col="blue")

popl <- ddply(trait, .(Note), summarise,
  wmean = mean(wmean))

#################################################################################
# KRN wmean for GenSel 4-2-2011
#################################################################################

setwd("/Users/yangjl/Documents/workingSpace/NAMGWAS/Pheno/nam_krn")
trait <- read.delim("nam_pheno_krn_040211.txt", header=TRUE)
dim(trait)
head(trait)

wmean <- trait[,c("Genotype", "wmean","FID")]
head(wmean)
dim(wmean)
wmean$Fix <- 1
wmean <- wmean[, c(1,2,4,3)]
#FID -> FID$
write.table(wmean, "pheno_KRN_wmean_040211.txt", sep="\t", row.names=FALSE, quote=FALSE)

#################################################################################
# KRN wmean for Regional GWAS 4-2-2011
#################################################################################
###wmean -----------------------------############
mylm <- lm(wmean~Note, data=trait)
coeff <- mylm$coefficients

p <- as.character(unique(trait$Note))

for (i in 2:25){
	myp <- p[i]
	trait[trait$Note==myp,]$wmean <- trait[trait$Note==myp,]$wmean - coeff[i]
}

## mean -----------------------------############
mylmean <- lm(mean~Note, data=trait)
coeff <- mylmean$coefficients

p <- as.character(unique(trait$Note))

for (i in 2:25){
	myp <- p[i]
	trait[trait$Note==myp,]$mean <- trait[trait$Note==myp,]$mean - coeff[i]
}

mylm2 <- lm(mean~Note, data=trait)

### TR -----------------------------############
mylTR <- lm(TR~Note, data=trait)
coeff <- mylTR$coefficients
p <- as.character(unique(trait$Note))
for (i in 2:25){
	myp <- p[i]
	trait[trait$Note==myp,]$TR <- trait[trait$Note==myp,]$TR - coeff[i]
}

mylm2 <- lm(TR~Note, data=trait)

### CIR -> PSSR -----------------------------############
mylPSSR <- lm(PSSR~Note, data=trait)
coeff <- mylPSSR$coefficients
p <- as.character(unique(trait$Note))
for (i in 2:25){
	myp <- p[i]
	trait[trait$Note==myp,]$PSSR <- trait[trait$Note==myp,]$PSSR - coeff[i]
}

mylm2 <- lm(PSSR~Note, data=trait)

write.table(trait, "nam_pheno_krn-p_040211.txt", row.names=FALSE, quote=TRUE)
#################################################################################
save.image("pheno_NAM_KRN.Rdata")




