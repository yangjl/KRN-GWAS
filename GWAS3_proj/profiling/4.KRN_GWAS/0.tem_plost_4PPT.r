# Jinliang Yang
# Purpose: plot results
# start: 2.9.2012

########## plot the Missing Rate
par(mfrow=c(1,2))

tem2 <- read.table("~/Desktop/zmhapmap1_v2_012622.lmiss", header=T)
hist(tem2$F_MISS, main="HapMap1 SNP", xlab="Missing Rate", xlim=c(0, 1), ylim=c(0, 1e6))
abline(v=0.6, col="red")
hist(tem$mr, main="NAM RNA-seq SNP", xlab="Missing Rate", xlim=c(0,1))
abline(v=0.6, col="red")


# plot the MAF
maf1 <- read.table("~/Desktop/zmhapmap1_v2_012622.frq", header=T)
maf2 <- read.table("~/Desktop/NAM_RNAseq_012712.frq", header=T)
nrow(maf1[maf1$MAF<0.1,])
nrow(maf2[maf2$MAF<0.1,])

hist(maf1$MAF, main="HapMap1 SNP", xlab="MAF", xlim=c(0,0.5), ylim=c(0, 3e5))
abline(v=0.1, col="red")
hist(maf2$MAF, main="NAM RNA-seq SNP", xlab="MAF", xlim=c(0, 0.5), ylim=c(0,3e5))
abline(v=0.1, col="red")

