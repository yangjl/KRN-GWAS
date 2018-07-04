# Jinliang Yang
# Purpose: Calculate the LD block
# start: 8.25.2011

# location: 129.186.85.9
setwd("/mnt/02/yangjl/GWAS/plot")
load("data4stackplot.RData")
library("gplots")

######################################################################################################
# PLINK results to find the SNPID
plink_pval <- read.table("wmean_merged_GWAS.qassoc", header=TRUE)
plink_pval$snpid <- paste(plink_pval$CHR, plink_pval$BP, sep="_")

SNP18 <- merge(union, plink_pval, by="snpid")
#SNP18 <- SNP18$SNP
#write.table(SNP18, "snp18.list", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

SNP18 <- SNP18[,c("snpid","SNP", "CHR.x","BP.x")]
names(SNP18) <- c("snpid","SNP", "CHR","BP")

######################################################################################################
# COPY to snp18.list to .7
gcta_mac  --bfile HapApex_PLINK_half  --ld snp18.list  --ld-window 5000  --ld-sig 0.05  --out snp18

######################################################################################################
# COPY results and change back to .9

snptable <- read.table("snp18.rsq.ld", header=TRUE)
r <- read.table("snp18.r.ld",sep=" ",fill=TRUE, header=FALSE)
snp <- read.table("snp18.snp.ld", sep=" ", fill=TRUE, header=FALSE)
ld.list <- list()
for(i in 1:18){
	tsnp <- snptable$target_SNP[i]
	
	temr <- r[i,]
	temr <- t(temr)
	temr <- temr[,1]
	temr <- temr[!is.na(temr)]
	
	temsnp <- snp[i,]
	temsnp <- t(temsnp)
	temsnp <- temsnp[,1]
	temsnp <- temsnp[temsnp!="" & !is.na(temsnp)]

	snp1 <- data.frame(snpid=temsnp, r=temr)

	snpld <- merge(snp1, plink_pval, by.x="snpid", by.y="SNP")
	ld.list[[as.character(tsnp)]] <- snpld	
}



#### Select the 10 of the most signifcant one in the LD block
#### and the ones tightly linked with the target SNP r2 > 0.7
#### snp+- 25k
gene_v1 <- read.table("~/GWAS/Anno/ZmB73_4a.53_filtered_genes.first.transcript.gff", header=T)
gene_v1 <- subset(gene_v1, Feature=="gene")



#inputs: SNP18, snptable

LENGTH= 25000
gene18.list <- list()

for(i in 1:18){
	snpid <- SNP18$snpid[i]
	snpnm <- as.character(SNP18$SNP[i])

	mygene <- subset(gene_v1, (Chr==SNP18$CHR[i] & Start <= SNP18$BP[i]+LENGTH & Start >= SNP18$BP[i]-LENGTH)
					| (Chr==SNP18$CHR[i] & Start >= SNP18$BP[i]-LENGTH & End <= SNP18$BP[i]+LENGTH))
	
	myblock <- ld.list[[snpnm]]
	myblock$absr <- abs(myblock$r)
	
	# select 10 of the most significant p-values in the LD block
	linksnp1 <- head(myblock[order(myblock$P, decreasing=F),],10)
	
	# select the tightly linked SNPs with rsq > 0.7 or top10
	myr <- subset(snptable, target_SNP==snpnm)$max_rsq
	if(myr <= 0.7){
		linksnp2 <- head(myblock[order(myblock$absr),],10)
	}else{
		linksnp2 <- myblock[myblock$r >=0.7,]
	}
	
	linksnp <- rbind(linksnp1, linksnp2)
	
	for(j in 1:nrow(linksnp)){
		temgene <- subset(gene_v1, (Chr==linksnp$CHR[j] & Start <= linksnp$BP[j]+LENGTH & Start >= linksnp$BP[j]-LENGTH)
						 | (Chr==linksnp$CHR[j] & Start >= linksnp$BP[j]-LENGTH & End <= linksnp$BP[j]+LENGTH))
		mygene <- rbind(mygene, temgene)
	}
	mygene <- mygene[!duplicated(mygene$Gene),]
	dim(mygene)
	gene18.list[[snpnm]] <- mygene
	
}

######### Evidence supported genes #############
esgene <- data.frame()
for (i in 1:nrow(gene)){
	mygene <- subset(gene_v1, (Chr==gene$CHR[i] & Start <= gene$BP[i]+20000 & Start >= gene$BP[i]-20000)
					 | (Chr==gene$CHR[i] & Start >= gene$BP[i]-20000 & End <= gene$BP[i]+20000))
	esgene <- rbind(esgene, mygene)
}


#######################################################################################################################






snpld <- ld.list[[5]]
head(snpld[order(snpld$r, decreasing=TRUE),],30)
head(snpld[order(snpld$P, decreasing=F),],30)

palette(rich.colors(32)) # colors: 1 to 32

plot(snpld$BP, y=-log10(snpld$P), col=1 + 31*abs(snp1$r), pch=19, cex=0.4)
