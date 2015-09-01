#Jinliang yang
#Purpose: Use PLINK to trim the correlated SNPs
#updated: 6.7.2012

# Formatting Function:###################################################
# Modified form the GWAS_Hapmap1_SNP_prunning.r
######################################################################

#####---------RUN the PERL program to reformat file
#---------./SNPmatrix_reformat.pl NAM.AGPv2.SNPs+INDELs.pooled-20120217.txt NAM.AGPv2.SNPs+INDELs.pooled

	
setwd("/Users/yangjl/Documents/GWAS2_KRN/SNP/NAM_RNAseq")
tab5rows <- read.delim("NAM.AGPv2.SNPs+INDELs.pooled.snp", nrow=5, header=FALSE)
classes <- sapply(tab5rows, class)
RNAseq <- read.delim("NAM.AGPv2.SNPs+INDELs.pooled.snp", header = FALSE, colClasses = classes)


names(RNAseq) <- c("CHR", "rs","allele", "POSITION","B73","B97","CML103", "CML228","CML247", 
				"CML277","CML322","CML333","CML52","CML69","HP301",
				"IL14H","Ki11","Ki3","Ky21","M162W","M37W",
				"Mo17","Mo18W","Ms71","NC350","NC358","Oh43",
				"Oh7B","P39","Tx303","Tzi8")
	
names(RNAseq) <- c("chr", "rs","allele", "pos","B73", "Z001","Z002", "Z003","Z004", 
				"Z005","Z006","Z007","Z008","Z009","Z010",
				"Z011","Z012","Z013","Z014","Z015","Z016",
				"Z017","Z018","Z019","Z020","Z021","Z022",
				"Z023","Z024","Z025","Z026")

#remove the missing chr ==Mt, Pt and UNKNOWN
RNAseq <- subset(RNAseq, !is.na(chr));
dim(RNAseq)
#[1] 4364987      31

#######################################################################################################################

indel <- read.delim("NAM.AGPv2.SNPs+INDELs.pooled.indel", header=F)
names(indel) <- c("chr", "rs","allele", "pos","B73", "Z001","Z002", "Z003","Z004", 
                  "Z005","Z006","Z007","Z008","Z009","Z010",
                  "Z011","Z012","Z013","Z014","Z015","Z016",
                  "Z017","Z018","Z019","Z020","Z021","Z022",
                  "Z023","Z024","Z025","Z026")
indel <- indel[!is.na(indel$chr),]

indel$pos2 <- c(0, indel$pos[-nrow(indel)])
indel$dis <- indel$pos - indel$pos2
indel <- indel[indel$dis!=1,]
indel <- indel[,1:31]
# 561149


##########################
RNAseq <- rbind(RNAseq, indel)
RNAseq <- RNAseq[order(RNAseq$chr, RNAseq$pos), ]
dim(RNAseq)
#[1] 4926136      31

dsnp <- RNAseq;


#-------->output the .tmap
nm <- names(RNAseq)[5:31]
tmap<- data.frame(family= nm, individual=c(27, 1:26), paternal=0, maternal=0, sex=1, pheno=0)
write.table(tmap, "NAM.AGPv2.SNPs.pooled.tfam", quote=FALSE, row.names=FALSE, col.names=FALSE)

#-------->output the .tped
### Haploid to Diploid:
for (i in 5:31){
  RNAseq[,i] <- paste(RNAseq[,i], RNAseq[,i], sep=" ")
}
RNAseq$geno <- 0;
RNAseq$id <- paste(RNAseq$chr, RNAseq$pos, sep="_")
RNAseq <- RNAseq[, c(1,33,32,4,5:31)]

write.table(RNAseq, "NAM.AGPv2.SNPs.pooled.tped", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE);
#################################
plink --tped NAM.AGPv2.SNPs.pooled.tped --tfam NAM.AGPv2.SNPs.pooled.tfam --missing-genotype N --make-bed --out NAM.AGPv2.SNPs.pooled














