#Jinliang yang
#Purpose: Use PLINK to trim the correlated SNPs
#start: 1.23.2012
#updated: June.21.2014

#Note: Running on server 129.186.85.9

# Formatting Function:###################################################
# 1. get rid of the duplicates 
# 2. generate the PLINK map file
# 3. generate the PLINK ped file of parents
# 4. generate parental SNP for imputation
# 5. filter the B73 inconsistent calls
######################################################################

#chr2ped <- function(chr=1){
	
setwd("~/DBcenter/VariationDB/hmp1/")
tab5rows <- read.delim("maizeHapMapV1_B73RefGenV2_20110309_ALL.hmp.txt", header=T, nrow=5)
classes <- sapply(tab5rows, class)

hmp1 <- read.delim("maizeHapMapV1_B73RefGenV2_20110309_ALL.hmp.txt", header=TRUE, colClasses = classes)

nrow(hmp1); #[1] 1624651
	
hmp1<- hmp1[, c(1:5, 12:38)]
names(hmp1) <- c("rs","alleles","chr","pos","strand","B73","B97","CML103",
				"CML228","CML247","CML277","CML322","CML333","CML52","CML69","HP301",
				"IL14H","KI11","KI3","KY21","M162W","M37W","MO17","MO18W","MS71","NC350","NC358",
				"OH43","OH7B","P39","TX303","TZI8")
	
names(hmp1) <- c("rs","alleles","chr","pos","strand","B73","Z001","Z002",
				"Z003","Z004","Z005","Z006","Z007","Z008","Z009","Z010",
				"Z011","Z012","Z013","Z014","Z015","Z016","Z017","Z018","Z019","Z020","Z021",
				"Z022","Z023","Z024","Z025","Z026")

### remove the duplicated rows ###
hmp1$id <- paste(hmp1$chr, hmp1$pos, sep="_");
hmp1 <- hmp1[!duplicated(hmp1$id), ]

hmp1$alleles <- as.character(hmp1$alleles)
hmp1$ref <- substr(hmp1$alleles, 1, 1);
hmp1$alt <- substr(hmp1$alleles, nchar(hmp1$alleles), nchar(hmp1$alleles));

hmp1[hmp1$ref=="-", ]$alleles <- "-/+";
hmp1[hmp1$alt=="-", ]$alleles <- "+/-";
hmp1$ref <- substr(hmp1$alleles, 1, 1);
hmp1$alt <- substr(hmp1$alleles, nchar(hmp1$alleles), nchar(hmp1$alleles));

backup <- hmp1
### remove the one that B73 not consistent with Ref ###
hmp1 <- subset(hmp1, B73 == "N" | B73 == ref);
hmp1$B73 <- hmp1$ref;

### change the -/+ to A/T  ###
alleles <- hmp1$alleles;
hmp1 <- hmp1[, -5];

hmp1[hmp1 == "-"] <- "A";
hmp1[hmp1 == "+"] <- "T";
hmp1$alleles <- alleles;

dsnp <- hmp1;


dim(hmp1)
#1] 1573250      34
#get rid of the loci contain multiple alleles


#-------->output the .tmap
nm <- names(hmp1)[5:31]
hmp1.tmap<- data.frame(family= nm, individual=c(27, 1:26), paternal=0, maternal=0, sex=1, pheno=0)
write.table(hmp1.tmap, "zmhmp1_v2.tfam", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

#-------->output the .tped
### Haploid to Diploid:
for (i in 5:31){
	hmp1[,i] <- paste(hmp1[,i], hmp1[,i], sep=" ")
}
hmp1$geno <- 0
hmp1 <- hmp1[, c(3,32,35,4, 5:31)]

write.table(hmp1, "zmhmp1_v2.tped", sep="\t", quote=FALSE, col.names=FALSE,row.names=FALSE)


################################################################################
# Then PLINK Command:
################################################################################
# 1. make the bed
plink --tped zmhmp1_v2.tped --tfam zmhmp1_v2.tfam --missing-genotype N --make-bed --out zmhmp1_v2_060712







