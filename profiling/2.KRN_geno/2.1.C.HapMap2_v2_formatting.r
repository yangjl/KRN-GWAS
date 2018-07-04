#Jinliang yang
#Purpose: formatting the HapMap2 V2 SNPs for imputation
#start: 6.6.2012

#Note: Running on server 129.186.85.7

# Formatting Function:###################################################
# 1. get rid of the duplicates 
# 2. generate the PLINK map file
# 3. generate the PLINK ped file of parents
# 4. generate parental SNP for imputation
# 5. filter the B73 inconsistent calls
######################################################################


setwd("/Users/yangjl/Documents/VariationDB/VCF/HapMap2/")
source("chr2ped.R")
chr2ped(chrfile="maizeHapMapV2_B73RefGenV2_201203028_chr2.hmp.txt", output="HapMap2_s1/maizeHapMap2_chr2_s1")

chr2ped(chrfile="maizeHapMapV2_B73RefGenV2_201203028_chr1.hmp.txt", output="HapMap2_s1/maizeHapMap2_chr1_s1")
chr2ped(chrfile="maizeHapMapV2_B73RefGenV2_201203028_chr3.hmp.txt", output="HapMap2_s1/maizeHapMap2_chr3_s1")
chr2ped(chrfile="maizeHapMapV2_B73RefGenV2_201203028_chr4.hmp.txt", output="HapMap2_s1/maizeHapMap2_chr4_s1")
chr2ped(chrfile="maizeHapMapV2_B73RefGenV2_201203028_chr5.hmp.txt", output="HapMap2_s1/maizeHapMap2_chr5_s1")
chr2ped(chrfile="maizeHapMapV2_B73RefGenV2_201203028_chr6.hmp.txt", output="HapMap2_s1/maizeHapMap2_chr6_s1")
chr2ped(chrfile="maizeHapMapV2_B73RefGenV2_201203028_chr7.hmp.txt", output="HapMap2_s1/maizeHapMap2_chr7_s1")
chr2ped(chrfile="maizeHapMapV2_B73RefGenV2_201203028_chr8.hmp.txt", output="HapMap2_s1/maizeHapMap2_chr8_s1")
chr2ped(chrfile="maizeHapMapV2_B73RefGenV2_201203028_chr9.hmp.txt", output="HapMap2_s1/maizeHapMap2_chr9_s1")
chr2ped(chrfile="maizeHapMapV2_B73RefGenV2_201203028_chr10.hmp.txt", output="HapMap2_s1/maizeHapMap2_chr10_s1")


################################################################################
# Then PLINK Command:
################################################################################
# 1. make the bed
plink --tped maizeHapMap2_chr1_s1.tped --tfam maizeHapMap2_chr1_s1.tfam --missing-genotype N --make-bed --out zmHapMap2v2_chr1_s1
plink --tped maizeHapMap2_chr2_s1.tped --tfam maizeHapMap2_chr2_s1.tfam --missing-genotype N --make-bed --out zmHapMap2v2_chr2_s1
plink --tped maizeHapMap2_chr3_s1.tped --tfam maizeHapMap2_chr3_s1.tfam --missing-genotype N --make-bed --out zmHapMap2v2_chr3_s1
plink --tped maizeHapMap2_chr4_s1.tped --tfam maizeHapMap2_chr4_s1.tfam --missing-genotype N --make-bed --out zmHapMap2v2_chr4_s1
plink --tped maizeHapMap2_chr5_s1.tped --tfam maizeHapMap2_chr5_s1.tfam --missing-genotype N --make-bed --out zmHapMap2v2_chr5_s1
plink --tped maizeHapMap2_chr6_s1.tped --tfam maizeHapMap2_chr6_s1.tfam --missing-genotype N --make-bed --out zmHapMap2v2_chr6_s1
plink --tped maizeHapMap2_chr7_s1.tped --tfam maizeHapMap2_chr7_s1.tfam --missing-genotype N --make-bed --out zmHapMap2v2_chr7_s1
plink --tped maizeHapMap2_chr8_s1.tped --tfam maizeHapMap2_chr8_s1.tfam --missing-genotype N --make-bed --out zmHapMap2v2_chr8_s1
plink --tped maizeHapMap2_chr9_s1.tped --tfam maizeHapMap2_chr9_s1.tfam --missing-genotype N --make-bed --out zmHapMap2v2_chr9_s1
plink --tped maizeHapMap2_chr10_s1.tped --tfam maizeHapMap2_chr10_s1.tfam --missing-genotype N --make-bed --out zmHapMap2v2_chr10_s1



# 2. calculate the missing rate:
plink --bfile zmHapMap2v2_chr1_s1 --missing --out zmHapMap2v2_chr1_s1
plink --bfile zmHapMap2v2_chr1_s1 --freq --out zmHapMap2v2_chr1_s1
plink --bfile zmHapMap2v2_chr2_s1 --missing --out zmHapMap2v2_chr2_s1
plink --bfile zmHapMap2v2_chr2_s1 --freq --out zmHapMap2v2_chr2_s1
plink --bfile zmHapMap2v2_chr3_s1 --missing --out zmHapMap2v2_chr3_s1
plink --bfile zmHapMap2v2_chr3_s1 --freq --out zmHapMap2v2_chr3_s1
plink --bfile zmHapMap2v2_chr4_s1 --missing --out zmHapMap2v2_chr4_s1
plink --bfile zmHapMap2v2_chr4_s1 --freq --out zmHapMap2v2_chr4_s1
plink --bfile zmHapMap2v2_chr5_s1 --missing --out zmHapMap2v2_chr5_s1
plink --bfile zmHapMap2v2_chr5_s1 --freq --out zmHapMap2v2_chr5_s1
plink --bfile zmHapMap2v2_chr6_s1 --missing --out zmHapMap2v2_chr6_s1
plink --bfile zmHapMap2v2_chr6_s1 --freq --out zmHapMap2v2_chr6_s1
plink --bfile zmHapMap2v2_chr7_s1 --missing --out zmHapMap2v2_chr7_s1
plink --bfile zmHapMap2v2_chr7_s1 --freq --out zmHapMap2v2_chr7_s1
plink --bfile zmHapMap2v2_chr8_s1 --missing --out zmHapMap2v2_chr8_s1
plink --bfile zmHapMap2v2_chr8_s1 --freq --out zmHapMap2v2_chr8_s1
plink --bfile zmHapMap2v2_chr9_s1 --missing --out zmHapMap2v2_chr9_s1
plink --bfile zmHapMap2v2_chr9_s1 --freq --out zmHapMap2v2_chr9_s1
plink --bfile zmHapMap2v2_chr10_s1 --missing --out zmHapMap2v2_chr10_s1
plink --bfile zmHapMap2v2_chr10_s1 --freq --out zmHapMap2v2_chr10_s1


###read in the individual and loci missing rate:
snpfilter <- function(lmissfile="./HapMap2_s1/zmHapMap2v2_s1.lmiss", maffile="./HapMap2_s1/zmHapMap2v2_s1.frq",
                      maf=0.01, mr=0.6) {
  lmiss <- read.table(lmissfile, header=TRUE);
  frq <- read.table(maffile, header=TRUE)
  
  lmiss1 <- merge(lmiss, frq, by="SNP")
  lmiss2 <- subset(lmiss1, MAF >= maf & F_MISS <= mr)
  return(lmiss2)
}

chr1 <- snpfilter(lmissfile="zmHapMap2v2_chr1_s1.lmiss", maffile="zmHapMap2v2_chr1_s1.frq", maf=0.01)
chr2 <- snpfilter(lmissfile="zmHapMap2v2_chr2_s1.lmiss", maffile="zmHapMap2v2_chr2_s1.frq", maf=0.01)
chr3 <- snpfilter(lmissfile="zmHapMap2v2_chr3_s1.lmiss", maffile="zmHapMap2v2_chr3_s1.frq", maf=0.01)
chr4 <- snpfilter(lmissfile="zmHapMap2v2_chr4_s1.lmiss", maffile="zmHapMap2v2_chr4_s1.frq", maf=0.01)
chr5 <- snpfilter(lmissfile="zmHapMap2v2_chr5_s1.lmiss", maffile="zmHapMap2v2_chr5_s1.frq", maf=0.01)
chr6 <- snpfilter(lmissfile="zmHapMap2v2_chr6_s1.lmiss", maffile="zmHapMap2v2_chr6_s1.frq", maf=0.01)
chr7 <- snpfilter(lmissfile="zmHapMap2v2_chr7_s1.lmiss", maffile="zmHapMap2v2_chr7_s1.frq", maf=0.01)
chr8 <- snpfilter(lmissfile="zmHapMap2v2_chr8_s1.lmiss", maffile="zmHapMap2v2_chr8_s1.frq", maf=0.01)
chr9 <- snpfilter(lmissfile="zmHapMap2v2_chr9_s1.lmiss", maffile="zmHapMap2v2_chr9_s1.frq", maf=0.01)
chr10 <- snpfilter(lmissfile="zmHapMap2v2_chr10_s1.lmiss", maffile="zmHapMap2v2_chr10_s1.frq", maf=0.01)

chr <- rbind(chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, chr9, chr10)
write.table(chr_maf5, "HapMap2v2_s2_snplist_23M.txt", sep="\t", row.names=FALSE, quote=FALSE)

####### prepare the density SNP set #########

chr2 <- read.table("maizeHapMap2_chr2_s1.dsnp", header=T, colClasses = classes)
chr3 <- read.table("maizeHapMap2_chr3_s1.dsnp", header=T, colClasses = classes)
chr4 <- read.table("maizeHapMap2_chr4_s1.dsnp", header=T, colClasses = classes)
chr5 <- read.table("maizeHapMap2_chr5_s1.dsnp", header=T, colClasses = classes)
chr6 <- read.table("maizeHapMap2_chr6_s1.dsnp", header=T, colClasses = classes)
chr7 <- read.table("maizeHapMap2_chr7_s1.dsnp", header=T, colClasses = classes)
chr8 <- read.table("maizeHapMap2_chr8_s1.dsnp", header=T, colClasses = classes)
chr9 <- read.table("maizeHapMap2_chr9_s1.dsnp", header=T, colClasses = classes)
chr10 <- read.table("maizeHapMap2_chr10_s1.dsnp", header=T, colClasses = classes)


snplist <- read.table("HapMap2v2_s2_snplist_23M.txt", header=TRUE);
chr2ped_v2 <- function (inputfile="maizeHapMap2_chr1_s1.dsnp", outputfile="maizeHapMap2_chr1_s2"){
  chr5rows<- read.table(inputfile, header=T, nrow=5)
  classes <- sapply(chr5rows, class)
  chr <- read.table(inputfile, header=T, colClasses = classes)
  
  #-------->output the .tmap
  nm <- names(chr)[5:31]
  tmap<- data.frame(family= nm, individual=c(27, 1:26), paternal=0, maternal=0, sex=0, pheno=0)
  write.table(tmap, paste(outputfile, "tfam", sep="."), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  
  
  chrtem <- subset(chr, rs %in% snplist$SNP)
  chrtem$id <- paste(chrtem$chr, chrtem$pos, sep="_")
  chrtem <- chrtem[!duplicated(chrtem$id),]
  #-------->output the .tped
  ### Haploid to Diploid:
  for (i in 5:31){
    chrtem[,i] <- paste(chrtem[,i], chrtem[,i], sep=" ")
  }
  chrtem$geno <- 0
  tped <- chrtem[, c(3,32,33,4, 5:31)]
  write.table(tped, paste(outputfile, "tped", sep="."), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  return(nrow(tped))
  
}

chr2ped_v2(inputfile="maizeHapMap2_chr1_s1.dsnp", outputfile="maizeHapMap2_chr1_s2") #[1] 2844416
chr2ped_v2(inputfile="maizeHapMap2_chr2_s1.dsnp", outputfile="maizeHapMap2_chr2_s2") #[1] 2098022
chr2ped_v2(inputfile="maizeHapMap2_chr3_s1.dsnp", outputfile="maizeHapMap2_chr3_s2") #[1] 2109517
chr2ped_v2(inputfile="maizeHapMap2_chr4_s1.dsnp", outputfile="maizeHapMap2_chr4_s2") #[1] 2140971
chr2ped_v2(inputfile="maizeHapMap2_chr5_s1.dsnp", outputfile="maizeHapMap2_chr5_s2") #[1] 1860437
chr2ped_v2(inputfile="maizeHapMap2_chr6_s1.dsnp", outputfile="maizeHapMap2_chr6_s2") #[1] 1505918
chr2ped_v2(inputfile="maizeHapMap2_chr7_s1.dsnp", outputfile="maizeHapMap2_chr7_s2") #[1] 1513640
chr2ped_v2(inputfile="maizeHapMap2_chr8_s1.dsnp", outputfile="maizeHapMap2_chr8_s2") #[1] 1577175
chr2ped_v2(inputfile="maizeHapMap2_chr9_s1.dsnp", outputfile="maizeHapMap2_chr9_s2") #[1] 1438409
chr2ped_v2(inputfile="maizeHapMap2_chr10_s1.dsnp", outputfile="maizeHapMap2_chr10_s2") #[1] 1349013

# stage2. make the bed
plink --tped maizeHapMap2_chr1_s2.tped --tfam maizeHapMap2_chr1_s2.tfam --missing-genotype N --make-bed --out zmHapMap2v2_chr1_s2
plink --tped maizeHapMap2_chr2_s2.tped --tfam maizeHapMap2_chr2_s2.tfam --missing-genotype N --make-bed --out zmHapMap2v2_chr2_s2
plink --tped maizeHapMap2_chr3_s2.tped --tfam maizeHapMap2_chr3_s2.tfam --missing-genotype N --make-bed --out zmHapMap2v2_chr3_s2
plink --tped maizeHapMap2_chr4_s2.tped --tfam maizeHapMap2_chr4_s2.tfam --missing-genotype N --make-bed --out zmHapMap2v2_chr4_s2
plink --tped maizeHapMap2_chr5_s2.tped --tfam maizeHapMap2_chr5_s2.tfam --missing-genotype N --make-bed --out zmHapMap2v2_chr5_s2
plink --tped maizeHapMap2_chr6_s2.tped --tfam maizeHapMap2_chr6_s2.tfam --missing-genotype N --make-bed --out zmHapMap2v2_chr6_s2
plink --tped maizeHapMap2_chr7_s2.tped --tfam maizeHapMap2_chr7_s2.tfam --missing-genotype N --make-bed --out zmHapMap2v2_chr7_s2
plink --tped maizeHapMap2_chr8_s2.tped --tfam maizeHapMap2_chr8_s2.tfam --missing-genotype N --make-bed --out zmHapMap2v2_chr8_s2
plink --tped maizeHapMap2_chr9_s2.tped --tfam maizeHapMap2_chr9_s2.tfam --missing-genotype N --make-bed --out zmHapMap2v2_chr9_s2
plink --tped maizeHapMap2_chr10_s2.tped --tfam maizeHapMap2_chr10_s2.tfam --missing-genotype N --make-bed --out zmHapMap2v2_chr10_s2


