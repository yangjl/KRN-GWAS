# Jinliang Yang
# 6/20/2014
# update: comparing the old and new version of analysis

ob <- load("~/Documents/KRN_GWAS_v3/GWAS3_proj/cache/cv_pvals.RData")

setwd("~/Documents/KRN_GWAS_v3/GWAS3_proj/data/")
elite1 <- read.csv("Elite_pval1.csv")
elite2 <- read.csv("Elite_pval2.csv")
#Binomial XBSA
xbsa2 <- read.csv("XBSA_binary.csv")
bsle2 <- read.csv("BSLE_binary.csv")


com2pvals <- function(newp=elite_pval, oldp=elite1){
  newp$SNP <- paste(newp$Chromosome, newp$Position, sep="_")
  twop <- merge(newp, oldp, by="SNP")
  print(cor.test(twop$qval.x, twop$qval.y))  
}

####
com2pvals(newp=elite_pval, oldp=elite1)
com2pvals(newp=xbsa_pval1, oldp=xbsa2)
com2pvals(newp=xbsa_pval2, oldp=xbsa2)

com2pvals(newp=BSLE_pval, oldp=bsle2)





