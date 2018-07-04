## Jinliang Yang
## updated: June 20, 2014

XBSA2dsf3_pheno4 <- function(){
  dsnp <- read.csv("~/Documents/KRN_GWAS_v3/GWAS3_proj/data/S.table11_usda_pi_geno.csv")
  ###
  col3 <- dsnp[, 1:3]
  names(col3) <- c("snpid", "chr", "pos")
  source("~/Documents/Rcodes/snpid_chrpos.R")
  col3 <- snpid_chrpos(df=col3, which_dir="snpid2pos")
  pigeno <- merge(col3, dsnp, by="snpid")
  source("~/Documents/Rcodes/dsnp2GenABEL.R")
  dsnp2GenABEL(dsf3=pigeno, geno_cols=4:ncol(pigeno), 
               output="~/Documents/KRN_GWAS_v3/GWAS3_proj/cache/pi10.dat")
  
  ### prepare phenotypic data for GenABEL:
  pheno <- read.csv("~/Documents/KRN_GWAS_v3/GWAS3_proj/data/S.table6.1_elite_pheno.csv")
  pheno <- subset(pheno, Population == "USDA Extreme KRN accessions")
  
  pheno$bKRN <- NA
  pheno[pheno$Note=="High KRN", ]$bKRN <- 1
  pheno[pheno$Note=="Random KRN", ]$bKRN <- 0;
  pheno[pheno$Note=="Low KRN",]$bKRN <- -1;
  
  pheno <- pheno[, c(1:2,6)]
  names(pheno) <- c("id", "KRN", "bKRN")
  
  pheno$sex <- 1
  pheno$age <- 90
  pheno <- pheno[, c("id", "sex", "age", "KRN", "bKRN")]
  write.table(pheno, "~/Documents/KRN_GWAS_v3/GWAS3_proj/cache/xbsa_pheno.dat")
}

### prepare geno and pheno for GWAS
XBSA2dsf3_pheno4()

#dsnp2GenABEL(dsf3=elitegeno, geno_cols=4:ncol(elitegeno), output="test2.raw")

##########################################################################################################

xbsa_gwas <- function(){
  library("GenABEL")
  convert.snp.text("~/Documents/KRN_GWAS_v3/GWAS3_proj/cache/pi10.dat",
                   "~/Documents/KRN_GWAS_v3/GWAS3_proj/cache/genos_xbsa10.raw")
  
  setwd("~/Documents/KRN_GWAS_v3/GWAS3_proj/cache")
  xbsa <- load.gwaa.data(phe = "xbsa_pheno.dat", gen = "genos_xbsa10.raw",force = T)
  descriptives.trait(xbsa)
  qc1 <- check.marker(xbsa, p.level = 0, het.fdr = 0.001, ibs.threshold = 0.99,
                      callrate=0.4, maf=0.05, perid.call=0.5)
  # In total, 125 (100%) markers passed all criteria
  # In total, 196 (100%) people passed all criteria
  sub1 <- xbsa[qc1$idok, qc1$snpok];
  res1 <- qtscore(KRN, data=sub1, trait.type = "gaussian")
  #http://www.genabel.org/sites/default/files/pdfs/GenABEL-tutorial.pdf
  pval1 <- results(res1)
  source("~/Documents/Rcodes/get_qval.R")
  pval1$qval <-  get_qval(pval=pval1, pcol="P1df", method="fdr")
  
  res2 <- qtscore(bKRN, data=sub1, trait.type = "gaussian")
  pval2 <- results(res2)
  pval2$qval <-  get_qval(pval=pval2, pcol="P1df", method="fdr")
  resls <- list()
  resls[['pval1']] <- pval1;
  resls[['pval2']] <- pval2;
  return(resls)
}

### conducted GWAS with FDR control for binomial and quantitative traits
xbsa_pval <- xbsa_gwas()
xbsa_pval1 <- xbsa_pval[[1]] ### KRN quantative trait
xbsa_pval2 <- xbsa_pval[[2]] ## binary trait with high=1, random=0 and low=-1

source("~/Documents/Rcodes/save.append.R")
save.append(list=c("xbsa_pval1", "xbsa_pval2"), 
            file="~/Documents/KRN_GWAS_v3/GWAS3_proj/cache/cv_pvals.RData",
            description=c("KRN quantative trait",
                          "binary trait with high=1, random=0 and low=-1"))

ob <- load("~/Documents/KRN_GWAS_v3/GWAS3_proj/cache/cv_pvals.RData")


