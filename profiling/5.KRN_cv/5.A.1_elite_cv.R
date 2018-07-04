## Jinliang Yang
## updated: June 20, 2014

elite2dsf3_pheno4 <- function(){
  dsnp <- read.csv("~/Documents/KRN_GWAS_v3/GWAS3_proj/data/S.table10_elite_geno.csv")
  ###
  col3 <- dsnp[, 1:3]
  names(col3) <- c("snpid", "chr", "pos")
  source("~/Documents/Rcodes/snpid_chrpos.R")
  col3 <- snpid_chrpos(df=col3, which_dir="snpid2pos")
  elitegeno <- merge(col3, dsnp, by.x="snpid", by.y="SNPID")
  source("~/Documents/Rcodes/dsnp2GenABEL.R")
  dsnp2GenABEL(dsf3=elitegeno, geno_cols=4:ncol(elitegeno), 
               output="~/Documents/KRN_GWAS_v3/GWAS3_proj/cache/elite10.dat")
  
  ### prepare phenotypic data for GenABEL:
  pheno <- read.csv("~/Documents/KRN_GWAS_v3/GWAS3_proj/data/S.table6.1_elite_pheno.csv")
  pheno <- subset(pheno, Population == "Elite Inbred lines")
  names(pheno)[1:2] <- c("id", "KRN")
  pheno <- pheno[, 1:2]
  pheno$sex <- 1
  pheno$age <- 90
  pheno <- pheno[, c("id", "sex", "age", "KRN")]
  write.table(pheno, "~/Documents/KRN_GWAS_v3/GWAS3_proj/cache/elite_pheno.dat")
}

### prepare genotype phenotype for GWAS
elite2dsf3_pheno4()

#dsnp2GenABEL(dsf3=elitegeno, geno_cols=4:ncol(elitegeno), output="test2.raw")

##########################################################################################################

elite_gwas <- function(){
  library("GenABEL")
  convert.snp.text("~/Documents/KRN_GWAS_v3/GWAS3_proj/cache/elite10.dat",
                   "~/Documents/KRN_GWAS_v3/GWAS3_proj/cache/genos_elite10.raw")
  
  setwd("~/Documents/KRN_GWAS_v3/GWAS3_proj/cache")
  elite <- load.gwaa.data(phe = "elite_pheno.dat", gen = "genos_elite10.raw",force = T)
  descriptives.trait(elite)
  qc1 <- check.marker(elite, p.level = 0, het.fdr = 0.001, ibs.threshold = 0.99,
                      callrate=0.4, maf=0.05, perid.call=0.5)
  # In total, 125 (100%) markers passed all criteria
  # In total, 196 (100%) people passed all criteria
  sub1 <- elite[qc1$idok, qc1$snpok];
  res <- qtscore(KRN, data=sub1, trait.type = "gaussian")
  #http://www.genabel.org/sites/default/files/pdfs/GenABEL-tutorial.pdf
  pval1 <- results(res)
  source("~/Documents/Rcodes/get_qval.R")
  pval1$qval <-  get_qval(pval=pval1, pcol="P1df", method="fdr")
  return(pval1)
}

### conducted GWAS and FDR control for pval
elite_pval <- elite_gwas()

source("~/Documents/Rcodes/save.append.R")
save.append(list="elite_pval", file="~/Documents/KRN_GWAS_v3/GWAS3_proj/cache/cv_pvals.RData",
            description="qval after FDR control")

###################
ob <- load("~/Documents/KRN_GWAS_v3/GWAS3_proj/cache/cv_pvals.RData")


