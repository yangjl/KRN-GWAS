## Jinliang Yang
## updated: June 20, 2014

BSLE2dsf3_pheno4 <- function(){
  dsnp <- read.csv("~/Documents/KRN_GWAS_v3/GWAS3_proj/data/S.table12_BSLE_geno.csv")
  ###
  col3 <- dsnp[, 1:3]
  names(col3) <- c("snpid", "chr", "pos")
  source("~/Documents/Rcodes/snpid_chrpos.R")
  col3 <- snpid_chrpos(df=col3, which_dir="snpid2pos")
  bslegeno <- merge(col3, dsnp, by.x="snpid", by.y="SNPID")
  source("~/Documents/Rcodes/dsnp2GenABEL.R")
  dsnp2GenABEL(dsf3=bslegeno, geno_cols=4:ncol(bslegeno), 
               output="~/Documents/KRN_GWAS_v3/GWAS3_proj/cache/bsle10.dat")
  
  ### prepare phenotypic data for GenABEL:
  pheno <- read.csv("~/Documents/KRN_GWAS_v3/GWAS3_proj/data/S.table6.1_elite_pheno.csv")
  pheno <- subset(pheno, Population == "BSLE" & Note %in% c("C30SE", "C30LE"))
  
  pheno[pheno$Note=="C30SE", ]$KRN <- 1
  pheno[pheno$Note=="C30LE", ]$KRN <- 0;
    
  pheno <- pheno[, 1:2]
  names(pheno) <- c("id", "KRN")
  pheno$sex <- 1
  pheno$age <- 90
  pheno <- pheno[, c("id", "sex", "age", "KRN")]
  write.table(pheno, "~/Documents/KRN_GWAS_v3/GWAS3_proj/cache/bsle_pheno.dat")
}

### prepare geno and pheno for GWAS
BSLE2dsf3_pheno4()

#dsnp2GenABEL(dsf3=elitegeno, geno_cols=4:ncol(elitegeno), output="test2.raw")

##########################################################################################################

BSLE_gwas <- function(){
  library("GenABEL")
  convert.snp.text("~/Documents/KRN_GWAS_v3/GWAS3_proj/cache/bsle10.dat",
                   "~/Documents/KRN_GWAS_v3/GWAS3_proj/cache/genos_bsle10.raw")
  
  setwd("~/Documents/KRN_GWAS_v3/GWAS3_proj/cache")
  bsle <- load.gwaa.data(phe = "bsle_pheno.dat", gen = "genos_bsle10.raw",force = T)
  descriptives.trait(bsle)
  qc1 <- check.marker(bsle, p.level = 0, het.fdr = 0.001, ibs.threshold = 0.99,
                      callrate=0.4, maf=0.05, perid.call=0.5)
  # In total, 125 (100%) markers passed all criteria
  # In total, 196 (100%) people passed all criteria
  sub1 <- bsle[qc1$idok, qc1$snpok];
  res1 <- qtscore(KRN, data=sub1, trait.type = "binomial")
  #http://www.genabel.org/sites/default/files/pdfs/GenABEL-tutorial.pdf
  pval1 <- results(res1)
  source("~/Documents/Rcodes/get_qval.R")
  pval1$qval <-  get_qval(pval=pval1, pcol="P1df", method="fdr")
  
  ### estimate h2
  library("hglm")
  gkin <- ibs(sub1,w="freq")
  h2dm <- polygenic_hglm(KRN, kin=gkin, data=sub1, trait="binomial")
  message("h2: " ,h2dm$esth2)
  
  return(pval1)
}


##### BSLE simulation
source("~/Documents/KRN_GWAS_v3/GWAS3_proj/lib/BSLE_simulation.R")
raf <- read.csv("~/Documents/GWAS2_KRN/GWAS2_proj/cache/RAF_BSLE.csv")
raf <- merge(raf, sig[, c("snpid", "qval")], by="snpid")
simres <- simdrift(raf=raf, reps=10000, pop1=4000, pop2=300, cycle=30)

### conducted GWAS with FDR control for binomial and quantitative traits
BSLE_pval <- BSLE_gwas()

source("~/Documents/Rcodes/save.append.R")
save.append(list="BSLE_pval", 
            file="~/Documents/KRN_GWAS_v3/GWAS3_proj/cache/cv_pvals.RData",
            description=c("without simulation test"))

#####################
ob <- load("~/Documents/KRN_GWAS_v3/GWAS3_proj/cache/cv_pvals.RData")


