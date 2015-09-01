### Jinliang Yang
### June 22, 2014
### create a master table of cross-validation results

kav231 <- read.csv("~/Documents/KRN_GWAS_v3/GWAS3_proj/reports/S.table_kavsum.csv")
plos <- read.csv("~/Documents/KRN_GWAS_v3/GWAS3_proj/reports/S.table_plos_sum.csv")
names(plos)[1:3] <- c("snpid", "chr", "pos")
ck <- read.csv("~/Documents/KRN_GWAS_v3/GWAS3_proj/reports/Stable_cksum.csv")

kav231$type <- "KAV"
plos$type <- "PLoS KAV"
ck$type <- "Control"

nms <- c("type", "snpid", "chr", "pos", "elite1", "elite2", "xbsa", "bsle", "genotyped","cvd")
allkav <- rbind(kav231[, nms], plos[, nms], ck[, nms])
genod <- subset(allkav, genotyped==1)
message(sprintf("Genotyped SNPs: ", length(unique(genod$snpid))))

######
## KAVs:
kav_elite1 <- read.csv("~/Documents/KRN_GWAS_v3/GWAS3_proj/data/Elite_pval1.csv")
kav_elite2 <- read.csv("~/Documents/KRN_GWAS_v3/GWAS3_proj/data/Elite_pval2.csv")
kav_xbsa <- read.csv("~/Documents/KRN_GWAS_v3/GWAS3_proj/data/XBSA_binary.csv")
kav_bsle <- read.csv("~/Documents/KRN_GWAS_v3/GWAS3_proj/data/BSLE_binary.csv")

## Control:
ck_elite1 <- read.csv("~/Documents/KRN_GWAS_v3/GWAS3_proj/reports/Stable_ck_elite1.csv")
ck_elite2 <- read.csv("~/Documents/KRN_GWAS_v3/GWAS3_proj/reports/Stable_ck_elite2.csv")
ck_xbsa <- read.csv("~/Documents/KRN_GWAS_v3/GWAS3_proj/reports/Stable_ck_xbsa.csv")
ck_bsle <- read.csv("~/Documents/KRN_GWAS_v3/GWAS3_proj/reports/Stable_ck_bsle.csv")

## PLoS:
plos_elite1 <- read.csv("~/Documents/KRN_GWAS_v3/GWAS3_proj/reports/S.table_plos_e1.csv")
plos_elite2 <- read.csv("~/Documents/KRN_GWAS_v3/GWAS3_proj/reports/S.table_plos_e2.csv")
plos_xbsa <- read.csv("~/Documents/KRN_GWAS_v3/GWAS3_proj/reports/S.table_plos_xbsa.csv")
plos_bsle <- read.csv("~/Documents/KRN_GWAS_v3/GWAS3_proj/reports/S.table_plos_bsle.csv")

###########
snp_rbind <- function(nms=c("snpid", "qval", "dir"), kav=kav_elite1, plos=plos_elite1, ck=ck_elite1){
  
  kav <- kav[, c("SNP", "qval", "dir")]
  names(kav) <- nms
  #kav$type <- "KAV";
  
  ck$dir <- 0; 
  ck <- ck[, nms]
  #ck$type <- "CK"
  
  plos<- plos[, c("plossnpid", "qval", "dir")]
  #plos$type <- "plos"
  names(plos) <- nms
  snpbind <- rbind(kav, ck, plos)
  snpbind <- snpbind[!duplicated(snpbind$snpid),];
  
  message(sprintf("Total SNPs combined: %s", length(unique(snpbind$snpid))))
  message(sprintf("KAV: [%s], CK: [%s], PLoS: [%s]", length(unique(kav$snpid)),
                  length(unique(ck$snpid)), length(unique(plos$snpid))) )
  return(snpbind)
}
############
elite1 <- snp_rbind(nms=c("snpid", "qval", "dir"), kav=kav_elite1, plos=plos_elite1, ck=ck_elite1)
elite2 <- snp_rbind(nms=c("snpid", "qval", "dir"), kav=kav_elite2, plos=plos_elite2, ck=ck_elite2)
xbsa <- snp_rbind(nms=c("snpid", "qval", "dir"), kav=kav_xbsa, plos=plos_xbsa, ck=ck_xbsa)
bsle <- snp_rbind(nms=c("snpid", "qval", "dir"), kav=kav_bsle, plos=plos_bsle, ck=ck_bsle)

###########
### 120
names(elite1) <- c("snpid", "elite1_qval", "elite1_doe")
cvd <- merge(genod, elite1, by="snpid", all.x=TRUE)

names(elite2) <- c("snpid", "elite2_qval", "elite2_doe")
cvd <- merge(cvd, elite2, by="snpid", all.x=TRUE)

names(xbsa) <- c("snpid", "xbsa_qval", "xbsa_doe")
cvd <- merge(cvd, xbsa, by="snpid", all.x=TRUE)

names(bsle) <- c("snpid", "bsle_qval", "bsle_doe")
cvd <- merge(cvd, bsle, by="snpid", all.x=TRUE)

write.table(cvd, "~/Documents/KRN_GWAS_v3/GWAS3_proj/reports/Stable13_cv_summary.csv",
            sep=",", row.names=FALSE, quote=FALSE)
