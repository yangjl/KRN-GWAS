### Jinliang Yang
### plosv2

ck_cv_sum <- function(){
  ob <- load("~/Documents/KRN_GWAS_v3/GWAS3_proj/cache/cv_pvals.RData")
  ck <- read.csv("~/Documents/KRN_GWAS_v3/GWAS3_proj/data/KAS_list_control.csv")
  
  ### small DSF function:
  pval2dsf <- function(pval=elite_pval){
    names(pval)[1:2] <- c("chr", "pos");
    pval$snpid <- paste(pval$chr, pval$pos, sep="_");
    pval$bin10 <- paste(pval$chr, round(pval$pos/10, 0), sep="_")
    return(pval)
  }
  ck$bin10 <- paste(ck$chr, round(ck$pos/10, 0), sep="_")
  elite_pval <- pval2dsf(pval=elite_pval)
  xbsa_pval1 <- pval2dsf(pval=xbsa_pval1)
  bsle_pval <- pval2dsf(pval=bsle_pval)
  
  ### excluding KAV231 from following pval tables 
  e1 <- elite_pval[elite_pval$bin10 %in% ck$bin10,]
  x1 <- xbsa_pval1[xbsa_pval1$bin10 %in% ck$bin10,]
  b1 <- bsle_pval[bsle_pval$bin10 %in% ck$bin10,]
  
  ckres <- ck;
  ckres$elite1 <- ckres$elite2 <- ckres$xbsa <- ckres$bsle <- -9
  ckres[ckres$snpid %in% e1$snpid, ]$elite1 <- 0
  try(ckres[ckres$snpid %in% subset(e1, qval< 0.05)$snpid, ]$elite1 <- 1)
  ckres[ckres$snpid %in% e1$snpid, ]$elite2 <- 0
  try(ckres[ckres$snpid %in% subset(e1, qval< 0.05)$snpid, ]$elite2 <- 1)
  ckres[ckres$snpid %in% x1$snpid, ]$xbsa <- 0
  try(ckres[ckres$snpid %in% subset(x1, qval< 0.05)$snpid, ]$xbsa <- 1)
  ckres[ckres$snpid %in% b1$snpid, ]$bsle <- 0
  try(ckres[ckres$snpid %in% subset(b1, qval< 0.05)$snpid, ]$bsle <- 1)
  
  ### summary for results output
  ckres$cvd <- ckres$genotyped <- -9
  ### sucessfully genotyped IDs
  gid <- subset(ckres, bsle>=0 | xbsa>=0 | elite1>=0 | elite2 >=0)$snpid
  ### sucessfully cv IDs
  cvid <- subset(ckres, bsle==1 | xbsa==1 | elite1 ==1 | elite2==1)$snpid
  ckres[ckres$snpid %in% gid, ]$genotyped <- 1
  ckres[ckres$snpid %in% cvid, ]$cvd <- 1;
  cvbin10 <- subset(ckres, bsle==1 | xbsa==1 | elite1 ==1 | elite2==1)$bin10
  
  ### summary for results output
  resls <- list()
  resls[[1]] <- ckres;
  resls[["e1"]] <- e1;
  resls[["e2"]] <- e1;
  resls[["xbsa"]] <- x1;
  resls[["bsle"]] <- b1;
  return(resls)
}

##############################
ckres <- ck_cv_sum()
write.table(ckres[[1]], "~/Documents/KRN_GWAS_v3/GWAS3_proj/reports/Stable_cksum.csv", 
            sep=",", row.names=FALSE, quote=FALSE)
write.table(ckres[['e1']], "~/Documents/KRN_GWAS_v3/GWAS3_proj/reports/Stable_ck_elite1.csv", 
            sep=",", row.names=FALSE, quote=FALSE)
write.table(ckres[['e2']], "~/Documents/KRN_GWAS_v3/GWAS3_proj/reports/Stable_ck_elite2.csv", 
            sep=",", row.names=FALSE, quote=FALSE)
write.table(ckres[['xbsa']], "~/Documents/KRN_GWAS_v3/GWAS3_proj/reports/Stable_ck_xbsa.csv", 
            sep=",", row.names=FALSE, quote=FALSE)
write.table(ckres[['bsle']], "~/Documents/KRN_GWAS_v3/GWAS3_proj/reports/Stable_ck_bsle.csv", 
            sep=",", row.names=FALSE, quote=FALSE)


##################################################################
ck_cv_sum <- function(){
  ob <- load("~/Documents/KRN_GWAS_v3/GWAS3_proj/cache/cv_pvals.RData")
  kas <- read.csv("~/Documents/KRN_GWAS_v3/GWAS3_proj/reports/S.table.kav231.csv")
  plos <- read.csv("~/Documents/KRN_GWAS_v3/GWAS3_proj/cache/plosv2_KAV.csv")
  
  ### small DSF function:
  pval2dsf <- function(pval=elite_pval){
    names(pval)[1:2] <- c("chr", "pos");
    pval$snpid <- paste(pval$chr, pval$pos, sep="_");
    pval$bin10 <- paste(pval$chr, round(pval$pos/10, 0), sep="_")
    return(pval)
  }
  plos$bin10 <- paste(plos$chr, round(plos$pos/10, 0), sep="_")
  elite_pval <- pval2dsf(pval=elite_pval)
  xbsa_pval1 <- pval2dsf(pval=xbsa_pval1)
  bsle_pval <- pval2dsf(pval=bsle_pval)
  
  ### excluding KAV231 from following pval tables 
  e1 <- elite_pval[!(elite_pval$bin10 %in% c(plos$bin10, kas$bin10)),]
  x1 <- xbsa_pval1[!(xbsa_pval1$bin10 %in% c(plos$bin10, kas$bin10)),]
  b1 <- bsle_pval[!(bsle_pval$bin10 %in% c(plos$bin10, kas$bin10)),]
  
  ckres <- data.frame(snpid=unique(c(e1$snpid, x1$snpid, b1$snpid)), type="ck")
  ckres$elite1 <- ckres$elite2 <- ckres$xbsa <- ckres$bsle <- -9
  ckres[ckres$snpid %in% e1$snpid, ]$elite1 <- 0
  try(ckres[ckres$snpid %in% subset(e1, qval< 0.05)$snpid, ]$elite1 <- 1)
  ckres[ckres$snpid %in% e1$snpid, ]$elite2 <- 0
  try(ckres[ckres$snpid %in% subset(e1, qval< 0.05)$snpid, ]$elite2 <- 1)
  ckres[ckres$snpid %in% x1$snpid, ]$xbsa <- 0
  try(ckres[ckres$snpid %in% subset(x1, qval< 0.05)$snpid, ]$xbsa <- 1)
  ckres[ckres$snpid %in% b1$snpid, ]$bsle <- 0
  try(ckres[ckres$snpid %in% subset(b1, qval< 0.05)$snpid, ]$bsle <- 1)
}


