### Jinliang Yang
### plosv2

####################################################################################
plos_cv <- function(kav=kav231, elite1=elite1, elite2=elite2, xbsa=xbsa2, bsle=bsle2){
  kavres <- kav;
  kavres <- kavres[!duplicated(kavres$plossnpid),]
  #note: -9 missing, 0=genotyped, 1=cross-validated
  kavres$elite1 <- kavres$elite2 <- kavres$xbsa <- kavres$bsle <- -9
  kavres[kavres$bin10 %in% elite1$bin10, ]$elite1 <- 0
  try(kavres[kavres$bin10 %in% subset(elite1, qval< 0.05 & dir > 0)$bin10, ]$elite1 <- 1)
  kavres[kavres$bin10 %in% elite2$bin10, ]$elite2 <- 0
  try(kavres[kavres$bin10 %in% subset(elite2, qval< 0.05 & dir >0)$bin10, ]$elite2 <- 1)
  kavres[kavres$bin10 %in% xbsa$bin10, ]$xbsa <- 0
  try(kavres[kavres$bin10 %in% subset(xbsa, qval< 0.05 & dir >0)$bin10, ]$xbsa <- 1)
  kavres[kavres$bin10 %in% bsle$bin10, ]$bsle <- 0
  try(kavres[kavres$bin10 %in% subset(bsle, qval< 0.05 & dir >0)$bin10, ]$bsle <- 1)
  
  ### summary for results output
  kavres$cvd <- kavres$genotyped <- -9
  ### sucessfully genotyped IDs
  gid <- subset(kavres, bsle>=0 | xbsa>=0 | elite1>=0 | elite2 >=0)$plossnpid
  ### sucessfully cv IDs
  cvid <- subset(kavres, bsle==1 | xbsa==1 | elite1 ==1 | elite2==1)$plossnpid
  kavres[kavres$plossnpid %in% gid, ]$genotyped <- 1
  kavres[kavres$plossnpid %in% cvid, ]$cvd <- 1;
  cvbin10 <- subset(kavres, bsle==1 | xbsa==1 | elite1 ==1 | elite2==1)$bin10
  message(sprintf("PLoS: genotyped | cross-valdiated | bin10: [%s | %s | %s] ", 
                  length(unique(gid)), length(unique(cvid)),length(unique(cvbin10)) ))
  return(kavres)
}


########################################################################################
plos_cv_sum <- function(){
  ob <- load("~/Documents/KRN_GWAS_v3/GWAS3_proj/cache/cv_pvals.RData")
  kav231<- read.csv("~/Documents/KRN_GWAS_v3/GWAS3_proj/reports/S.table.kav231.csv")
  plosv2 <- read.csv("~/Documents/KRN_GWAS_v3/GWAS3_proj/cache/plosv2_KAV.csv")
  
  ### small DSF function:
  pval2dsf <- function(pval=elite_pval){
    names(pval)[1:2] <- c("chr", "pos");
    pval$snpid <- paste(pval$chr, pval$pos, sep="_");
    pval$bin10 <- paste(pval$chr, round(pval$pos/10, 0), sep="_")
    return(pval)
  }
  plosv2$bin10 <- paste(plosv2$chr, round(plosv2$pos/10, 0), sep="_")
  elite_pval <- pval2dsf(pval=elite_pval)
  xbsa_pval1 <- pval2dsf(pval=xbsa_pval1)
  bsle_pval <- pval2dsf(pval=BSLE_pval)
  
  names(plosv2)[1:3] <- c("plossnpid", "ploschr", "plospos")
  plosv2 <- plosv2[, c("plossnpid", "ploschr", "plospos", "effect", "p.value", "bin10")]
  ### excluding KAV231 from following pval tables 
  e1 <- elite_pval[!(elite_pval$bin10 %in% kav231$bin10),]
  e1 <- merge(e1, plosv2, by="bin10")
  e1$dir <- -e1$effect * e1$effB
  e2 <- e1
  e2$qval <- p.adjust(e1$qval, method="fdr")
  ### use quantitative trait!
  x1 <- xbsa_pval1[!(xbsa_pval1$bin10 %in% kav231$bin10),]
  x1 <- merge(x1, plosv2, by="bin10")
  x1$dir <- -x1$effect * x1$effB
  b1 <- bsle_pval[!(bsle_pval$bin10 %in% kav231$bin10),]
  b1 <- merge(b1, plosv2, by="bin10")
  b1$dir <- -b1$effect * b1$effB
  
  ploscv <- plos_cv(kav=plosv2, elite1=e1, elite2=e2, xbsa=x1, bsle=b1)
  ### summary for results output
  resls <- list()
  resls[[1]] <- ploscv;
  resls[["e1"]] <- e1;
  resls[["e2"]] <- e2;
  resls[["xbsa"]] <- x1;
  resls[["bsle"]] <- b1;
  return(resls)
}


##################
pls <- plos_cv_sum()

write.table(pls[[1]], "~/Documents/KRN_GWAS_v3/GWAS3_proj/reports/S.table_plos_sum.csv",
            sep=",", row.names=FALSE, quote=FALSE)
write.table(pls[['e1']], "~/Documents/KRN_GWAS_v3/GWAS3_proj/reports/S.table_plos_e1.csv",
            sep=",", row.names=FALSE, quote=FALSE)
write.table(pls[['e2']], "~/Documents/KRN_GWAS_v3/GWAS3_proj/reports/S.table_plos_e2.csv",
            sep=",", row.names=FALSE, quote=FALSE)
write.table(pls[['xbsa']], "~/Documents/KRN_GWAS_v3/GWAS3_proj/reports/S.table_plos_xbsa.csv",
            sep=",", row.names=FALSE, quote=FALSE)
write.table(pls[['bsle']], "~/Documents/KRN_GWAS_v3/GWAS3_proj/reports/S.table_plos_bsle.csv",
            sep=",", row.names=FALSE, quote=FALSE)




M <- as.table(rbind(c(21, 33), c(14, 31)))
dimnames(M) <- list(type = c("MyKAV","Plos"),
                    val = c("Genotyped","CVD"))
(Xsq <- chisq.test(M))


