## Jinliang Yang
## June 21, 2014

########
KAV_format <- function(kavfile="~/Documents/KRN_GWAS_v3/GWAS3_proj/data/Stable1.csv"){
  kav231 <- read.csv(kavfile);
  names(kav231)[1:4] <- c("snpid", "order", "chr", "pos");
  kav231$bin10 <- paste(kav231$chr, round(kav231$pos/10, 0), sep="_")
  message(sprintf("# of Total KAVs of three methods: [ %s ]", nrow(kav231)))
  
  kav231$bin <- paste(kav231$chr, round(kav231$pos/100000,0), sep="_")
  khead <- kav231[, 1:4]
  khead <- khead[!duplicated(khead$snpid),]
  khead <- khead[order(khead$chr, khead$pos),]
  khead$order <- 1:nrow(khead)
  message(sprintf("# of unique KAVs: [ %s ]", nrow(khead)))
  
  kav <- merge(khead, kav231[, -2:-4], by="snpid")
  return(kav)
}

kav231 <- KAV_format(kavfile="~/Documents/KRN_GWAS_v3/GWAS3_proj/data/Stable1.csv")
#### for publication
write.table(kav231, "~/Documents/KRN_GWAS_v3/GWAS3_proj/reports/S.table.kav231.csv",
            sep=",", row.names=FALSE, quote=FALSE)

############# KAS only
filtering <- function(ref=kav231, qval=elite1, cutoff=.05, BSLE=TRUE){
  ##############
  #sig1 <- subset(qval, qval <= cutoff);
  #sig1 <- sig1[!duplicated(sig1$bin10),]
  ##############
  #ref[ref$h2 < h2cutoff,]$effect <- 0;
  message(sprintf("qval: unique bin10 | unique KAVs: [ %s | %s ]", 
                  length(unique(qval$bin10)), length(unique(qval$SNP)) ))
  mkav <- merge(ref, qval, by="bin10");
  message(sprintf("matached with KAV231: unique bin10 | unique KAVs: [ %s | %s ]", 
                  length(unique(mkav$bin10)), length(unique(mkav$snpid)) ))
  if(BSLE){
    f2 <- subset(mkav, qval <= cutoff & dir > 0 & sim_qval >0.05);
  }else{
    f2 <- subset(mkav, qval <= cutoff & dir > 0);
  }
  
  message(sprintf("Significant of KAV231: unique bin10 | unique KAVs: [ %s | %s ]", 
                  length(unique(f2$bin10)), length(unique(f2$snpid)) ))
  
}
#############
setwd("~/Documents/KRN_GWAS_v3/GWAS3_proj/data/")
elite1 <- read.csv("Elite_pval1.csv")
elite2 <- read.csv("Elite_pval2.csv")
#Binomial XBSA
xbsa2 <- read.csv("XBSA_binary.csv")
bsle2 <- read.csv("BSLE_binary.csv")

pval1.1 <- filtering(ref=kav231, qval=elite1, cutoff=.05, BSLE=FALSE)
pval1.2 <- filtering(ref=kav231, qval=elite2, cutoff=.05, BSLE=FALSE)
pval2 <- filtering(ref=kav231, qval=xbsa2, cutoff=.05, BSLE=FALSE)

bsle2$sim_qval <- p.adjust(bsle2$sim_pval, method="fdr")
#### further filtering
pval3.1 <- filtering(ref=kav231, qval=bsle2, cutoff=.05, BSLE=FALSE)
pval3.2 <- filtering(ref=kav231, qval=bsle2, cutoff=.05, BSLE=TRUE)


####################################################################################
kav_cv_sum <- function(kav=kav231, elite1=elite1, elite2=elite2, xbsa=xbsa2, bsle=bsle2){
  kavres <- kav;
  kavres <- kavres[!duplicated(kavres$snpid),]
  #note: -9 missing, 0=genotyped, 1=cross-validated
  kavres$elite1 <- kavres$elite2 <- kavres$xbsa <- kavres$bsle <- -9
  kavres[kavres$bin10 %in% elite1$bin10, ]$elite1 <- 0
  try(kavres[kavres$bin10 %in% subset(elite1, qval< 0.05 & dir > 0)$bin10, ]$elite1 <- 1)
  kavres[kavres$bin10 %in% elite2$bin10, ]$elite2 <- 0
  try(kavres[kavres$bin10 %in% subset(elite2, qval< 0.05 & dir >0)$bin10, ]$elite2 <- 1)
  kavres[kavres$bin10 %in% xbsa$bin10, ]$xbsa <- 0
  try(kavres[kavres$bin10 %in% subset(xbsa, qval< 0.05 & dir >0)$bin10, ]$xbsa <- 1)
  kavres[kavres$bin10 %in% bsle$bin10, ]$bsle <- 0
  try(kavres[kavres$bin10 %in% subset(bsle, qval< 0.05 & dir >0 & sim_qval > 0.05)$bin10, ]$bsle <- 1)
  
  ### summary for results output
  kavres$cvd <- kavres$genotyped <- -9
  ### sucessfully genotyped IDs
  gid <- subset(kavres, bsle>=0 | xbsa>=0 | elite1>=0 | elite2 >=0)$snpid
  ### sucessfully cv IDs
  cvid <- subset(kavres, bsle==1 | xbsa==1 | elite1 ==1 | elite2==1)$snpid
  kavres[kavres$snpid %in% gid, ]$genotyped <- 1
  kavres[kavres$snpid %in% cvid, ]$cvd <- 1;
  cvbin10 <- subset(kavres, bsle==1 | xbsa==1 | elite1 ==1 | elite2==1)$bin10
  message(sprintf("KAV231: genotyped | cross-valdiated | bin10: [%s | %s | %s] ", 
                  length(unique(gid)), length(unique(cvid)),length(unique(cvbin10)) ))
  return(kavres)
}

kavsum <- kav_cv_sum(kav=kav231, elite1=elite1, elite2=elite2, xbsa=xbsa2, bsle=bsle2)
write.table(kavsum, "~/Documents/KRN_GWAS_v3/GWAS3_proj/reports/S.table_kavsum.csv",
            sep=",", row.names=FALSE, quote=FALSE)




