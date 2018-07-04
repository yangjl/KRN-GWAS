### Jinliang Yang
### 6/20/2014
### get PLoS KAV version2


#chr2ped <- function(chr=1){

get_plosv2 <- function(){
  ### hmp1
  setwd("~/DBcenter/VariationDB/HapMap1/")
  tab5rows <- read.delim("maizeHapMapV1_B73RefGenV2_20110309_ALL.hmp.txt", header=T, nrow=5)
  classes <- sapply(tab5rows, class)
  
  hmpraw <- read.delim("maizeHapMapV1_B73RefGenV2_20110309_ALL.hmp.txt", header=TRUE, colClasses = classes)
  nrow(hmp1); #[1] 1624651
  hmp1<- hmpraw[, c(1:5, 11)]
  
  ###################
  plos <- read.csv("~/Documents/KRN_GWAS_v3/validation/pgen_kav.csv")
  plos <- subset(plos, trait=="ERN")[, 1:7]
  plos$pzeid <- paste("PZE", sprintf('%02.0f', plos$chrom), sprintf('%08.0f', plos$bp), sep="")
  plosv2 <- merge(plos[, c(1,3,6:8)], hmp1, by.x="pzeid", by.y="QCcode")
  
  #### DSF: density SNP format:
  plosv2$snpid <- paste(plosv2$chrom, plosv2$pos, sep="_")
  plosv2 <- plosv2[, c("snpid", "chrom", "pos", "alleles", "strand", "trait", "cM", "effect", 
                       "p.value", "rs.", "pzeid")]
  names(plosv2) <- c("snpid", "chr", "pos", "alleles", "strand", "trait", "cM", "effect", 
                     "p.value", "rs", "pzeid_v1")
  message(sprintf("merged snpid / input snpid: [%s / %s]", nrow(plos), nrow(plosv2)))
  return(plosv2)
}

plosv2 <- get_plosv2()
write.table(plosv2, "~/Documents/KRN_GWAS_v3/GWAS3_proj/cache/plosv2_KAV.csv",
            sep=",", row.names=FALSE, quote=FALSE)


