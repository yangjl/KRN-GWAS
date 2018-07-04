# Jinliang Yang
# July 14th, 2014

bins <- read.csv("~/Documents/KRN_GWAS_v3/GWAS3_proj/data/S.table5.KAV.bins.csv")
names(bins) <- c("snpid", "chr", "pos", "value", "method", "bin")
dim(bins)
#[1] 986   6

bins$bin <- paste(bins$chr, round(bins$pos/100000, 0), sep="_")
length(unique(bins$bin))
#764

table(bins$method)
#bayes single-variant       stepwise 
#442            257            300 
length(unique(subset(bins, method=="single-variant")$bin))
# 192
length(unique(subset(bins, method=="stepwise")$bin))
# 296
length(unique(subset(bins, method=="Bayesian-based")$bin))
# 343


############################
kav231 <- read.csv("~/Documents/KRN_GWAS_v3/GWAS3_proj/data/S.table7.prioritized_KAVs.csv")
names(kav231) <- c("order", "snpid", "chr", "pos", "value", "method", "bin", "effect", "h2")
length(unique(kav231$snpid)) #231
length(unique(kav231$bin)) #125

sel <- subset(bins, snpid %in% kav231$snpid)
length(unique(sel$snpid))

subset(kav231, !(snpid %in% sel$snpid))


##############
report_overlap <- function(bins=bins){
  
  message(sprintf("total variants [ %s ] and bins [ %s ]", 
                  length(unique(bins$snpid)), length(unique(bins$bin))))
  
  bin1 <- bin2 <- bin3 <- bin1a <- bin1b <- bin1c <- bin2a <- bin2b <- bin2c <- data.frame()
  for(bini in unique(bins$bin)){
    mybin <- subset(bins, bin %in% bini)
    if(length(unique(mybin$method)) == 3){
      bin3 <- rbind(bin3, mybin)
    }
    if(length(unique(mybin$method)) == 2){
      bin2 <- rbind(bin2, mybin)
      if(sum(unique(mybin$method) %in% c("single-variant", "Bayesian-based")) == 2 ){
        bin2a <- rbind(bin2a, mybin)
      }else if(sum(unique(mybin$method) %in% c("single-variant", "stepwise")) == 2){
        bin2b <- rbind(bin2b, mybin)
      }else if(sum(unique(mybin$method) %in% c("Bayesian-based", "stepwise")) == 2){
        bin2c <- rbind(bin2c, mybin)
      }
    }
    if(length(unique(mybin$method)) == 1){
      bin1 <- rbind(bin1, mybin)
      if(unique(mybin$method) == "single-variant"){
        bin1a <- rbind(bin1a, mybin)
      }else if(unique(mybin$method) == "Bayesian-based"){
        bin1b <- rbind(bin1b, mybin)
      }else if(unique(mybin$method) == "stepwise"){
        bin1c <- rbind(bin1c, mybin)
      }
    }   
  }
  message(sprintf("Detected by three: [ %s ] variants and bins [ %s ]", 
                  length(unique(bin3$snpid)), length(unique(bin3$bin))))
  message(sprintf("Detected by two: [ %s ] variants and bins [ %s ]", 
                  length(unique(bin2$snpid)), length(unique(bin2$bin))))
  message(sprintf("Detected by single-variant and bayes: [ %s ] variants and bins [ %s ]", 
                  length(unique(bin2a$snpid)), length(unique(bin2a$bin))))
  message(sprintf("Detected by single-variant and stepwise: [ %s ] variants and bins [ %s ]", 
                  length(unique(bin2b$snpid)), length(unique(bin2b$bin))))
  message(sprintf("Detected by bayes and stepwise: [ %s ] variants and bins [ %s ]", 
                  length(unique(bin2c$snpid)), length(unique(bin2c$bin))))
  message(sprintf("Detected by one: [ %s ] variants and bins [ %s ]", 
                  length(unique(bin1$snpid)), length(unique(bin1$bin))))
  message(sprintf("Detected by single-variant ONLY: [ %s ] variants and bins [ %s ]", 
                  length(unique(bin1a$snpid)), length(unique(bin1a$bin))))
  message(sprintf("Detected by bayes ONLY: [ %s ] variants and bins [ %s ]", 
                  length(unique(bin1b$snpid)), length(unique(bin1b$bin))))
  message(sprintf("Detected by stepwise ONLY: [ %s ] variants and bins [ %s ]", 
                  length(unique(bin1c$snpid)), length(unique(bin1c$bin))))
  
  a1 <- subset(bins, method=="single-variant");
  a2 <- subset(bins, method=="Bayesian-based")
  a3 <- subset(bins, method=="stepwise")
  message(sprintf("Detected by single-variant ALL: [ %s ] variants and bins [ %s ]", 
                  length(unique(a1$snpid)), length(unique(a1$bin))))
  message(sprintf("Detected by Bayesian-based ALL: [ %s ] variants and bins [ %s ]", 
                  length(unique(a2$snpid)), length(unique(a2$bin))))
  message(sprintf("Detected by stepwise ALL: [ %s ] variants and bins [ %s ]", 
                  length(unique(a3$snpid)), length(unique(a3$bin))))
  lsbin <- list()
  lsbin[['bin1']] <- bin1
  lsbin[['bin2']] <- bin2
  lsbin[['bin3']] <- bin3
  return(lsbin)
}

########
res <- report_overlap(bins=bins)
