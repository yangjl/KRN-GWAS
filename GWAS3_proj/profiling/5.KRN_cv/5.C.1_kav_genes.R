### Jinliang Yang
### 6/22/2014
### find KAV-linked genes


kav231 <- read.csv("~/Documents/KRN_GWAS_v3/GWAS3_proj/reports/S.table.kav231.csv")


fgsv2 <- read.delim("~/db/AGPv2/ZmB73_5b_FGS.gff", header=FALSE)
names(fgsv2) <- c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute" )
gene <- subset(fgsv2, feature == "gene")
### 39656

get_kav_linked_genes <- function(kav=kav231, fgs=gene, bin=500000){
  kav <- kav[!duplicated(kav$snpid),]
  
  res <- data.frame()
  for(i in 1:nrow(kav)){
    tem <- subset(fgs, chr == kav$chr[i] & start >= kav$pos[i]-bin & end <= kav$pos[i]+bin);
    tem$KAV <- kav$snpid[i]
    res <- rbind(res, tem)
    print(i)
  }
  #res <- res[!duplicated(res$attribute),]
  message(sprintf("[ %s ] of KAV-linked genes", nrow(res)))
  return(res)
}

kavgenes <- get_kav_linked_genes(kav=kav231, fgs=gene, bin=500000)
#[2690] of KAV-linked genes
kavgenes$attribute <- gsub("ID=", "", kavgenes$attribute)
kavgenes$attribute <- gsub(";Name=.*", "", kavgenes$attribute)

get_KRN_genes <- function(pid=50, pcov=50 ){
  setwd("~/Documents/KRN_GWAS_v3/validation/")
  auxin <- read.csv("Auxin_KEGG_CoGe.csv");
  auxin$set <- "auxin"
  cyto <- read.csv("Cytokinin_KEGG_CoGe.csv");
  cyto$set <- "cytokinin"
  other <- read.csv("Classical_gene2_CoGe.csv");
  other$set <- "other"
  kg <- rbind(auxin, cyto, other)
  
  kg$Coverage <- as.character(as.numeric(gsub("%", "", kg$Coverage)));
  kg$PercID <- as.character(as.numeric(gsub("%", "", kg$PercID)))
  kg <- subset(kg, PercID>pid & Coverage>pcov)
  kg <- kg[!duplicated(kg$ClosestFGS), ]
  return(kg)
}

kg <- get_KRN_genes(pid=50, pcov=50)
krngenes <- merge(kg, kavgenes, by.x="ClosestFGS", by.y="attribute")
write.table(krngenes, "~/Documents/KRN_GWAS_v3/GWAS3_proj/reports/S.table14_krngenes.csv",
            row.names=FALSE, quote=FALSE)
