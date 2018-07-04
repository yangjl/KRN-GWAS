### Jinliang Yang


stable5 <- read.csv("~/Documents/KRN_GWAS_v3/GWAS3_proj/reports/S.table5.KAV.bins.csv")
dim(stable5)
head(stable5)
names(stable5) <- c("snpid", "chr", "pos", "value", "method", "bin")
length(unique(stable5$bin))


#####
stable9 <- read.csv("~/Documents/KRN_GWAS_v3/GWAS3_proj/reports/S.table9.cv.samples.csv")
dim(stable9)
table(stable9$Population)

###
Stable13 <- read.csv("~/Documents/KRN_GWAS_v3/GWAS3_proj/reports/Stable13_cv_summary.csv")
head(Stable13)
dim(Stable13) #120
length(unique(Stable13$snpid))

elite1 <- subset(Stable13, type=="KAV" & elite1 >= 0)
dim(elite1)
dim(subset(elite1, elite1>0))

elite2 <- subset(Stable13, type=="KAV" & elite2 >= 0)
dim(elite2)
dim(subset(elite2, elite2>0))

xbsa <- subset(Stable13, type=="KAV" & xbsa >= 0)
dim(xbsa)
dim(subset(xbsa, xbsa>0))

bsle <- subset(Stable13, type=="KAV" & bsle >= 0)
dim(bsle)
dim(subset(bsle, bsle>0))

kav <- subset(Stable13, type=="KAV")
dim(subset(kav, elite1 >= 0 | elite2 >= 0 | xbsa >=0 | bsle >=0))

dim(subset(kav, elite1 >= 1 | elite2 >= 1 | xbsa >=1 | bsle >=1))

ck <- subset(Stable13, type=="Control")
ck

##############
plos <- subset(Stable13, type=="PLoS KAV")
dim(subset(plos, elite1 >= 0 | elite2 >= 0 | xbsa >=0 | bsle >=0))
dim(subset(plos, elite1 >= 1 | elite2 >= 1 | xbsa >=1 | bsle >=1))
