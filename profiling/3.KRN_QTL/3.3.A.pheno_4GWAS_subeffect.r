# Jinliang Yang
# last update: 2.27.2011
# 3.14.2011
# purpose: subtract the QTL effect by chromosome:

# read.me
# plot the QTL effects heat map of CIM results
# Summary the QTL from CIM results
# geno code: IBM AA 0 Mo17; aa 1 B73
#			 NAM AA 2 nonB73; aa 0 B73
# QTL Additive Effect: should be 2x <<<<<---- in QTL summary, I transfer to -2x
# LOD = 0.217LR


# set the working directory
setwd("/home/yangjl/Research/QTL/QTLcart/NAMKRN_joint")


###############################################################################
#                      MAP File                                               #
###############################################################################
map <- read.delim("~/Research/QTL/QTLcart/NAM_map_20080419.Eddy.imputed.txt", header=TRUE)
dim(map)
#1106 2
map <- map[order(map$ch, map$genetpos, map$imputepos),]



###############################################################################
#                      Geno and Pheno File                                    #
###############################################################################

tagSNP <- read.table("markergenotypes062508.txt", comment.char="#", header=TRUE)
dim(tagSNP)
#[1] 4699 1111
tagSNP <-tagSNP[,c(-2,-3,-4,-5)]
names(tagSNP) <- c("Geno", as.character(map$marker))
dim(tagSNP)
#[1] 4699 1107



trait <- read.table("~/Research/QTL/QTLcart/nam_pheno_krn-p_040211.txt", header=TRUE)
dim(trait)

cross <- merge(trait, tagSNP, by.x="Genotype", by.y="Geno")
dim(cross)
#[1] 4142 1121 #1:15 col are traits, 16:1121 col are geno
###############################################################################
#                      CIM6 QTL Effect                                        #
###############################################################################
# subtracting the other QTL effects and return the interval boundary:

effect <- read.delim("~/Research/QTL/QTLcart/Joint_CIM6_summary.txt", header=TRUE)
map <- map[!is.na(map$pos),]

getpheno <- function(qtl=1){
	
	mycross <- cross
	myeffect <- effect[-qtl,]
	
	for (i in 1:nrow(myeffect)){
		mrk <- as.character(myeffect$SigMarker[i])
		eff <- myeffect$Additive[i]
		
		idx <- which(cross[,mrk]==2)
		mycross[idx,]$wmean <- mycross[idx,]$wmean + eff
	}
#mycross[,2:15]
	out.nm <- paste("pheno_wmean_qtl",qtl, sep="")
#mycross$IID <- gsub("E", "", mycross$IID)
	write.table(mycross[,2:15], out.nm, sep="\t", row.names=FALSE, quote=FALSE, col.names=FALSE)
	
	mychr <- as.numeric(effect[qtl, "Chr"])
	mymap <- map[map$ch== mychr,]
	left <- effect[qtl, "Left"]
	left.idx <- which.min(abs(mymap$genetpos - left))
	left.pos <- mymap[left.idx,]$pos
	
	right <- effect[qtl, "Right"]
	right.idx <- which.min(abs(mymap$genetpos - right))
	right.pos <- mymap[right.idx,]$pos
	
	return(c(mychr, left.pos, right.pos))
	
}


qtl1 <- getpheno(qtl=1)
qtl1
#[1]        1 16530875 26283293
qtl2 <- getpheno(qtl=2)
qtl2

qtl3 <- getpheno(qtl=3)
qtl4 <- getpheno(qtl=4)
qtl5 <- getpheno(qtl=5)
qtl6 <- getpheno(qtl=6)
qtl7 <- getpheno(qtl=7)
qtl8 <- getpheno(qtl=8)
qtl9 <- getpheno(qtl=9)
qtl10 <- getpheno(qtl=10)
qtl11 <- getpheno(qtl=11)
qtl12 <- getpheno(qtl=12)
qtl13 <- getpheno(qtl=13)
qtl14 <- getpheno(qtl=14)
qtl15 <- getpheno(qtl=15)
qtl16 <- getpheno(qtl=16)



summary(lm(wmean~as.factor(PZA00658.21), data=chr1))

summary(lm(wmean~as.factor(PZA00694.6), data=cross))
















