# Data Formatting for QTLcart Use
# start: 02/18/2011
# update: 04/08/2011

# read.me
# plot the QTL effects heat map of CIM results
# Summary the QTL from CIM results
# geno code: IBM AA 0 Mo17; aa 1 B73
#			 NAM AA 2 nonB73; aa 0 B73
# QTL Additive Effect: should be 2x
# LOD = 0.217LR



#################################################################
setwd("/Users/yangjl/Documents/KRN_GWAS/QTL/QTLcart/NAMKRN_sep")

#STEP ONE: get the QTL matrix

trait <- read.table("NAM_name.txt", header=FALSE)
mynote <- as.character(trait$V1)

#------CIM mapping results of wmean
outqtl <- read.table("myB97_krn.eqt", header=FALSE)
names(outqtl) <- c("order", "Chr", "Marker", "Position", "LOD", "Additive", "Dominance", "R2",
"TR2", "S")
outqtl$pop <- "B97"
outqtl$trait <- "wmean"
outqtl <- outqtl[!duplicated(outqtl$Marker),]

for (i in c(2:11, 13:25)){ # 12 Ki11 has no QTL
	b <- mynote[i]
	stem <- paste(b, "krn", sep="_")
	
	mystem.eqt<- paste("my", stem, ".eqt", sep="")
	
	
	qtl <- read.table(mystem.eqt, header=FALSE)
	names(qtl) <- c("order", "Chr", "Marker", "Position", "LOD", "Additive", "Dominance", "R2",
					"TR2", "S")
	qtl$pop <- b
	qtl$trait <- "wmean"
		
	outqtl <- rbind(outqtl, qtl)
	
}

dim(outqtl)
#[1] 181  12


#-------------Mo17
setwd("/Users/yangjl/Documents/KRN_GWAS/QTL/QTLcart/IBM_KRN")

ibm <- read.table("myibm_ril_krn.eqt6", header=FALSE)
names(ibm) <- c("order", "Chr", "Marker", "Position", "LOD", "Additive", "Dominance", "R2",
"TR2", "S")
ibm$pop <- "Mo17"
ibm$trait <- "wmean"


#--------------JOINT
setwd("/Users/yangjl/Documents/KRN_GWAS/QTL/QTLcart/NAMKRN_joint")

joint <- read.table("mynam_krn.eqt6", header=FALSE)
names(joint) <- c("order", "Chr", "Marker", "Position", "LOD", "Additive", "Dominance", "R2",
"TR2", "S")
joint$pop <- "joint"
joint$trait <- "wmean"


outqtl <- rbind(outqtl, ibm)
outqtl$Additive <- 2*outqtl$Additive
dim(outqtl)
#[1] 194  12

############-------------------------------------------------------------#########

setwd("/Users/yangjl/Documents/KRN_GWAS/QTL/QTLcart/")

##################################################################################
#                  plot Joint QTL LOD
##################################################################################
setwd("/Users/yangjl/Documents/KRN_GWAS/QTL/QTLcart/NAMKRN_joint")


cim6 <- read.table("mynam_krn.z6", header=TRUE)
cim6 <- cim6[,1:8]
cim6 <- cim6[, c(1,3,4)]
cim6$position <- 100* as.numeric(cim6$position)
names(cim6)[3] <- "lod"
cim6$c <- as.numeric(cim6$c)
cim6$lod <- 0.217*as.numeric(cim6$lod)


####-------determining the threshold
thred2 <- read.table("mynam_krn.z6e", header=FALSE)
cutoff2 <- quantile(thred2$V2, 0.95)*0.217




############-------------------------------------------------------------#########
##################################################################################
#                  plot the heatmap
##################################################################################
###caculate the actual position
data2 <- cim6

data2$UID <- paste(as.character(data2$c), as.character(data2$position), sep="_")

chr <- data2[, 1] #chr
pos <- newpos <- data2[, 2] #position


tablechr <- table(chr)
allchr <- as.vector(tablechr)

gap <- 20
CMindex <- cumsum(allchr) # cumulatived marker number in each chr
for (i in 1:10) {
	u <- CMindex[i]
	l <- CMindex[i] - allchr[i] + 1
	chr <- l:u
	d <- diff(pos[chr])
	newpos[chr] <- c(gap, d)
}
CM <- cumsum(as.numeric(newpos))    #calculate the actually position

data2$pos <- CM


outqtl$UID <- paste(as.character(outqtl$Chr), as.character(outqtl$Position), sep="_")

data3 <- merge(data2, outqtl, by="UID")

save.image("heatmap.RData")






#####################Plotting the CIM6#######################
####-------start to plot
ylim =c(0, 174)
colors <- rep(c("blue", "green"), 5)
mymhtplot(cim6, colors=colors, usepos = TRUE, gap=20, pch=19, cex=0.5, logscale=FALSE, cutoffs=cutoff2,
main="QTL Position and Effect", ylim=ylim)
axis(2, at=c(0, 10,20,30,40,50,60), labels=c(0, 10,20,30,40,50,60), las=2, cex=0.8)


data3$color <- round(data3$Additive, 1)*10 +14
data3 <- data3[order(data3$color),];
#color = colorRampPalette(c("blue3","darkgoldenrod", "firebrick2"))(26)
color = colorRampPalette(c("blue3","white", "firebrick2"))(26)
#color= heat.colors(21)

for (i in c(1:25)){
  b <- mynote[i]
	suby <- 174 - 4*i + 3
	abline(h= suby, col = "gray60")
	
	if (i!=12){
		subdata <- data3[data3$pop == b,]
		subdata$y <- suby
		order <- subdata$color
		points(as.numeric(subdata$pos), subdata$y, pch=19, col=color[order], cex=1.3)
		axis(2, at = suby, labels = b, las=2, cex=0.3)
		
	} else{
		axis(2, at = suby, labels = b, las=2, cex=0.3)
	}
}

par(fig=c(0.7,1,0.2,0.5), new=TRUE)
image(matrix(1:26, ncol=1, byrow=FALSE), col=color, xlab="", axes=FALSE)
axis(1, at=c(0, 0.5, 1), labels=c("1.3", "0", "-1.3"), las=1, cex=0.5)

#legend(1500,30, c("1.5","0","-2.6"), pch=19, col = c("navy", "white", "firebrick3"), 
#text.col =  c("navy", "white", "firebrick3"), bg = 'gray90')

dev.off()














################################################################################
# mymhtplot function
################################################################################
mymhtplot <- function (data, usepos = FALSE, logscale = TRUE, base = 10, cutoffs = c(4, 
6, 8), colors = NULL, labels = NULL, gap = NULL, ...) 
{
    data2 <- data[!apply(is.na(data), 1, any), ]
    chr <- data2[, 1] #chr
    pos <- newpos <- data2[, 2] #position
    p <- data2[, 3] #p-value
    tablechr <- table(chr)
    allchr <- as.vector(tablechr)
    n.chr <- length(allchr)
    colorlist <- colors()
    if (is.null(colors)) 
	colors <- sample(colorlist, n.chr)
    if (is.null(labels)) 
	labels <- names(tablechr)
    if (is.null(gap)) 
	gap <- 0
    CMindex <- cumsum(allchr) # cumulatived marker number in each chr
    for (i in 1:n.chr) {
        u <- CMindex[i]
        l <- CMindex[i] - allchr[i] + 1
        chr <- l:u
        if (usepos) 
		d <- diff(pos[chr])
        else d <- rep(1, allchr[i] - 1)
        newpos[chr] <- c(gap, d)
    }
    CM <- cumsum(as.numeric(newpos))    #calculate the actually position
	
    args <- list(...)
    if ("ylim" %in% names(args)) 
	dp <- seq(ylim[1], ylim[2], length = sum(allchr))
	else dp <- seq(min(p), max(p), length = sum(allchr))
#else dp <- seq(-2.5, 3.5, length = sum(allchr))
    if (logscale) 
	y <- -log(dp, base)
    else y <- dp
    par(xaxt = "n", yaxt = "n")
    plot(CM, y, type = "n", xlab = "", ylab = "", axes = FALSE, 
		 ...)
    axis(1, tick = FALSE)
    axis(2, tick = FALSE)
    par(xaxt = "s", yaxt = "s")
    for (i in 1:n.chr) {
        u <- CMindex[i]
        l <- CMindex[i] - allchr[i] + 1
		mid <- 
        chr <- l:u
        cat("Plotting points ", l, "-", u, "\n")
        if (logscale) 
		y <- -log(p[chr], base)
        else y <- p[chr]
        points(CM[chr], y, col = colors[i], ...)
		
		yl <- CMindex[i] - allchr[i]/2 + 1;
        axis(1, at = CM[yl], labels = paste("chr", labels[i], sep=""), 
			 ...)
		abline(v= ifelse(i == 1, CM[1], CM[l]), col = "gray60")
    }
    abline(h = cutoffs, col="red")
#axis(3, at = cutoffs, lwd = 0)
    mtext(ifelse(logscale, paste("-log", base, "(Observed value)", 
								 sep = ""), ""), 2, line = 2.5, las = 0)
    if ("xlab" %in% names(args)) 
	xlabel <- xlab
    else xlabel <- ifelse(is.null(names(chr)), "Chromosome", 
						  names(chr))
    mtext(xlabel, 1, line = 2.5, las = 0)
}
