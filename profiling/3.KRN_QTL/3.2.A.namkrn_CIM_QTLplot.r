# Data Formatting for QTLcart Use
# start: 02/18/2011
# update: 04/08/2011

# read.me
# plot the results output from QTLcart
# geno code: IBM AA 0 Mo17; aa 1 B73
#			 NAM AA 2 nonB73; aa 0 B73
# QTL Additive Effect: should be 2x
# LOD = 0.217LR

setwd("/home/yangjl/Research/QTL/QTLcart/")
pdf("QTLcart_040811.pdf")
##################################################################################
#                  plot Joint QTL LOD
##################################################################################
setwd("/home/yangjl/Research/QTL/QTLcart/NAMKRN_joint")

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


####-------start to plot
par (oma=c(0,0,0,0))
ylim=c(0,75)
colors <- rep(c("blue", "green"), 5)
mymhtplot(cim6, colors=colors, usepos = TRUE, gap=20, pch=19, cex=0.5, logscale=FALSE, cutoffs=cutoff2,
main="Joint (Composite Interval Mapping)", ylim=ylim)
axis(2, at=c(0, 10,20,30,40,50,60,70), labels=c(0, 10,20,30,40,50,60,70), las=2)

##################################################################################
#                  plot individual QTL by popl.-Mo17
##################################################################################
setwd("/home/yangjl/Research/QTL/QTLcart/IBM_KRN")

####-------formatting the output file
cim3 <- read.table("myibm_ril_krn.z3", header=TRUE)
cim3 <- cim3[,1:8]
cim3 <- cim3[, c(1,3,4)]
cim3$position <- 100* as.numeric(cim3$position)
names(cim3)[3] <- "lod"
cim3$c <- as.numeric(cim3$c)
cim3$lod <- 0.217*as.numeric(cim3$lod)


####-------determining the threshold
thred <- read.table("myibm_ril_krn.z3e", header=FALSE)
cutoff1 <- quantile(thred$V2, 0.95)*0.217


####-------start to plot
par (mfrow=c(2,1), oma=c(0,0,0,0))
ylim=c(0,12)
colors <- rep(c("blue", "green"), 5)
mymhtplot(cim3, colors=colors, usepos = TRUE, gap=20, pch=19, cex=0.5, logscale=FALSE, cutoffs=cutoff1,
main="Mo17 (Interval Mapping)", ylim=ylim)
axis(2, at=c(1,3,5,7,9,11), labels=c(1,3,5,7,9,11), las=2)

#####################Plotting the CIM6#######################
####-------formatting the output file
cim6 <- read.table("myibm_ril_krn.z6", header=TRUE)
cim6 <- cim6[,1:8]
cim6 <- cim6[, c(1,3,4)]
cim6$position <- 100* as.numeric(cim6$position)
names(cim6)[3] <- "lod"
cim6$c <- as.numeric(cim6$c)
cim6$lod <- 0.217*as.numeric(cim6$lod)


####-------determining the threshold
thred2 <- read.table("myibm_ril_krn.z6e", header=FALSE)
cutoff2 <- quantile(thred2$V2, 0.95)*0.217


####-------start to plot
mymhtplot(cim6, colors=colors, usepos = TRUE, gap=20, pch=19, cex=0.5, logscale=FALSE, cutoffs=cutoff2,
main="Mo17 (Composite Interval Mapping)", ylim=ylim)
axis(2, at=c(1,3,5,7,9,11), labels=c(1,3,5,7,9,11), las=2)


##################################################################################
#                  plot individual QTL by popl.-NAM parents
##################################################################################
setwd("/home/yangjl/Research/QTL/QTLcart/NAMKRN_sep")

mynote <- read.delim("NAM_name.txt", header=FALSE)

for (i in c(1:11,13:25)){ #12, Ki11 has no QTL
	b <- as.character(mynote$V1[i])
	stem <- paste(b, "krn", sep="_")
	
	stem.z6<- paste(stem, ".z6", sep="")
	mystem.z6 <- paste("my", stem.z6, sep="")

	stem.z6e<- paste(stem, ".z6e", sep="")
	mystem.z6e <- paste("my", stem.z6e, sep="")
	
####-------start to plot

	
#####################Plotting the CIM6#######################
####-------formatting the output file
	cim6 <- read.table(mystem.z6, header=TRUE)
	cim6 <- cim6[,1:8]
	cim6 <- cim6[, c(1,3,4)]
	cim6$position <- 100* as.numeric(cim6$position)
	names(cim6)[3] <- "lod"
	cim6$c <- as.numeric(cim6$c)
	cim6$lod <- 0.217*as.numeric(cim6$lod)
	
	
####-------determining the threshold
	thred2 <- read.table(mystem.z6e, header=FALSE)
	cutoff2 <- quantile(thred2$V2, 0.95)*0.217
	
	main2=paste(b, "(Composite Interval Mapping)", sep=" ")
####-------start to plot
	par (mfrow=c(2,1), oma=c(0,0,0,0))
	ylim=c(0,12)
	
	mymhtplot(cim6, colors=colors, usepos = TRUE, gap=20, pch=19, cex=0.5, logscale=FALSE, cutoffs=cutoff2,
			  main=main2, ylim=ylim)
	axis(2, at=c(1,3,5,7,9,11), labels=c(1,3,5,7,9,11), las=2)
}

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
        chr <- l:u
        cat("Plotting points ", l, "-", u, "\n")
        if (logscale) 
		y <- -log(p[chr], base)
        else y <- p[chr]
        points(CM[chr], y, col = colors[i], ...)
        axis(1, at = ifelse(i == 1, CM[1], CM[l]), labels = labels[i], 
			 ...)
    }
    abline(h = cutoffs, col="red")
#axis(3, at = cutoffs, lwd = 0)
    mtext(ifelse(logscale, paste("-log", base, "(Observed value)", 
								 sep = ""), "LOD"), 2, line = 2.5, las = 0)
    if ("xlab" %in% names(args)) 
	xlabel <- xlab
    else xlabel <- ifelse(is.null(names(chr)), "Chromosome", 
						  names(chr))
    mtext(xlabel, 1, line = 2.5, las = 0)
}
