#Purpose: QTL analysis using QTLcart
#last updated: 2/17/2011

setwd("/home/yangjl/Research/QTL/QTLcart")
geno <- read.csv("geno.ibm.snp_1016.csv")
bckup <- geno
pheno <- read.csv("pheno.ibm.ril_15.csv")

tem <- geno[1:2,]
tem <- t(tem)
write.table(tem, "ibm_map.inp", quote=FALSE)


geno[is.na(geno)] <- "--"
pheno[is.na(pheno)] <- "."

d <- merge(geno, pheno, by="ID")

td <- t(d)
write.table(td, "ibm_trait.inp", quote=FALSE)

###--------
# Plot the results

cim3 <- read.table("mycim3", header=TRUE)
cim3 <- cim3[, 1:8]

subcim3 <- cim3[,c(1,3,4)]
subcim3$position <- 100*subcim3$position


colors <- rep(c("blue", "green"), 5)
mymhtplot(subcim3, colors=colors, gap=10, pch=19, logscale=FALSE, cutoffs=20)


###--------CIM6
cim6 <- read.table("mycim6", header=TRUE)
cim6 <- cim6[,1:8]
subcim6 <- cim6[, c(1,3,4)]
subcim6$position <- 100*subcim6$position

thred <- read.table("myqtlcart.z6e", header=FALSE)
quantile(thred$V2, 0.95)
colors <- rep(c("blue", "green"), 5)
mymhtplot(subcim6, colors=colors, gap=10, pch=19, logscale=FALSE, cutoffs=10.7)

####results plotting
qtlcart.z6c
qtlcart.z6e

joint <- read.table("mycim_namkrn_joint_2wm", header=TRUE, comment.char = "#")
joint <- joint[,1:8]
sub <- joint[, c(1,3,4)]
sub$position <- 100* as.numeric(sub$position)
sub$c <- as.numeric(sub$c)
sub$H0.H1 <- as.numeric(sub$H0.H1)


thred <- read.table("myqtlcart.z6e", header=FALSE)
quantile(thred$V2, 0.95)
colors <- rep(c("blue", "green"), 5)
mymhtplot(sub, colors=colors, gap=10, pch=19, logscale=TURE, cutoffs=3)





################################################################################
# mymhtplot function
################################################################################
mymhtplot <- function (data, usepos = FALSE, logscale = TRUE, base = 10, cutoffs = c(4, 
6, 8), colors = NULL, labels = NULL, gap = NULL, ...) 
{
    data2 <- data[!apply(is.na(data), 1, any), ]
    chr <- data2[, 1]
    pos <- newpos <- data2[, 2]
    p <- data2[, 3]
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
    CMindex <- cumsum(allchr)
    for (i in 1:n.chr) {
        u <- CMindex[i]
        l <- CMindex[i] - allchr[i] + 1
        chr <- l:u
        if (usepos) 
		d <- diff(pos[chr])
        else d <- rep(1, allchr[i] - 1)
        newpos[chr] <- c(gap, d)
    }
    CM <- cumsum(as.numeric(newpos))
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
    abline(h = cutoffs)
    axis(2, at = cutoffs, lwd = 0)
    mtext(ifelse(logscale, paste("-log", base, "(Observed value)", 
								 sep = ""), "Observed value"), 2, line = 2.5, las = 0)
    if ("xlab" %in% names(args)) 
	xlabel <- xlab
    else xlabel <- ifelse(is.null(names(chr)), "Chromosome", 
						  names(chr))
    mtext(xlabel, 1, line = 2.5, las = 0)
}
