# Jinliang Yang
# Purpose: quick plot for GWAS results
# start: 2.11.2012
# add the RNA-seq LM regression data




######################################################################################################
###--- KRN v2 different chromosome
qmht <- function(infile="KRN_meta_chr1.mrkRes1", gap=10000000, window=FALSE, col=col, cex=0.3, cutoff=-9,
				gwas="GenSel", sub="ALL"){
# OPTIONS: 
# window, use ModelFreq or Windowfreq for plotting, if FALSE using ModelFreq; if TRUE using WindowFreq;
# col: color for each chromosome
# cutoff: Setting the cutoff line, < 0 will not be plotted
# sub: subset the data, sub=ALL/ADD/GENO_2NF

#############
# location: 129.186.85.7
source("/Users/yangjl/Documents/Codes/reusableRcode.r")
##########################################################################################
	if(gwas=="GenSel"){
		bayes <- read.table(infile, header=TRUE)
		bayes$CHR = -9;
		bayes$BP = bayes$ID;
	
		outb <- data.frame();
		for(chr in 1:10){
			idx <- grep(paste(chr,"_", sep=""), bayes$ID);
			bayes2 <- bayes[idx,]
			bayes2$CHR = chr;
			bayes2$BP <- sub(paste(chr,"_", sep=""), "", bayes2$BP);
		
			outb <- rbind(outb, bayes2);
		}
	
		outb$CHR <- as.numeric(as.character(outb$CHR));
		outb$BP <- as.numeric(as.character(outb$BP));
		outb <- newpos(outb, GAP=gap);	
		
		if(window==FALSE){
			plot(x=-10, y=-10, ylim=c(0,1), xlim=c(0, max(outb$pos)), xlab="pos (bp)", ylab="Model Freq");
			for(i in 1:10){
				p <- subset(outb, CHR==i);
				points(x=p$pos, y=p$ModelFreq, col=col[i], cex=cex);
			}
		} else{
			plot(x=-10, y=-10, ylim=c(0,1), xlim=c(0, max(outb$pos)), xlab="pos (bp)", ylab="Window Freq");
			for(i in 1:10){
				p <- subset(outb, CHR==i);
				points(x=p$pos, y=p$WindowFreq, col=col[i], cex=cex);
			}
		}
	} else if(gwas=="PLINK"){
		pval <- read.table(infile, header=TRUE);
		pval <- pval[!is.na(pval$P),]
		
		if(sub!="ALL"){
			pval <- subset(pval, TEST==sub);
		}
				
		pval$log <- -log10(pval$P);
		
		if(nrow(pval[pval$P==0,])>0){
			tem <- pval$P;
			pval[pval$P==0, ]$log <- -log10(min(tem[tem>0]));
		}
		
		outb <- newpos(pval, GAP=gap);
		
		plot(x=-10, y=-10, ylim=c(0,max(p$log)), xlim=c(0, max(outb$pos)), xlab="pos (bp)", ylab="pval(-log10)");
			for(i in 1:10){
				p <- subset(outb, CHR==i);
				points(x=p$pos, y=p$log, col=col[i], cex=cex);
		}
		
	} else {
		print("###---study type error!---#");
	}
	

	if(cutoff>0){
		abline(h=cutoff, col="red", lty=2);
		
	}
	return(outb);
}

##################################################################
bayes <- qmht(infile="KRN_2nd_run.mrkRes1", gap=10000000, col=rep(c("blue","darkgreen"), 5), 
		gwas="GenSel", window=FALSE )

pval  <- qmht(infile="~/Documents/GWAS2_KRN/Method/PLINK/run_2nd/krn_2run_nocovar.qassoc",
		 gap=10000000, col=rep(c("blue","darkgreen"), 5), window=TRUE, gwas="PLINK", 
		 cutoff=-9, sub="ALL")

pval2  <- qmht(infile="~/Documents/GWAS2_KRN/Method/PLINK/run_2nd/krn_2run_dom.assoc.linear",
		 gap=10000000, col=rep(c("blue","darkgreen"), 5), window=TRUE, gwas="PLINK", 
		 cutoff=-9, sub="ADD")






pval1 <- read.table("~/Documents/GWAS2_KRN/Method/PLINK/run_2nd/krn_2run_nocovar.qassoc", header=T)
pval2 <- read.table("krn_2run_nocovar.qassoc", header=T)

pval1 <- pval1[!is.na(pval1$P),]
pval1 <- pval1[order(pval1$P, decreasing=F),]










