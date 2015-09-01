# Jinliang Yang
# Purpose: Gene Co-Occurrence Test and permutation test
# start: 8.25.2011

# location: 129.186.85.9
setwd("/mnt/02/yangjl/GWAS/plot")

######################################################################################################
# PLINK results to find the SNPID
plink_pval <- read.table("wmean_merged_GWAS.qassoc", header=TRUE)
plink_pval$snpid <- paste(plink_pval$CHR, plink_pval$BP, sep="_")
plink_pval <- plink_pval[,c("SNP", "CHR", "BP", "P", "snpid")]

gene_v1 <- read.table("~/GWAS/Anno/ZmB73_4a.53_filtered_genes.first.transcript.gff", header=T)
gene_v1 <- subset(gene_v1, Feature=="gene")

######### Evidence supported genes #############
#### Select 10 of the most signifcant ones in the LD block
#### and the ones tightly linked with the target SNP r2 > 0.7
#### snp+- 25k
#### requirement: should have title "Gene", which is the gene model names
evidence_gene <- read.csv("Excel_CoGeBlast_maize_filtered.csv")
evidence_gene <- merge(gene_v1,evidence_gene,  by.x="Gene", by.y="Closest.Feature")
evidence_gene <- evidence_gene[,c(1:7, 11,15)]

#snp18 <- read.table("snp18.rsq.ld", header=TRUE, sep="\t")
#ldblock <- read.table("snp1000t.rsq.ld", header=TRUE)

#######################################################################################################################
# Gene Co-Occurrence Test
#######################################################################################################################
# Data formatting for the Gene Co-Occurrence Test:
# calculate each permutation:
# input=output from the GWAS software:
gco.test <- function(nperm=1, perm=TRUE, output="gco.perm.txt", input="snp1000t", 
	pvalnum=10, rsq_cutoff=0.7, LENGTH=1000000, knownGeneList=evidence_gene, geneModel=gene_v1) {
		
		ldblock <- read.table(paste(input, "rsq.ld",sep="."), header=TRUE)
		
		if(perm){
			inputr <- paste(input, "r.ld",sep=".");
			inputsnp <- paste(input, "snp.ld", sep=".");
			startrow=18*(nperm-1)+1;
			endrow=18*nperm;
			temrfile <- paste("temrfile", nperm, sep="")
			temsnpfile <- paste("temsnpfile", nperm, sep="")
		
			stdin1 <- paste("sed -n '", startrow,",", endrow, "p' ", inputr, " > ", temrfile, sep="")
			stdin2 <- paste("sed -n '", startrow,",", endrow, "p' ", inputsnp, " > ", temsnpfile, sep="")
		
			system(stdin1, intern=FALSE)
			system(stdin2, intern=FALSE)
		} else{
			temrfile <- paste(input, "r.ld", sep=".");
			temsnpfile <- paste(input, "snp.ld", sep=".");
		}

		hitgene <- data.frame();

		for(i in 1:18){
			perm1.r <- read.table(temrfile, sep=" ",fill=TRUE, header=FALSE, nrows=1, skip=i-1)
			perm1.snp <- read.table(temsnpfile, sep=" ", fill=TRUE, header=FALSE, nrows=1, skip=i-1)
	
			temr <- t(perm1.r)
			temr <- temr[,1]
			temr <- temr[!is.na(temr)]
	
			temsnp <- t(perm1.snp)
			temsnp <- temsnp[,1]
			temsnp <- temsnp[temsnp!="" & !is.na(temsnp)]
	
			snp1 <- data.frame(SNP=temsnp, r=temr)
			snpld <- merge(snp1, plink_pval , by="SNP")
			snpld$absr <- abs(snpld$r)
	
			lsnp <- as.character(ldblock[18*(nperm-1)+i,"L_snp"])
			rsnp <- as.character(ldblock[18*(nperm-1)+i,"R_snp"])
			snpld <- snpld[order(snpld$BP),]
			
			if(lsnp != head(snpld,1)$SNP & rsnp != tail(snpld, 1)$SNP) {
				stop("target SNP not in the LD block")
			} else {
				#### Target SNP linked genes
				targetsnp <- merge(ldblock[18*(nperm-1)+i, ], plink_pval, by.x="target_SNP", by.y="SNP")
				mygene <- subset(gene_v1, (Chr==targetsnp$CHR & Start <= targetsnp$BP+LENGTH & Start >= targetsnp$BP-LENGTH)
						 | (Chr==targetsnp$CHR & Start >= targetsnp$BP-LENGTH & End <= targetsnp$BP+LENGTH))
		
		
				#### LD SNP linked genes
				if(ldblock[18*(nperm-1)+i,]$max_rsq >= rsq_cutoff){
					linksnp <- snpld[snpld$absr >=rsq_cutoff,]
				}else{
				#### top10 p-value SNP linked genes
					linksnp <- head(snpld[order(snpld$P, decreasing=F),],pvalnum)
				}
		
					###########################################################################
					# Get the SNP linked genes
					for(j in 1:nrow(linksnp)){
						temgene <- subset(gene_v1, (Chr==linksnp$CHR[j] & Start <= linksnp$BP[j]+LENGTH & Start >= linksnp$BP[j]-LENGTH)
							  | (Chr==linksnp$CHR[j] & Start >= linksnp$BP[j]-LENGTH & End <= linksnp$BP[j]+LENGTH))
						mygene <- rbind(mygene, temgene)
					}
				
					mygene <- mygene[!duplicated(mygene$Gene),]
					if(nrow(mygene) > 0) mygene$targetsnp <- targetsnp$snpid
					hitgene_tem <- merge(knownGeneList, mygene, by="Gene")
				if(!perm & nrow(mygene)>0) write.table(mygene, output, sep="\t", col.names=F, row.names=FALSE, quote=F, append=TRUE)

			} #### the end of else ### Target SNP linked genes
			
			hitgene <- rbind(hitgene, hitgene_tem)
		} # End of the for loop
		
		############################################################################
		##### for the original test:
		
		if(!perm){
			return(hitgene);	
		} else { #### for the permutation test:
			stdin3 <- paste("rm", temrfile, temsnpfile, sep=" ")
			system(stdin3, intern=FALSE);
			
			ntimes <- nrow(hitgene)
			report <- data.frame(hitgene=ntimes, Rep=nperm)
			write.table(report, output, sep="\t", col.names=F, row.names=F, quote=F, append=TRUE)
			return(ntimes);
		}
		
}# End of the function

save.image("GCO.perm.RData")
#####################################################################################

library(snow, lib="~/bin/Rlib/")
library(snowfall, lib="~/bin/Rlib/")

sfInit(parallel=TRUE, cpus=20)
sfExportAll()
sfLapply(1:1000, gco.test, output="gco.perm1000", input="snp1000t")
sfRemoveAll()
sfStop()

########################### SERVER.7 ############################
load("GCO.perm.RData")
cob_gene <- read.csv("COB_genotype.csv")
dim(cob_gene)
cob_gene <- merge(gene_v1, cob_gene,  by="Gene")
cob_gene <- cob_gene[,1:6]

gco.test(nperm=1, perm=FALSE, output="gcotest_cob.txt", input="snp18", 
pvalnum=40, rsq_cutoff=0.8, LENGTH=1000000, knownGeneList=cob_gene, geneModel=gene_v1) 


a <- matrix(c(2344,32540,10,67), nrow=2)
chisq.test(a)

b <- matrix(c(4,2340, 3, 30193), nrow=2)
chisq.test(b)

fisher.exact(b)


#################################################################
list1 <- read.table("gcotest_cob.txt", sep="\t", header=FALSE)
cyc2 <- read.delim("maizecyc_v2_pwy.txt", header=TRUE)

#idx <- grep("cytokinin",cyc2$pathway_name, ignore.case=TRUE)

#cytokinins degradation = PWY-2841
pwy2841 <- subset(cyc2, pathway_id=="PWY-2841")
u1 <- merge(pwy2841, list1, by.x="gramene_gene_id", by.y="V1")

#cytokinins o-glucoside biosynthesis = PWY-2902
pwy2902 <- subset(cyc2, pathway_id=="PWY-2902")
u2 <- merge(pwy2902, list1, by.x="gramene_gene_id", by.y="V1")
dim(u2)
d <- matrix(c(18,257,2326,29939), nrow=2)
fisher.test(d)

tem2 <- read.table("gco.perm1000", sep="\t", header=FALSE)













