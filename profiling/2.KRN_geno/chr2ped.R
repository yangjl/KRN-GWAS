#Jinliang yang
#Purpose: formatting the HapMap2 V2 SNPs for imputation
#start: 6.6.2012

chr2ped <- function(chrfile="maizeHapMapV2_B73RefGenV2_201203028_chr2.hmp.txt", output="maizeHapMap2_chr2_s1"){
  #### this function read in the maizeHapMap2 file by chr and then do the following:
  # 1. remove the B73 != Ref SNP
  # 2. remove the duplicated SNP
  # 3. change the coding of -/+ to A/T;
  # 4. change the Y coding to alt SNP
  # 5. Output the following three files 1) output.dsnp; 2) output.tmap; 3) output.tped
  
  begTime <- Sys.time();
  # read in the first 5 rows of the inputfile
  tab5rows <- read.delim(chrfile, header = TRUE, nrows = 5)
  classes <- sapply(tab5rows, class)
  chr <- read.delim(chrfile, header = TRUE, colClasses = classes);
  print("###---Finished reading the file ---###");
  
  # NAM files
  nam <- read.table("~/Documents/VariationDB/nam_parents.txt", header=TRUE)
  namp <- paste(toupper(nam$parent), ".MZ", sep="");
  
  #
  chr <- chr[, c("rs.", "alleles", "chrom", "pos", namp)];
  names(chr) <- c("rs", "alleles", "chr", "pos", as.character(nam$pop));
  tot <- nrow(chr); 
  
  ### remove the duplicated rows ###
  chr <- chr[!duplicated(chr$rs), ]
  chr$ref <- substr(chr$alleles, 1, 1);
  chr$alt <- substr(chr$alleles, 3, 3);
  flt1 <- nrow(chr);
  
  
  ### remove the one that B73 not consistent with Ref ###
  chr <- subset(chr, B73 == "N" | B73 == ref);
  chr$B73 <- chr$ref;
  flt2 <- nrow(chr);
  
  ### change the -/+ to A/T  ###
  alleles <- chr$alleles;
  chr[chr == "-"] <- "A";
  chr[chr == "+"] <- "T";
  chr$alleles <- alleles;
  
  
  
  ### re-coding the "Y", "R" ... to alt allele ###
  for(i in 6:31){
    chr[chr[,i]!="A" & chr[,i]!="T" & chr[,i] !="C" & chr[,i] !="G" & chr[,i]!="N", i] <- 
      chr[chr[,i]!="A" & chr[,i]!="T" & chr[,i] !="C" & chr[,i] !="G" & chr[,i]!="N", ]$alt;
  }
  chr <- chr[, 1:31];
  
  ######################## OUTPUT for different purpose ##################################
  #for imputation use
  print(cat("total SNP", tot, "\n",
            "duplicated", tot-flt1, "\n",
            "B73 != Ref", flt1-flt2, "\n",
            "remaining SNP", flt2, "\n", sep=" "));
  write.table(chr[, 1:31], paste(output, "dsnp", sep="."), sep="\t", quote=FALSE, row.names=FALSE)
  
  #-------->output the .tmap
  nm <- names(chr)[5:31]
  chr.tmap<- data.frame(family= nm, individual=c(27, 1:26), paternal=0, maternal=0, sex=0, pheno=0)
  write.table(chr.tmap, paste(output, "tfam", sep="."), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  
  #-------->output the .tped
  ### Haploid to Diploid:
  for (i in 5:31){
    chr[,i] <- paste(chr[,i], chr[,i], sep=" ")
  }
  chr$geno <- 0;
  chr <- chr[, c(3,1,32,4,5:31)]
  
  write.table(chr, paste(output, "tped", sep="."), sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE); 
  print(Sys.time() - begTime);
}

