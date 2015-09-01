# Data Formatting for QTLcart Use
# update: 02/18/2011
# update: 02/09/2012

# set the working directory
setwd("~/Documents/GWAS2_KRN/QTL/")
######First Joint Mapping###################

###############################################################################
#                      MAP File                                               #
###############################################################################

map <- read.table("markers061208.txt", header=FALSE)
map <- map[order(map$V4),]
map <- map[,c(1,3)]
write.table(map, "nam_map_021911.inp", sep="\t", quote=FALSE, row.names=FALSE)
###------------
# Manually changed the file:
 

###############################################################################
#                      Geno and Pheno File                                    #
###############################################################################

# imputed nam tagging SNP map
tagSNP <- read.table("markergenotypes062508.txt", comment.char="#", header=TRUE)
dim(tagSNP)
#[1] 4699 1111
tagSNP <-tagSNP[,c(-2,-3,-4,-5)]
names(tagSNP) <- c("Geno", as.character(map$V1))
dim(tagSNP)
#[1] 4699 1107

trait <- read.table("nam_pheno_krn-p_040211.txt", header=TRUE)
dim(trait)
#[1] 4780   15

trait$FID <- as.character(trait$FID)
for (i in 1:25){
	
	p <- unique(trait$FID)[i]
	trait[trait$FID == p,]$FID <- as.numeric(i)
}
	

joint <- merge(trait, tagSNP, by.x="Genotype", by.y="Geno")
dim(joint)
#[1] 4541 1121

###############################################################################
### start the output => Joint
###############################################################################
setwd("/home/yangjl/Research/QTL/QTLcart/NAMKRN_joint")
outputfile = "nam_krn_040211.inp"
cat("# 123456789     -filetype cross.inp", 
"# Documentation at the end",
"-Cross RI1",
"-traits 5",
"-otraits 4",
paste("-SampleSize", nrow(joint), sep="\t"),
"-case   no",
"-TranslationTable",
"AA      2       2",
"Aa      1       -",
"aa      0       0",
"A-      12      1.5",
"a-      10      0.5",
"--      -1      1",
"-start individuals markers",file=outputfile, sep = "\n")


mygeno <- joint[,c(1, 16:1121)]
write.table(mygeno, outputfile,sep="\t",append=TRUE, quote=FALSE, col.names=FALSE, row.names=FALSE)

cat("-stop individuals markers",
"-start individual traits 5 pop mean wmean TR PSSR named", append=TRUE,sep="\n", file = outputfile)

myt1 <- joint[, c(1,2, 6,7,10,14)]
write.table(myt1, outputfile,sep="\t",append=TRUE, quote=FALSE, col.names=FALSE,  row.names=FALSE)

cat("-stop individuals traits", "-start individual otraits 4 TP PSSP +Note Year named", sep="\n",
append=TRUE, file = outputfile)

#otrait
myot <- joint[,c("Genotype", "TP", "PSSP", "Note", "Year")]
write.table(myot, outputfile, append=TRUE, sep="\t",quote=FALSE, col.names=FALSE, row.names=FALSE)

cat("-stop individuals otraits", "-quit", sep="\n", append=TRUE, file=outputfile)







###############################################################################
#                      Seperate Analysis                                      #
###############################################################################
setwd("/home/yangjl/Research/QTL/QTLcart/NAMKRN_sep")

mynote <- unique(joint$Note)
for(i in 1:25){
	b <- as.character(mynote[i])
	myfile <- joint[joint$Note==b,]
	print(paste(b, dim(myfile), "\n", sep="----"))
	
	outputfile = paste(b, "_krn_040211.inp", sep="")
	cat("# 123456789     -filetype cross.inp", 
		"# Documentation at the end",
		"-Cross RI1",
		"-traits 4",
		"-otraits 4",
		paste("-SampleSize", nrow(myfile), sep="\t"),
		"-case   no",
		"-TranslationTable",
		"AA      2       2",
		"Aa      1       -",
		"aa      0       0",
		"A-      12      1.5",
		"a-      10      0.5",
		"--      -1      1",
		"-start individuals markers",file=outputfile, sep = "\n")
	

	mygeno <- myfile[,c(1, 16:1121)]
	write.table(mygeno, outputfile,sep="\t",append=TRUE, quote=FALSE, col.names=FALSE, row.names=FALSE)
	
	cat("-stop individuals markers",
		"-start individual traits 4 mean wmean TR PSSR named", append=TRUE,sep="\n", file = outputfile)
	
	myt1 <- myfile[, c(1,6,7,10,14)]
	write.table(myt1, outputfile,sep="\t",append=TRUE, quote=FALSE, col.names=FALSE,  row.names=FALSE)
		
	cat("-stop individuals traits", "-start individual otraits 4 TP PSSP Note Year named", sep="\n",
		append=TRUE, file = outputfile)
	
	#otrait
	myot <- myfile[,c("Genotype", "TP", "PSSP", "Note", "Year")]
	write.table(myot, outputfile, append=TRUE, sep="\t",quote=FALSE, col.names=FALSE, row.names=FALSE)
	
	cat("-stop individuals otraits", "-quit", sep="\n", append=TRUE, file=outputfile)
	
}

###############################################################################
#               Preparing the QTLcart Commands:
###############################################################################
setwd("/home/yangjl/Research/QTL/QTLcart/NAMKRN_sep")


#### create the map in the directory
Rmap -i nam_map_021911.inp -X nam_map -A -V

mynote <- as.character(unique(joint$Note))
for(i in 1:25){
	b <- mynote[i]
	stem <- paste(b, "krn", sep="_")
	
	stem.sh <- paste(stem, ".sh", sep="")
	stem.inp <- paste(stem, "_040211.inp", sep="")
	stem.cro <- paste(stem, ".cro", sep="")
	stem.log <- paste(stem, ".log", sep="")
	stem.map <- paste(stem, ".map",sep="")
	
	cat("cp", "nam_map.map", stem.map,"\n", file=stem.sh, sep=" ")
	
#Rcross -i B97_krn_022111.inp -o B97_krn.cro -e B97_krn.log -m nam_map.map -A -V	
	cat("Rcross",
		"-i", stem.inp,
		"-o", stem.cro,
		"-e", stem.log,
		"-m", stem.map,
		"-A -V","\n",
		file=stem.sh, append=TRUE,sep = " ")
	
#SRmapqtl  -X B97_krn -M 2 -t 2 -F 0.05 -B 0.05 -A -V

	cat("SRmapqtl",
		"-X", stem,
		"-M 2 -t 2 -F 0.05 -B 0.05",
		"-A -V",
		"\n",
		file=stem.sh, append=TRUE, sep = " ")
	
# Interval mapping
#Zmapqtl -i B97_krn.cro -o B97_krn.z3 -e B97_krn.log -m nam_map.map -M 3 -t 2 -d 1 -r 0 -A -V
#Zmapqtl -X B97_krn -M 3 -t 2 -d 1 -r 1000 -A 
	cat("Zmapqtl",
		"-i", stem.cro,
		"-o", paste(stem, ".z3", sep=""),
		"-e", stem.log,
		"-m", stem.map,
		"-M 3 -t 2 -d 1 -r 0 -A -V",
		"\n",
		"Zmapqtl",
		"-X", stem,
		"-M 3 -t 2 -d 1 -r 1000 -A -V",
		"\n",
		file=stem.sh, append=TRUE, sep = " ")
	
# Composite Interval Mapping:
#Zmapqtl -X B97_krn -M 6 -t 2 -d 1 -A 
#Zmapqtl -X B97_krn -M 6 -t 2 -d 1 -r 1000 -A -V
	cat("Zmapqtl",
		"-i", stem.cro,
		"-o", paste(stem, ".z6", sep=""),
		"-e", stem.log,
		"-m", stem.map,
		"-M 6 -t 2 -d 1 -r 0 -A -V",
		"\n",
		"Zmapqtl",
		"-X", stem,
		"-M 6 -t 2 -d 1 -r 1000 -A -V",
		"\n",
		file=stem.sh, append=TRUE, sep = " ")

# get the QTL results
#Eqtl -z B97_krn.z3 -o B97_krn.eqt -m B97_krn.map -H 10 -a 0.05 3 -L 1 -I PZ -A -V
	cat("Eqtl",
		"-z", paste(stem, ".z6", sep=""),
		"-o", paste(stem, ".eqt", sep=""),
		"-m", stem.map,
		"-H 10 -a 0.05 -L 1 -I PZ -A -V",
		"\n",
		
		file=stem.sh, append=TRUE, sep = " ")
}


# format to the cross file
Rcross -i B97_krn_022111.inp -o B97_krn.cro -e B97_krn.log -m nam_map.map -A -V
#LRmapqtl -i B97_krn.cro -o B97_krn.lr -m nam_map.map -r 1000 -t 2 -A -V
SRmapqtl  -i B97_krn.cro -o B97_krn.sr -e B97_krn.log -m nam_map.map -M 2 -t 2 -F 0.05 -B 0.05 -A -Vll

# Interval mapping
Zmapqtl -i B97_krn.cro -o B97_krn.z3 -e B97_krn.log -m nam_map.map -M 3 -t 2 -d 1 -r 0 -A
Zmapqtl -X B97_krn -M 3 -t 2 -d 1 -r 1000 -A 

# Composite Interval Mapping:
Zmapqtl -X B97_krn -M 6 -t 2 -d 1 -A 
Zmapqtl -X B97_krn -M 6 -t 2 -d 1 -r 1000 -A -V

# get the QTL results
Eqtl -z B97_krn.z3 -o B97_krn.eqt -m B97_krn.map -H 10 -a 0.05 3 -L 1 -I PZ -A -V

screen
sh B97_krn.sh; sh CML103_krn.sh
sh CML228_krn.sh; sh CML247_krn.sh; sh CML277_krn.sh; sh CML322_krn.sh
sh CML333_krn.sh; sh CML52_krn.sh; sh CML69_krn.sh
sh Hp301_krn.sh; sh Il14H_krn.sh

sh Ki11_krn.sh;	sh Ki3_krn.sh; sh Ky21_krn.sh; sh M162W_krn.sh  
sh M37W_krn.sh; sh Mo18W_krn.sh; sh Ms71_krn.sh; sh NC350_krn.sh; sh NC358_krn.sh
sh Oh43_krn.sh; sh Oh7B_krn.sh; sh P39_krn.sh; sh Tx303_krn.sh; sh Tzi8_krn.sh  
###############################################################################
#               Pasering QTL to a matrix
###############################################################################
setwd("/home/yangjl/Research/QTL/QTLcart/NAMKRN_sep")
#./myqtlpaser.pl B97_krn.z3 > myB97_krn.z3
#./myqtlpaser.pl B97_krn.z6 > myB97_krn.z6
#./myqtlpaser.pl B97_krn.z3e > myB97_krn.z3e
#./myqtlpaser.pl B97_krn.z6e > myB97_krn.z6e

trait <- read.table("nam_pheno_krn-p_040211.txt", header=TRUE)
mynote <- as.character(unique(trait$Note))
for(i in 1:25){
	b <- mynote[i]
	stem <- paste(b, "krn", sep="_")
	
	stem.z3<- paste(stem, ".z3", sep="")
	mystem.z3 <- paste("my", stem.z3, sep="")
	stem.z6<- paste(stem, ".z6", sep="")
	mystem.z6 <- paste("my", stem.z6, sep="")
	stem.z3e<- paste(stem, ".z3e", sep="")
	mystem.z3e <- paste("my", stem.z3e, sep="")
	stem.z6e<- paste(stem, ".z6e", sep="")
	mystem.z6e <- paste("my", stem.z6e, sep="")
	
	cat("./myqtlpaser.pl", stem.z3, ">", mystem.z3, "\n", append=TRUE, file="get_mytable.sh", sep=" ")
	cat("./myqtlpaser.pl", stem.z6, ">", mystem.z6, "\n", append=TRUE, file="get_mytable.sh", sep=" ")
	cat("./myqtlpaser.pl", stem.z3e, ">", mystem.z3e, "\n", append=TRUE, file="get_mytable.sh", sep=" ")
	cat("./myqtlpaser.pl", stem.z6e, ">", mystem.z6e, "\n", append=TRUE, file="get_mytable.sh", sep=" ")
}









