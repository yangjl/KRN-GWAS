# Jinliang Yang
# last update: 4.8.2011 | 2.9.2011
# purpose: summary the CIM6 QTL results:

# read.me
# plot the QTL effects heat map of CIM results
# Summary the QTL from CIM results
# geno code: IBM AA 0 Mo17; aa 1 B73
#			 NAM AA 2 nonB73; aa 0 B73
# QTL Additive Effect: should be 2x
# LOD = 0.217LR

# location: server.7

setwd("/Users/yangjl/Documents/GWAS2_KRN/QTL/NAMKRN_joint/")

########################################################################
# PERL to extract the useful info
########################################################################
./myqtlpaser.pl nam_krn.eqt6 > mynam_krn.eqt6
./myqtlpaser.pl nam_krn.z6 > mynam_krn.z6



################## INPUT physical pos DATA #############################
#--------------read in the marker names and postions from JZmapqtl results:
mrk <- read.table("mynam_krn.jz63", header=FALSE)
names(mrk) <- c("Chr", "Num", "Marker", "Pos", "LR", "Add")
mrk$UID <- paste(mrk$Chr, mrk$Num, sep="_")
mrk <- mrk[,c(1:4,7)]
dim(mrk)
#[1] 1864    5

#--read in the physical map (missing ones been imputed)
phymap <- read.delim("../NAM_map_AGPv2.imputed.txt", header=TRUE)
dim(phymap)
phymap <- phymap[, c(1:4,8,9)]
#[1] 1106    6

myphymap <- merge(mrk, phymap, by.x="Marker", by.y="marker")
dim(myphymap)
#[1] 1864   10

#############
qtlz6 <- read.table("mynam_krn.z6", header=TRUE)
qtlz6$UID <- paste(qtlz6$c, qtlz6$m, sep="_")

qtlz6 <- qtlz6[order(qtlz6$UID),]
myphymap <- myphymap[order(myphymap$UID),]

qtlz6 <- cbind(qtlz6[, c("UID", "H0.H1")], myphymap, by="UID")
names(qtlz6)[1:2] <- c("UID1", "LOD");
qtlz6$LOD <- 0.217*qtlz6$LOD
qtlz6 <- qtlz6[order(qtlz6$ch, qtlz6$genetpos),]
qtlz6 <- qtlz6[!duplicated(qtlz6$Marker),]
qtlz6$Pos <- 100 * qtlz6$Pos
dim(qtlz6)
# 886 13

################## QTL DATA ###########################################
#------------ Joint CIM6 mapping QTL results
qtl <- read.table("mynam_krn.eqt6", header=FALSE)
qtl <- qtl[,2:6]
names(qtl) <- c("Chr", "Marker", "Position", "LOD", "Additive")
qtl$UID <- paste(qtl$Chr, qtl$Marker, sep="_")
dim(qtl)
#82 6

#------ FIND the marker names
myqtl <- merge(qtl, myphymap, by="UID")
myqtl <- myqtl[!duplicated(myqtl$UID),]

myqtl <- myqtl[, c("UID","Chr.x", "Position", "LOD", "Additive", "Marker.y", "genetpos", "imputepos")]
names(myqtl) <- c("UID","Chr", "Position", "LOD", "Additive", "Marker", "genetpos", "imputepos")
myqtl <- myqtl[order(myqtl$Chr, myqtl$Position),]

#############################
#Determining the QTL interval:
#############################

jbackup <- myqtl

jcim6 <- myqtl
jcim6$diff <- c(0, diff(as.numeric(jcim6$Position)))

chrvector <- unique(jcim6$Chr)
for(i in chrvector){	
jcim6[jcim6$Chr==i,]$diff[1] <- 0	
}

jcim6$interval <- 0
intval=0
for (i in 1:nrow(jcim6)){
	if (jcim6$diff[i] > 13 | jcim6$diff[i]==0) intval <- intval +1
	jcim6$interval[i] <- intval	
}

#######-------------->>>> get the significant interval
out <- data.frame(Chr=1, LOD=1,Effect=1, SigMarker="A",SigGenetic=1, SigPhy=1, 
		LeftMarker="A",LeftGenetic=1,LeftPhy=1,RightMarker="A", RightGenetic=1,RightPhy=1)

outqtl <- data.frame();
for (i in 1:max(jcim6$interval)){
	sub <- jcim6[jcim6$interval==i,]
	
	idx <- which.max(sub$LOD)
	myout <- out
	myout$Chr <- unique(sub$Chr)
	myout$LOD <- max(sub$LOD)
	myout$Effect <- -2*sub[idx,]$Additive
	
	myout$SigMarker <- as.character(sub[idx,]$Marker)
	myout$SigGenetic <- as.character(sub[idx,]$Position)
	myout$SigPhy <- as.character(sub[idx,]$imputepos)
	
	idx0 <- which(qtlz6$Marker == as.character(sub[idx,]$Marker));
	
	##########################################
	LOD0 <- qtlz6[idx0,]$LOD;
	if(qtlz6[idx0,]$Chr == qtlz6[idx0-1,]$Chr){
		idx1 <- idx0-1;
		while ((LOD0 - qtlz6[idx1,]$LOD < 1.5) & (qtlz6[idx1,]$Chr == qtlz6[idx1-1,]$Chr)){
			idx1 = idx1 -1;
		}
	} else{
		idx1 <- idx0;
	}
	
	myout$LeftMarker <- as.character(qtlz6[idx1,]$Marker);
	myout$LeftGenetic <- as.character(qtlz6[idx1,]$Pos);
	myout$LeftPhy <- as.character(qtlz6[idx1,]$imputepos);

	##########################################
	if(qtlz6[idx0,]$Chr == qtlz6[idx0+1,]$Chr) {
		idx2 <- idx0+1;
		while (LOD0 - qtlz6[idx2,]$LOD < 1.5 & qtlz6[idx2,]$Chr == qtlz6[idx2+1,]$Chr){
			idx2 = idx2 + 1;
		} 
	} else {
		idx2 = idx0
	};
	
	myout$RightMarker <- as.character(qtlz6[idx2,]$Marker);
	myout$RightGenetic <- as.character(qtlz6[idx2,]$Pos);
	myout$RightPhy <- as.character(qtlz6[idx2,]$imputepos);
	
	outqtl <- rbind(outqtl,myout) 
}


outqtl$trait <- "wmean"
outqtl$method <- "CIM6"
outqtl$pop <- "joint"

write.table(outqtl, "../Joint_CIM6_summary_v2.txt", sep="\t", row.names=FALSE, quote=FALSE)


#############################################################################################################################
########-------------------------------------Seperate QTL summary----------------------------------------------##############
setwd("/Users/yangjl/Documents/KRN_GWAS/QTL/QTLcart/NAMKRN_sep")

########################################################################
# PERL to extract the useful info
########################################################################
nam <- read.table("NAM_name.txt", header=FALSE);
outfile = "get_mytable.sh"
for (i in 1:25){
	line1 <- paste("./myqtlpaser.pl", 
				   paste(nam$V1[i], "_krn.z6", sep=""),
				   ">",
				   paste("my", nam$V1[i], "_krn.z6", sep=""),
				   sep = " ")
	line2 <- paste("./myqtlpaser.pl", 
				   paste(nam$V1[i], "_krn.z6e", sep=""),
				   ">",
				   paste("my", nam$V1[i], "_krn.z6e", sep=""),
				   sep = " ")
	line3 <- paste("./myqtlpaser.pl", 
				   paste(nam$V1[i], "_krn.eqt", sep=""),
				   ">",
				   paste("my", nam$V1[i], "_krn.eqt", sep=""),
				   sep = " ")
	
	cat(line1, line2, line3, sep="\n", file=outfile, append=TRUE)
	}

######----->RUN SHELL PROGRAM
sh get_mytable.sh


################## INPUT physical pos DATA #############################
#--------------read in the marker names and postions from JZmapqtl results:
mrk <- read.table("../NAMKRN_joint/mynam_krn.jz63", header=FALSE)
names(mrk) <- c("Chr", "Num", "Marker", "Pos", "LR", "Add")
mrk$UID <- paste(mrk$Chr, mrk$Num, sep="_")
mrk <- mrk[,c(1:4,7)]
dim(mrk)
#[1] 1864    5

#--read in the physical map (missing ones been imputed)
#phymap <- read.delim("../NAM_map_20080419.Eddy.imputed.txt", header=TRUE)
dim(phymap)
#[1] 1106    6

myphymap <- merge(mrk, phymap, by.x="Marker", by.y="marker")
dim(myphymap)
#[1] 1864   10

################## QTL DATA ###########################################
#------------ Seperate QTL mapping results
nam <- read.table("NAM_name.txt", header=FALSE);
for (k in c(1:11,13:25)){

	p <- paste(nam$V1[k],"_krn.eqt", sep="")
	myp <- paste("my", p, sep="")

	qtl <- read.table(myp, header=FALSE)
	qtl <- qtl[,2:6]
	names(qtl) <- c("Chr", "Marker", "Position", "LOD", "Additive")
	qtl$UID <- paste(qtl$Chr, qtl$Marker, sep="_")
	#dim(qtl)
	#7 6
	if (nrow(qtl) >= 1){
	#------ FIND the marker names
		myqtl <- merge(qtl, myphymap, by="UID")
		myqtl <- myqtl[!duplicated(myqtl$UID),]

		myqtl <- myqtl[, c("UID","Chr.x", "Position", "LOD", "Additive", "Marker.y", "genetpos", "imputepos")]
		names(myqtl) <- c("UID","Chr", "Position", "LOD", "Additive", "Marker", "genetpos", "imputepos")
		myqtl <- myqtl[order(myqtl$Chr, myqtl$Position),]
		
		#############
		p.z6 <- paste("my",nam$V1[k], "_krn.z6", sep="")
		qtlz6 <- read.table(p.z6, header=TRUE)
		qtlz6$UID <- paste(qtlz6$c, qtlz6$m, sep="_")
		
		qtlz6 <- qtlz6[order(qtlz6$UID),]
		myphymap <- myphymap[order(myphymap$UID),]
		
		qtlz6 <- cbind(qtlz6[, c("UID", "H0.H1")], myphymap, by="UID")
		names(qtlz6)[1:2] <- c("UID1", "LOD");
		qtlz6$LOD <- 0.217*qtlz6$LOD
		qtlz6 <- qtlz6[order(qtlz6$ch, qtlz6$genetpos),]
		qtlz6 <- qtlz6[!duplicated(qtlz6$Marker),]
		qtlz6$Pos <- 100 * qtlz6$Pos
		#dim(qtlz6)

		#############################
		#Determining the QTL interval:
		#############################

		jcim6 <- myqtl
		jcim6$diff <- c(0, diff(as.numeric(jcim6$Position)))

		chrvector <- unique(jcim6$Chr)
		for(i in chrvector){	
			jcim6[jcim6$Chr==i,]$diff[1] <- 0	
		}

		jcim6$interval <- 0
		intval=0
		for (i in 1:nrow(jcim6)){
			if (jcim6$diff[i] > 13 | jcim6$diff[i]==0) intval <- intval +1
			jcim6$interval[i] <- intval	
		}


	#######-------------->>>> get the significant interval
		out <- data.frame(Chr=1, LOD=1,Effect=1, SigMarker="A",SigGenetic=1, SigPhy=1, 
						  LeftMarker="A",LeftGenetic=1,LeftPhy=1,RightMarker="A", RightGenetic=1,RightPhy=1)

		outqtl <- data.frame();
		for (i in 1:max(jcim6$interval)){
			sub <- jcim6[jcim6$interval==i,]
			
			idx <- which.max(sub$LOD)
			myout <- out
			myout$Chr <- unique(sub$Chr)
			myout$LOD <- max(sub$LOD)
			myout$Effect <- -2*sub[idx,]$Additive
			
			myout$SigMarker <- as.character(sub[idx,]$Marker)
			myout$SigGenetic <- as.character(sub[idx,]$Position)
			myout$SigPhy <- as.character(sub[idx,]$imputepos)
			
			idx0 <- which(qtlz6$Marker == as.character(sub[idx,]$Marker));
			
##########################################
			LOD0 <- qtlz6[idx0,]$LOD;
			if(qtlz6[idx0,]$Chr == qtlz6[idx0-1,]$Chr){
				idx1 <- idx0-1;
				while ((LOD0 - qtlz6[idx1,]$LOD < 1.5) & (qtlz6[idx1,]$Chr == qtlz6[idx1-1,]$Chr)){
					idx1 = idx1 -1;
				}
			} else{
				idx1 <- idx0;
			}
			
			myout$LeftMarker <- as.character(qtlz6[idx1,]$Marker);
			myout$LeftGenetic <- as.character(qtlz6[idx1,]$Pos);
			myout$LeftPhy <- as.character(qtlz6[idx1,]$imputepos);
			
##########################################
			if(qtlz6[idx0,]$Chr == qtlz6[idx0+1,]$Chr) {
				idx2 <- idx0+1;
				while (LOD0 - qtlz6[idx2,]$LOD < 1.5 & qtlz6[idx2,]$Chr == qtlz6[idx2+1,]$Chr){
					idx2 = idx2 + 1;
				} 
			} else {
				idx2 = idx0
			};
			
			myout$RightMarker <- as.character(qtlz6[idx2,]$Marker);
			myout$RightGenetic <- as.character(qtlz6[idx2,]$Pos);
			myout$RightPhy <- as.character(qtlz6[idx2,]$imputepos);
			
			outqtl <- rbind(outqtl,myout) 
		}
		
		

		outqtl$trait <- "wmean"
		outqtl$method <- "CIM6"
		outqtl$pop <- nam$V1[k]

		write.table(outqtl, "../Sep_CIM6_summary_v2.txt", sep="\t", 
					row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE)
		
	} else{
	print (paste("NO QTL", nam$V1[k], sep="\n"))
	}

}# the first layer of the loop #######
###############
myqtlsum <- read.delim("../Sep_CIM6_summary_v2.txt", header=FALSE)
names(myqtlsum) <- c(names(out),"trait", "method", "pop")
write.table(myqtlsum, "../Sep_CIM6_summary_v2.txt", sep="\t", row.names=FALSE, quote=FALSE)
########---------------------------------End of Seperate QTL summary-------------------------------------------##############
#############################################################################################################################




########################################################################
# Impute the missing physical position:
phymap <- read.csv("NAM_map_AGPv2.csv")
phymap$imputepos <- 0
phymap <- phymap[order(phymap$ch, phymap$genetpos, phymap$v2),]

for (i in 1:nrow(phymap)){
	mymrk <- phymap[i,]
	if(!is.na(mymrk$v2)){
		phymap[i,]$imputepos <-mymrk$v2
	}
	
	else {
		t1 <- i-1
		while (is.na(phymap[t1,]$v2)){
			t1 <- t1-1
		}
		up <- phymap[t1,]
		
		t2 <- i+1
		while(is.na(phymap[t2,]$v2)){
			t2 <- t2+1
		}
		down <- phymap[t2,]
		
		if(up$ch == mymrk$ch & down$ch==mymrk$ch){
			if(down$genetpos-up$genetpos == 0)
			{phymap[i,]$imputepos = round(mean(up$v2+down$v2),0)}
			else{
# Major imputation formular:######
				phymap[i,]$imputepos <- round((down$v2-up$v2)*(mymrk$genetpos-up$genetpos)/(down$genetpos-up$genetpos)+up$v2,0)
			}
		}#if end#
		
		if(up$ch==mymrk$ch & down$ch!=mymrk$ch){
			phymap[i,]$imputepos <- up$v2
		}
		
		if(up$ch!=mymrk$ch & down$ch==mymrk$ch){
			phymap[i,]$imputepos <- down$v2
		}
		
	}#else end#
}# forend#
write.table(phymap, "NAM_map_AGPv2.imputed.txt", sep="\t", row.names=FALSE, quote=FALSE)
########################################################################













