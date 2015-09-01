# Jinliang Yang
# Purpose: to create a summarized plot stacking all the information
# 8.23.2011

# location: 129.186.85.9
####################################################################################
library("grid")
#setwd("/mnt/02/yangjl/GWAS/plot")
ob <- load("~/Documents/KRN_GWAS_v3/cache_large/data4stackplot_v2.RData")	
# use predefined gap in the previous code: gap=10000000

# the basic plot information:
# sets colors based on colors argument.
colors <- rep(c("slateblue", "cyan4"), 5)
#xscale:
cl <- read.csv("~/Documents/KRN_GWAS_v3/GWAS3_proj/data/chr_length_B73v2.csv")
cl <- newpos(cl, GAP=gap)
xmax <- max(cl$pos)

ticks <- cl$pos[1]/2
chrlines <- cl$pos[1]+gap/2
for(i in 2:10){
	ticks <- c(ticks, cl$pos[i-1] + (cl$pos[i]-cl$pos[i-1])/2)
	chrlines <- c(chrlines, cl$pos[i]+gap/2)	
}

#kas82 <- kas82[, c('snpid','Select','chr','pos','ModelFreq')]
#names(kas82)[3:4] <- c('CHR', 'BP')
#kas82 <- newpos(kas82, GAP=gap)


##############################################################################
grid.newpage();

#######################################
####-----1st layer layout
#######################################
layout1 = grid.layout(3,3,
	heights=unit(c(3,1,5), c('lines', 'null', 'lines')),
	widths=unit(c(5,1,3), c('lines','null','lines'))
)
pushViewport(viewport(layout=layout1))

##### yaxis label:
vp21=viewport(layout.pos.row=2,layout.pos.col=1,
	name="vp21",
	clip="off");

pushViewport(vp21);
grid.text(label="Model Frequency",x=unit(2, 'lines'), y=unit(0.8, "native"), just="center", rot=90)
grid.text(label=expression(-log[10](italic(p))),x=unit(2, 'lines'), y=unit(0.4, "native"), just="center", rot=90)
grid.text(label="LOD",x=unit(2, 'lines'), y=unit(0.1, "native"), just="center", rot=90)
upViewport()

##### right side of yaxis:
vp23=viewport(layout.pos.row=2,layout.pos.col=3,
name="vp23",
clip="off");

pushViewport(vp23);
grid.text(label="NAM Bayesian GWAS",x=unit(1, 'lines'), y=unit(0.8, "native"), gp=gpar(cex=0.8), just="center", rot=90)
grid.text(label="NAM Regular GWAS",x=unit(1, 'lines'), y=unit(0.4, "native"), gp=gpar(cex=0.8),just="center", rot=90)
grid.text(label="NAM QTL(JCIM)",x=unit(1, 'lines'), y=unit(0.1, "native"), gp=gpar(cex=0.8),just="center", rot=90)
upViewport()


#######################################
####-----2nd layer layout
#######################################
### GENE UNION SNP lines --------------------------------------------------
topvp <- dataViewport(
xscale=c(0, xmax), yscale=unit(c(0,1),'npc'),
extension=c(0,.05),
layout.pos.row=2,layout.pos.col=2,
name="topvp",
clip="off");

pushViewport(topvp)
grid.rect();


layout2 = grid.layout(3,1, 
	heights=unit(c(2,3,2),c('null','null','null','null','null','null'))
	)
pushViewport(viewport(layout=layout2))

### QTL --------------------------------------------------
qtlVp=dataViewport(
	xscale=c(0, xmax), yscale=c(-3,80),
	extension=c(0,.05),
	layout.pos.row=3,layout.pos.col=1,
	name="qtl",
	clip="off");

pushViewport(qtlVp);
grid.xaxis(at =ticks, label=c("chr1", "chr2", "chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10"));
grid.yaxis(at=c(10,30,50,70), label=c(10,30,50,70), gp=gpar(cex=0.8));
grid.rect();

grid.lines(x=unit(c(0, xmax), "native"), y=unit(c(qtlcutoff, qtlcutoff), "native"), gp=gpar(col="red", lty="dashed"))
#cim6big <- cim6[cim6$lod >= qtlcutoff,]
for (i in unique(cim6$CHR)){
		grid.points(x=cim6[cim6$CHR==i,]$pos, y=cim6[cim6$CHR==i,]$lod,pch=19,gp=gpar(col=colors[i],alpha=0.5,pch=2, cex=0.3));	
}

for (i in 1:nrow(qtl)){
	grid.lines(x=unit(c(qtl$LeftPos[i], qtl$RightPos[i]), "native"), y=unit(c(0.7,0.7), "npc"), gp=gpar(lwd=4, lineend="butt"));
#grid.points(x=unit(qtl$SigPos[i], "native"), y=unit(0.7,"npc"), gp=gpar(cex=0.3, col="red", fill="red"))
}

upViewport()

### GWAS pval------------------------------------------------------
gwasVp=dataViewport(
	xscale=c(0, xmax), yscale=c(0,max(gwas$logp)+5),
	extension=c(0,.05),
	layout.pos.row=2,layout.pos.col=1,
	name="gwas",
	clip="off");

pushViewport(gwasVp);
#grid.xaxis();
grid.yaxis(at=c(10,30,50,70,90), label=c(10,30,50,70,90), gp=gpar(cex=0.8));
#grid.rect();

for (i in 1:10){
	grid.points(x=gwas[gwas$CHR==i,]$pos, y=gwas[gwas$CHR==i,]$logp, pch=19, gp=gpar(col=colors[i],alpha=0.3, cex=0.2, pch=2));	
}


#tem1 <- head(gwas[order(gwas$logp, decreasing = T),], 150);
#grid.points(x=tem1$pos, y=tem1$logp, pch=19, gp=gpar(col="red", cex=0.3));



#############################################################
kas82$Select <- as.character(kas82$Select)
grid.points(x=kas82[kas82$Select=="Control",]$pos, 
	y=kas82[kas82$Select=="Control",]$logp, pch=19, gp=gpar(col="yellow", cex=0.2));

grid.points(x=kas82[kas82$Select!="Control",]$pos, 
	y=kas82[kas82$Select!="Control",]$logp, pch=19, gp=gpar(col="red", cex=0.2));
############################################################
upViewport()

### Bayesin Modelfrq------------------------------------------------------	
bayesVp =dataViewport(xscale=c(0, xmax), yscale=c(-0.01,0.8),
	extension=c(0,.05), layout.pos.row=1, layout.pos.col=1,
	name='bayes', clip='off'
);
pushViewport(bayesVp);
#grid.xaxis();
grid.yaxis(at=c(0.1,0.3,0.5,0.7), label=c(.1,0.3,0.5,0.7), gp=gpar(cex=0.8));
#grid.rect();
for(i in 1:10) {
	grid.points(x=bayes[bayes$CHR==i,]$pos, y=bayes[bayes$CHR==i,]$ModelFreq, pch=19,gp=gpar(pch=20,
	col=colors[i],alpha=0.3, cex=0.3));	
}
grid.lines(x=unit(c(0, 1), 'npc'), y=c(0,0))

#tem2 <- head(bayes[order(bayes$ModelFreq, decreasing = T),], 150);
#grid.points(x=tem2$pos, y=tem2$ModelFreq, pch=19, gp=gpar(col="red", cex=0.3));



kas82$Select <- as.character(kas82$Select)
grid.points(x=kas82[kas82$Select=="Control",]$pos, 
	y=kas82[kas82$Select=="Control",]$ModelFreq, pch=19, gp=gpar(col="yellow", cex=0.2));

grid.points(x=kas82[kas82$Select!="Control",]$pos, 
	y=kas82[kas82$Select!="Control",]$ModelFreq, pch=19, gp=gpar(col="red", cex=0.2));

for (i in 1:nrow(rnaseq[1:50,])){
  grid.lines(x=unit(c(rnaseq$pos[i], rnaseq$pos[i]), "native"), y=unit(c(0.1,0.9), "npc"), gp=gpar(lwd=1, lty=2));
  #grid.points(x=unit(qtl$SigPos[i], "native"), y=unit(0.7,"npc"), gp=gpar(cex=0.3, col="red", fill="red"))
}


upViewport()



