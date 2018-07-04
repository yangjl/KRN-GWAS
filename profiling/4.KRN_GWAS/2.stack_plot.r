# Jinliang Yang
# Purpose: to create a summarized plot stacking all the information
# 8.23.2011

# location: 129.186.85.9
####################################################################################
library("grid")
setwd("/mnt/02/yangjl/GWAS/plot")
load("data4stackplot.RData")	
# use predefined gap in the previous code: gap=10000000

# the basic plot information:
# sets colors based on colors argument.
colors <- rep(c("slateblue", "cyan4"), 5)
#xscale:
cl <- read.csv("chr_length_B73v1.csv")
cl <- newpos(cl, GAP=gap)
xmax <- max(cl$pos)

ticks <- cl$pos[1]/2
chrlines <- cl$pos[1]+gap/2
for(i in 2:10){
	ticks <- c(ticks, cl$pos[i-1] + (cl$pos[i]-cl$pos[i-1])/2)
	chrlines <- c(chrlines, cl$pos[i]+gap/2)	
}

bayes82 <- bayes82[, c('snpid','Select','chr','pos','ModelFreq')]
names(bayes82)[3:4] <- c('CHR', 'BP')
bayes82 <- newpos(bayes82, GAP=gap)


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
grid.text(label=expression(-log[10](italic(p))),x=unit(2, 'lines'), y=unit(0.85, "native"), rot=90)
grid.text(label="Model Frequency",x=unit(2, 'lines'), y=unit(0.6, "native"), just="center", rot=90)
grid.text(label=expression(-log[10](italic(p))),x=unit(2, 'lines'), y=unit(0.35, "native"), just="center", rot=90)
grid.text(label="LOD",x=unit(2, 'lines'), y=unit(0.1, "native"), just="center", rot=90)
upViewport()

##### right side of yaxis:
vp23=viewport(layout.pos.row=2,layout.pos.col=3,
name="vp23",
clip="off");

pushViewport(vp23);
grid.text(label="Extreme PI",x=unit(1, 'lines'), y=unit(0.95, "native"), gp=gpar(cex=0.8), rot=90)
grid.text(label="Off-PVP",x=unit(1, 'lines'), y=unit(0.85, "native"), gp=gpar(cex=0.8),rot=90)
grid.text(label="BSLE",x=unit(1, 'lines'), y=unit(0.75, "native"),gp=gpar(cex=0.8), rot=90)

grid.text(label="NAM Bayesian GWAS",x=unit(1, 'lines'), y=unit(0.6, "native"), gp=gpar(cex=0.8), just="center", rot=90)
grid.text(label="NAM Regular GWAS",x=unit(1, 'lines'), y=unit(0.3, "native"), gp=gpar(cex=0.8),just="center", rot=90)
grid.text(label="NAM QTL(JCIM)",x=unit(1, 'lines'), y=unit(0.1, "native"), gp=gpar(cex=0.8),just="center", rot=90)
upViewport()

##### top
vp12=dataViewport(xscale=c(0, xmax), yscale=c(0,1),
	extension=c(0,.05),
	layout.pos.row=1,layout.pos.col=2,
	name="vp12",
	clip="off");

pushViewport(vp12);
geneaty <- rep(c(1,2), times=6)
for(i in 1:nrow(gene)){
	grid.text(label=gene$Gene[i], x=unit(gene$pos[i], "native"), y=unit(geneaty[i], "lines"))
}
upViewport()

### GENE UNION SNP lines --------------------------------------------------
topvp <- dataViewport(
	xscale=c(0, xmax), yscale=unit(c(0,1),'npc'),
	extension=c(0,.05),
	layout.pos.row=2,layout.pos.col=2,
	name="topvp",
	clip="off");

pushViewport(topvp)
grid.rect();


for(i in 1:nrow(union)){
	grid.lines(x=unit(c(union$pos[i], union$pos[i]), "native"), y=unit(c(0,1), "npc"), gp=gpar(lty="dashed", col="bisque3"))
	
}
for(i in 1:9){
	grid.lines(x=unit(c(chrlines[i], chrlines[i]), "native"), y=unit(c(0,1), "npc"), gp=gpar())
}			   


#######################################
####-----2nd layer layout
#######################################
layout2 = grid.layout(6,1, 
	heights=unit(c(1,1,1,2,3,2),c('null','null','null','null','null','null'))
	)
pushViewport(viewport(layout=layout2))

### QTL --------------------------------------------------
qtlVp=dataViewport(
	xscale=c(0, xmax), yscale=c(-3,80),
	extension=c(0,.05),
	layout.pos.row=6,layout.pos.col=1,
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
	layout.pos.row=5,layout.pos.col=1,
	name="gwas",
	clip="off");

pushViewport(gwasVp);
#grid.xaxis();
grid.yaxis(at=c(10,30,50,70), label=c(10,30,50,70), gp=gpar(cex=0.8));
#grid.rect();

for (i in 1:10){
	grid.points(x=gwas[gwas$CHR==i,]$pos, y=gwas[gwas$CHR==i,]$logp, pch=19, gp=gpar(col=colors[i],alpha=0.3, cex=0.2, pch=2));	
}

gwas82$Select <- as.character(gwas82$Select)
grid.points(x=gwas82[gwas82$Select=="Control",]$pos, 
	y=gwas82[gwas82$Select=="Control",]$logp, pch=19, gp=gpar(col="yellow", cex=0.2));

grid.points(x=gwas82[gwas82$Select!="Control",]$pos, 
	y=gwas82[gwas82$Select!="Control",]$logp, pch=19, gp=gpar(col="red", cex=0.2));

upViewport()

### Bayesin Modelfrq------------------------------------------------------	
bayesVp =dataViewport(xscale=c(0, xmax), yscale=c(-0.01,0.8),
	extension=c(0,.05), layout.pos.row=4, layout.pos.col=1,
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

bayes82$Select <- as.character(bayes82$Select)
grid.points(x=bayes82[bayes82$Select=="Control",]$pos, 
	y=bayes82[bayes82$Select=="Control",]$ModelFreq, pch=19, gp=gpar(col="yellow", cex=0.2));

grid.points(x=bayes82[bayes82$Select!="Control",]$pos, 
	y=bayes82[bayes82$Select!="Control",]$ModelFreq, pch=19, gp=gpar(col="red", cex=0.2));

upViewport()

### BSLE ------------------------------------------------------	
bsleVp =dataViewport(xscale=c(0, xmax), yscale=c(-1,max(bsle$logp)+2),
	extension=c(0,.05), layout.pos.row=3, layout.pos.col=1,
	name='bsle', clip='off'
);
pushViewport(bsleVp);

for(i in 1:10) {
	grid.points(x=bsle[bsle$CHR==i,]$pos, y=bsle[bsle$CHR==i,]$logp, pch=19, gp=gpar(col=colors[i],alpha=0.5, cex=0.4));	
}
grid.yaxis(at=c(2,6,10), label=c(2,6,10), gp=gpar(cex=0.8))
grid.lines(x=unit(c(0, 1), 'npc'), y=c(0,0))

bsle$Select <- as.character(bsle$Select)
grid.points(x=bsle[bsle$Pbonf <= 0.05 & bsle$Select!="Control",]$pos, 
	y=bsle[bsle$Pbonf <= 0.05 & bsle$Select!="Control",]$logp, pch=19, gp=gpar(col="red", cex=0.4));	
grid.points(x=bsle[bsle$Pbonf <= 0.05 & bsle$Select=="Control",]$pos, 
	y=bsle[bsle$Pbonf <= 0.05 & bsle$Select=="Control",]$logp, pch=19, gp=gpar(col="yellow", cex=0.4));	

upViewport()

### OFF-PVP ------------------------------------------------------	
pvpVp =dataViewport(xscale=c(0, xmax), yscale=c(-1,6),
	extension=c(0,.05), layout.pos.row=2, layout.pos.col=1,
	name='pvp', clip='off'
);

pushViewport(pvpVp);
#grid.xaxis();
grid.yaxis(at=c(2,4), label=c(2,4), gp=gpar(cex=0.8));
#grid.rect();

for(i in unique(offpvp$CHR)) {
	grid.points(x=offpvp[offpvp$CHR==i,]$pos, y=offpvp[offpvp$CHR==i,]$logp, pch=19, gp=gpar(col=colors[i],alpha=0.5, cex=0.5));	
}

offpvp$Select.x <- as.character(offpvp$Select.x)

grid.points(x=offpvp[offpvp$Pbonf <= 0.05 & offpvp$Select.x!="Control",]$pos, 
	y=offpvp[offpvp$Pbonf <= 0.05 & offpvp$Select.x!="Control",]$logp, pch=19, gp=gpar(col="red", cex=0.4));	
grid.points(x=offpvp[offpvp$Pbonf <= 0.05 & offpvp$Select.x=="Control",]$pos, 
	y=offpvp[offpvp$Pbonf <= 0.05 & offpvp$Select.x=="Control",]$logp, pch=19, gp=gpar(col="yellow", cex=0.4));	

grid.lines(x=unit(c(0, 1), 'npc'), y=c(0,0))

upViewport()


### GERMPLASM ------------------------------------------------------	
germVp =dataViewport(xscale=c(0, xmax), yscale=c(-1,5),
	extension=c(0,.05), layout.pos.row=1, layout.pos.col=1,
	name='germ', clip='off'
);
pushViewport(germVp);
#grid.xaxis();
grid.yaxis(at=c(2,4), label=c(2,4), gp=gpar(cex=0.8));
#grid.rect();
for(i in unique(germ$CHR)) {
	grid.points(x=germ[germ$CHR==i,]$pos, y=germ[germ$CHR==i,]$logp, pch=19,gp=gpar(col=colors[i],alpha=0.5, cex=0.5));	
}

germ$Select <- as.character(germ$Select)

grid.points(x=germ[germ$Pbonf <= 0.05 & germ$Select!="Control",]$pos, 
	y=germ[germ$Pbonf <= 0.05 & germ$Select!="Control",]$logp, pch=19, gp=gpar(col="red", cex=0.4));	
grid.points(x=germ[germ$Pbonf <= 0.05 & germ$Select=="Control",]$pos, 
	y=germ[germ$Pbonf <= 0.05 & germ$Select=="Control",]$logp, pch=19, gp=gpar(col="yellow", cex=0.4));	

grid.lines(x=unit(c(0, 1), 'npc'), y=c(0,0))

upViewport()


