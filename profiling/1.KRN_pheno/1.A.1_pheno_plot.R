# Jinliang Yang
# 7.10. 2014
# phenotypic data of KRN

krn <- read.csv("~/Documents/KRN_GWAS_v3/GWAS3_proj/data/S.table9_krn_6230.csv")
dim(krn)

mean(krn$KRN)
range(krn$KRN)



krn_plot <- function(krnall=krn){
  plot(density(subset(krnall, pop=="NAMRIL")$KRN ), col="red", lwd=2, ylim=c(0,.5), main="Density Plot of KRN", xlab="KRN")
  lines(density(subset(krnall, pop=="Diallel")$KRN), col="bisque4", lwd=2)
  lines(density(subset(krnall, pop=="BxRIL")$KRN), col="cadetblue4", lwd=2)
  lines(density(subset(krnall, pop=="MxRIL")$KRN), col="chocolate4", lwd=2)
  abline(v=17.1, col="blue", lty=2, lwd=2)
  temp <- legend("topright", legend =c(" ", " ", " ", " "), title="Population", text.width=strwidth("NAMRIL"),
                 lty=1, lwd=3, col=c("red", "bisque4", "cadetblue4", "chocolate4"), xjust=1, yjust=1)
  text(temp$rect$left + temp$rect$w, temp$text$y,
       c("NAMRIL", "Diallel", "BxRIL", "MxRIL"), pos=2)
}

names(krn) <- c("geno", "pop", "FID", "KRN")
krn_plot(krnall=krn)
