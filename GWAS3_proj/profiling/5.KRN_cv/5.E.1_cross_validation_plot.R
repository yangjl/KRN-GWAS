# Jinliang Yang
# Functions to draw the cross-validation plot



xy_rescale <- function(tem=myinput){
  library("plotrix")
  ##############################################################
  #rescale the log10 values:
  tem0 <- subset(tem, log10 <= 2)
  tem0$log10_rs <- rescale(tem0$log10, c(0, 10))
  tem1 <- subset(tem, log10 > 2 & log10 <=20)
  tem1$log10_rs <- rescale(tem1$log10, c(10, 20))
  tem2 <- subset(tem, log10 > 20 & log10 <= 50)
  tem2$log10_rs <- rescale(tem2$log10, c(20, 30))
  tem3 <- subset(tem, log10 > 50)
  tem3$log10_rs <- rescale(tem3$log10, c(30, 40))
  tem <- rbind(tem0, tem1, tem2, tem3)
  
  #rescale the log10 values:
  tem0 <- subset(tem, ModelFreq <= 0.002)
  tem0$mf_rs <- rescale(tem0$ModelFreq, c(0, 1))
  tem1 <- subset(tem, ModelFreq > 0.002 & ModelFreq <= 0.02)
  tem1$mf_rs <- rescale(tem1$ModelFreq, c(1, 2))
  tem2 <- subset(tem, ModelFreq > 0.02 & ModelFreq <= 0.2)
  tem2$mf_rs <- rescale(tem2$ModelFreq, c(2, 3))
  tem3 <- subset(tem, ModelFreq > 0.2)
  tem3$mf_rs <- rescale(tem3$ModelFreq, c(3, 4))
  
  tem <- rbind(tem0, tem1, tem2, tem3)
  #############################################
  return(tem)
}

###########
cvplot <- function(tem=tem, ck=ck){
  par(mfrow=c(1,1))
  plot(x=-100, y=-100, xlim=c(0, 43), ylim=c(-0.5, 4.5), yaxt="n",xaxt="n",
       main="", 
       xlab="Single-variant GWAS (-log10(pval))", 
       ylab="Bayesian-based GWAS (ModelFreq)")
  
  #tem9 <- subset(tem, method!="Control")
  #t1 <- subset(tem9, ModelFreq>=0.1);
  if(!is.null(ck)){
    
    points(ck$log10_rs, ck$mf_rs, col="grey", cex=2, pch=16);
    subcvd <- subset(ck, cvd == 1)
    if(nrow(subcvd) > 0 ){
      points(subcvd$log10_rs, subcvd$mf_rs, col="red", cex=1.2, pch=17)
    }
  }
  points(tem$log10_rs, tem$mf_rs, col="#0000ff80", cex=2, pch=16);
  cvd <- subset(tem, elite1 == 1 | elite2 == 1 | xbsa1 == 1 | bsle1 == 1)
  points(cvd$log10_rs, cvd$mf_rs, col="red",cex=1.2, pch=17)
  
  axis(2, at= c(0, 1, 2, 3, 4), labels= c(0, 0.002,0.02, 0.2, 1), las=1, cex=0.8)
  axis(1, at= c(0, 10, 20, 30, 40), labels= c(0, 2, 20, 50, 150), las=1, cex=0.8)
  abline(h=2, lty=2)
  abline(v=20, lty=2)
  
  
  legend("topright", pch=c(16, 16, 17), col=c("#0000ff80", "grey", "red"),
         legend=c("Informative KAVs","Control Variants", "Cross-Validated"))
  
}




