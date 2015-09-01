######################################################################
simdrift <- function(raf=raf, reps=1000, pop1=4000, pop2=300, cycle=30){
  raf$pval <- 1;
  for(i in 1:nrow(raf)){
    myfrq <- raf$C0[i];
    diff <- raf$C30LE[i] - raf$C30SE[i];
    
    test1 <- replicate(reps, drift(pop1=pop1, pop2=pop2, frq=myfrq, cycle=cycle));
    test2 <- replicate(reps, drift(pop1=pop1, pop2=pop2, frq=myfrq, cycle=cycle));
    
    raf$pval[i] <- (sum( abs(test1-test2) > abs(diff) & (test1-test2)*diff >0 )+1)/reps;  
  }
  return(raf)
}

####################
drift <- function(pop1=4000, pop2=300, frq=0.1, cycle=30) {
  #popl = 200 # evulated popl. size
  #frq = 0.1 # c0 allele frq
  for(i in 1:cycle){
    a <- sample(c(1,0), 2*pop1, replace = TRUE, prob=c(frq, 1-frq))
    frq <- sum(a)/(2*pop1)
    #print(frq)
    b <- sample(c(1,0), 2*pop2, replace=TRUE, prob=c(frq, 1-frq));
    frq <- sum(b)/(2*pop2)    
  }
  return(frq)
}