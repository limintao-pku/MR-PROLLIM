mywald_test<-function(eff,vcov){
  if(anyNA(eff)){return(c(NA,NA))}
  chisq<-t(eff)%*%solve(vcov)%*%cbind(eff)
  p<-pchisq(chisq,df=length(eff),lower.tail=F)
  return(c(chisq,p))
}

