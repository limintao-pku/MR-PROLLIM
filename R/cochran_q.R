cochran_q<-function(eff,se){
  loc<-which(!is.na(eff))
  eff<-eff[loc]
  se<-se[loc]
  w<-1/se^2/sum(1/se^2)
  eff_m<-sum(w*eff)
  q=sum((eff-eff_m)^2/se^2)
  return(q)
}

