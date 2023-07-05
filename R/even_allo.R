even_allo<-function(n,m){
  e<-floor(n/m)
  out<-rep(e,m)
  if(n%%m>0){
    out[1:(n%%m)]<-e+1
  }
  return(out)
}

