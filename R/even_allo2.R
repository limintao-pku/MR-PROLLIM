even_allo2<-function(n,cores){
  x<-even_allo(n,cores)
  out<-list()
  for(i in 1:cores){
    out[[i]]<-list(x[i],i)
  }
  return(out)
}

