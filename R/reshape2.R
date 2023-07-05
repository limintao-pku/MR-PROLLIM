reshape2<-function(matr,n){
  out<-matrix(NA,nrow=nrow(matr)*n,ncol=ncol(matr))
  for(i in 1:n){
    out[seq(i,nrow(out),by=n),]<-matr
  }
  return(out)
}

