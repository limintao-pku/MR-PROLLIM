my_eigen<-function(x){
  z<-.Internal(La_rs(x,F))
  ord<-rev(seq_along(z$values))
  list(values=z$values[ord],vectors=z$vectors[,ord])
}

