mvdnorm3<-function(x,u,sigma){
  if(is.matrix(x)){x<-t(x)}
  if(is.matrix(u)){u<-t(u)}
  
  c<-cbind(x-u)
  x<-(t(solve(sigma))%*%c)*c
  y<-.colSums(x,nrow(x),ncol(x))
  exp(y*(-0.5))/sqrt( ((2*pi)^nrow(c))*det(sigma) )
}

