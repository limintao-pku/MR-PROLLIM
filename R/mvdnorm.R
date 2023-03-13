mvdnorm<-function(x1,x2,u1,u2,s1,s2,r){
  p<-1/(2*pi*s1*s2*sqrt(1-r^2))*exp( -1/(2-2*r^2)*
                                       ( (x1-u1)^2/s1/s1-2*r*(x1-u1)*(x2-u2)/s1/s2+(x2-u2)^2/s2/s2 ) )
  return(p)
}
