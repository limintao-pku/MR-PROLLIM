f_log_linear<-function(beta,x,g,c,grad=T,der_indiv=F){
  d<-cbind(g,c,1)
  px<-exp(d%*%beta)
  px[px>=1]<-1
  out<-sum(-log(px^x*(1-px)^(1-x)))
  if(grad){
    if(der_indiv){
      return(-c((x/px+(x-1)/(1-px))*px)*d)
    }
    grad_out<-t((x/px+(x-1)/(1-px))*px)%*%d
    attr(out,'gradient')=-grad_out
  }
  return(out)         
}

