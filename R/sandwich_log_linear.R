sandwich_log_linear<-function(loc,beta_opt,x,g,c,return_indiv=F,d=.Machine$double.eps^(1/3),name=NULL){
  if(anyNA(beta_opt)){return(matrix(NA,ncol=length(loc),nrow=length(loc)))}
  n<-length(x)
  n2<-length(beta_opt)
  B<-matrix(NA,nrow=n2,ncol=n2)
  for(i in 1:n2){
    beta1<-beta2<-beta_opt
    beta1[i]<-beta1[i]+d
    beta2[i]<-beta2[i]-d
    
    f1<-f_log_linear(beta=beta1,x=x,g=g,c=c)
    f1<-attr(f1,"gradient")
    f2<-f_log_linear(beta=beta2,x=x,g=g,c=c)
    f2<-attr(f2,"gradient")
    B[,i]<-(f1-f2)/2/d/n
  }
  
  M_indiv<-cbind(f_log_linear(beta_opt,x,g,c,grad=T,der_indiv=T))
  fit<-est_indiv(B=B,M_indiv=M_indiv,loc=loc,eff=beta_opt[loc],name=name)
  if(return_indiv){return(fit)}
  out<-est_vcov(fit)
  colnames(out)<-rownames(out)<-NULL
  return(out)
}

