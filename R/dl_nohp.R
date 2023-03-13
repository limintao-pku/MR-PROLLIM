dl_nohp<-function(est_mkp,vcov_mkp=NULL,mkp_indiv,return_vcov=F,name,nonNA_loc,length_all){
  m1<-est_mkp[1]
  m2<-est_mkp[2]
  k1<-est_mkp[3]
  k2<-est_mkp[4]
  p_td<-est_mkp[5]
  
  fb101<-(exp(m1)-1)/(exp(k1)-1)/p_td
  fb102<-(exp(m2)-1)/(exp(k2)-1)/p_td
  
  o1<-c(exp(m1)/p_td/(exp(k1)-1), 0, -(exp(m1)-1)/p_td*exp(k1)/((exp(k1)-1)^2), 0,
        -(exp(m1)-1)/(exp(k1)-1)/p_td/p_td)
  o2<-c(0, exp(m2)/p_td/(exp(k2)-1), 0, -(exp(m2)-1)/p_td*exp(k2)/((exp(k2)-1)^2),
        -(exp(m2)-1)/(exp(k2)-1)/p_td/p_td)
  if(any(c(o1,o2)%in%c(-Inf,Inf,NA,NaN))){
    warning(paste(name,": NAs or Infs detected in dl_noph."))
    return(list(est=NA,se=NA,est_indiv=rep(NA,length_all)))
  }
  if(is.null(vcov_mkp)){
    vcov_mkp<-est_vcov(mkp_indiv)
  }
  vcov<-rbind(o1,o2)%*%vcov_mkp%*%cbind(o1,o2)
  
  if(return_vcov){
    return(list(eff=c(fb101,fb102),vcov=vcov))
  }
  
  fb1_indiv<-scale(mkp_indiv,scale=F)%*%cbind(o1,o2)+matrix(c(fb101,fb102),nrow=nrow(mkp_indiv),ncol=2,byrow=T)
  out<-meta_eff2(c(fb101,fb102),vcov,fb1_indiv)
  est_indiv<-rep(NA,length_all)
  est_indiv[nonNA_loc]<-out$eff_m_indiv
  
  return(list(est=out$eff_m,se=out$se_m,est_indiv=est_indiv))
}
