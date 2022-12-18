est_proc_bi_p1.2_dl<-function(dl_data1_unsuit,fb1,se_fb1,fb1_indiv,fb1_w,parallel_trace){
  calcu_b4<-function(est_mkp,est_indiv,fb1,se_fb1,fb1_indiv,fb1_w,nonNA_loc,length_all,name){
    mycov<-function(v,m,w){
      out<-rep(NA,ncol(m))
      for(i in 1:ncol(m)){
        s<-cov(v,m[,i],use="na.or.complete")
        s<-s*sum(!is.na(v)&!is.na(m[,i]))
        if(is.na(s)){s<-0}
        s<-s/sum(!is.na(v))/sum(!is.na(m[,i]))
        out[i]<-s
      }
      out<-t(w)%*%out
      return(out)
    }
    
    m1<-est_mkp[1]
    m2<-est_mkp[2]
    k1<-est_mkp[3]
    k2<-est_mkp[4]
    p_td<-est_mkp[5]
    deno1<-(exp(k1)-1)*p_td*fb1+1
    deno2<-(exp(k2)-1)*p_td*fb1+1
    h1<-suppressWarnings(m1-log(deno1))
    h2<-suppressWarnings(m2-log(deno2))
    o1<-c(1,0,-p_td*fb1*exp(k1)/deno1,0,-fb1*(exp(k1)-1)/deno1,-p_td*(exp(k1)-1)/deno1)
    o2<-c(0,1,0,-p_td*fb1*exp(k2)/deno2,-fb1*(exp(k2)-1)/deno2,-p_td*(exp(k2)-1)/deno2)
    vcov_mkp<-est_vcov(est_indiv)
    vcov<-matrix(NA,6,6)
    vcov[1:5,1:5]<-vcov_mkp
    
    v6<-rep(NA,6)
    for(i in 1:5){
      s<-rep(NA,length_all)
      s[nonNA_loc]<-est_indiv[,i]
      v6[i]<-mycov(s,fb1_indiv,fb1_w)
    }
    v6[6]<-se_fb1^2
    vcov[,6]<-vcov[6,]<-v6
    vcov_h<-rbind(o1,o2)%*%vcov%*%cbind(o1,o2)
    
    if(anyNA(vcov_h)|anyNA(c(h1,h2))){
      warning(paste(name,": NAs detected in the indirect hp test. This SNP is removed."))
      list(eff=c(NA,NA),vcov=matrix(NA,2,2),wald_p=NA)
    }
    
    wald<-mywald_test(c(h1,h2),vcov_h)
    rownames(vcov_h)<-colnames(vcov_h)<-NULL
    return(list(eff=c(h1,h2),vcov=vcov_h,wald_p=wald[2]))
  }
  stopifnot(length(dl_data1_unsuit)!=0)
  
  out<-list()
  for(i in 1:length(dl_data1_unsuit)){
    s<-dl_data1_unsuit[[i]]
    fit<-calcu_b4(c(s$m_hat[1,],s$k_hat[1,]),s$mkp_ind,fb1=fb1,se_fb1=se_fb1,fb1_indiv=fb1_indiv,fb1_w=fb1_w,s$nonNA_loc,nrow(fb1_indiv),names(dl_data1_unsuit)[i])
    out<-c(out,list(fit))
    if(parallel_trace){my_moni("SNP",i,length(dl_data1_unsuit))}
  }
  names(out)<-names(dl_data1_unsuit)
  
  return(out)
}
