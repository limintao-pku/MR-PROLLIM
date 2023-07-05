outlier_test<-function(eff,se,vcov=NULL,p_cut=0.05,se_cut_k=1.5,dt=T,q_test=T){
  #Ture denotes normal SNPs
  out<-!is.na(eff)
  out_se<-rep(T,length(out))
  
  if(!is.null(se_cut_k)&sum(out)>1){
    out_se<-se<=quantile(se,c(0.75),na.rm=T)+se_cut_k*IQR(se,na.rm=T)
    if(dt){cat(sum(!out_se,na.rm=T),"SNP(s) is/are considered as outlier(s) due to large standard error.\r\n")}
    out_se[is.na(out_se)]<-F
  }
  out<-out&out_se
  
  if(q_test){
    equal_wald<-function(eff,vcov){
      cochran_q<-function(eff,se){
        loc<-which(!is.na(eff))
        eff<-eff[loc]
        se<-se[loc]
        w<-1/se^2/sum(1/se^2)
        eff_m<-sum(w*eff)
        q=sum((eff-eff_m)^2/se^2)
        return(q)
      }
      if(length(eff)==1){
        return(c(chisq=NA,p.value=1))
      }
      if(is.vector(vcov)){
        q<-cochran_q(eff,sqrt(vcov))
        p<-pchisq(q,length(eff)-1,lower.tail=F)
        return(c(chisq=q,p.value=p))
      }
      R<-cbind(rep(1,length(eff)-1),-diag(length(eff)-1))
      chisq_stat<-t(R%*%eff)%*%solve(R%*%vcov%*%t(R))%*%(R%*%eff)
      p<-pchisq(chisq_stat,length(eff)-1,lower.tail=F)
      return(c(chisq=chisq_stat,p.value=p))
    }
    ind<-function(x,loc1,loc2){
      if(is.vector(x)){
        return(x[loc1])
      }
      return(x[loc1,loc2])
    }
    
    if(is.null(vcov)){
      vcov<-se^2
    }
    
    if(sum(out)==0){stop("No valid inputs for outlier_test.")}
    if(sum(out)==1){message("The number of valid inputs is <=1. outlier_test will not work.")}
    
    while(T){
      if(equal_wald(eff[out],ind(vcov,out,out))[2]>p_cut){break}
      if(sum(out)<=2){
        message("Only two SNPs remained, and they showed heterogeneity.")
        break
      }
      loc<-which(out)
      p<-c()
      for(i in 1:length(loc)){
        out2<-out
        out2[loc[i]]<-F
        p[i]<-equal_wald(eff[out2],ind(vcov,out2,out2))[2]
      }
      out[loc[which.max(p)]]<-F
    }
  }
  return(out)
}
