equal_wald<-function(eff,vcov){
    if(is.vector(vcov)){
      q<-cochran_q(eff,sqrt(vcov))
      p<-pchisq(q,length(eff)-1,lower.tail=F)
      return(c(chisq=q,p.value=p))
    }
    if(length(eff)==1){
      return(c(chisq=NA,p.value=1))
    }
    R<-cbind(rep(1,length(eff)-1),-diag(length(eff)-1))
    chisq_stat<-t(R%*%eff)%*%solve(R%*%vcov%*%t(R))%*%(R%*%eff)
    p<-pchisq(chisq_stat,length(eff)-1,lower.tail=F)
    return(c(chisq=chisq_stat,p.value=p))
  }

