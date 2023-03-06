est_k_prior<-function(k_matr,s_matr,s_list=NULL,ss_list=NULL,p_cut,u1_sp=NULL,p0_sp=NULL,start=c(0,0,-1,-1,1),p0_start=NULL,mc.cores=1,
                      dt=T,PSOCK=F,cpp=T,nlminb_control=list(rel.tol=1e-12,sing.tol=1e-12,step.min=0.8,eval.max=300,iter.max=300),...){
  f_est_k_prior<-function(beta,k_matr,s_list,ss_list,p_cut,u1_sp=NULL,p0_sp=NULL,cpp=T){
    if(is.null(p0_sp)){
      p0<-1/(1+1/exp(beta[1]))
    }else{
      stopifnot(p0_sp==0)
      p0<-p0_sp
    }
    
    if(is.null(u1_sp)){
      u1<-beta[2]
      u2<-2*u1
    }else{
      u1<-u1_sp
      u2<-2*u1
    }
    
    s12<-exp(beta[3])
    s22<-exp(beta[4])
    rho<-1/(1+1/exp(beta[5]))
    sigma<-matrix(c(s12,sqrt(s12*s22)*rho,sqrt(s12*s22)*rho,s22),2)
    
    mychisq<-qchisq(1-p_cut,2)
    if(cpp){
      if(!is.null(p0_sp)){p0_sp<- -1}else{p0_sp<-2}
      p<-f_est_k_prior_loop(k_matr,s_list,ss_list,mychisq,p0_sp,p0,p_cut,u1,u2,sigma)
    }else{
      p<-rep(NA,nrow(k_matr))
      if(is.null(p0_sp)){
        for(i in 1:nrow(k_matr)){
          p1<-my_pmvnEll(chisq=mychisq,sigma=s_list[[i]]+sigma,mu=c(u1,u2),e=ss_list[[i]],x0=c(0,0),lower.tail=F)
          #p1<-1-mean(mvdnorm3(list_r[[i]],c(u1,u2),s_list[[i]]+sigma))*area[i]
          p00<-p0*p_cut/(p0*p_cut+(1-p0)*p1)
          p[i]<-(p00*mvdnorm3(k_matr[i,],c(0,0),s_list[[i]])/p_cut)+
            ((1-p00)*mvdnorm3(k_matr[i,],c(u1,u2),s_list[[i]]+sigma)/p1)
        }
      }else{
        for(i in 1:nrow(k_matr)){
          p1<-my_pmvnEll(chisq=mychisq,sigma=s_list[[i]]+sigma,mu=c(u1,u2),e=ss_list[[i]],x0=c(0,0),lower.tail=F)
          #p1<-1-mean(mvdnorm3(list_r[[i]],c(u1,u2),s_list[[i]]+sigma))*area[i]
          p[i]<-(mvdnorm3(k_matr[i,],c(u1,u2),s_list[[i]]+sigma)/p1)
        }
      }
    }
    
    -sum(log(p))
  }
  cut_num<-function(x,n){
    if(length(x)<n){n<-length(x)}
    out<-list()
    n1<-length(x)
    x1<-cut(1:n1,breaks=seq(0,n1,length.out=n+1),labels=1:n)
    for(i in 1:n){
      out[[i]]<-x[which(x1==i)]
    }
    return(out)
  }
  my_trans<-function(beta,u1_sp=NULL,p0_sp=NULL){
    if(is.null(p0_sp)){
      p0<-1/(1+1/exp(beta[1]))
    }else{
      p0<-p0_sp
    }
    if(is.null(u1_sp)){
      u1<-beta[2]
      u2<-2*u1
    }else{
      u1<-u1_sp
      u2<-2*u1
    }
    s12<-exp(beta[3])
    s22<-exp(beta[4])
    rho<-1/(1+1/exp(beta[5]))
    c(p0=p0,u1=u1,s12=s12,s22=s22,rho=rho)
  }
  
  if(is.null(s_list)){
    s_list<-list()
    ss_list<-list()
    for(i in 1:nrow(s_matr)){
      s_list[[i]]<-matrix(c(s_matr[i,1],s_matr[i,2],s_matr[i,2],s_matr[i,3]),2)
      ss_list[[i]]<-solve(s_list[[i]])
    }
  }
  
  if(dt){cat("Start nlminb optimizing for parameters of k_prior.\r\n")}
  
  if(is.null(p0_start)){
    if(T){
      fit<-nlminb2nlm(nlminb2(f=f_est_k_prior,p=start,k_matr=k_matr,s_list=s_list,ss_list=ss_list,p_cut=p_cut,u1_sp=u1_sp,p0_sp=p0_sp,cpp=cpp,
                              control=nlminb_control))
      if(fit$code!=0){message("Abnormal nlminb code in est_k_prior.")}
    }else{
      fit<-suppressWarnings(nlm(f=f_est_k_prior,p=start,k_matr=k_matr,s_list=s_list,ss_list=ss_list,p_cut=p_cut,u1_sp=u1_sp,p0_sp=p0_sp,cpp=cpp,...))
    }
  }else{
    stopifnot(length(p0_start)>=1)
    
    my_task<-function(p0_start,start){
      fit<-list()
      myrc<-NULL
      for(i in 1:length(p0_start)){
        start[1]<-log(p0_start[i]/(1-p0_start[i]))
        fit[[i]]<-withCallingHandlers({tryCatch({nlminb2nlm(nlminb2(f=f_est_k_prior,p=start,k_matr=k_matr,s_list=s_list,
                                                                    ss_list=ss_list,p_cut=p_cut,u1_sp=u1_sp,p0_sp=p0_sp,cpp=cpp,
                                                                    control=nlminb_control))},error=function(e){list(minimum=Inf,error=e)})},warning=function(w){myrc<<-c(myrc,w$message);invokeRestart("muffleWarning")})
      }
      return(fit)
    }

    mynum<-cut_num(p0_start,mc.cores)
    mycheck<-"pass"
    exec_base_func<-function(x){
      suppressWarnings(library(MRprollim,quietly=T))
    }
    myfit<-withCallingHandlers({my_parallel(X=mynum,FUN=my_task,start=start,mc.cores=mc.cores,PSOCK=PSOCK,dt=dt,
                                            print_message=F,export_parent=T,exec_base_func=exec_base_func,seed=NULL)},warning=function(w){mycheck<<-w})
    if((!identical(mycheck,"pass"))&mc.cores!=1){
      stop(as.character(mycheck))
    }
    mybind_list_parallel<-function(list_parallel){
      out1<-list()
      for(i in 1:length(list_parallel)){
        out1<-c(out1,list_parallel[[i]])
      }
      return(out1)
    }
    myfit<-mybind_list_parallel(myfit)
    minimum<-unlist(lapply(myfit,FUN=function(x){x$minimum}))
    if(all(is.infinite(minimum))){stop("Infinities or errors detected for all values in p0_start.")}
    fit<-myfit[[which.min(minimum)]]
    if(fit$code!=0){message("Abnormal nlminb code in est_k_prior.")}
  }
  
  par<-as.list(my_trans(fit$estimate,u1_sp=u1_sp,p0_sp=p0_sp))
  
  return(list(par=par,nlminb_out=fit))
}
