mynlminb<-function(f,p,...,control=list(rel.tol=1e-12,sing.tol=1e-12,step.min=0.8,eval.max=300,iter.max=300),check=T,allow_error=T,check_delta=c(1e-2,1e-3,1e-4),return_na=T,name=NULL){
  if(check){
    e<-f(p,...)
    if(e%in%c(Inf,-Inf,NA,NaN)){
      stop(paste0("Inappropriate initial values for ",paste(name,collapse=", "),"."))
    }
    fit<-tryCatch({nlminb2nlm(nlminb2(f=f,p=p,...,control=control,gradient=T))},error=function(e){e})
    if("error"%in%class(fit)){
      if(T){
        warning("Errors occured in nlminb. Retrying with nlm.")
        fit<-suppressWarnings(nlm(f=f,p=p,...,iterlim=300,gradtol=1e-8,steptol=1e-6,stepmax=2))
        if(anyNA(fit$gradient)){fit$code<-"NA_gradient"}
        fit$message<-paste("nlm code",fit$code)
        if(fit$code%in%c(1,2)){fit$code<-0}
      }
    }
    if(!fit$code%in%c(0)){
      if(!stringr::str_detect(fit$message,"nlm code")){
        warning("Abnormal nlminb code detected. Retrying with nlm.")
        fit<-suppressWarnings(nlm(f=f,p=p,...,iterlim=300,gradtol=1e-8,steptol=1e-6,stepmax=2))
        if(anyNA(fit$gradient)){
          fit$code<-"NA_gradient"
        }
        fit$message<-paste("nlm code",fit$code)
        if(fit$code%in%c(1,2)){fit$code<-0}
      }
      if(!fit$code%in%c(0)){
        if(!allow_error){stop(paste(paste(name,collapse=", "),":","Abnormal nlm code:",fit$message))}else{
          warning(paste(paste(name,collapse=", "),":","Abnormal optimization code:",fit$message))
        }
        if(return_na){fit$estimate<-rep(NA,length(fit$estimate))}
        return(fit)
      }
    }
    for(i in 1:length(check_delta)){
      delta<-max(fit$estimate[1]*check_delta[i],check_delta[length(check_delta)])
      est1<-est2<-fit$estimate
      est1[1]<-est1[1]+delta
      est2[1]<-est2[1]-delta
      l1<-f(est1,...)
      l2<-f(est2,...)
      grad1<-attr(l1,"gradient")[1]
      grad2<-attr(l2,"gradient")[1]
      if(!anyNA(c(grad1,grad2))){break}
    }
    if(anyNA(c(grad1,grad2))){
      if(!allow_error){stop(paste(paste(name,collapse=", "),":","The estimate for the first parameter may be inappropriate"))}else{
        warning(paste(paste(name,collapse=", "),":","The estimate for the first parameter may be inappropriate"))
      }
      if(return_na){fit$estimate<-rep(NA,length(fit$estimate))}
    }else{
      if(!grad1*grad2<0){
        if(!allow_error){stop(paste(paste(name,collapse=", "),":","The estimate for the first parameter may be inappropriate"))}else{
          warning(paste(paste(name,collapse=", "),":","The estimate for the first parameter may be inappropriate"))
        }
        if(return_na){fit$estimate<-rep(NA,length(fit$estimate))}
      }
    }
    return(fit)
  }else{
    fit<-nlminb2nlm(nlminb2(f=f,p=p,...,control=control,gradient=T))
    return(fit)
  }
}

