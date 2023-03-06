mynlminb<-function(f,p,...,control,check_delta_k=c(1e-3,1e-4),name=NULL){
  e<-f(p,...)
  if(e%in%c(Inf,-Inf,NA,NaN)){
    stop(paste0(name,": inappropriate initial values."))
  }
  
  nlminb_scale<-control$nlminb_control$scale
  if(is.null(nlminb_scale)){nlminb_scale<-1}
  control$nlminb_control[which(names(control$nlminb_control)=="scale")]<-NULL
  nlminb_c_list<-control$nlminb_control
  
  nlm_c_list_org<-list(hessian=F,fscale=1,print.level=0,ndigit=12,
                       gradtol=1e-10,stepmax=2,steptol=1e-10,iterlim=300,check.analyticals=F)
  nlm_c_list<-match.list(control$nlm_control,nlm_c_list_org)
  
  fit<-tryCatch({nlminb2nlm(nlminb2(f=f,p=p,...,gradient=T,scale=nlminb_scale,control=nlminb_c_list))},error=function(e){e})
  s1_fit<-F
  if(!"error"%in%class(fit)){
    fit_message<-fit$message
    if(fit$code==0|(stringr::str_detect(fit_message,"(singular convergence)|(false convergence)"))){
      s1_fit<-T
    }else{
      warning(paste0(name,": nlminb returns an abnormal message: ",fit_message))
    }
  }else{
    fit_message<-"nlminb gives an error"
    warning(paste0(name,": nlminb gives an error."))
  }
  
  if(s1_fit){s2_start<-fit$estimate}else{s2_start<-p}
  k<-list(rep(1,length(p)),rep(0.999,length(p)),rep(1.001,length(p)))
  fit2_list<-list()
  
  for(i in 1:length(nlm_c_list$steptol)){
    for(j in 1:3){
      if(i==1&j==1){check.analyticals<-nlm_c_list$check.analyticals}else{check.analyticals<-F}
      fit2<-tryCatch({suppressWarnings(nlm(f=f,p=s2_start*k[[j]],...,hessian=nlm_c_list$hessian,
                                           fscale=nlm_c_list$fscale,print.level=nlm_c_list$print.level,ndigit=nlm_c_list$ndigit,
                                           gradtol=nlm_c_list$gradtol,stepmax=nlm_c_list$stepmax,steptol=nlm_c_list$steptol[i],
                                           iterlim=nlm_c_list$iterlim,check.analyticals=check.analyticals))},
                     error=function(e){e})
      if(!"error"%in%class(fit2)){fit2_list<-c(fit2_list,list(fit2))}
      if(tryCatch({fit2$code%in%c(1)},error=function(e){F})){break}
    }
    if(tryCatch({fit2$code%in%c(1)},error=function(e){F})){break}
  }
  
  if((!tryCatch({fit2$code%in%c(1)},error=function(e){F}))&length(fit2_list)>0){
    fit2_gr<-unlist(lapply(fit2_list,FUN=function(x){sum(abs(x$gradient))}))
    fit2_gr[is.na(fit2_gr)]<-Inf
    fit2<-fit2_list[[which.min(fit2_gr)]]
  }
  
  if("error"%in%class(fit2)){
    warning(paste0(name,": ","nlm gives an error; this SNP will be removed."))
    return(list(estimate=rep(NA,length(p)),message=paste0("nlminb message: ",fit_message)))
  }
  
  if(anyNA(fit2$gradient)){
    fit2$code<-"NA_gradient"
  }
  fit2$message<-paste0("nlminb message: ",fit_message,"; ","nlm code: ",fit2$code)
  
  if(s1_fit){s2_code_acc<-c(1,2,3)}else{s2_code_acc<-c(1,2)}
  if(fit2$code%in%s2_code_acc){fit2$code<-0}else{fit2$code<-1}
  
  if(fit2$code%in%c(0)){
    for(i in 1:length(check_delta_k)){
      delta<-max(fit2$estimate[1]*check_delta_k[i],check_delta_k[length(check_delta_k)])
      est1<-est2<-fit2$estimate
      est1[1]<-est1[1]+delta
      est2[1]<-est2[1]-delta
      l1<-f(est1,...)
      l2<-f(est2,...)
      grad1<-attr(l1,"gradient")[1]
      grad2<-attr(l2,"gradient")[1]
      if(!anyNA(c(grad1,grad2))){break}
    }
    if(anyNA(c(grad1,grad2))|grad1*grad2>=0){
      warning(paste0(name,": ","the estimate for the first parameter may be inappropriate; this SNP will be removed."))
      fit2$estimate<-rep(NA,length(fit2$estimate))
    }
  }else{
    warning(paste0(name,": ","nlm returns an abnormal message: ",fit2$message,"; this SNP will be removed."))
    fit2$estimate<-rep(NA,length(fit2$estimate))
  }
  return(fit2)
}
