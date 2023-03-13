summary.mrp<-function(est_out,sd,robust_u1=T,s1_cut_k=0.01,p12_hat_cut=0.5){
  trans_proc_p3.2<-function(beta_norm,p1_sp=NULL,p2_sp=NULL,r_sp=NULL,model_u2=T,est_type="b"){
    if(anyNA(beta_norm)){return(c(b1=NA,p1=NA,p2=NA,u1=NA,s1=NA,r=NA,u2=NA))}
    mod_f<-function(x,r=F){
      if(is.na(x)){return(NA)}
      if(x>1){return(1)}
      if(r){
        if(x< -1){return(-1)}
      }else{
        if(x<0){return(0)}
      }
      x
    }
    names(beta_norm)<-NULL
    if(est_type=="b"){
      b1<-suppressWarnings(-log(1-beta_norm[1]))
      if(is.na(b1)){
        #warning("The estimated b1_td is > 1.")
        b1<-Inf
      }
    }else{
      b1<-beta_norm[1]
    }
    if(is.null(p1_sp)){p1<-mod_f(beta_norm[2])}else{p1<-p1_sp}
    if(is.null(p2_sp)){p2<-mod_f(beta_norm[3])}else{p2<-p2_sp}
    u1<-beta_norm[4]
    s1<-beta_norm[5]
    if(s1<0){s1<-0}
    if(is.null(r_sp)){r<-tanh(beta_norm[6])}else{r<-r_sp}
    if(model_u2){u2<-beta_norm[7]}else{u2<-0}
    return(c(b1=b1,p1=p1,p2=p2,u1=u1,s1=s1,r=r,u2=u2))
  }
  get_model<-function(p1_sp=NULL,p2_sp=NULL,r_sp=NULL,model_u2=T,Egger_info){
    if(identical(p1_sp,NULL)&identical(p2_sp,NULL)){
      return("Full model")
    }
    if(identical(p1_sp,NULL)&identical(p2_sp,1)){
      return("Reduced model")
    }
    if(identical(p1_sp,NULL)&identical(p2_sp,0)){
      return("Zero-correlation model")
    }
    if(identical(p1_sp,1)&identical(p2_sp,0)&identical(Egger_info,"Egger")){
      return("Egger model")
    }
    if(identical(p1_sp,1)&identical(p2_sp,0)&identical(Egger_info,"intercept")){
      return("Intercept model")
    }
    return("Unknown model")
  }
  trans_b1<-function(x,fb1=F){
    if(is.null(x)){return(x)}
    if(is.na(x)){return(x)}
    if(fb1){
      out<- suppressWarnings(-log(1-x))
      if(is.na(out)){
        #warning("The 95% confidence interval gives an fb1 that can not be transformed to b1. Inf is used instead.")
        out<-Inf
      }
    }else{out<-x}
    return(out)
  }
  
  if(!identical(class(est_out),"MR-PROLLIM output")){
    stop("est_out should be an MR-PROLLIM output.")
  }
  
  stopifnot(sd>0)
  
  #p12 and mmqr
  out<-NULL
  if(identical(est_out$model,"dll")){m<-T}else{m<-F}
  if(m|!is.null(est_out$parameter$t_b1)){sd<-1}
  
  if(!is.null(est_out[["IVW_suit"]])){
    n<-length(est_out[["IVW_suit"]]$SNP_name)
    s<-est_out[["IVW_suit"]]
    if(n>0){
      eff<-s$eff_m/sd
      se<-s$se_m/sd
      z<-eff/se
      p.value<-2*pnorm(-abs(z))
      lower<-eff-se*qnorm(0.975)
      upper<-eff+se*qnorm(0.975)
      d<-data.frame(eff=trans_b1(eff,m),se_norm=se,z_norm=z,p.value=p.value,lower=trans_b1(lower,m),upper=trans_b1(upper,m),n_SNP=n)
      rownames(d)<-"IVW_suit"
      out<-rbind(out,d)
    }
  }
  if(!is.null(est_out[["IVW_all"]])){
    n<-length(est_out[["IVW_all"]]$SNP_name)
    s<-est_out[["IVW_all"]]$result
    if(n>0){
      for(j in 1:nrow(s)){
        eff<-s[j,1]/sd
        se<-s[j,2]/sd
        z<-eff/se
        p.value<-2*pnorm(-abs(z))
        lower<-eff-se*qnorm(0.975)
        upper<-eff+se*qnorm(0.975)
        d<-data.frame(eff=trans_b1(eff,m),se_norm=se,z_norm=z,p.value=p.value,lower=trans_b1(lower,m),upper=trans_b1(upper,m),n_SNP=n)
        rownames(d)<-rownames(s)[j]
        out<-rbind(out,d)
      }
    }
  }
  get_data<-function(v,n_IV,name,sd,m){
    v[1]<-v[1]/sd
    v[2]<-v[2]/sd
    z<-v[1]/v[2]
    p.value<-2*pnorm(-abs(z))
    lower<-v[1]-v[2]*qnorm(0.975)
    upper<-v[1]+v[2]*qnorm(0.975)
    out<-(data.frame(eff=trans_b1(v[1],m),se_norm=v[2],z_norm=z,p.value=p.value,
                     lower=trans_b1(lower,m),upper=trans_b1(upper,m),n_SNP=n_IV))
    rownames(out)<-name
    return(out)
  }
  if(!is.null(est_out[["median"]])){
    s<-est_out[["median"]]
    out<-rbind(out,
               get_data(s$simple_median,length(s$SNP_name),"simple_median",sd,m),
               get_data(s$weighted_median,length(s$SNP_name),"weighted_median",sd,m))
  }
  if(!is.null(est_out[["mode"]])){
    s<-est_out[["mode"]]
    out<-rbind(out,
               get_data(s$simple_mode,length(s$SNP_name),"simple_mode",sd,m),
               get_data(s$weighted_mode,length(s$SNP_name),"weighted_mode",sd,m))
  }
  if(!is.null(est_out[["IVW_outlier_removal"]])){
    s<-est_out[["IVW_outlier_removal"]]
    out<-rbind(out,
               get_data(c(s$eff_m,s$se_m),length(s$SNP_name),"IVW_outlier_removal",sd,m))
  }
  if(!is.null(est_out[["add_random"]])){
    s<-est_out[["add_random"]]
    out<-rbind(out,
               get_data(c(s$eff_m,s$se_m),length(s$SNP_name),"add_random",sd,m))
  }
  if(!is.null(est_out[["IVW_simple"]])){
    s<-est_out[["IVW_simple"]]
    out<-rbind(out,
               get_data(c(s$eff_m,s$se_m),length(s$SNP_name),"IVW_simple",sd,m))
  }
  
  if(!is.null(out)){return(out)}
  
  #p3
  est_list<-est_out$estimate
  par<-est_out$parameter
  se_asy<-sqrt(diag(est_list$vcov))
  names(se_asy)<-colnames(est_list$vcov)
  se_all<-rep(0,7)
  names(se_all)<-names(est_list$beta_norm)
  se_asy<-match.list(se_asy,se_all,T)
  if(robust_u1){
    se_u1<-est_list$se_u1_sandwich
    if(is.na(se_u1)){se_u1<-0}
    se_asy[4]<-se_u1
  }
  se_boot<-as.numeric(sqrt(diag(est_list$vcov_boot)))
  if(identical(se_boot,numeric(0))){se_boot<-NA}
  beta<-trans_proc_p3.2(est_list$beta_norm,p1_sp=par[["p1_sp"]],p2_sp=par[["p2_sp"]],r_sp=par[["r_sp"]],model_u2=par[["model_u2"]],est_type=par[["est_type"]])
  out_message<-NULL
  if(median(sqrt(est_out$data$m_sigma[,1]))*s1_cut_k>beta["s1"]){out_message<-c(out_message,"s1_hat may be too small.")}
  model<-get_model(p1_sp=par[["p1_sp"]],p2_sp=par[["p2_sp"]],r_sp=par[["r_sp"]],model_u2=par[["model_u2"]],Egger_info=par[["Egger_info"]])
  p.value_asy<-2*pnorm(-abs(est_list$beta_norm/se_asy))
  p.value_asy[se_asy==0]<-0
  p.value_boot<-2*pnorm(-abs(est_list$beta_norm/se_boot))
  p.value_boot[se_boot==0]<-0
  lower_asy<-trans_proc_p3.2(est_list$beta_norm-se_asy*qnorm(0.975),p1_sp=par[["p1_sp"]],p2_sp=par[["p2_sp"]],r_sp=par[["r_sp"]],model_u2=par[["model_u2"]],est_type=par[["est_type"]])
  upper_asy<-trans_proc_p3.2(est_list$beta_norm+se_asy*qnorm(0.975),p1_sp=par[["p1_sp"]],p2_sp=par[["p2_sp"]],r_sp=par[["r_sp"]],model_u2=par[["model_u2"]],est_type=par[["est_type"]])
  lower_boot<-trans_proc_p3.2(est_list$beta_norm-se_boot*qnorm(0.975),p1_sp=par[["p1_sp"]],p2_sp=par[["p2_sp"]],r_sp=par[["r_sp"]],model_u2=par[["model_u2"]],est_type=par[["est_type"]])
  upper_boot<-trans_proc_p3.2(est_list$beta_norm+se_boot*qnorm(0.975),p1_sp=par[["p1_sp"]],p2_sp=par[["p2_sp"]],r_sp=par[["r_sp"]],model_u2=par[["model_u2"]],est_type=par[["est_type"]])
  
  if(!is.null(est_out$data$boot_data)){
    empirical_ci<-t(apply(est_out$data$boot_data,2,FUN=function(x){quantile(x,c(0.025,0.975))}))
    empirical_ci[,1]<-trans_proc_p3.2(empirical_ci[,1],p1_sp=par[["p1_sp"]],p2_sp=par[["p2_sp"]],r_sp=par[["r_sp"]],model_u2=par[["model_u2"]],est_type=par[["est_type"]])
    empirical_ci[,2]<-trans_proc_p3.2(empirical_ci[,2],p1_sp=par[["p1_sp"]],p2_sp=par[["p2_sp"]],r_sp=par[["r_sp"]],model_u2=par[["model_u2"]],est_type=par[["est_type"]])
  }else{
    empirical_ci=NULL
  }
  sum_data<-cbind(est=beta,se_asy_norm=se_asy,p.value_asy=p.value_asy,lower_asy=lower_asy,upper_asy=upper_asy,se_boot_norm=se_boot,p.value_boot=p.value_boot,lower_boot=lower_boot,upper_boot=upper_boot)
  rownames(sum_data)<-names(beta)
  loc<-which(stringr::str_detect(rownames(sum_data),"b1|u1"))
  for(i in 1:length(loc)){
    sum_data[loc[i],1:2]<-sum_data[loc[i],1:2]/sd
    sum_data[loc[i],4:6]<-sum_data[loc[i],4:6]/sd
    sum_data[loc[i],8:9]<-sum_data[loc[i],8:9]/sd
    if(!is.null(empirical_ci)){
      empirical_ci[loc[i],]<-empirical_ci[loc[i],]/sd
    }
  }
  
  if(beta[2]*beta[3]>p12_hat_cut){
    out_message<-c(out_message,"A large p1_hat*p2_hat is detected.")
    p_u1<-est_list$p_u1_sandwich
    if(!is.na(p_u1)){
      out_message<-c(out_message,paste("Sandwich P value for u1 = 0 according to the current model is",p_u1))
      if(!identical(model,"Reduced model")){
        out_message<-c(out_message,"A more precise P value for u1 = 0 can be obtained from the reduced model.")
      }
    }
  }
  
  if(identical(par[["est_type"]],"b")){
    if(sum_data[1,1]>10|any(is.infinite(sum_data[1,4:5]))){
      out_message<-c(out_message,"A large b1_hat or infinity is detected. Recoding the binary exposure may help.")
    }
    if(sum_data[1,1]>10&!par[["t_b1"]]){
      out_message<-c(out_message,"A large b1_hat is detected. Set control_p3$t_b1 = TRUE is recommended.")
    }
  }
  
  sum_data[c(2,3,5,6),c(3,7)]<-NA
  if(all(is.na(sum_data[,6]))){
    sum_data<-sum_data[,1:5]
  }
  
  if(is.null(empirical_ci)){
    return(list(sum_data=sum_data,model=model,message=out_message))
  }
  return(list(sum_data=sum_data,empirical_ci=empirical_ci,model=model,message=out_message))
}
