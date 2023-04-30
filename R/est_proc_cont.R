est_proc_cont<-function(x,y,g,c=NULL,c_inherit=T,dt=T,mc.cores=1,PSOCK=F,parallel_trace=F,
                        est_type=c("p3","p2","p1","me_mo_q_re"),
                        twosample_data=NULL,start=NULL,snp_exp_check=T,p_snp="ask",
                        cd=T,max_unique=10,cd_g_code=T,
                        n_min_limit=100,
                        control_limit_c=NULL,
                        control_p12=NULL,data_dhp=NULL,data_nohp=NULL,
                        control_p3=NULL,data_p3_k=NULL,data_p3=NULL,data_p3_opt=NULL,
                        control_est_k_prior=NULL,control_est_k_post=NULL,control_global_search=NULL,
                        control_me_mo_q_re=NULL,data_me_mo_q_re=NULL,...
){
  c_cat<-F#legacy
  stopifnot(is.matrix(g))
  est_type<-match.arg(est_type)
  
  #functions
  ncol2<-function(x){
    if(is.null(x)){return(0)}else{return(ncol(x))}
  }
  update_twosample_data<-function(data,loc,c_inherit){
    if(is.null(data)){return(data)}
    if(c_inherit){
      return(list(g=myselect(data$g,loc),c=data$c))
    }else{
      return(list(g=myselect(data$g,loc),c=data$c[loc]))
    }
  }
  update_dum_loc<-function(dum_loc,loc,c_inherit){
    if(is.null(dum_loc)){return(dum_loc)}
    if(identical(dum_loc,"auto")){return(dum_loc)}
    if(c_inherit){
      return(dum_loc)
    }else{
      return(dum_loc[loc])
    }
  }
  
  #check ...
  check_genoud_par<-function(MemoryMatrix=-pi,boundary.enforcement=-pi,
                             BFGS=-pi,hessian=-pi,unif.seed=-pi,int.seed=-pi,share.type=-pi,instance.number=-pi,output.path=-pi,
                             output.append=-pi,project.path=-pi,P1=0,P2=-pi,P3=-pi,P4=-pi,P5=-pi,P6=-pi,P7=100, 
                             P8=-pi,P9=-pi,P9mix=-pi,BFGSfn=-pi,BFGShelp=-pi,optim.method=-pi,debug=-pi,return_list=F){
    if(return_list){
      org<-list(MemoryMatrix=-pi,boundary.enforcement=-pi,
                BFGS=-pi,hessian=-pi,unif.seed=-pi,int.seed=-pi,share.type=-pi,instance.number=-pi,output.path=-pi,
                output.append=-pi,project.path=-pi,P1=-pi,P2=-pi,P3=-pi,P4=-pi,P5=-pi,P6=-pi,P7=-pi, 
                P8=-pi,P9=-pi,P9mix=-pi,BFGSfn=-pi,BFGShelp=-pi,optim.method=-pi,debug=-pi)
      now<-list(MemoryMatrix=MemoryMatrix,boundary.enforcement=boundary.enforcement,
                BFGS=BFGS,hessian=hessian,unif.seed=unif.seed,int.seed=int.seed,share.type=share.type,instance.number=instance.number,output.path=output.path,
                output.append=output.append,project.path=project.path,P1=P1,P2=P2,P3=P3,P4=P4,P5=P5,P6=P6,P7=P7, 
                P8=P8,P9=P9,P9mix=P9mix,BFGSfn=BFGSfn,BFGShelp=BFGShelp,optim.method=optim.method,debug=debug)
      loc<-NA
      for(i in 1:length(now)){
        loc[i]<-!identical(now[[i]],org[[i]])
      }
      return(now[loc])
    }
    return("pass")
  }
  check_genoud_par(...,return_list=F)
  
  #check data
  if(cd){
    z<-check_data(x=x,y=y,g=g,c=c,c_inherit=c_inherit,type="c",u_limit=max_unique,twosample_data=twosample_data,cd_g_code=cd_g_code)
    g<-z[[1]]
    twosample_data<-z[[2]]
    rm(z);gc()
    if(dt){cat("Appropriate input data: true.\r\n")}
  }
  
  #par p1/p2
  control_p12_org<-list(suit=NULL,return_suit=F,suit_p_cut=0.1,suit_adj_m="fdr",
                        return_data_dhp=F,
                        p2_hp_p=0.05,p2_adj_m="none",
                        stage1_simplification=T,ind_hp=T,
                        outlier_detect=T,outlier_p=0.05,se_cut_k=1.5,
                        start_ind_hp=NULL,hp_p=0.1,adj_m="fdr",
                        start_nohp=NULL,
                        c_type=c("linear","median&mode"),boot_n=10000,mode_phi=1,
                        nlminb_control=list(nlminb_control=list(scale=1,eval.max=300,iter.max=300),
                                            nlm_control=list(gradtol=1e-10,stepmax=2,steptol=10^c(-10,-12,-8))))
  #c_type is combination type
  control_p12<-match.list(control_p12,control_p12_org)
  if(T){
    suit=control_p12$suit
    return_suit=control_p12$return_suit
    outlier_detect=control_p12$outlier_detect;outlier_p=control_p12$outlier_p
    se_cut_k=control_p12$se_cut_k
    start_ind_hp=control_p12$start_ind_hp;start_nohp=control_p12$start_nohp
    hp_p=control_p12$hp_p
  }
  check_start(g,start,"start")
  check_start(g,start_ind_hp,"start_ind_hp")
  check_start(g,start_nohp,"start_nohp")
  check_start(g,suit,"suit")
  
  #trinary g
  control_limit_c_org<-list(limit_c=T,dum_loc="auto",limit_c_twosample=T,dum_loc_twosample="auto",quantile=c(0,1),outlier=T)
  control_limit_c<-match.list(control_limit_c,control_limit_c_org)
  check_dum_loc(c,control_limit_c$dum_loc,"control_limit_c$dum_loc")
  check_dum_loc(twosample_data$c,control_limit_c$dum_loc_twosample,"control_limit_c$dum_loc_twosample")
  if(est_type!="me_mo_q_re"){
    j<-apply(g,2,FUN=function(x){
      y<-unique2(x)
      length(y)==3&min(table(x))>=n_min_limit&identical(as.numeric(sort(y)),c(0,1,2))
    })
    if(sum(j)<ncol(g)){
      message(sum(!j)," SNPs are not trinary and removed.")
      g<-myselect(g,which(j))
      start<-start[j]
      if(!c_inherit){
        c<-c[j]
      }
      start_ind_hp<-start_ind_hp[j]
      start_nohp<-start_nohp[j]
      suit<-suit[j]
      control_limit_c$dum_loc<-update_dum_loc(control_limit_c$dum_loc,j,c_inherit)
    }
    if(!is.null(twosample_data)){
      j0<-apply(twosample_data$g,2,FUN=function(x){
        y<-unique2(x)
        length(y)==3&min(table(x))>=n_min_limit&identical(as.numeric(sort(y)),c(0,1,2))
      })
      if(!identical(j,j0)){stop("Removing non-trinary SNPs leads to mismatched g and twosample_data$g.")}
      twosample_data<-update_twosample_data(twosample_data,which(j0),c_inherit)
      control_limit_c$dum_loc_twosample<-update_dum_loc(control_limit_c$dum_loc_twosample,which(j0),c_inherit)
    }
    stopifnot(ncol(g)!=0)
    gc()
  }
  
  #limit control variables
  limit_c<-function(matr,dum_loc=NULL,quantile=c(0.025,0.975),outlier=F){
    dum_loc<-as.integer(dum_loc)
    limit<-function(x,q,outlier){
      if(outlier){
        y<-quantile(x,c(0.25,0.75))
        y[2]<-y[2]+1.5*IQR(x)
        y[1]<-y[1]-1.5*IQR(x)
        x[x>y[2]]<-y[2]
        x[x<y[1]]<-y[1]
      }
      y<-quantile(x,q)
      x[x>y[2]]<-y[2]
      x[x<y[1]]<-y[1]
      x
    }
    for(i in 1:ncol(matr)){
      if(i%in%dum_loc){next}
      matr[,i]<-limit(matr[,i],quantile,outlier)
    }
    return(matr)
  }
  if(!is.null(c)&control_limit_c$limit_c){
    if(identical(control_limit_c$dum_loc,"auto")){
      control_limit_c$dum_loc<-lapply(c,FUN=function(x){
        z<-apply(x,2,FUN=function(y){(length(unique(y))<=2)})
        out<-which(z)
        if(length(out)==0){out<-NULL}
        out
      })
    }
    if(dt){cat("Start modifying the extreme values in c_list.\r\n")}
    for(i in 1:length(c)){
      c[[i]]<-limit_c(c[[i]],dum_loc=control_limit_c$dum_loc[[i]],quantile=control_limit_c$quantile,outlier=control_limit_c$outlier)
    }
  }
  if(!is.null(twosample_data$c)&control_limit_c$limit_c_twosample){
    if(identical(control_limit_c$dum_loc_twosample,"auto")){
      control_limit_c$dum_loc_twosample<-lapply(twosample_data$c,FUN=function(x){
        z<-apply(x,2,FUN=function(y){(length(unique(y))<=2)})
        out<-which(z)
        if(length(out)==0){out<-NULL}
        out
      })
    }
    if(length(control_limit_c$dum_loc_twosample)!=length(twosample_data$c)){stop("length(control_limit_c$dum_loc_twosample) != length(twosample_data$c).")}
    if(dt){cat("Start modifying the extreme values in twosample_data$c_list.\r\n")}
    for(i in 1:length(twosample_data$c)){
      twosample_data$c[[i]]<-limit_c(twosample_data$c[[i]],dum_loc=control_limit_c$dum_loc_twosample[[i]],quantile=control_limit_c$quantile,outlier=control_limit_c$outlier)
    }
  }
  
  #scaling
  if(dt){cat("x and continuous c (if appropriate) are automatically scaled.\r\n")}
  x<-as.numeric(scale(x))
  if(!is.null(c)){
    c<-lapply(c,FUN=function(x){
      loc<-apply(x,2,FUN=function(x){length(unique(x))>2})
      x_out<-as.matrix(scale(x))
      x_out[,!loc]<-x[,!loc]
      return(x_out)
    })
  }
  if(!is.null(twosample_data$c)){
    twosample_data$c<-lapply(twosample_data$c,FUN=function(x){
      loc<-apply(x,2,FUN=function(x){length(unique(x))>2})
      x_out<-as.matrix(scale(x))
      x_out[,!loc]<-x[,!loc]
      return(x_out)
    })
  }
  
  #snp_exp_check
  if(est_type!="p3"&snp_exp_check){
    if(dt){cat("Start snp_exp_check.\r\n")}
    if(is.null(twosample_data)){
      pre_snp_check<-get_sig_snp(x=x,g=g,c=c,c_inherit=c_inherit,start=NULL,type="c",
                                 p_cut=c(5e-8),return_dt=T,
                                 cd=F,trinary_only=F,
                                 control_limit_c=list(limit_c=F),scale=F,
                                 mc.cores=mc.cores,PSOCK=PSOCK,parallel_trace=parallel_trace,dt=dt,
                                 nlm=F)
    }else{
      pre_snp_check<-get_sig_snp(x=x,g=twosample_data$g,c=twosample_data$c,c_inherit=c_inherit,start=NULL,type="c",
                                 p_cut=c(5e-8),return_dt=T,
                                 cd=F,trinary_only=F,
                                 control_limit_c=list(limit_c=F),scale=F,
                                 mc.cores=mc.cores,PSOCK=PSOCK,parallel_trace=parallel_trace,dt=dt,
                                 nlm=F)
    }
    wald_p<-pre_snp_check$p
    
    if(is.character(p_snp)){
      print(formatC(quantile(wald_p,c(0,0.25,0.5,0.75,1),na.rm=T),format="e"),quote=FALSE)
      while(T){
        u_input<-eval(parse(text=readline("Please select a cutoff P value:")))
        u_input2<-readline(paste(sum(wald_p<u_input,na.rm=T),"SNPs remained. Continue? (y/n):"))
        if(identical(u_input2,"y")){break}
      }
      loc_pval<-which(wald_p<u_input)
    }else{
      loc_pval<-which(wald_p<p_snp)
      if(dt){cat(length(loc_pval),"SNPs satisfy p_snp.\r\n")}
    }
    
    g<-myselect(g,loc_pval)
    suit<-suit[loc_pval]
    start<-start[loc_pval]
    start_ind_hp<-start_ind_hp[loc_pval]
    start_nohp<-start_nohp[loc_pval]
    if(!is.null(twosample_data)){
      twosample_data<-update_twosample_data(twosample_data,loc_pval,c_inherit)
    }
    if(!c_inherit){
      c<-c[loc_pval]
    }
    stopifnot(ncol(g)!=0)
    gc()
  }
  
  #functions
  ivw_random<-function(eff,se,nlminb_control=list()){
    nlminb_control<-nlminb_control$nlminb_control
    nlminb_scale<-nlminb_control$scale
    if(is.null(nlminb_scale)){nlminb_scale<-1}
    nlminb_control[which(names(nlminb_control)=="scale")]<-NULL
    
    f_ivw_random<-function(beta,eff,se,est_vcov=F){
      if(est_vcov){
        s2<-beta[2]^2
        p<-dnorm(eff,beta[1],sqrt(s2+se^2))
        return(-log(p))
      }
      s2<-exp(beta[2])
      p<-dnorm(eff,beta[1],sqrt(s2+se^2))
      return(sum(-log(p)))
    }
    
    loc<-!(is.na(eff)|is.na(se))
    myeff<-eff[loc]
    myse<-se[loc]
    if(length(myse)<=3){
      warning("The final number of SNPs is <=3. Additive random effects model is not used.")
      return(list(eff=NA,se=NA,tau2=NA))
    }
    fit<-nlminb2nlm(nlminb2(f=f_ivw_random,p=c(0,-1),eff=myeff,se=myse,scale=nlminb_scale,control=nlminb_control))
    if(fit$code%in%c(1)){
      warning("Abnormal nlminb code in ivw_random.")
    }
    
    vcov<-est_variance_ml(f=f_ivw_random,beta_opt=c(fit$estimate[1],sqrt(exp(fit$estimate[2]))),n=length(myeff),loc=1:2,beta_name=c("u","s"),name="ivw_random",eff=myeff,se=myse,est_vcov=T,hessian=T,sandwich=F)
    return(list(eff=fit$estimate[1],se=sqrt(vcov[1,1]),tau2=exp(fit$estimate[2])))
  }
  mbe<-function(eff,w=NULL,phi=c(1,0.5,0.25)){
    if(is.null(w)){w<-rep(1/length(eff),length(eff))}
    s<-0.9*(min(sd(eff),mad(eff)))/length(eff)^(1/5)
    if(is.na(s)|s==0){return(NA)}
    beta<-NULL
    for(cur_phi in phi){
      h<-s*cur_phi
      densityIV<-density(eff,weights=w,bw=h)
      beta[length(beta)+1]<-mean(densityIV$x[densityIV$y==max(densityIV$y)])
    }
    return(beta)
  }
  median_w<-function(eff,w){
    loc<-which(!is.na(eff))
    eff<-eff[loc]
    w<-w[loc]
    names(eff)<-names(w)<-NULL
    
    eff.order<-eff[order(eff)]
    w.order<-w[order(eff)]
    w.sum<-cumsum(w.order)-0.5*w.order
    w.sum<-w.sum/sum(w)
    below<-max(which(w.sum<0.5))
    w.est<-eff.order[below]+(eff.order[below+1]-eff.order[below])/(w.sum[below+1]-w.sum[below])*(0.5-w.sum[below])
    return(w.est) 
  }
  mycombine2<-function(eff1,eff2,se1,se2,type=c("linear","median&mode"),boot_n=10000,mode_phi=1,dt=T,parallel_trace=T){
    vcov<-matrix(0,length(c(eff1,eff2)),length(c(eff1,eff2)))
    diag(vcov)<-c(se1,se2)^2
    out<-NULL
    
    if("linear"%in%type){
      fit<-meta_eff2(c(eff1,eff2),vcov)
      out<-rbind(out,data.frame(eff_m=fit$eff_m,se_m=fit$se_m))
      rownames(out)<-"linear_combination"
    }
    
    if("median&mode"%in%type){
      boot_data<-mvrnorm(boot_n,c(eff1,eff2),vcov)
      boot_s_me<-boot_w_me<-boot_s_mo<-boot_w_mo<-rep(NA,boot_n)
      n1<-length(eff1)
      n2<-n1+1
      n3<-length(eff1)+length(eff2)
      vcov0<-vcov[n2:n3,n2:n3]
      
      m2<-meta_eff2(eff2,vcov0)
      w<-1/c(diag(vcov)[1:n1],m2$se_m^2)
      w<-w/sum(w)
      eff<-c(eff1,m2$eff_m)
      s_me<-median(eff)
      w_me<-median_w(eff,w)
      s_mo<-mbe(eff,NULL,mode_phi)
      w_mo<-mbe(eff,w,mode_phi)
      
      if(dt){cat("Start bootstrapping.\r\n")}
      for(i in 1:boot_n){
        m2<-meta_eff2(boot_data[i,n2:n3],vcov0)
        eff<-c(boot_data[i,1:n1],m2$eff_m)
        boot_s_me[i]<-median(eff)
        boot_w_me[i]<-median_w(eff,w)
        boot_s_mo[i]<-mbe(eff,NULL,mode_phi)
        boot_w_mo[i]<-mbe(eff,w,mode_phi)
        if(parallel_trace){my_moni("Bootstrap repeat",i,boot_n)}
      }
      mydata<-data.frame(eff_m=c(s_me,w_me,s_mo,w_mo),
                         se_m=c(sd(boot_s_me),sd(boot_w_me),sd(boot_s_mo),sd(boot_w_mo))
      )
      rownames(mydata)<-c("simple_median","weighted_median","simple_mode","weighted_mode")
      out<-rbind(out,mydata)
    }
    return(out)
  }
  crt_parboot<-function(eff,se,n){
    out<-suppressWarnings(matrix(rnorm(length(eff)*n,mean=eff,sd=se),nrow=length(eff),ncol=n))
    return(out)
  }
  
  #####main estimating procedure#####
  if(est_type=="me_mo_q_re"){
    control_mmqr_org<-list(mode_phi=1,outlier_p=0.05,se_cut_k=NULL,boot_n=10000,inspect_data_mmqr=F,
                           nlminb_control=list(nlminb_control=list(scale=1,eval.max=300,iter.max=300),
                                               nlm_control=list(gradtol=1e-10,stepmax=2,steptol=10^c(-10,-12,-8))))
    control_mmqr<-match.list(control_me_mo_q_re,control_mmqr_org)
    
    if(is.null(data_me_mo_q_re)){
      if(dt){cat("Start estimating with the assumption of no hp.\r\n")}
      hp1<-est_proc_cont_p1.3(x=x,y=y,g_nohp=g,c=c,c_inherit=c_inherit,start=start,mc.cores=mc.cores,PSOCK=PSOCK,parallel_trace=parallel_trace,dt=dt,
                              twosample_data=twosample_data,nlminb_control=control_mmqr$nlminb_control)
      if(identical(class(hp1),"myerror")){return(hp1)}
      eff<-hp1$eff
      se<-hp1$se
    }else{
      hp1<-data_me_mo_q_re
      eff<-data_me_mo_q_re$eff
      se<-data_me_mo_q_re$se
    }
    
    if(control_mmqr$inspect_data_mmqr){return(hp1)}
    
    myeff<-eff[1,][!is.na(eff[1,])]
    myse<-se[1,][!is.na(se[1,])]
    stopifnot(length(myeff)==length(myse))
    if(length(myeff)==0){
      warning("Estimates are all NAs.")
      return(hp1)
    }
    
    w<-1/(myse^2)/sum(1/(myse^2))
    
    s_me<-median(myeff)
    w_me<-median_w(eff=myeff,w=w)
    
    s_mo<-mbe(myeff,NULL,control_mmqr$mode_phi)
    w_mo<-mbe(myeff,w,control_mmqr$mode_phi)
    
    ivw_re<-ivw_random(myeff,myse,control_mmqr$nlminb_control)
    SNP_names_re<-names(myeff)
    s_ivw<-meta_eff(myeff,myse)
    
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
        return(list(chisq=q,p.value=p,SNP_num=length(eff)))
      }
      R<-cbind(rep(1,length(eff)-1),-diag(length(eff)-1))
      chisq_stat<-t(R%*%eff)%*%solve(R%*%vcov%*%t(R))%*%(R%*%eff)
      p<-pchisq(chisq_stat,length(eff)-1,lower.tail=F)
      return(list(chisq=chisq_stat,p.value=p,SNP_num=length(eff)))
    }
    h_test<-equal_wald(myeff,myse^2)
    non_outlier<-outlier_test(myeff,myse,p_cut=control_mmqr$outlier_p,se_cut_k=control_mmqr$se_cut_k,dt=dt)
    loc<-!is.na(hp1$eff[1,])
    loc[loc]<-non_outlier
    if(dt){cat("Q-based outlier removal:",sum(!loc),"SNP(s) is/are judged to be outlier(s).\r\n")}
    eff_m<-meta_eff(myeff[non_outlier],myse[non_outlier])$est[1]
    se_m<-meta_eff(myeff[non_outlier],myse[non_outlier])$est[2]
    IVs_m<-colnames(hp1$eff)[loc]
    
    boot_s_me<-boot_w_me<-boot_s_mo<-boot_w_mo<-NULL
    if(T){
      if(dt){cat("Start bootstrapping.\r\n")}
      boot_data<-crt_parboot(eff=myeff,se=myse,n=control_mmqr$boot_n)
      boot_s_me<-boot_w_me<-boot_s_mo<-boot_w_mo<-rep(NA,control_mmqr$boot_n)
      for(i in 1:control_mmqr$boot_n){
        boot_s_me[i]<-median(boot_data[,i],na.rm=T)
        boot_w_me[i]<-median_w(eff=boot_data[,i],w=w)
        boot_s_mo[i]<-mbe(boot_data[,i],NULL,control_mmqr$mode_phi)
        boot_w_mo[i]<-mbe(boot_data[,i],w,control_mmqr$mode_phi)
        if(parallel_trace){my_moni("Bootstrap repeat",i,control_mmqr$boot_n)}
      }
    }
    
    myout<-list(median=list(simple_median=c(eff=s_me,boot_se=sd(boot_s_me,na.rm=T)),
                            weighted_median=c(eff=w_me,boot_se=sd(boot_w_me,na.rm=T)),
                            boot_data=list(simple_median=boot_s_me,weighted_median=boot_w_me),
                            SNP_name=SNP_names_re),
                mode=list(simple_mode=c(eff=s_mo,boot_se=sd(boot_s_mo,na.rm=T)),
                          weighted_mode=c(eff=w_mo,boot_se=sd(boot_w_mo,na.rm=T)),
                          boot_data=list(simple_mode=boot_s_mo,weighted_mode=boot_w_mo),
                          SNP_name=SNP_names_re),
                IVW_outlier_removal=list(eff_m=eff_m,se_m=se_m,SNP_name=IVs_m),
                add_random=list(eff_m=ivw_re$eff,se_m=ivw_re$se,SNP_name=SNP_names_re,tau2=ivw_re$tau2),
                IVW_simple=list(eff_m=s_ivw$est[1],se_m=s_ivw$est[2],SNP_name=SNP_names_re),
                data_me_mo_q_re=hp1,
                heterogeneity_test=h_test)
    class(myout)<-"MR-PROLLIM output"
    if(dt){cat("MR-PROLLIM classical methods (Procedure me_mo_q_re) finished.\r\n")}
    return(myout)
  }
  
  if(est_type%in%c("p1","p2")){
    #suitable SNPs
    if(!is.null(data_dhp)){suit<-data_dhp$suit}
    if(is.null(suit)){
      if(dt){cat("Start screening suitable SNPs for direct HP test.\r\n")}
      suit<-select_suit(x=x,g=g,c=c,c_inherit=c_inherit,twosample_data=twosample_data,mc.cores=mc.cores,
                        dt=dt,parallel_trace=F,PSOCK=PSOCK,p_cut=control_p12$suit_p_cut,adj_m=control_p12$suit_adj_m)
      if(identical(class(suit),"myerror")){return(suit)}
    }
    stopifnot(length(suit)==ncol(g))
    
    if(sum(suit)==0){stop("No suitable SNP for direct HP test.")}
    if(dt){cat(sum(suit),"SNP(s) is/are judged to be suitable for the direct HP test.\r\n")}
    if(return_suit){return(suit)}
    
    g_suit<-myselect(g,which(suit))
  }
  
  if(est_type=="p1"){
    if(is.null(data_dhp)){
      if(dt){cat("Start direct HP test.\r\n")}
      if(c_inherit){
        hp1<-est_proc_cont_p1(x=x,y=y,g_suit=g_suit,c=c,c_inherit=T,start=start[which(suit)],mc.cores=mc.cores,PSOCK=PSOCK,parallel_trace=parallel_trace,dt=dt,
                              twosample_data=update_twosample_data(twosample_data,which(suit),T),nlminb_control=control_p12$nlminb_control)
      }else{
        hp1<-est_proc_cont_p1(x=x,y=y,g_suit=g_suit,c=c[which(suit)],c_inherit=F,start=start[which(suit)],mc.cores=mc.cores,PSOCK=PSOCK,parallel_trace=parallel_trace,dt=dt,
                              twosample_data=update_twosample_data(twosample_data,which(suit),F),nlminb_control=control_p12$nlminb_control)
      }
    }else{
      hp1<-data_dhp$hp1
    }
    if(identical(class(hp1),"myerror")){return(hp1)}
    
    if(control_p12$return_data_dhp){return(list(suit=suit,hp1=hp1))}
    
    #Conduct outlier detection
    outlier_name<-NULL
    loc_no<-which(!is.na(hp1$eff[1,]))
    if(outlier_detect){
      if(dt){cat("Start outlier detection.\r\n")}
      non_outlier<-outlier_test(hp1$eff[1,],hp1$se[1,],vcov=NULL,p_cut=outlier_p,se_cut_k=se_cut_k,dt=dt,q_test=T)
      non_outlier2<-outlier_test(hp1$eff[1,],hp1$se[1,],vcov=NULL,p_cut=outlier_p,se_cut_k=se_cut_k,dt=F,q_test=F)
      if(dt){cat(sum(!non_outlier),"SNP(s) is/are judged to be outlier(s). Those with large SEs are moved to the cohort of unsuitable SNPs.\r\n")}
      loc_no<-which(non_outlier)
      outlier_name<-colnames(hp1$eff)[which(!non_outlier)]
      suit[suit==T]<-non_outlier2
    }
    
    if(length(loc_no)==0){
      warning("All suitable SNPs are judged to be outliers. Already extracted data are returned.")
      myout<-list( cochran_q_suit=list(Q=NA,Q_p=NA),
                   IVW_suit=list(eff_m=NA,se_m=NA,SNP_name=NA),
                   IVW_all=list(result=NULL,SNP_name=NULL),
                   d_hp_test=hp1,
                   SNP_outlier=outlier_name)
      class(myout)<-"MR-PROLLIM output"
      return(myout)
    }
    
    Q<-cochran_q(hp1$eff[1,loc_no],hp1$se[1,loc_no])
    Q_p<-pchisq(Q,length(loc_no)-1,lower.tail=F)
    
    eff_m<-meta_eff(hp1$eff[1,loc_no],hp1$se[1,loc_no])$est[1]
    se_m<-meta_eff(hp1$eff[1,loc_no],hp1$se[1,loc_no])$est[2]
    IVs_m<-colnames(hp1$eff)[loc_no]
    
    if(!control_p12$ind_hp){
      myout<-list( cochran_q_suit=list(Q=Q,Q_p=Q_p),
                   IVW_suit=list(eff_m=eff_m,se_m=se_m,SNP_name=IVs_m),
                   IVW_all=list(result=NULL,SNP_name=NULL),
                   d_hp_test=hp1,
                   SNP_outlier=outlier_name)
      class(myout)<-"MR-PROLLIM output"
      return(myout)
    }
    if(sum(!suit)==0){
      cat("No SNP is considered unsuitable for direct HP tests. Already extracted data are returned.\r\n")
      myout<-list( cochran_q_suit=list(Q=Q,Q_p=Q_p),
                   IVW_suit=list(eff_m=eff_m,se_m=se_m,SNP_name=IVs_m),
                   IVW_all=list(result=NULL,SNP_name=NULL),
                   d_hp_test=hp1,
                   SNP_outlier=outlier_name)
      class(myout)<-"MR-PROLLIM output"
      return(myout)
    }
    
    g_unsuit<-myselect(g,which(!suit))
    
    if(dt){cat("Start indirect HP test.\r\n")}
    if(c_inherit){
      hp2<-est_proc_cont_p1.2(x=x,y=y,g_unsuit=g_unsuit,c=c,
                              c_inherit=T,start=start_ind_hp[which(!suit)],
                              b1_a1_list=hp1$w1[loc_no],w_a1_list=hp1$w2[loc_no],
                              b1_sp=eff_m,
                              b1_indiv_matr=myselect(hp1$b1_indiv,loc_no),se_indiv=hp1$se[1,loc_no],
                              se_m=se_m,mc.cores=mc.cores,parallel_trace=parallel_trace,dt=dt,PSOCK=PSOCK,
                              twosample_data=update_twosample_data(twosample_data,which(!suit),T),nlminb_control=control_p12$nlminb_control)
    }else{
      hp2<-est_proc_cont_p1.2(x=x,y=y,g_unsuit=g_unsuit,c=c[which(!suit)],
                              c_inherit=F,start=start_ind_hp[which(!suit)],
                              b1_a1_list=hp1$w1[loc_no],w_a1_list=hp1$w2[loc_no],
                              b1_sp=eff_m,
                              b1_indiv_matr=myselect(hp1$b1_indiv,loc_no),se_indiv=hp1$se[1,loc_no],
                              se_m=se_m,mc.cores=mc.cores,parallel_trace=parallel_trace,dt=dt,PSOCK=PSOCK,
                              twosample_data=update_twosample_data(twosample_data,which(!suit),F),nlminb_control=control_p12$nlminb_control)
    }
    if(identical(class(hp2),"myerror")){return(hp2)}
    
    nohp<-p.adjust(hp2$wald_test[2,],control_p12$adj_m)>hp_p
    nohp[is.na(nohp)]<-F
    if(sum(nohp)==0){
      cat("All SNPs unsuitable for direct HP tests are detected with HP. Already extracted data are returned.\r\n")
      myout<-list( cochran_q_suit=list(Q=Q,Q_p=Q_p),
                   IVW_suit=list(eff_m=eff_m,se_m=se_m,SNP_name=IVs_m),
                   IVW_all=list(result=NULL,SNP_name=NULL),
                   d_hp_test=hp1,
                   ind_hp_test=hp2,
                   SNP_outlier=outlier_name)
      class(myout)<-"MR-PROLLIM output"
      return(myout)
    }
    if(dt){cat(sum(nohp),"SNP(s) is/are considered to be without HP.\r\n")}
    
    g_nohp<-myselect(g_unsuit,which(nohp))
    
    loc<-which(colnames(g)%in%colnames(g_nohp))
    
    if(is.null(data_nohp)){
      if(dt){cat("Start estimating using SNPs assumed to be without HP.\r\n")}
      if(c_inherit){
        nohp1<-est_proc_cont_p1.3(x=x,y=y,g_nohp=g_nohp,c=c,c_inherit=T,start=start_nohp[loc],mc.cores=mc.cores,PSOCK=PSOCK,parallel_trace=parallel_trace,dt=dt,
                                  twosample_data=update_twosample_data(twosample_data,loc,T),nlminb_control=control_p12$nlminb_control)
      }else{
        nohp1<-est_proc_cont_p1.3(x=x,y=y,g_nohp=g_nohp,c=c[loc],c_inherit=F,start=start_nohp[loc],mc.cores=mc.cores,PSOCK=PSOCK,parallel_trace=parallel_trace,dt=dt,
                                  twosample_data=update_twosample_data(twosample_data,loc,F),nlminb_control=control_p12$nlminb_control)
      }
      if(identical(class(nohp1),"myerror")){return(nohp1)}
    }else{
      loc_nohp<-which(colnames(data_nohp$eff)%in%colnames(g_nohp))
      if(length(loc_nohp)<ncol(g_nohp)*0.9){warning("SNP loss rate is > 10% in the given data_nohp. Please check whether data_nohp is correct.")}
      nohp1<-list(eff=myselect(data_nohp$eff,loc_nohp),se=myselect(data_nohp$se,loc_nohp))
    }
    loc2<-which(!is.na(nohp1$eff[1,]))
    stopifnot(length(loc2)>0)
    fit_c<-mycombine2(eff1=nohp1$eff[1,loc2],eff2=hp1$eff[1,loc_no],se1=nohp1$se[1,loc2],se2=hp1$se[1,loc_no],control_p12$c_type,control_p12$boot_n,control_p12$mode_phi,dt=dt,parallel_trace=parallel_trace)
    IVs_m_all<-c(names(nohp1$eff[1,loc2]),names(hp1$eff[1,loc_no]))
    
    myout<-list( cochran_q_suit=list(Q=Q,Q_p=Q_p),
                 IVW_suit=list(eff_m=eff_m,se_m=se_m,SNP_name=IVs_m),
                 IVW_all=list(result=fit_c,SNP_name=IVs_m_all),
                 d_hp_test=hp1,
                 ind_hp_test=hp2,
                 no_hp=nohp1,
                 IV_outlier=outlier_name)
    class(myout)<-"MR-PROLLIM output"
    if(dt){cat("Fixed-effects MR-PROLLIM (Procedure 1) finished.\r\n")}
    return(myout)
  }
  
  if(est_type=="p2"){
    if(is.null(data_dhp)){
      if(dt){cat("Start direct HP test.\r\n")}
      if(c_inherit){
        hp1<-est_proc_cont_p1(x=x,y=y,g_suit=g_suit,c=c,c_inherit=T,start=start[which(suit)],mc.cores=mc.cores,PSOCK=PSOCK,parallel_trace=parallel_trace,dt=dt,
                              twosample_data=update_twosample_data(twosample_data,which(suit),T),nlminb_control=control_p12$nlminb_control)
      }else{
        hp1<-est_proc_cont_p1(x=x,y=y,g_suit=g_suit,c=c[which(suit)],c_inherit=F,start=start[which(suit)],mc.cores=mc.cores,PSOCK=PSOCK,parallel_trace=parallel_trace,dt=dt,
                              twosample_data=update_twosample_data(twosample_data,which(suit),F),nlminb_control=control_p12$nlminb_control)
      }
    }else{
      hp1<-data_dhp$hp1
    }
    if(identical(class(hp1),"myerror")){return(hp1)}
    
    if(control_p12$return_data_dhp){return(list(suit=suit,hp1=hp1))}
    
    hp_test<-function(eff,se,p_cut,adj_m){
      which(p.adjust(pnorm(-abs(eff/se))*2,adj_m)>p_cut)
    }
    withouthp<-hp_test(eff=hp1$eff[2,],se=hp1$se[2,],p_cut=control_p12$p2_hp_p,adj_m=control_p12$p2_adj_m)
    
    if(length(withouthp)==0){
      warning("All suitable SNPs are considered to be with HP. Already extracted data are returned.")
      myout<-list( cochran_q_suit=list(Q=NA,Q_p=NA),
                   IVW_suit=list(eff_m=NA,se_m=NA,SNP_name=NA),
                   IVW_all=list(result=NULL,SNP_name=NULL),
                   d_hp_test=hp1)
      class(myout)<-"MR-PROLLIM output"
      return(myout)
    }
    
    if(dt){cat(length(withouthp),"suitable SNPs are considered to be without HP.\r\n")}
    
    loc<-which(colnames(g)%in%(colnames(g_suit)[withouthp]))
    nohp0<-NULL
    if(control_p12$stage1_simplification){
      if(dt){cat("Start stage1_simplification.\r\n")}
      if(c_inherit){
        nohp0<-est_proc_cont_p1.3(x=x,y=y,g_nohp=myselect(g,loc),c=c,c_inherit=T,start=start_nohp[loc],mc.cores=mc.cores,parallel_trace=parallel_trace,dt=dt,PSOCK=PSOCK,
                                  twosample_data=update_twosample_data(twosample_data,loc,T),nlminb_control=control_p12$nlminb_control,return_b1_indiv=T)
      }else{
        nohp0<-est_proc_cont_p1.3(x=x,y=y,g_nohp=myselect(g,loc),c=c[loc],c_inherit=F,start=start_nohp[loc],mc.cores=mc.cores,parallel_trace=parallel_trace,dt=dt,PSOCK=PSOCK,
                                  twosample_data=update_twosample_data(twosample_data,loc,F),nlminb_control=control_p12$nlminb_control,return_b1_indiv=T)
      }
      if(identical(class(nohp0),"myerror")){return(nohp0)}
      if(anyNA(nohp0$eff[1,])){
        if(dt){cat(sum(is.na(nohp0$eff[1,])),"SNPs are removed due to NAs in stage1_simplification.\r\n")}
      }
      withouthp<-which(!is.na(nohp0$eff[1,]))
      stopifnot(length(withouthp)>0)
      eff<-nohp0$eff[1,withouthp]
      se<-nohp0$se[1,withouthp]
      w1<-nohp0$w1[withouthp]
      w2<-nohp0$w2[withouthp]
      b1_indiv<-myselect(nohp0$b1_indiv,withouthp)
    }else{
      eff<-hp1$eff[1,withouthp]
      se<-hp1$se[1,withouthp]
      w1<-hp1$w1[withouthp]
      w2<-hp1$w2[withouthp]
      b1_indiv<-myselect(hp1$b1_indiv,withouthp)
    }
    
    #Conduct outlier detection
    outlier_name<-NULL
    non_outlier<-rep(T,length(eff))
    if(outlier_detect){
      if(dt){cat("Start outlier detection.\r\n")}
      non_outlier<-outlier_test(eff,se,vcov=NULL,p_cut=outlier_p,se_cut_k=se_cut_k,dt=dt,q_test=T)
      non_outlier2<-outlier_test(eff,se,vcov=NULL,p_cut=outlier_p,se_cut_k=se_cut_k,dt=F,q_test=F)
      if(dt){cat(sum(!non_outlier),"SNP(s) is/are judged to be outlier(s). Those with large SEs are moved to the cohort of unsuitable SNPs.\r\n")}
      suit[colnames(g)%in%(names(eff)[!non_outlier2])]<-F
      outlier_name<-names(eff)[!non_outlier]
    }
    
    eff<-eff[non_outlier]
    se<-se[non_outlier]
    w1<-w1[non_outlier]
    w2<-w2[non_outlier]
    b1_indiv<-myselect(b1_indiv,which(non_outlier))
    gc()
    
    if(length(eff)==0){
      warning("All suitable SNPs without HP are judged to be outliers. Already extracted data are returned.")
      myout<-list( cochran_q_suit=list(Q=NA,Q_p=NA),
                   IVW_suit=list(eff_m=NA,se_m=NA,SNP_name=NA),
                   IVW_all=list(result=NULL,SNP_name=NULL),
                   d_hp_test=hp1,
                   fs_no_hp=nohp0,
                   SNP_outlier=outlier_name)
      class(myout)<-"MR-PROLLIM output"
      return(myout)
    }
    
    Q<-cochran_q(eff,se)
    Q_p<-pchisq(Q,length(eff)-1,lower.tail=F)
    
    eff_m<-meta_eff(eff,se)$est[1]
    se_m<-meta_eff(eff,se)$est[2]
    IVs_m<-names(eff)
    #fs_data<-list(eff=eff,se=se)
    
    if(!control_p12$ind_hp){
      myout<-list( cochran_q_suit=list(Q=Q,Q_p=Q_p),
                   IVW_suit=list(eff_m=eff_m,se_m=se_m,SNP_name=IVs_m),
                   IVW_all=list(result=NULL,SNP_name=NULL),
                   d_hp_test=hp1,
                   fs_no_hp=nohp0,
                   SNP_outlier=outlier_name)
      class(myout)<-"MR-PROLLIM output"
      return(myout)
    }
    if(sum(!suit)==0){
      cat("No SNP is considered unsuitable for direct HP tests. Already extracted data are returned.\r\n")
      myout<-list( cochran_q_suit=list(Q=Q,Q_p=Q_p),
                   IVW_suit=list(eff_m=eff_m,se_m=se_m,SNP_name=IVs_m),
                   IVW_all=list(result=NULL,SNP_name=NULL),
                   d_hp_test=hp1,
                   fs_no_hp=nohp0,
                   SNP_outlier=outlier_name)
      class(myout)<-"MR-PROLLIM output"
      return(myout)
    }
    
    g_unsuit<-myselect(g,which(!suit))
    
    if(dt){cat("Start indirect HP test.\r\n")}
    if(c_inherit){
      hp2<-est_proc_cont_p1.2(x=x,y=y,g_unsuit=g_unsuit,c=c,
                              c_inherit=T,start=start_ind_hp[which(!suit)],
                              b1_a1_list=w1,w_a1_list=w2,
                              b1_sp=eff_m,
                              b1_indiv_matr=b1_indiv,se_indiv=se,
                              se_m=se_m,mc.cores=mc.cores,parallel_trace=parallel_trace,dt=dt,PSOCK=PSOCK,
                              twosample_data=update_twosample_data(twosample_data,which(!suit),T),nlminb_control=control_p12$nlminb_control)
    }else{
      hp2<-est_proc_cont_p1.2(x=x,y=y,g_unsuit=g_unsuit,c=c[which(!suit)],
                              c_inherit=F,start=start_ind_hp[which(!suit)],
                              b1_a1_list=w1,w_a1_list=w2,
                              b1_sp=eff_m,
                              b1_indiv_matr=b1_indiv,se_indiv=se,
                              se_m=se_m,mc.cores=mc.cores,parallel_trace=parallel_trace,dt=dt,PSOCK=PSOCK,
                              twosample_data=update_twosample_data(twosample_data,which(!suit),F),nlminb_control=control_p12$nlminb_control)
    }
    if(identical(class(hp2),"myerror")){return(hp2)}
    
    nohp<-p.adjust(hp2$wald_test[2,],control_p12$adj_m)>hp_p
    nohp[is.na(nohp)]<-F
    if(sum(nohp)==0){
      cat("All SNPs unsuitable for direct HP tests are detected with HP. Already extracted data are returned.\r\n")
      myout<-list( cochran_q_suit=list(Q=Q,Q_p=Q_p),
                   IVW_suit=list(eff_m=eff_m,se_m=se_m,SNP_name=IVs_m),
                   IVW_all=list(result=NULL,SNP_name=NULL),
                   d_hp_test=hp1,
                   fs_no_hp=nohp0,
                   ind_hp_test=hp2,
                   SNP_outlier=outlier_name)
      class(myout)<-"MR-PROLLIM output"
      return(myout)
    }
    if(dt){cat(sum(nohp),"SNP(s) is/are considered to be without HP.\r\n")}
    
    g_nohp<-myselect(g_unsuit,which(nohp))
    
    loc<-which(colnames(g)%in%colnames(g_nohp))
    
    if(is.null(data_nohp)){
      if(dt){cat("Start estimating using SNPs assumed to be without HP.\r\n")}
      if(c_inherit){
        nohp1<-est_proc_cont_p1.3(x=x,y=y,g_nohp=g_nohp,c=c,c_inherit=T,start=start_nohp[loc],mc.cores=mc.cores,PSOCK=PSOCK,parallel_trace=parallel_trace,dt=dt,
                                  twosample_data=update_twosample_data(twosample_data,loc,T),nlminb_control=control_p12$nlminb_control)
      }else{
        nohp1<-est_proc_cont_p1.3(x=x,y=y,g_nohp=g_nohp,c=c[loc],c_inherit=F,start=start_nohp[loc],mc.cores=mc.cores,PSOCK=PSOCK,parallel_trace=parallel_trace,dt=dt,
                                  twosample_data=update_twosample_data(twosample_data,loc,F),nlminb_control=control_p12$nlminb_control)
      }
      if(identical(class(nohp1),"myerror")){return(nohp1)}
    }else{
      loc_nohp<-which(colnames(data_nohp$eff)%in%colnames(g_nohp))
      if(length(loc_nohp)<ncol(g_nohp)*0.9){warning("SNP loss rate is > 10% in the given data_nohp. Please check whether data_nohp is correct.")}
      nohp1<-list(eff=myselect(data_nohp$eff,loc_nohp),se=myselect(data_nohp$se,loc_nohp))
    }
    loc2<-which(!is.na(nohp1$eff[1,]))
    stopifnot(length(loc2)>0)
    fit_c<-mycombine2(eff1=nohp1$eff[1,loc2],eff2=eff,se1=nohp1$se[1,loc2],se2=se,control_p12$c_type,control_p12$boot_n,control_p12$mode_phi,dt=dt,parallel_trace=parallel_trace)
    IVs_m_all<-c(names(nohp1$eff[1,loc2]),names(eff))
    
    myout<-list( cochran_q_suit=list(Q=Q,Q_p=Q_p),
                 IVW_suit=list(eff_m=eff_m,se_m=se_m,SNP_name=IVs_m),
                 IVW_all=list(result=fit_c,SNP_name=IVs_m_all),
                 d_hp_test=hp1,
                 fs_no_hp=nohp0,
                 ind_hp_test=hp2,
                 no_hp=nohp1,
                 IV_outlier=outlier_name)
    class(myout)<-"MR-PROLLIM output"
    if(control_p12$stage1_simplification){
      pl<-"with fs)"
    }else{
      pl<-"without fs)"
    }
    if(dt){cat("Fixed-effects MR-PROLLIM (Procedure 2,",pl,"finished.\r\n")}
    return(myout)
  }
  
  if(est_type=="p3"){
    control_est_k_prior_org<-list(p_snp="previous",start=c(0,0,-1,-1,1),p0_start=c(seq(0.1,0.9,by=0.2),0.99,0.01),
                                  auto_s=T,p0_sp=NULL,p0_cut=1e-8,u1_sp=NULL,
                                  nlminb_control=list(rel.tol=1e-10,sing.tol=1e-10,step.min=1,eval.max=300,iter.max=300))
    control_est_k_post_org<-list(n_post=3000,n0=10000,p_cover=0.9999,f=NULL)
    control_est_k_prior<-match.list(control_est_k_prior,control_est_k_prior_org)
    control_est_k_post<-match.list(control_est_k_post,control_est_k_post_org)
    
    control_p3_org<-list(n_snp_limit=50,n_min_limit=n_min_limit,
                         inspect_data_p3_k=F,inspect_data_p3=F,inspect_data_p3_opt=F,
                         nome=F,p_snp=p_snp,
                         s_filter=F,s_filter_ask=F,s_k_r_limit=c(0.2,0.2),s_k_a_limit="auto",auto_k_limit=0.2,
                         beta_start=c(0,0,0,0,-1,2,0),auto_s1=T,ask=F,
                         p1_sp=NULL,p2_sp=NULL,r_sp=NULL,model_u2=T,Egger="auto",
                         model_select=T,check_fit_upper=0.999,check_fit_lower=0.001,s1_cut_k=0.01,
                         sandwich=T,hessian=T,adj_rs=F,se_adj_boot_n=10000,
                         boot_se=F,n_boot=3000,n_rep_max=3,
                         nlminb_control0=list(nlminb_control=list(scale=1,eval.max=300,iter.max=300),
                                              nlm_control=list(gradtol=1e-10,stepmax=2,steptol=10^c(-10,-12,-8))),
                         nlminb_control=list(nlminb_control=list(rel.tol=1e-12,sing.tol=1e-12,step.min=1,eval.max=300,iter.max=300),
                                             nlm_control=list(gradtol=1e-8,stepmax=2,steptol=1e-8)))
    control_p3<-match.list(control_p3,control_p3_org)
    nlm_c_list_org<-list(hessian=F,fscale=1,print.level=0,ndigit=12,
                         gradtol=1e-8,stepmax=2,steptol=1e-8,iterlim=300,check.analyticals=T)
    nlm_c_list<-match.list(control_p3$nlminb_control$nlm_control,nlm_c_list_org)
    nlminb_scale<-control_p3$nlminb_control$nlminb_control$scale
    if(is.null(nlminb_scale)){nlminb_scale<-1}
    control_p3$nlminb_control$nlminb_control[which(names(control_p3$nlminb_control$nlminb_control)=="scale")]<-NULL
    nlminb_c_list<-control_p3$nlminb_control$nlminb_control
    
    control_sec<-list(nlminb_control=control_p3$nlminb_control0)
    mod_egger<-function(x){
      if(identical(x,"T")){return(T)}
      if(identical(x,"F")){return(F)}
      if(identical(x,"auto")){return("auto")}
      if(is.logical(x)){return(x)}
      stop("Egger in control_p3 should be T, F, or 'auto'.")
    }
    control_p3$Egger<-mod_egger(control_p3$Egger)
    p1_sp<-control_p3$p1_sp;p2_sp<-control_p3$p2_sp;r_sp<-control_p3$r_sp;model_u2<-control_p3$model_u2
    beta_start<-control_p3$beta_start
    boot_se<-control_p3$boot_se
    #n_boot_max<-control_p3$n_boot_max;n_boot<-control_p3$n_boot
    ask<-control_p3$ask
    n_boot<-control_p3$n_boot
    n_rep_max<-control_p3$n_rep_max
    s1_cut_k<-control_p3$s1_cut_k
    
    control_global_search_org<-list(global_search=T,global_search_EI=T,gs_type="genoud",
                                    genoud_control=list(pop.size=min(max(1000,mc.cores*300),10000),pop.size.EI=100,max.generations=60,wait.generations=10,
                                                        hard.generation.limit=TRUE,
                                                        solution.tolerance=1e-5,
                                                        gradient.check=F,print.level=0,
                                                        BFGSburnin=10,balance=F,
                                                        optim_control=list()),
                                    GenSA_control=list(maxit=3000,maxit.EI=300,max.time=NULL,verbose=F,trace.mat=F,seed=-100377),
                                    lower=c(-10,-10,-10,-10,-30,-10,-5),upper=c(10,10,10,10,10,10,5),auto_s1_k=5)
    control_gs<-match.list(control_global_search,control_global_search_org) 
    genoud_control<-match.list(control_gs$genoud_control,control_global_search_org$genoud_control)
    stopifnot(identical(names(control_gs$GenSA_control)[1],"maxit"))
    stopifnot(identical(names(control_gs$GenSA_control)[2],"maxit.EI"))
    update_GenSA<-function(list,exc){
      list<-list[-exc]
      names(list)[1]<-"maxit"
      list
    }
    update_opt_c<-function(optim_control,beta_loc){
      if(length(optim_control)>0){
        for(i in 1:length(optim_control)){
          if(length(optim_control[[i]])==7){
            optim_control[[i]]<-optim_control[[i]][beta_loc]
          }
        }
      }
      optim_control
    }
    join_par<-function(expr1,expr2){
      c1<-paste0(capture.output(expr1),collapse="")
      c2<-paste0(capture.output(expr2),collapse="")
      if(identical(c2,"expression()")){
        return(paste0(stringr::str_remove_all(c1,"expression\\(|\\)\\)$"),")"))
      }
      cc<-paste0(stringr::str_remove_all(c1,"expression\\(|\\)\\)$"),", ",
                 stringr::str_remove_all(c2,"expression\\(|\\)$"),")")
      cc
    }
    
    #functions required
    trans_par_norm<-function(beta,p1_sp=NULL,p2_sp=NULL,r_sp=NULL,model_u2=T,est_type="b",t_b1=F){
      if(est_type=="b"&!t_b1){b1<-1-exp(-beta[1])}else{b1<-beta[1]}
      if(is.null(p1_sp)){p1<-1/(1+1/exp(beta[2]))}else{p1<-p1_sp}
      if(is.null(p2_sp)){p2<-1/(1+1/exp(beta[3]))}else{p2<-p2_sp}
      u1<-beta[4]
      s1<-sqrt(exp(beta[5]))
      if(is.null(r_sp)){r<-1/(1+1/exp(beta[6]));r_t<-atanh(r)}else{r<-r_sp;r_t<-atanh(r)}
      if(model_u2){u2<-beta[7]}else{u2<-0}
      
      if(est_type=="b"){
        return(c(b1_t=b1,p1=p1,p2=p2,u1=u1,s1=s1,r_t=r_t,u2=u2))
      }
      return(c(b1=b1,p1=p1,p2=p2,u1=u1,s1=s1,r_t=r_t,u2=u2))
    }
    crt_data2<-function(y,x,g,c,start=NULL,name=NULL,twosample_data=NULL,
                        stage1=F,data_k_hat=NULL,data_k_sigma=NULL,nlminb_control=list()){
      #c is a matrix.
      #start is a vector or NULL.
      #list(matrix) c in twosample_data.
      #vector g in twosample_data.
      ncol2<-function(x){
        if(is.null(x)){return(0)}else{return(ncol(x))}
      }
      
      loc_m<-which(!is.na(g))
      y<-y[loc_m]
      if(is.null(twosample_data)){x<-x[loc_m]}
      g<-g[loc_m]
      
      j<-F
      if(is.null(c)&is.null(twosample_data)){
        o1<-y[g==1]
        o2<-y[g==2]
        o3<-y[g==0]
        o4<-x[g==1]
        o5<-x[g==2]
        o6<-x[g==0]
        
        if(stage1){
          j<-unlist(lapply(list(o1,o2,o3),FUN=function(x){length(unique(x))==1}))
          if(T%in%j){
            return(list(k_hat=c(NA,NA),k_sigma=matrix(NA,2,2)))
          }
        }
        
        if(T){
          out3<-(mean(o4))-(mean(o6))#k1
          out4<-(mean(o5))-(mean(o6))#k2
          
          m2<-var(o6)/length(o6)
          vcov_k<-rbind(
            c( var(o4)/length(o4)+m2, m2),
            c( m2, var(o5)/length(o5)+m2)
          )
          
          if(stage1){
            return(list(k_hat=c(out3,out4),k_sigma=vcov_k))
          }
          
          out1<-log(mean(o1))-log(mean(o3))#m1
          out2<-log(mean(o2))-log(mean(o3))#m2
          
          m1<-1/(mean(o3)^2)*var(o3)/length(o3)
          vcov_m<-rbind(
            c( 1/(mean(o1)^2)*var(o1)/length(o1)+m1, m1),
            c( m1, 1/(mean(o2)^2)*var(o2)/length(o2)+m1)
          )
          
          m3<-1/mean(o3)*cov(o3,o6)/length(o3)
          vcov_mk<-rbind(
            c(1/mean(o1)*cov(o1,o4)/length(o1)+m3,m3),
            c(m3,1/mean(o2)*cov(o2,o5)/length(o2)+m3)
          )
          
          m_hat<-rbind(c(out1,out2))
          m_sigma<-rbind(c(vcov_m[1,1],vcov_m[1,2],vcov_m[2,2]))
          k_hat<-rbind(c(out3,out4))
          k_sigma<-rbind(c(vcov_k[1,1],vcov_k[1,2],vcov_k[2,2]))
          
          mk_vcov<-matrix(NA,4,4)
          mk_vcov[1:2,1:2]<-vcov_m
          mk_vcov[3:4,3:4]<-vcov_k
          mk_vcov[1:2,3:4]<-vcov_mk
          mk_vcov[3:4,1:2]<-t(vcov_mk)
          #j<-!cover_p(c(mean(o1),mean(o2),mean(o3)),sqrt(c(var(o1)/length(o1),var(o2)/length(o2),var(o3)/length(o3))),p=p_limit)
        }
      }
      if((!is.null(c))|(!is.null(twosample_data))){
        if(is.null(twosample_data)){
          c<-myselect(c,loc_m,"r")
          g_dum<-crt_dum(g)
          n1<-ncol(g_dum)
          
          if(stage1){
            d_matr<-cbind(g_dum,c,1)
            a10<-solve(t(d_matr)%*%d_matr)%*%t(d_matr)%*%x
            a1<-a10[1:n1]
            a1_indiv_c<-est_indiv(B=-t(d_matr)%*%d_matr/length(x),M_indiv=c(x-d_matr%*%a10)*d_matr,1:length(a1),rep(0,length(a1)))
            k_sigma0<-est_vcov(a1_indiv_c)
            if(anyNA(k_sigma0)){a1<-c(NA,NA)}
            return(list(k_hat=a1,k_sigma=k_sigma0))
          }else{
            a1<-data_k_hat
          }
          
          if(is.null(start)){start_m<-c(rep(0,n1),rep(0,ncol(c)),-1)}else{
            start_m<-start[1:length(c(rep(0,n1),rep(0,ncol(c)),-1))]
          }
          fit1<-mynlminb(f=f_cont_c_mr,p=start_m,y=y,g=g_dum,g_dum=g_dum,c=c,a1=a1,b4_prior=rep(1,n1),b1_sp=0,name=name,control=nlminb_control)
          vcov<-sandwich_cont_c_mr(c( 1:n1,(n1+2):(2*n1+1) ),beta_opt=fit1$estimate,y=y,g=g_dum,g_dum=g_dum,c=c,
                                   a1=a1,x=x,b4_prior=rep(1,n1),b1_sp=0,b1_sp_indiv=rep(0,length(y)),direct_M=T,name=name)
          if(anyNA(vcov)){fit1$estimate<-rep(NA,length(fit1$estimate))}
          m_hat<-rbind(fit1$estimate[1:n1])
          m_sigma<-rbind(c(vcov[3,3],vcov[3,4],vcov[4,4]))
          k_hat<-rbind(a1)
          k_sigma<-rbind(c(vcov[1,1],vcov[1,2],vcov[2,2]))
          
          mk_vcov<-vcov
          mk_vcov[,3:4]<-vcov[,1:2]
          mk_vcov[,1:2]<-vcov[,3:4]
          o<-mk_vcov[1:2,]
          mk_vcov[1:2,]<-mk_vcov[3:4,]
          mk_vcov[3:4,]<-o
        }else{
          if(stage1){
            g_x<-c(twosample_data$g)
            c_x<-twosample_data$c[[1]]
            if(identical(as.numeric(length(unique(c_x[,1]))),1)){c_x<-NULL}
            loc_m2<-which(!is.na(g_x))
            g_x<-g_x[loc_m2]
            c_x<-myselect(c_x,loc_m2,"r")
            x<-x[loc_m2]
            
            g_dum2<-crt_dum(g_x)
            d_matr2<-cbind(g_dum2,c_x,1)
            b<-(solve(t(d_matr2)%*%d_matr2)%*%t(d_matr2)%*%x)
            a1_2<-b[1:ncol(g_dum2)]
            
            a1_indiv_c<-est_indiv(B=-t(d_matr2)%*%d_matr2/length(x),M_indiv=c(x-d_matr2%*%b)*d_matr2,1:length(a1_2),rep(0,length(a1_2)))
            k_sigma0<-est_vcov(a1_indiv_c)
            if(anyNA(k_sigma0)){a1_2<-c(NA,NA)}
            return(list(k_hat=a1_2,k_sigma=k_sigma0))
          }else{
            a1_2<-data_k_hat
            k_vcov<-data_k_sigma
          }
          
          c<-myselect(c,loc_m,"r")
          g_dum<-crt_dum(g)
          n1<-ncol(g_dum)
          if(is.null(start)){start_m<-c(rep(0,n1),rep(0,ncol2(c)),-1)}else{
            start_m<-start[1:length(c(rep(0,n1),rep(0,ncol2(c)),-1))]
          }
          fit1<-mynlminb(f=f_cont_c_mr,p=start_m,y=y,g=g_dum,g_dum=g_dum,c=c,a1=a1_2,b4_prior=rep(1,n1),b1_sp=0,name=name,control=nlminb_control)
          fit_indiv<-sandwich_cont_c_mr(c( 1:n1,(n1+2):(2*n1+1) ),beta_opt=fit1$estimate,y=y,g=g_dum,g_dum=g_dum,c=c,a1=a1_2,a1_indiv=matrix(a1_2,ncol=length(a1_2),nrow=length(y),byrow=T),
                                        b4_prior=rep(1,n1),b1_sp=0,b1_sp_indiv=rep(0,length(y)),
                                        return_indiv=T,return_w=T,direct_M=T,name=name)
          
          #a1_indiv<-est_indiv(B=-t(d_matr2)%*%d_matr2/length(x),M_indiv=c(x-d_matr2%*%b)*d_matr2,1:length(a1_2),a1_2)
          vcov<-est_vcov(fit_indiv[[1]])+fit_indiv[[2]]%*%k_vcov%*%t(fit_indiv[[2]])
          if(anyNA(vcov)){fit1$estimate<-rep(NA,length(fit1$estimate))}
          m_hat<-rbind(fit1$estimate[1:n1])
          m_sigma<-rbind(c(vcov[3,3],vcov[3,4],vcov[4,4]))
          k_hat<-rbind(a1_2)
          k_sigma<-rbind(c(vcov[1,1],vcov[1,2],vcov[2,2]))
          
          mk_vcov<-vcov
          mk_vcov[,3:4]<-vcov[,1:2]
          mk_vcov[,1:2]<-vcov[,3:4]
          o<-mk_vcov[1:2,]
          mk_vcov[1:2,]<-mk_vcov[3:4,]
          mk_vcov[3:4,]<-o
          #stop("need to check mk_vcov 0.")
        }
      }
      
      if(sum(j)!=0){
        m_hat[1,]<-c(NA,NA)
        m_sigma[1,]<-c(NA,NA,NA)
        k_hat[1,]<-c(NA,NA)
        k_sigma[1,]<-c(NA,NA,NA)
        mk_vcov<-matrix(NA,4,4)
      }
      return(list(m_hat=m_hat,m_sigma=m_sigma,k_hat=k_hat,k_sigma=k_sigma,mk_sigma=mk_vcov))
    }
    diagnose_fit<-function(fit,p1_sp=NULL,p2_sp=NULL,r_sp=NULL,model_u2=T,upper=0.999,lower=0.001,myEgger,dt=T){
      loc<-c(is.null(p1_sp),is.null(p2_sp),is.null(r_sp))
      #c(T,T,T) full model
      #c(T,T,F) full model with r_sp=1/0
      #c(T,F,T) reduced model
      #c(T,F,F) reduced model with r_sp=1/0
      #c(F,F,T) E/I model
      #c(F,F,F) E/I model with r_sp=1/0
      m<-trans_proc_p3(fit$estimate,p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2)
      x<-m[c("p1","p2","r")]
      
      if(dt){print(m)}
      
      if(anyNA(x)){return("fail")}
      if(identical(loc,c(F,F,F))){return("pass")}
      
      if((x[1]<lower|x[1]>upper)&loc[1]==T){
        #if(x[3]<lower&loc[3]==T){r_sp<-0;if(dt){message("Model selection: set r_sp = 0.")}}
        #if(x[3]>upper&loc[3]==T){r_sp<-1;if(dt){message("Model selection: set r_sp = 1.")}}
        #if(x[1]>upper){message("Estimate of p1 is too close to 1.");Egger_info<<-"p1 -> 1";if(identical(myEgger,"auto")){myEgger<-T}}else{Egger_info<<-"p1 -> 0";if(identical(myEgger,"auto")){myEgger<-F}}
        #if(myEgger){o<-"Egger model."}else{o<-"intercept model."}
        #if(dt){message(paste("Model selection: to",o))}
        if(dt){message("Estimate of p1 is too close to 0 or 1.")}
        return(list(p1_sp=1,p2_sp=0,r_sp=NULL,model_u2=T,myEgger=myEgger))
      }
      if((x[2]<lower|x[2]>upper)&loc[2]==T){
        if(x[3]<lower&loc[3]==T){r_sp<-0;if(dt){message("Model selection: set r_sp = 0.")}}
        if(x[3]>upper&loc[3]==T){r_sp<-1;if(dt){message("Model selection: set r_sp = 1.")}}
        if(dt){message("Model selection: to reduced model.")}
        return(list(p1_sp=NULL,p2_sp=1,r_sp=r_sp,model_u2=model_u2,myEgger=myEgger))
      }
      if(x[3]<lower&loc[3]==T){
        if(dt){message("Model selection: set r_sp = 0.")}
        return(list(p1_sp=p1_sp,p2_sp=p2_sp,r_sp=0,model_u2=model_u2,myEgger=myEgger))
      }
      if(x[3]>upper&loc[3]==T){
        if(dt){message("Model selection: set r_sp = 1.")}
        return(list(p1_sp=p1_sp,p2_sp=p2_sp,r_sp=1,model_u2=model_u2,myEgger=myEgger))
      }
      return("pass")
    }
    check_fit<-function(fit1,p1_sp=NULL,p2_sp=NULL,r_sp=NULL,model_u2=T,lower=0.001,upper=0.999){
      if(fit1$code%in%c(0,1)){
        m<-trans_proc_p3(fit1$estimate,p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2)
        m2<-m[c("p1","p2","r")][c(is.null(p1_sp),is.null(p2_sp),is.null(r_sp))]
        j<-m2>upper|m2<lower
        if(sum(j)!=0){return("abnormal value detected by check_fit")}
      }
      return(fit1$code)
    }
    nrow2<-function(x){
      if(is.null(x)){return(0)}else{
        return(nrow(x))
      }
    }
    mybind_list<-function(list){
      out<-NULL
      for(i in 1:length(list)){
        out<-rbind(out,list[[i]])
      }
      return(out)
    }
    format_fit1<-function(fit1){
      if(!is.list(fit1)){return(fit1)}
      list(minimum=fit1$value,estimate=fit1$par,gradient=fit1$gradients,
           code=0,iterations=fit1$generations,message=list(type="GenSA or genoud",peakgeneration=fit1$peakgeneration,operators=fit1$operators))
    }
    withWarnings<-function(expr){
      mywarnings<-NULL 
      value<-withCallingHandlers(expr,warning=function(w){
        mywarnings<<-c(mywarnings,conditionMessage(w))
        invokeRestart("muffleWarning")})
      out<-list(value=value,warning=mywarnings)
      class(out)<-"mywarning"
      out
    }
    mywald_test<-function(eff,vcov){
      if(anyNA(eff)){return(c(NA,NA))}
      chisq<-t(eff)%*%solve(vcov)%*%cbind(eff)
      p<-pchisq(chisq,df=length(eff),lower.tail=F)
      return(c(chisq,p))
    }
    crt_beta_loc<-function(beta,p1_sp,p2_sp,r_sp,model_u2){
      exc<-NULL
      if(!is.null(p1_sp)){
        exc<-c(exc,"p1")
      }
      if(!is.null(p2_sp)){
        exc<-c(exc,"p2")
        if(identical(p2_sp,0)){
          exc<-c(exc,"u1")
        }
      }
      if(!is.null(r_sp)){
        exc<-c(exc,"r")
      }
      if(!model_u2){
        exc<-c(exc,"u2")
      }
      loc<-which(!c("b1","p1","p2","u1","s1","r","u2")%in%exc)
      beta<-beta[loc]
      return(list(beta,loc))
    }
    ext_estimate<-function(est,beta_start,beta_loc){
      beta_start[beta_loc]<-est
      beta_start
    }
    
    #functions for postsample
    if(T){
      #This is the new version
      get_rotate_sample<-function(s_s,mychisq,n=10000,center=c(0,0)){
        if(mychisq==0){return(list(sample=matrix(center,nrow=n,ncol=2,byrow=T),area=0))}
        A=s_s[1,1]/mychisq
        B=s_s[2,2]/mychisq
        C=s_s[1,2]/mychisq
        
        if(abs(A-B)/abs(A)<1e-3){
          r<-pi/4
          c2<-C/(sin(r)*cos(r))
          c3<-A+B
          a<-sqrt(2/(c3+c2))
          b<-sqrt(2/(c3-c2))
        }else{
          f<-function(x,v){-sin(x)*cos(x)/cos(2*x)-v}
          r<-uniroot(f=f,interval=c(-pi/4,pi/4),v=C/(B-A),extendInt ="no")$root
          b<-sqrt((-cos(2*r))/(A*sin(r)^2-B*cos(r)^2))
          a<-sqrt((-cos(2*r))/(-A*cos(r)^2+B*sin(r)^2))
        }
        
        mysample<-NULL
        while(T){
          mysample1<-cbind(runif(round(n*1.5),-a,a),runif(round(n*1.5),-b,b))
          id<-mysample1[,1]^2/a^2+mysample1[,2]^2/b^2<=1
          mysample<-rbind(mysample,mysample1[id,])
          if(nrow(mysample)>=n){break}
        }
        mysample<-mysample[1:n,]
        myrotate<-function(x,y,r,center){
          cbind(x*cos(r)-y*sin(r)+center[1],y*cos(r)+x*sin(r)+center[2])
        }
        return(list(sample=myrotate(mysample[,1],mysample[,2],r,center),area=pi*a*b))
      }
      
      if(F){
        est_rotate_par<-function(s_s,mychisq){
          A=s_s[1,1]/mychisq
          B=s_s[2,2]/mychisq
          C=s_s[1,2]/mychisq
          
          f<-function(x,v){-sin(x)*cos(x)/cos(2*x)-v}
          r<-uniroot(f=f,interval=c(-pi/4,pi/4),v=C/(B-A),extendInt ="no")$root
          b<-sqrt((-cos(2*r))/(A*sin(r)^2-B*cos(r)^2))
          a<-sqrt((-cos(2*r))/(-A*cos(r)^2+B*sin(r)^2))
          return(list(a=a,b=b,r=r))
        }
      }
      
      est_rotate_par<-function(s_s,mychisq){
        A=s_s[1,1]/mychisq
        B=s_s[2,2]/mychisq
        C=s_s[1,2]/mychisq
        
        if(abs(A-B)/abs(A)<1e-3){
          r<-pi/4
          c2<-C/(sin(r)*cos(r))
          c3<-A+B
          a<-sqrt(2/(c3+c2))
          b<-sqrt(2/(c3-c2))
        }else{
          f<-function(x,v){-sin(x)*cos(x)/cos(2*x)-v}
          r<-uniroot(f=f,interval=c(-pi/4,pi/4),v=C/(B-A),extendInt ="no")$root
          b<-sqrt((-cos(2*r))/(A*sin(r)^2-B*cos(r)^2))
          a<-sqrt((-cos(2*r))/(-A*cos(r)^2+B*sin(r)^2))
        }
        return(list(a=a,b=b,r=r))
      }
      myrotate<-function(x,y,r,center){
        cbind(x*cos(r)-y*sin(r)+center[1],y*cos(r)+x*sin(r)+center[2])
      }
      get_rotate_sample<-function(par,n=10000,center=c(0,0)){
        a<-par$a
        b<-par$b
        r<-par$r
        while(T){
          mysample<-cbind(runif(n,-a,a),runif(n,-b,b))
          id<-mysample[,1]^2/a^2+mysample[,2]^2/b^2<=1
          mysample<-mysample[id,]
          if(nrow(mysample)!=0){break}
        }
        return(myrotate(mysample[,1],mysample[,2],r,center))
      }
      get_ellipse_p2<-function(mychisq,mysigma_ellipse,myu,mysigma){
        #my_ellipse is located at c(0,0)
        #myu and mysigma are parameters for background MVN 
        if(mychisq==0){return(0)}
        my_pmvnEll(mychisq,mysigma,myu,solve(mysigma_ellipse),x0=c(0,0),T)
      }
    }
    get_postsample<-function(k_hat,s_matr,p_cover=0.9999,n=3000,n_f1=10000,p_p0,p_u1,p_sigma,name){
      mychisq<-qchisq(p_cover,2)
      s_s<-solve(s_matr[1:2,1:2])
      r_par<-est_rotate_par(s_s,mychisq)
      
      p1<-tryCatch({get_ellipse_p2(mychisq,s_matr[1:2,1:2],c(p_u1,2*p_u1)-k_hat[1:2],p_sigma)},error=function(e){"e"})
      
      if(is.null(p1)|is.na(p1)|identical(p1,"e")){p1<-0}
      if(p1<0|p1>1){p1<-0}
      if(p1<1e-9){
        warning(paste0(name,": normal posterior distribution is used (possible outlier SNP, p1 < 1e-9)."))
        return(mvrnorm(n,k_hat,s_matr))
      }
      
      f10<-function(n,r_par,k_hat,p_u1,p_sigma){
        sample_n<-get_rotate_sample(r_par,round(n*1.28*5),k_hat[1:2])
        p_n<-mvdnorm3(sample_n,c(p_u1,2*p_u1),p_sigma)
        max(p_n)*1.2
      }
      f1<-function(n,r_par,k_hat,p_u1,p_sigma,s_p1){
        out<-NULL
        while(T){
          sample_n<-get_rotate_sample(r_par,round(n*1.28*5),k_hat[1:2])
          p_n<-mvdnorm3(sample_n,c(p_u1,2*p_u1),p_sigma)
          p_n<-p_n/s_p1
          p_n[p_n>=1]<-1
          id<-as.logical(rbinom(length(p_n),1,p_n))
          out<-rbind(out,sample_n[id,])
          if(nrow(out)>=n){break}
        }
        return(out[1:n,])
      }
      s_p1<-f10(n_f1,r_par,k_hat,p_u1,p_sigma)
      
      sample1<-f1(n_f1,r_par,k_hat,p_u1,p_sigma,s_p1)
      
      p2<-mvdnorm3(k_hat,sample1,s_matr)
      p31<-mean(p2)
      p32<-mvdnorm3(k_hat,rbind(c(0,0)),s_matr)
      p10.0<-(1-p_p0)*p1/((1-p_p0)*p1+p_p0)
      p10<-p10.0*p31/(p10.0*p31+(1-p10.0)*p32)
      if(p10<0.05){
        warning(paste0(name,": normal posterior distribution is used (possible outlier SNP, p10 < 0.05)."))
        return(mvrnorm(n,k_hat,s_matr))
      }
      
      f20<-function(n1,r_par,k_hat,p_u1,p_sigma,s_matr,s_p1){
        sample1<-f1(n1,r_par,k_hat,p_u1,p_sigma,s_p1)
        p<-mvdnorm3(k_hat,sample1,s_matr)
        max(p)*1.2
      }
      f2<-function(n1,n2,r_par,k_hat,p_u1,p_sigma,s_matr,s_p1,s_p2){
        out<-NULL
        while(T){
          sample1<-f1(n1,r_par,k_hat,p_u1,p_sigma,s_p1)
          p<-mvdnorm3(k_hat,sample1,s_matr)
          p<-p/s_p2
          p[p>=1]<-1
          id<-as.logical(rbinom(length(p),1,p))
          out<-rbind(out,rbind(sample1[id,]))
          if(nrow(out)>=n2){break}
        }
        return(out[1:n2,])
      }
      s_p2<-f20(n_f1,r_par,k_hat,p_u1,p_sigma,s_matr,s_p1)
      
      sample2<-f2(round(n*p10),round(n*p10),r_par,k_hat,p_u1,p_sigma,s_matr,s_p1,s_p2)
      
      out<-rbind(sample2,matrix(0,nrow=n-round(n*p10),ncol=2))
      return(out[sample(1:n,n,replace=F),])
    }
    crt_signk<-function(k_hat,p1_sp,p2_sp,model_u2,Egger){
      stopifnot(!anyNA(k_hat))
      if(identical(Egger,"auto")){Egger<-T}
      j<-p1_sp==1&p2_sp==0&model_u2&Egger
      if(length(j)==0){j<-F}
      if(j){
        egger<-T
      }else{
        egger<-F
      }
      if(egger){
        sign_k1<-as.integer(2*(as.numeric(k_hat[,1]>=0)-0.5))
        sign_k2<-as.integer(2*(as.numeric(k_hat[,2]>=0)-0.5))
      }else{
        sign_k1<-sign_k2<-rep(1L,nrow(k_hat))
      }
      return(list(sign_k1=sign_k1,sign_k2=sign_k2))
    }
    
    if(is.null(data_p3_opt)){
      if(is.null(data_p3)){
        #suit<-apply(g,2,FUN=function(x){identical(as.numeric(sort(unique(x))),c(0,1,2))})
        #if(sum(suit)==0){stop("No SNP has exactly 3 unique values, or SNPs are not encoded as 0,1,2.")}
        #g_suit<-myselect(g,which(suit))
        #g_name<-colnames(g);rm(g)
        #suit2<-apply(g_suit,2,FUN=function(x){sum(table(x)<control_p3$n_case_lower)==0})
        #suit[suit==T]<-suit2
        if(ncol(g)<control_p3$n_snp_limit){
          stop("The number of SNPs < n_snp_limit.")
        }
        
        #g_suit<-myselect(g_suit,which(suit2))
        #start_suit<-start[which(suit)]
        #twosample_data<-update_twosample_data(twosample_data,which(suit),c_inherit)
        #if(c_inherit){
        #  c_suit<-c
        #}else{
        #  c_suit<-c[which(suit)]
        #}
        g_suit<-g;rm(g)
        c_suit<-c;rm(c)
        start_suit<-start;rm(start)
        gc()
        if(dt){cat(ncol(g_suit),"SNPs are suitable for Procedure 3.\r\n")}
        
        if(is.null(data_p3_k)){
          if(dt){cat("Start preparing stage_1 data required for Procedure 3.\r\n")}
          
          my_task<-function(my_loc){
            t0<-Sys.time()
            my_proc<-0
            pid<-my_loc[[2]]
            my_loc<-my_loc[[1]]
            mywarn<-paste("Child process",pid,"done.")
            g_suit<-myselect(g_suit,my_loc)
            start_suit<-start_suit[my_loc]
            
            withCallingHandlers(
              {
                n1<-ncol(g_suit)
                k_hat_list<-k_sigma_list<-list()
                
                for(i in 1:n1){
                  if(c_inherit){c_m<-c_suit[[1]]}else{
                    c_m<-c_suit[[my_loc[i]]]
                  }
                  if(identical(as.numeric(length(unique(c_m[,1]))),1)){c_m<-NULL}
                  out<-crt_data2(y=y,x=x,g=g_suit[,i],c=c_m,start=start_suit[[i]],name=colnames(g_suit)[i],
                                 twosample_data=update_twosample_data(twosample_data,my_loc[i],c_inherit),
                                 stage1=T,data_k_hat=NULL,data_k_sigma=NULL,nlminb_control=control_sec$nlminb_control)
                  k_hat_list[[i]]<-out[[1]]
                  k_sigma_list[[i]]<-out[[2]]
                  
                  if(parallel_trace&mc.cores==1){my_moni("SNP",i=i,n=n1)}
                  if(mc.cores!=1&parallel_trace){
                    my_proc<-my_moni2(paste0("Child process ",pid,":"),i,n1,my_proc,time=T,t0=t0)
                  }
                }
                names(k_hat_list)<-names(k_sigma_list)<-colnames(g_suit)
              }
              ,warning=function(w){mywarn<<-c(mywarn,w$message);invokeRestart("muffleWarning")}
            )
            if(T){
              if(length(mywarn)>1|parallel_trace){message_parallel(mywarn)}
            }
            return(list(k_hat_list,k_sigma_list))
          }
          
          mynum<-cut_num(1:ncol(g_suit),mc.cores)
          mycheck<-"pass"
          add_obj_list<-list(var=c("g_suit","start_suit","c_inherit","c_suit","twosample_data","parallel_trace","mc.cores","control_sec"),
                             env=environment())
          exec_base_func<-function(x){
            suppressWarnings(library(MRprollim,quietly=T))
          }
          myfit<-withCallingHandlers({my_parallel(X=mynum,FUN=my_task,mc.cores=mc.cores,PSOCK=PSOCK,dt=dt,
                                                  print_message=parallel_trace,add_obj_list=add_obj_list,exec_base_func=exec_base_func,export_parent_func=T,seed=NULL)},warning=function(w){mycheck<<-w})
          if((!identical(mycheck,"pass"))&mc.cores!=1){
            warning("An error occurred. Output of my_parallel with errors is returned.")
            message(mycheck)
            class(myfit)<-"myerror"
            return(myfit)
          }
          gc()
          mybind_list_parallel<-function(list_parallel){
            out1<-out2<-NULL
            for(i in 1:length(list_parallel)){
              out1<-c(out1,list_parallel[[i]][[1]])
              out2<-c(out2,list_parallel[[i]][[2]])
            }
            return(list(out1,out2))
          }
          myfit<-mybind_list_parallel(myfit)
          k_hat_list<-myfit[[1]]
          k_sigma_list<-myfit[[2]]
        }else{
          k_hat_list<-data_p3_k$k_hat_list
          k_sigma_list<-data_p3_k$k_sigma_list
        }
        
        if(control_p3$inspect_data_p3_k){
          return(list(k_hat_list=k_hat_list,k_sigma_list=k_sigma_list))
        }
        
        wald_p<-rep(NA,length(k_hat_list))
        for(i in 1:length(k_hat_list)){
          wald_p[i]<-mywald_test(k_hat_list[[i]][1:2],k_sigma_list[[i]][1:2,1:2])[2]
        }
        
        #user input
        u_input<-NULL
        if(is.character(control_p3$p_snp)){
          print(formatC(quantile(wald_p,c(0,0.25,0.5,0.75,1),na.rm=T),format="e"),quote=FALSE)
          while(T){
            u_input<-eval(parse(text=readline("Please select a cutoff P value:")))
            u_input2<-readline(paste(sum(wald_p<u_input,na.rm=T),"SNPs remained. Continue? (y/n):"))
            if(identical(u_input2,"y")){break}
          }
          loc_pval<-which(wald_p<u_input)
        }else{
          loc_pval<-which(wald_p<control_p3$p_snp)
          if(dt){cat(length(loc_pval),"SNPs satisfy p_snp.\r\n")}
        }
        
        if(length(loc_pval)<control_p3$n_snp_limit){
          warning("The number of remaining SNPs < n_snp_limit. Already extracted data are returned.")
          return(list(k_hat_list=k_hat_list,k_sigma_list=k_sigma_list))
        }
        
        g_suit<-myselect(g_suit,loc_pval)
        start_suit<-start_suit[loc_pval]
        k_hat_list<-k_hat_list[loc_pval]
        k_sigma_list<-k_sigma_list[loc_pval]
        wald_p<-wald_p[loc_pval]
        twosample_data<-update_twosample_data(twosample_data,loc_pval,c_inherit)
        if(!c_inherit){
          c_suit<-c_suit[loc_pval]
        }
        gc()
        
        if(control_p3$s_filter){
          u_input2<-"y"
          while(T){
            if(identical(control_p3$s_k_a_limit,"auto")){
              k_hat<-matrix(unlist(k_hat_list),ncol=2,byrow=T)
              control_p3$s_k_a_limit<-apply(k_hat,2,FUN=sd)*control_p3$auto_k_limit
            }
            
            fi<-rep(NA,length(k_hat_list))
            for(i in 1:length(fi)){
              se1<-sqrt(k_sigma_list[[i]][1,1])
              se2<-sqrt(k_sigma_list[[i]][2,2])
              fi[i]<-(se1*qnorm(0.975)/abs(k_hat_list[[i]][1])<control_p3$s_k_r_limit[1]|se1<control_p3$s_k_a_limit[1])&
                (se2*qnorm(0.975)/abs(k_hat_list[[i]][2])<control_p3$s_k_r_limit[2]|se2<control_p3$s_k_a_limit[2])
            }
            
            loc<-fi
            loc1<-which(loc)
            if(control_p3$s_filter_ask){
              cat(sum(!loc),"SNPs will be removed due to large standard errors of k_hat.\r\n")
              cat(sum(loc),"SNPs will remain for the following analysis.\r\n")
              u_input2<-readline("Continue? (y/n):")
              if(identical(u_input2,"n")){
                #auto<-"auto"
                control_p3$s_k_a_limit<-eval(parse(text=readline("s_k_a_limit:")))
                control_p3$auto_k_limit<-eval(parse(text=readline("auto_k_limit:")))
                control_p3$s_k_r_limit<-eval(parse(text=readline("s_k_r_limit:")))
              }
            }
            if(identical(u_input2,"y")){break}
          }
          g_suit<-myselect(g_suit,loc1)
          start_suit<-start_suit[loc1]
          k_hat_list<-k_hat_list[loc1]
          k_sigma_list<-k_sigma_list[loc1]
          wald_p<-wald_p[loc1]
          twosample_data<-update_twosample_data(twosample_data,loc1,c_inherit)
          if(!c_inherit){
            c_suit<-c_suit[loc1]
          }
          gc()
          
          if(dt){cat(sum(!loc),"SNPs are removed due to large standard errors of k_hat.\r\n")}
          if(sum(loc)<control_p3$n_snp_limit){
            warning("The number of remaining SNPs < n_snp_limit. Remaining data are returned.")
            return(list(k_hat_list=k_hat_list,k_sigma_list=k_sigma_list))
          }
        }
        
        if(dt){cat("Start preparing stage_2 data required for Procedure 3.\r\n")}
        
        my_task<-function(my_loc){
          t0<-Sys.time()
          my_proc<-0
          pid<-my_loc[[2]]
          my_loc<-my_loc[[1]]
          mywarn<-paste("Child process",pid,"done.")
          g_suit<-myselect(g_suit,my_loc)
          start_suit<-start_suit[my_loc]
          k_hat_list1<-k_hat_list[my_loc]
          k_sigma_list1<-k_sigma_list[my_loc]
          
          withCallingHandlers(
            {
              n1<-ncol(g_suit)
              m_hat<-k_hat<-matrix(NA,ncol=2,nrow=n1)
              sigma1_matr<-sigma2_matr<-matrix(NA,ncol=3,nrow=n1)
              mk_sigma<-list()
              
              for(i in 1:n1){
                if(c_inherit){c_m<-c_suit[[1]]}else{
                  c_m<-c_suit[[my_loc[i]]]
                }
                if(identical(as.numeric(length(unique(c_m[,1]))),1)){c_m<-NULL}
                out<-crt_data2(y=y,x=x,g=g_suit[,i],c=c_m,start=start_suit[[i]],name=colnames(g_suit)[i],
                               twosample_data=update_twosample_data(twosample_data,my_loc[i],c_inherit),
                               stage1=F,data_k_hat=k_hat_list1[[i]],data_k_sigma=k_sigma_list1[[i]],nlminb_control=control_sec$nlminb_control)
                m_hat[i,]<-out$m_hat
                sigma1_matr[i,]<-out$m_sigma
                k_hat[i,]<-out$k_hat
                sigma2_matr[i,]<-out$k_sigma
                mk_sigma[[i]]<-out$mk_sigma
                
                if(parallel_trace&mc.cores==1){my_moni("SNP",i=i,n=n1)}
                if(mc.cores!=1&parallel_trace){
                  my_proc<-my_moni2(paste0("Child process ",pid,":"),i,n1,my_proc,time=T,t0=t0)
                }
              }
              rownames(m_hat)<-rownames(k_hat)<-rownames(sigma1_matr)<-rownames(sigma2_matr)<-names(mk_sigma)<-colnames(g_suit)
            }
            ,warning=function(w){mywarn<<-c(mywarn,w$message);invokeRestart("muffleWarning")}
          )
          if(T){
            if(length(mywarn)>1|parallel_trace){message_parallel(mywarn)}
          }
          return(list(m_hat,sigma1_matr,k_hat,sigma2_matr,mk_sigma))
        }
        mynum<-cut_num(1:ncol(g_suit),mc.cores)
        add_obj_list<-list(var=c("g_suit","start_suit","k_hat_list","k_sigma_list","c_inherit","c_suit","twosample_data","parallel_trace","mc.cores","control_sec"),
                           env=environment())
        exec_base_func<-function(x){
          suppressWarnings(library(MRprollim,quietly=T))
        }
        mycheck<-"pass"
        myfit<-withCallingHandlers({my_parallel(X=mynum,FUN=my_task,mc.cores=mc.cores,PSOCK=PSOCK,dt=dt,
                                                print_message=parallel_trace,add_obj_list=add_obj_list,exec_base_func=exec_base_func,export_parent_func=T,seed=NULL)},warning=function(w){mycheck<<-w})
        
        if((!identical(mycheck,"pass"))&mc.cores!=1){
          warning("An error occurred. Output of my_parallel with errors is returned.")
          message(mycheck)
          class(myfit)<-"myerror"
          return(myfit)
        }
        mybind_list_parallel<-function(list_parallel){
          out1<-out2<-out3<-out4<-out5<-NULL
          for(i in 1:length(list_parallel)){
            out1<-rbind(out1,list_parallel[[i]][[1]])
            out2<-rbind(out2,list_parallel[[i]][[2]])
            out3<-rbind(out3,list_parallel[[i]][[3]])
            out4<-rbind(out4,list_parallel[[i]][[4]])
            out5<-c(out5,list_parallel[[i]][[5]])
          }
          return(list(out1,out2,out3,out4,out5))
        }
        myfit<-mybind_list_parallel(myfit)
        m_hat<-myfit[[1]]
        sigma1_matr<-myfit[[2]]
        k_hat<-myfit[[3]]
        sigma2_matr<-myfit[[4]]
        mk_sigma_list<-myfit[[5]]
        
        loc<-which(!is.na(m_hat[,1]))
        m_hat<-myselect(m_hat,loc,"r")
        k_hat<-myselect(k_hat,loc,"r")
        sigma1_matr<-myselect(sigma1_matr,loc,"r")
        sigma2_matr<-myselect(sigma2_matr,loc,"r")
        mk_sigma_list<-mk_sigma_list[loc]
        wald_p<-wald_p[loc]
        gc()
        if(nrow(m_hat)<control_p3$n_snp_limit){
          warning("The number of SNPs < n_snp_limit. Some of the SNPs failed to give the data required for Procedure 3. Already extracted data are returned.")
          return(list(m_hat=m_hat,m_sigma=sigma1_matr,k_hat=k_hat,k_sigma=sigma2_matr,mk_sigma_list=mk_sigma_list,wald_p=wald_p,u_input=u_input))
        }
      }else{
        m_hat<-data_p3$m_hat
        sigma1_matr<-data_p3$m_sigma
        k_hat<-data_p3$k_hat
        sigma2_matr<-data_p3$k_sigma
        mk_sigma_list<-data_p3$mk_sigma_list
        wald_p<-data_p3$wald_p
        u_input<-data_p3$u_input
      }
      
      if(control_p3$inspect_data_p3){
        return(list(m_hat=m_hat,m_sigma=sigma1_matr,k_hat=k_hat,k_sigma=sigma2_matr,mk_sigma_list=mk_sigma_list,wald_p=wald_p,u_input=u_input))
      }
    }
    
    if(control_p3$nome==F){
      if(is.null(data_p3_opt)){
        #est parameters
        if(is.null(control_est_k_post$f)){
          if(identical(control_est_k_prior$p_snp,"max")){
            k_prior_p<-max(wald_p)
          }else if(identical(control_est_k_prior$p_snp,"previous")){
            if(is.character(control_p3$p_snp)){
              k_prior_p<-u_input
            }else{
              k_prior_p<-control_p3$p_snp
            }
          }else{
            k_prior_p<-control_est_k_prior$p_snp
          }
          
          if(any(k_prior_p<wald_p)){
            message("SNPs with wald_p > control_est_k_prior$p_snp detected. max(wald_p) is used for est_k_prior.")
            k_prior_p<-max(wald_p)
          }
          if(dt){cat("Cutoff P value used for est_k_prior:",k_prior_p,"\r\n")}
          
          if(k_prior_p<control_est_k_prior$p0_cut){
            control_est_k_prior$p0_sp<-0
            cat("p0_sp = 0 is used for est_k_prior.\r\n")
          }
          
          if(control_est_k_prior$auto_s){
            auto_s_out<-apply(k_hat,2,FUN=var)
            stopifnot(!anyNA(auto_s_out))
            control_est_k_prior$start[3:4]<-floor(log(auto_s_out))
          }
          fit_k<-tryCatch({est_k_prior(k_hat,sigma2_matr,NULL,NULL,p_cut=k_prior_p,u1_sp=control_est_k_prior$u1_sp,p0_sp=control_est_k_prior$p0_sp,start=control_est_k_prior$start,
                                       p0_start=control_est_k_prior$p0_start,mc.cores=mc.cores,dt=dt,PSOCK=PSOCK,
                                       nlminb_control=control_est_k_prior$nlminb_control)},error=function(e){e})
          
          if("error"%in%class(fit_k)){
            warning(fit_k)
            if(k_prior_p<1e-9){
              warning("Cutoff P value may be too small for est_k_prior. Normal posterior will be used.")
              fit_k<-list(par=list(p0=0,u1=0,s12=var(k_hat[,1]),s22=var(k_hat[,2]),rho=1),nlm_out=fit_k)
            }else{
              message("Data already extracted are returned.")
              return(list(fit_k_prior=fit_k,m_hat=m_hat,m_sigma=sigma1_matr,k_hat=k_hat,k_sigma=sigma2_matr,mk_sigma_list=mk_sigma_list,wald_p=wald_p,u_input=u_input))
            }
          }
          #if(dt){print(unlist(fit_k$par))}
        }
        
        v2s_matr<-function(x){
          matrix(c(x[1],x[2],x[2],x[3]),2)
        }
        if(dt){cat("Start sampling k posterior distribution.\r\n")}
        gc()
        
        my_task<-function(my_loc){
          t0<-Sys.time()
          my_proc<-0
          pid<-my_loc[[2]]
          my_loc<-my_loc[[1]]
          mywarn<-paste("Child process",pid,"done.")
          
          k_hat_m<-myselect(k_hat,my_loc,"r")
          sigma2_matr_m<-myselect(sigma2_matr,my_loc,"r")
          n1<-nrow(k_hat_m)
          post_sample_k1<-post_sample_k2<-matrix(0,nrow=n1,ncol=control_est_k_post$n_post)
          p_sigma<-matrix(c(fit_k$par$s12,sqrt(fit_k$par$s12*fit_k$par$s22)*fit_k$par$rho,
                            sqrt(fit_k$par$s12*fit_k$par$s22)*fit_k$par$rho,fit_k$par$s22),2)
          my_name<-rownames(k_hat_m)
          withCallingHandlers(
            {
              for(i in 1:n1){
                if(is.null(control_est_k_post$f)){
                  post<-get_postsample(k_hat_m[i,],v2s_matr(sigma2_matr_m[i,]),p_cover=control_est_k_post$p_cover,n=control_est_k_post$n_post,n_f1=control_est_k_post$n0,p_p0=fit_k$par$p0,p_u1=fit_k$par$u1,
                                       p_sigma=p_sigma,name=my_name[i])
                }else{
                  post<-control_est_k_post$f(k_hat_m[i,],v2s_matr(sigma2_matr_m[i,]))
                }
                post_sample_k1[i,]<-post[,1]
                post_sample_k2[i,]<-post[,2]
                if(mc.cores==1&parallel_trace){my_moni("SNP",i=i,n=n1)}
                if(mc.cores!=1&parallel_trace){
                  my_proc<-my_moni2(paste0("Child process ",pid,":"),i,n1,my_proc,time=T,t0=t0)
                }
              }
            },
            warning=function(w){mywarn<<-c(mywarn,w$message);invokeRestart("muffleWarning")})
          if(T){
            #if(length(mywarn)>1|parallel_trace){message_parallel(mywarn)}
            if(parallel_trace){message_parallel(mywarn[1])}
          }
          return(list(post_sample_k1,post_sample_k2,mywarn))
        }
        
        mynum<-cut_num(1:nrow(m_hat),mc.cores)
        add_obj_list<-list(var=c("k_hat","sigma2_matr","fit_k","control_est_k_post","parallel_trace","mc.cores"),
                           env=environment())
        exec_base_func<-function(x){
          suppressWarnings(library(MRprollim,quietly=T))
        }
        mycheck<-"pass"
        myfit<-withCallingHandlers({my_parallel(X=mynum,FUN=my_task,mc.cores=mc.cores,PSOCK=PSOCK,dt=dt,
                                                print_message=parallel_trace,add_obj_list=add_obj_list,exec_base_func=exec_base_func,export_parent_func=T,seed=round(runif(1,1,.Machine$integer.max)))},warning=function(w){mycheck<<-w})
        
        if((!identical(mycheck,"pass"))&mc.cores!=1){
          warning("An error occurred. Output of my_parallel with errors is returned.")
          message(mycheck)
          class(myfit)<-"myerror"
          return(myfit)
        }
        mybind_list_parallel<-function(list_parallel){
          out1<-out2<-out3<-NULL
          for(i in 1:length(list_parallel)){
            out1<-rbind(out1,list_parallel[[i]][[1]])
            out2<-rbind(out2,list_parallel[[i]][[2]])
            out3<-c(out3,list_parallel[[i]][[3]])
          }
          return(list(out1,out2,out3))
        }
        myfit<-mybind_list_parallel(myfit)
        post_sample_k1<-myfit[[1]]
        post_sample_k2<-myfit[[2]]
        post_messages<-myfit[[3]]
        if(length(post_messages)>mc.cores){message("Some SNPs may be outliers. Normal posterior is used. See post_messages.")}
        
        prepare_p3.1<-function(k_hat,post_sample_k1,post_sample_k2,mk_sigma_list){
          out1<-out2<-matrix(NA,nrow=nrow(post_sample_k1),ncol=ncol(post_sample_k1))
          out3<-list()
          for(i in 1:nrow(k_hat)){
            s<-mk_sigma_list[[i]][1:2,3:4]%*%solve(mk_sigma_list[[i]][3:4,3:4])
            x<-cbind(k_hat[i,1]-post_sample_k1[i,],k_hat[i,2]-post_sample_k2[i,])%*%t(s)
            out1[i,]<-x[,1]
            out2[i,]<-x[,2]
            out3[[i]]<-mk_sigma_list[[i]][1:2,1:2]-s%*%mk_sigma_list[[i]][3:4,1:2]
          }
          return(list(f1_matr=out1,f2_matr=out2,sigma_prime_list=out3))
        }
        data_opt_p3.1<-prepare_p3.1(k_hat,post_sample_k1,post_sample_k2,mk_sigma_list)
      }else{
        m_hat<-data_p3_opt$m_hat
        sigma1_matr<-data_p3_opt$m_sigma
        k_hat<-data_p3_opt$k_hat
        sigma2_matr<-data_p3_opt$k_sigma
        mk_sigma_list<-data_p3_opt$mk_sigma_list
        post_sample_k1<-data_p3_opt$post_sample_k1;post_sample_k2<-data_p3_opt$post_sample_k2
        data_opt_p3.1<-data_p3_opt$data_opt_p3.1
        wald_p<-data_p3_opt$wald_p;u_input<-data_p3_opt$u_input;post_messages<-data_p3_opt$post_messages
        fit_k<-data_p3_opt$fit_k
      }
      
      if(control_p3$inspect_data_p3_opt){
        return(list(m_hat=m_hat,
                    m_sigma=sigma1_matr,
                    k_hat=k_hat,
                    k_sigma=sigma2_matr,
                    mk_sigma_list=mk_sigma_list,
                    post_sample_k1=post_sample_k1,post_sample_k2=post_sample_k2,
                    data_opt_p3.1=data_opt_p3.1,
                    wald_p=wald_p,u_input=u_input,post_messages=post_messages,
                    fit_k=fit_k))
      }
      
      myEgger0<-control_p3$Egger
      signk<-crt_signk(k_hat,p1_sp,p2_sp,model_u2,myEgger0)
      if(identical(myEgger0,"auto")){
        Egger_info<-"NULL"
      }else{
        if(myEgger0){Egger_info<-"Egger"}else{
          Egger_info<-"intercept"
        }
      }
      if(control_p3$auto_s1){
        beta_start[5]<-floor(log(var(m_hat[,1])))
        control_gs$lower[5]<-beta_start[5]-control_gs$auto_s1_k
        control_gs$upper[5]<-beta_start[5]+control_gs$auto_s1_k
      }
      
      #check initial value
      if(dt){cat("Check initial value.\r\n")}
      while(T){
        if(identical(p1_sp,1)&identical(p2_sp,0)&is.null(r_sp)&identical(model_u2,T)&identical(myEgger0,"auto")){
          signk<-crt_signk(k_hat,p1_sp,p2_sp,model_u2,F)
          initial_test1<-suppressWarnings(est_proc_cont_p3.1_f(beta_start,m_matrix=m_hat,sign_k1=signk$sign_k1,sign_k2=signk$sign_k2,
                                                               f1_matr=data_opt_p3.1$f1_matr,f2_matr=data_opt_p3.1$f2_matr,sigma_prime_list=data_opt_p3.1$sigma_prime_list, 
                                                               post_sample_k1=post_sample_k1,post_sample_k2=post_sample_k2,
                                                               p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2))
          signk<-crt_signk(k_hat,p1_sp,p2_sp,model_u2,T)
          initial_test2<-suppressWarnings(est_proc_cont_p3.1_f(beta_start,m_matrix=m_hat,sign_k1=signk$sign_k1,sign_k2=signk$sign_k2,
                                                               f1_matr=data_opt_p3.1$f1_matr,f2_matr=data_opt_p3.1$f2_matr,sigma_prime_list=data_opt_p3.1$sigma_prime_list, 
                                                               post_sample_k1=post_sample_k1,post_sample_k2=post_sample_k2,
                                                               p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2))
          if(initial_test1%in%c(NA,NaN,Inf,-Inf)|initial_test2%in%c(NA,NaN,Inf,-Inf)){initial_test<-NA}else{initial_test<-0}
        }else{
          initial_test<-suppressWarnings(est_proc_cont_p3.1_f(beta_start,m_matrix=m_hat,sign_k1=signk$sign_k1,sign_k2=signk$sign_k2,
                                                              f1_matr=data_opt_p3.1$f1_matr,f2_matr=data_opt_p3.1$f2_matr,sigma_prime_list=data_opt_p3.1$sigma_prime_list, 
                                                              post_sample_k1=post_sample_k1,post_sample_k2=post_sample_k2,
                                                              p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2))
        }
        if(ask){
          if(initial_test%in%c(NA,NaN,Inf,-Inf)){
            beta_start<-eval(parse(text=readline("Original beta_start for Procedure 3 is not appropriate.\r\nPlease provide a new one (vector):")))
          }else{break}
        }else{
          if(initial_test%in%c(NA,NaN,Inf,-Inf)){
            warning("beta_start is not appropriate. data_input is returned.")
            return(list(m_hat=m_hat,m_sigma=sigma1_matr,k_hat=k_hat,k_sigma=sigma2_matr,
                        post_sample_k1=post_sample_k1,post_sample_k2=post_sample_k2,mk_sigma_list=mk_sigma_list,wald_p=wald_p,u_input=u_input,post_messages=post_messages))
          }else{break}
        }
      }
      
      #optimizing
      est_proc_cont_p3.1_f0<-function(...){
        out<-MRprollim:::est_proc_cont_p3.1_f(...)
        if(is.na(out)){return(Inf)}
        out
      }
      if(mc.cores>1&control_gs$global_search&control_gs$gs_type=="genoud"){
        eval(parse(text=paste0(".mrp_p3_cl<-parallel::makePSOCKcluster(",mc.cores,")")),envir=.GlobalEnv)
        add_par_list<-check_genoud_par(...,return_list=T)
        assign(".mrp_p3_f1",withWarnings,envir=.GlobalEnv)
        eval(expression(environment(.mrp_p3_f1)<-.GlobalEnv),envir=.GlobalEnv)
        assign(".mrp_p3_f2",est_proc_cont_p3.1_f0,envir=.GlobalEnv)
        eval(expression(environment(.mrp_p3_f2)<-.GlobalEnv),envir=.GlobalEnv)
      }else{
        mycl<-F
      }
      if(dt){cat("Start optimizing.\r\n")}
      auto_c<-1
      fit_gs<-fit_nlminb<-list()
      beta_start0<-beta_start
      stopifnot(length(beta_start0)==7)
      while(T){
        crtb<-crt_beta_loc(beta_start0,p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2)
        beta_start<-crtb[[1]]
        beta_loc<-crtb[[2]]
        if(auto_c!=1){
          if(identical(p1_sp,1)&identical(p2_sp,0)&is.null(r_sp)&identical(model_u2,T)&identical(myEgger0,"auto")){
            signk<-crt_signk(k_hat,p1_sp,p2_sp,model_u2,F)
            initial_test1<-suppressWarnings(est_proc_cont_p3.1_f(beta_start,m_matrix=m_hat,sign_k1=signk$sign_k1,sign_k2=signk$sign_k2,
                                                                 f1_matr=data_opt_p3.1$f1_matr,f2_matr=data_opt_p3.1$f2_matr,sigma_prime_list=data_opt_p3.1$sigma_prime_list, 
                                                                 post_sample_k1=post_sample_k1,post_sample_k2=post_sample_k2,
                                                                 p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2,beta_loc=beta_loc))
            signk<-crt_signk(k_hat,p1_sp,p2_sp,model_u2,T)
            initial_test2<-suppressWarnings(est_proc_cont_p3.1_f(beta_start,m_matrix=m_hat,sign_k1=signk$sign_k1,sign_k2=signk$sign_k2,
                                                                 f1_matr=data_opt_p3.1$f1_matr,f2_matr=data_opt_p3.1$f2_matr,sigma_prime_list=data_opt_p3.1$sigma_prime_list, 
                                                                 post_sample_k1=post_sample_k1,post_sample_k2=post_sample_k2,
                                                                 p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2,beta_loc=beta_loc))
            if(initial_test1%in%c(NA,NaN,Inf,-Inf)|initial_test2%in%c(NA,NaN,Inf,-Inf)){initial_test<-NA}else{initial_test<-0}
          }else{
            initial_test<-suppressWarnings(est_proc_cont_p3.1_f(beta_start,m_matrix=m_hat,sign_k1=signk$sign_k1,sign_k2=signk$sign_k2,
                                                                f1_matr=data_opt_p3.1$f1_matr,f2_matr=data_opt_p3.1$f2_matr,sigma_prime_list=data_opt_p3.1$sigma_prime_list, 
                                                                post_sample_k1=post_sample_k1,post_sample_k2=post_sample_k2,
                                                                p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2,beta_loc=beta_loc))
          }
        }
        if(initial_test%in%c(NA,NaN,Inf,-Inf)){
          fit$code<-"inappropriate initial values for the next model"
          warning("Inappropriate initial values for the next model.")
          return(list(fit=fit,fit_nlminb=fit_nlminb,fit_gs=fit_gs,m_hat=m_hat,m_sigma=sigma1_matr,k_hat=k_hat,k_sigma=sigma2_matr,
                      post_sample_k1=post_sample_k1,post_sample_k2=post_sample_k2,
                      mk_sigma_list=mk_sigma_list,wald_p=wald_p,u_input=u_input,post_messages=post_messages))
        }else{
          if(identical(p1_sp,1)&identical(p2_sp,0)&is.null(r_sp)&identical(model_u2,T)&identical(myEgger0,"auto")){
            signk<-crt_signk(k_hat,p1_sp,p2_sp,model_u2,F)
            fit_int<-tryCatch({nlminb2nlm(nlminb2(f=est_proc_cont_p3.1_f,p=beta_start,m_matrix=m_hat,sign_k1=signk$sign_k1,sign_k2=signk$sign_k2,
                                                  f1_matr=data_opt_p3.1$f1_matr,f2_matr=data_opt_p3.1$f2_matr,sigma_prime_list=data_opt_p3.1$sigma_prime_list, 
                                                  post_sample_k1=post_sample_k1,post_sample_k2=post_sample_k2,
                                                  p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2,beta_loc=beta_loc,
                                                  control=nlminb_c_list))},
                              error=function(e){return(e)})
            fit_nlminb<-c(fit_nlminb,list(fit_int))
            if(control_gs$global_search&control_gs$global_search_EI){
              if(control_gs$gs_type=="genoud"){
                if(mc.cores>1){
                  .mrp_p3_par<-list(beta_start=beta_start,beta_loc=beta_loc,
                                    print.level=genoud_control$print.level,max.generations=genoud_control$max.generations,pop.size.EI=genoud_control$pop.size.EI,
                                    Domains=cbind(control_gs$lower[beta_loc],control_gs$upper[beta_loc]),
                                    solution.tolerance=genoud_control$solution.tolerance,
                                    gradient.check=genoud_control$gradient.check,
                                    control=update_opt_c(genoud_control$optim_control,beta_loc),
                                    BFGSburnin=genoud_control$BFGSburnin,
                                    balance=genoud_control$balance,
                                    wait.generations=genoud_control$wait.generations,hard.generation.limit=genoud_control$hard.generation.limit,
                                    m_hat=m_hat,
                                    signk=signk,
                                    data_opt_p3.1=data_opt_p3.1,
                                    post_sample_k1=post_sample_k1,post_sample_k2=post_sample_k2,
                                    p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2)
                  assign(".mrp_p3_par",.mrp_p3_par,envir=.GlobalEnv)
                  eval(expression(parallel::clusterExport(.mrp_p3_cl,varlist=".mrp_p3_par",envir=.GlobalEnv)),envir=.GlobalEnv)
                  exp0<-join_par(expression(rgenoud::genoud(
                    fn=.mrp_p3_f2,nvars=length(.mrp_p3_par$beta_start),
                    starting.values=.mrp_p3_par$beta_start,
                    print.level=.mrp_p3_par$print.level,max.generations=.mrp_p3_par$max.generations,pop.size=.mrp_p3_par$pop.size.EI,
                    Domains=.mrp_p3_par$Domains,cluster=.mrp_p3_cl,
                    solution.tolerance=.mrp_p3_par$solution.tolerance,
                    gradient.check=.mrp_p3_par$gradient.check,
                    control=.mrp_p3_par$control,
                    BFGSburnin=.mrp_p3_par$BFGSburnin,
                    balance=.mrp_p3_par$balance,
                    wait.generations=.mrp_p3_par$wait.generations,hard.generation.limit=.mrp_p3_par$hard.generation.limit,
                    m_matrix=.mrp_p3_par$m_hat,sign_k1=.mrp_p3_par$signk$sign_k1,sign_k2=.mrp_p3_par$signk$sign_k2,
                    f1_matr=.mrp_p3_par$data_opt_p3.1$f1_matr,f2_matr=.mrp_p3_par$data_opt_p3.1$f2_matr,sigma_prime_list=.mrp_p3_par$data_opt_p3.1$sigma_prime_list, 
                    post_sample_k1=.mrp_p3_par$post_sample_k1,post_sample_k2=.mrp_p3_par$post_sample_k2,
                    p1_sp=.mrp_p3_par$p1_sp,p2_sp=.mrp_p3_par$p2_sp,r_sp=.mrp_p3_par$r_sp,model_u2=.mrp_p3_par$model_u2,beta_loc=.mrp_p3_par$beta_loc
                  )),as.expression(add_par_list))
                  eval(parse(text=paste0(".mrp_p3_fit<-tryCatch({.mrp_p3_f1(",exp0,")},error=function(e){e})")),envir =.GlobalEnv)
                  fit1_int<-get(".mrp_p3_fit",envir=.GlobalEnv)
                }else{
                  fit1_int<-tryCatch({withWarnings(rgenoud::genoud(
                    fn=est_proc_cont_p3.1_f0,nvars=length(beta_start),
                    starting.values=beta_start,
                    print.level=genoud_control$print.level,max.generations=genoud_control$max.generations,pop.size=genoud_control$pop.size.EI,
                    Domains=cbind(control_gs$lower[beta_loc],control_gs$upper[beta_loc]),cluster=mycl,
                    solution.tolerance=genoud_control$solution.tolerance,
                    gradient.check=genoud_control$gradient.check,
                    control=update_opt_c(genoud_control$optim_control,beta_loc),
                    BFGSburnin=genoud_control$BFGSburnin,
                    balance=genoud_control$balance,
                    wait.generations=genoud_control$wait.generations,hard.generation.limit=genoud_control$hard.generation.limit,
                    m_matrix=m_hat,sign_k1=signk$sign_k1,sign_k2=signk$sign_k2,
                    f1_matr=data_opt_p3.1$f1_matr,f2_matr=data_opt_p3.1$f2_matr,sigma_prime_list=data_opt_p3.1$sigma_prime_list, 
                    post_sample_k1=post_sample_k1,post_sample_k2=post_sample_k2,
                    p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2,beta_loc=beta_loc,...
                  ))},error=function(e){return(e)})
                }
              }else{
                fit1_int<-tryCatch({GenSA::GenSA(fn=est_proc_cont_p3.1_f0,par=beta_start,
                                                 lower=control_gs$lower[beta_loc],
                                                 upper=control_gs$upper[beta_loc],
                                                 control=update_GenSA(control_gs$GenSA_control,1),
                                                 m_matrix=m_hat,sign_k1=signk$sign_k1,sign_k2=signk$sign_k2,
                                                 f1_matr=data_opt_p3.1$f1_matr,f2_matr=data_opt_p3.1$f2_matr,sigma_prime_list=data_opt_p3.1$sigma_prime_list, 
                                                 post_sample_k1=post_sample_k1,post_sample_k2=post_sample_k2,
                                                 p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2,beta_loc=beta_loc
                )},error=function(e){return(e)})
              }
              
              fit_gs<-c(fit_gs,list(fit1_int))
              
              if((!"error"%in%class(fit1_int))&(!"error"%in%class(fit_int))){
                if("mywarning"%in%class(fit1_int)){
                  warning_out<-fit1_int$warning[!stringr::str_detect(fit1_int$warning,"(Out of Boundary individual)|(NaNs produced)")]
                  if(length(warning_out)>0){
                    message("Warnings in genoud:\n",paste0(warning_out,collapse="\n"))
                  }
                  fit1_int<-fit1_int$value
                }
                if(fit1_int$value<fit_int$minimum){
                  fit_int<-format_fit1(fit1_int)
                }
              }else{
                message("Errors detected in fit_nlminb or fit_gs.")
              }
            }
            
            signk<-crt_signk(k_hat,p1_sp,p2_sp,model_u2,T)
            fit_egger<-tryCatch({nlminb2nlm(nlminb2(f=est_proc_cont_p3.1_f,p=beta_start,m_matrix=m_hat,sign_k1=signk$sign_k1,sign_k2=signk$sign_k2,
                                                    f1_matr=data_opt_p3.1$f1_matr,f2_matr=data_opt_p3.1$f2_matr,sigma_prime_list=data_opt_p3.1$sigma_prime_list, 
                                                    post_sample_k1=post_sample_k1,post_sample_k2=post_sample_k2,
                                                    p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2,beta_loc=beta_loc,
                                                    control=nlminb_c_list))},
                                error=function(e){return(e)})
            fit_nlminb<-c(fit_nlminb,list(fit_egger))
            if(control_gs$global_search&control_gs$global_search_EI){
              if(control_gs$gs_type=="genoud"){
                if(mc.cores>1){
                  .mrp_p3_par<-list(beta_start=beta_start,beta_loc=beta_loc,
                                    print.level=genoud_control$print.level,max.generations=genoud_control$max.generations,pop.size.EI=genoud_control$pop.size.EI,
                                    Domains=cbind(control_gs$lower[beta_loc],control_gs$upper[beta_loc]),
                                    solution.tolerance=genoud_control$solution.tolerance,
                                    gradient.check=genoud_control$gradient.check,
                                    control=update_opt_c(genoud_control$optim_control,beta_loc),
                                    BFGSburnin=genoud_control$BFGSburnin,
                                    balance=genoud_control$balance,
                                    wait.generations=genoud_control$wait.generations,hard.generation.limit=genoud_control$hard.generation.limit,
                                    m_hat=m_hat,
                                    signk=signk,
                                    data_opt_p3.1=data_opt_p3.1,
                                    post_sample_k1=post_sample_k1,post_sample_k2=post_sample_k2,
                                    p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2)
                  assign(".mrp_p3_par",.mrp_p3_par,envir=.GlobalEnv)
                  eval(expression(parallel::clusterExport(.mrp_p3_cl,varlist=".mrp_p3_par",envir=.GlobalEnv)),envir=.GlobalEnv)
                  exp0<-join_par(expression(rgenoud::genoud(
                    fn=.mrp_p3_f2,nvars=length(.mrp_p3_par$beta_start),
                    starting.values=.mrp_p3_par$beta_start,
                    print.level=.mrp_p3_par$print.level,max.generations=.mrp_p3_par$max.generations,pop.size=.mrp_p3_par$pop.size.EI,
                    Domains=.mrp_p3_par$Domains,cluster=.mrp_p3_cl,
                    solution.tolerance=.mrp_p3_par$solution.tolerance,
                    gradient.check=.mrp_p3_par$gradient.check,
                    control=.mrp_p3_par$control,
                    BFGSburnin=.mrp_p3_par$BFGSburnin,
                    balance=.mrp_p3_par$balance,
                    wait.generations=.mrp_p3_par$wait.generations,hard.generation.limit=.mrp_p3_par$hard.generation.limit,
                    m_matrix=.mrp_p3_par$m_hat,sign_k1=.mrp_p3_par$signk$sign_k1,sign_k2=.mrp_p3_par$signk$sign_k2,
                    f1_matr=.mrp_p3_par$data_opt_p3.1$f1_matr,f2_matr=.mrp_p3_par$data_opt_p3.1$f2_matr,sigma_prime_list=.mrp_p3_par$data_opt_p3.1$sigma_prime_list, 
                    post_sample_k1=.mrp_p3_par$post_sample_k1,post_sample_k2=.mrp_p3_par$post_sample_k2,
                    p1_sp=.mrp_p3_par$p1_sp,p2_sp=.mrp_p3_par$p2_sp,r_sp=.mrp_p3_par$r_sp,model_u2=.mrp_p3_par$model_u2,beta_loc=.mrp_p3_par$beta_loc
                  )),as.expression(add_par_list))
                  eval(parse(text=paste0(".mrp_p3_fit<-tryCatch({.mrp_p3_f1(",exp0,")},error=function(e){e})")),envir =.GlobalEnv)
                  fit1_egger<-get(".mrp_p3_fit",envir=.GlobalEnv)
                }else{
                  fit1_egger<-tryCatch({withWarnings(rgenoud::genoud(
                    fn=est_proc_cont_p3.1_f0,nvars=length(beta_start),
                    starting.values=beta_start,
                    print.level=genoud_control$print.level,max.generations=genoud_control$max.generations,pop.size=genoud_control$pop.size.EI,
                    Domains=cbind(control_gs$lower[beta_loc],control_gs$upper[beta_loc]),cluster=mycl,
                    solution.tolerance=genoud_control$solution.tolerance,
                    gradient.check=genoud_control$gradient.check,
                    control=update_opt_c(genoud_control$optim_control,beta_loc),
                    BFGSburnin=genoud_control$BFGSburnin,
                    balance=genoud_control$balance,
                    wait.generations=genoud_control$wait.generations,hard.generation.limit=genoud_control$hard.generation.limit,
                    m_matrix=m_hat,sign_k1=signk$sign_k1,sign_k2=signk$sign_k2,
                    f1_matr=data_opt_p3.1$f1_matr,f2_matr=data_opt_p3.1$f2_matr,sigma_prime_list=data_opt_p3.1$sigma_prime_list, 
                    post_sample_k1=post_sample_k1,post_sample_k2=post_sample_k2,
                    p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2,beta_loc=beta_loc,...
                  ))},error=function(e){return(e)})
                }
              }else{
                fit1_egger<-tryCatch({GenSA::GenSA(fn=est_proc_cont_p3.1_f0,par=beta_start,
                                                   lower=control_gs$lower[beta_loc],
                                                   upper=control_gs$upper[beta_loc],
                                                   control=update_GenSA(control_gs$GenSA_control,1),
                                                   m_matrix=m_hat,sign_k1=signk$sign_k1,sign_k2=signk$sign_k2,
                                                   f1_matr=data_opt_p3.1$f1_matr,f2_matr=data_opt_p3.1$f2_matr,sigma_prime_list=data_opt_p3.1$sigma_prime_list, 
                                                   post_sample_k1=post_sample_k1,post_sample_k2=post_sample_k2,
                                                   p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2,beta_loc=beta_loc
                )},error=function(e){return(e)})
              }
              
              fit_gs<-c(fit_gs,list(fit1_egger))
              
              if((!"error"%in%class(fit1_egger))&(!"error"%in%class(fit_egger))){
                if("mywarning"%in%class(fit1_egger)){
                  warning_out<-fit1_egger$warning[!stringr::str_detect(fit1_egger$warning,"(Out of Boundary individual)|(NaNs produced)")]
                  if(length(warning_out)>0){
                    message("Warnings in genoud:\n",paste0(warning_out,collapse="\n"))
                  }
                  fit1_egger<-fit1_egger$value
                }
                if(fit1_egger$value<fit_egger$minimum){
                  fit_egger<-format_fit1(fit1_egger)
                }
              }else{
                message("Errors detected in fit_nlminb or fit_gs.")
              }
            }
            
            if("error"%in%class(fit_int)|"error"%in%class(fit_egger)){
              warning("Errors detected in fit_int or fit_egger. Data_input and nlminb_output are returned.")
              return(list(fit_nlminb=fit_nlminb,fit_gs=fit_gs,m_hat=m_hat,m_sigma=sigma1_matr,k_hat=k_hat,k_sigma=sigma2_matr,
                          post_sample_k1=post_sample_k1,post_sample_k2=post_sample_k2,
                          mk_sigma_list=mk_sigma_list,wald_p=wald_p,u_input=u_input,post_messages=post_messages))
            }
            if(dt){message("Log likelihood for intercept model: ",-fit_int$minimum)}
            if(dt){message("Log likelihood for Egger model: ",-fit_egger$minimum)}
            if(fit_egger$minimum<fit_int$minimum){
              fit<-fit_egger;Egger_info<-"Egger";myEgger0<-T
              signk<-crt_signk(k_hat,p1_sp,p2_sp,model_u2,T)
            }else{
              fit<-fit_int;Egger_info<-"intercept";myEgger0<-F
              signk<-crt_signk(k_hat,p1_sp,p2_sp,model_u2,F)
            }
            if(dt){message("Model selection: to ",Egger_info," model.")}
          }else{
            fit<-tryCatch({nlminb2nlm(nlminb2(f=est_proc_cont_p3.1_f,p=beta_start,m_matrix=m_hat,sign_k1=signk$sign_k1,sign_k2=signk$sign_k2,
                                              f1_matr=data_opt_p3.1$f1_matr,f2_matr=data_opt_p3.1$f2_matr,sigma_prime_list=data_opt_p3.1$sigma_prime_list, 
                                              post_sample_k1=post_sample_k1,post_sample_k2=post_sample_k2,
                                              p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2,beta_loc=beta_loc,
                                              control=nlminb_c_list))},
                          error=function(e){return(e)})
            fit_nlminb<-c(fit_nlminb,list(fit))
            if(control_gs$global_search){
              if(identical(p1_sp,1)&identical(p2_sp,0)){
                pop.size<-genoud_control$pop.size.EI
                maxit.exc<-1
              }else{
                pop.size<-genoud_control$pop.size
                maxit.exc<-2
              }
              
              if(control_gs$gs_type=="genoud"){
                if(mc.cores>1){
                  .mrp_p3_par<-list(beta_start=beta_start,beta_loc=beta_loc,
                                    print.level=genoud_control$print.level,max.generations=genoud_control$max.generations,pop.size=pop.size,
                                    Domains=cbind(control_gs$lower[beta_loc],control_gs$upper[beta_loc]),
                                    solution.tolerance=genoud_control$solution.tolerance,
                                    gradient.check=genoud_control$gradient.check,
                                    control=update_opt_c(genoud_control$optim_control,beta_loc),
                                    BFGSburnin=genoud_control$BFGSburnin,
                                    balance=genoud_control$balance,
                                    wait.generations=genoud_control$wait.generations,hard.generation.limit=genoud_control$hard.generation.limit,
                                    m_hat=m_hat,
                                    signk=signk,
                                    data_opt_p3.1=data_opt_p3.1,
                                    post_sample_k1=post_sample_k1,post_sample_k2=post_sample_k2,
                                    p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2)
                  assign(".mrp_p3_par",.mrp_p3_par,envir=.GlobalEnv)
                  eval(expression(parallel::clusterExport(.mrp_p3_cl,varlist=".mrp_p3_par",envir=.GlobalEnv)),envir=.GlobalEnv)
                  exp0<-join_par(expression(rgenoud::genoud(
                    fn=.mrp_p3_f2,nvars=length(.mrp_p3_par$beta_start),
                    starting.values=.mrp_p3_par$beta_start,
                    print.level=.mrp_p3_par$print.level,max.generations=.mrp_p3_par$max.generations,pop.size=.mrp_p3_par$pop.size,
                    Domains=.mrp_p3_par$Domains,cluster=.mrp_p3_cl,
                    solution.tolerance=.mrp_p3_par$solution.tolerance,
                    gradient.check=.mrp_p3_par$gradient.check,
                    control=.mrp_p3_par$control,
                    BFGSburnin=.mrp_p3_par$BFGSburnin,
                    balance=.mrp_p3_par$balance,
                    wait.generations=.mrp_p3_par$wait.generations,hard.generation.limit=.mrp_p3_par$hard.generation.limit,
                    m_matrix=.mrp_p3_par$m_hat,sign_k1=.mrp_p3_par$signk$sign_k1,sign_k2=.mrp_p3_par$signk$sign_k2,
                    f1_matr=.mrp_p3_par$data_opt_p3.1$f1_matr,f2_matr=.mrp_p3_par$data_opt_p3.1$f2_matr,sigma_prime_list=.mrp_p3_par$data_opt_p3.1$sigma_prime_list, 
                    post_sample_k1=.mrp_p3_par$post_sample_k1,post_sample_k2=.mrp_p3_par$post_sample_k2,
                    p1_sp=.mrp_p3_par$p1_sp,p2_sp=.mrp_p3_par$p2_sp,r_sp=.mrp_p3_par$r_sp,model_u2=.mrp_p3_par$model_u2,beta_loc=.mrp_p3_par$beta_loc
                  )),as.expression(add_par_list))
                  eval(parse(text=paste0(".mrp_p3_fit<-tryCatch({.mrp_p3_f1(",exp0,")},error=function(e){e})")),envir =.GlobalEnv)
                  fit1<-get(".mrp_p3_fit",envir=.GlobalEnv)
                }else{
                  fit1<-tryCatch({withWarnings(rgenoud::genoud(
                    fn=est_proc_cont_p3.1_f0,nvars=length(beta_start),
                    starting.values=beta_start,
                    print.level=genoud_control$print.level,max.generations=genoud_control$max.generations,pop.size=pop.size,
                    Domains=cbind(control_gs$lower[beta_loc],control_gs$upper[beta_loc]),cluster=mycl,
                    solution.tolerance=genoud_control$solution.tolerance,
                    gradient.check=genoud_control$gradient.check,
                    control=update_opt_c(genoud_control$optim_control,beta_loc),
                    BFGSburnin=genoud_control$BFGSburnin,
                    balance=genoud_control$balance,
                    wait.generations=genoud_control$wait.generations,hard.generation.limit=genoud_control$hard.generation.limit,
                    m_matrix=m_hat,sign_k1=signk$sign_k1,sign_k2=signk$sign_k2,
                    f1_matr=data_opt_p3.1$f1_matr,f2_matr=data_opt_p3.1$f2_matr,sigma_prime_list=data_opt_p3.1$sigma_prime_list, 
                    post_sample_k1=post_sample_k1,post_sample_k2=post_sample_k2,
                    p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2,beta_loc=beta_loc,...
                  ))},error=function(e){return(e)})
                }
              }else{
                fit1<-tryCatch({GenSA::GenSA(fn=est_proc_cont_p3.1_f0,par=beta_start,
                                             lower=control_gs$lower[beta_loc],
                                             upper=control_gs$upper[beta_loc],
                                             control=update_GenSA(control_gs$GenSA_control,maxit.exc),
                                             m_matrix=m_hat,sign_k1=signk$sign_k1,sign_k2=signk$sign_k2,
                                             f1_matr=data_opt_p3.1$f1_matr,f2_matr=data_opt_p3.1$f2_matr,sigma_prime_list=data_opt_p3.1$sigma_prime_list, 
                                             post_sample_k1=post_sample_k1,post_sample_k2=post_sample_k2,
                                             p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2,beta_loc=beta_loc
                )},error=function(e){return(e)})
              }
              
              fit_gs<-c(fit_gs,list(fit1))
              
              if((!"error"%in%class(fit1))&(!"error"%in%class(fit))){
                if("mywarning"%in%class(fit1)){
                  warning_out<-fit1$warning[!stringr::str_detect(fit1$warning,"(Out of Boundary individual)|(NaNs produced)")]
                  if(length(warning_out)>0){
                    message("Warnings in genoud:\n",paste0(warning_out,collapse="\n"))
                  }
                  fit1<-fit1$value
                }
                if(fit1$value<fit$minimum){
                  fit<-format_fit1(fit1)
                }
              }else{
                message("Errors detected in fit_nlminb or fit_gs.")
              }
            }
          }
        }
        if("error"%in%class(fit)){
          warning("Errors detected. Already extracted data are returned.")
          return(list(fit=fit,fit_nlminb=fit_nlminb,fit_gs=fit_gs,m_hat=m_hat,m_sigma=sigma1_matr,k_hat=k_hat,k_sigma=sigma2_matr,
                      post_sample_k1=post_sample_k1,post_sample_k2=post_sample_k2,
                      mk_sigma_list=mk_sigma_list,wald_p=wald_p,u_input=u_input,post_messages=post_messages))
        }
        
        fit$estimate0<-fit$estimate
        fit$estimate<-ext_estimate(fit$estimate,beta_start0,beta_loc)
        
        if(control_p3$model_select==F){break}
        diag_out<-diagnose_fit(fit,p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2,lower=control_p3$check_fit_lower,upper=control_p3$check_fit_upper,myEgger=myEgger0,dt=dt)
        if(identical(diag_out,"pass")){
          if((median(sqrt(sigma1_matr[,1]))*s1_cut_k>sqrt(exp(fit$estimate[5])))&!identical(p1_sp,1)){
            message("Model selection: to Egger or intercept model. (s1 may be too small: ",sqrt(exp(fit$estimate[5])),")")
            diag_out<-list(p1_sp=1,p2_sp=0,r_sp=NULL,model_u2=T,myEgger=myEgger0)
          }else{
            break
          }
        }
        
        if(identical(diag_out,"fail")){fit$code<-"NAs detected by diagnose_fit";break}
        p1_sp<-diag_out$p1_sp
        p2_sp<-diag_out$p2_sp
        r_sp<-diag_out$r_sp
        model_u2<-diag_out$model_u2
        myEgger0<-diag_out$myEgger
        signk<-crt_signk(k_hat,p1_sp,p2_sp,model_u2,myEgger0)
        auto_c<-auto_c+1
      }
      if(mc.cores>1&control_gs$global_search&control_gs$gs_type=="genoud"){
        eval(expression({parallel::stopCluster(.mrp_p3_cl);rm(.mrp_p3_cl,.mrp_p3_f1,.mrp_p3_f2,.mrp_p3_par,.mrp_p3_fit)}),envir=.GlobalEnv)
        gc()
      }
      
      #use nlminb to check the gradients
      if(tryCatch({identical(fit$message$type,"GenSA or genoud")},error=function(e){F})){
        if(dt){cat("Genoud or GenSA finds estimates with a higher likelihood.\r\n")}
        crtb<-crt_beta_loc(fit$estimate,p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2)
        beta_start<-crtb[[1]]
        beta_loc<-crtb[[2]]
        fit0<-tryCatch({nlminb2nlm(nlminb2(f=est_proc_cont_p3.1_f,p=beta_start,m_matrix=m_hat,sign_k1=signk$sign_k1,sign_k2=signk$sign_k2,
                                           f1_matr=data_opt_p3.1$f1_matr,f2_matr=data_opt_p3.1$f2_matr,sigma_prime_list=data_opt_p3.1$sigma_prime_list, 
                                           post_sample_k1=post_sample_k1,post_sample_k2=post_sample_k2,
                                           p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2,beta_loc=beta_loc,
                                           control=nlminb_c_list))},
                       error=function(e){return(e)})
        fit_nlminb<-c(fit_nlminb,list(fit0))
        if("error"%in%class(fit0)){
          fit$code<-"nlminb outputs errors with the solutions given by genoud or GenSA"
        }else{
          fit<-fit0
          fit$estimate0<-fit$estimate
          fit$estimate<-ext_estimate(fit$estimate,beta_start0,beta_loc)
        }
      }
      
      #in case if model_select=F
      fit$code<-check_fit(fit,p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2,lower=control_p3$check_fit_lower,upper=control_p3$check_fit_upper)
      if(fit$code%in%c(0,1)){
        crtb<-crt_beta_loc(fit$estimate,p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2)
        beta_start<-crtb[[1]]
        beta_loc<-crtb[[2]]
        fit0<-tryCatch({suppressWarnings(nlm(f=est_proc_cont_p3.1_f,p=beta_start,m_matrix=m_hat,sign_k1=signk$sign_k1,sign_k2=signk$sign_k2,
                                             f1_matr=data_opt_p3.1$f1_matr,f2_matr=data_opt_p3.1$f2_matr,sigma_prime_list=data_opt_p3.1$sigma_prime_list, 
                                             post_sample_k1=post_sample_k1,post_sample_k2=post_sample_k2,
                                             p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2,beta_loc=beta_loc,
                                             hessian=nlm_c_list$hessian,
                                             fscale=nlm_c_list$fscale,print.level=nlm_c_list$print.level,ndigit=nlm_c_list$ndigit,
                                             gradtol=nlm_c_list$gradtol,stepmax=nlm_c_list$stepmax,steptol=nlm_c_list$steptol[1],
                                             iterlim=nlm_c_list$iterlim,check.analyticals=nlm_c_list$check.analyticals))},
                       error=function(e){return(e)})
        fit_nlminb<-c(fit_nlminb,list(fit0))
        if("error"%in%class(fit0)){
          fit$code<-"nlm outputs errors for the current solution"
        }else{
          if(anyNA(fit0$gradient)){fit0$code<-"NA_gradient"}
          fit0$message<-paste0("nlm code: ",fit0$code)
          if(fit0$code%in%c(1,2,3)){fit0$code<-0}else{fit0$code<-1}
          fit0$estimate0<-fit0$estimate
          fit0$estimate<-ext_estimate(fit0$estimate,beta_start0,beta_loc)
          fit0$code<-check_fit(fit0,p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2,lower=control_p3$check_fit_lower,upper=control_p3$check_fit_upper)
          fit<-fit0
        }
      }
      if(!fit$code%in%c(0)){
        message("Abnormal optimization code; already extracted data are returned.")
        return(list(fit=fit,fit_nlminb=fit_nlminb,fit_gs=fit_gs,m_hat=m_hat,m_sigma=sigma1_matr,k_hat=k_hat,k_sigma=sigma2_matr,
                    post_sample_k1=post_sample_k1,post_sample_k2=post_sample_k2,
                    mk_sigma_list=mk_sigma_list,wald_p=wald_p,u_input=u_input,post_messages=post_messages))
      }
      
      #bootstrap
      vcov_boot<-out_boot<-NULL
      if(boot_se){
        if(dt){cat("Start bootstrapping.\r\n")}
        n1<-nrow(m_hat)
        
        my_task<-function(n_rep){
          t0<-Sys.time()
          pid<-n_rep[[2]]
          n_tar<-n_rep[[3]]
          n_rep<-n_rep[[1]]
          my_proc<-0
          out_boot<-NULL
          
          crtb<-crt_beta_loc(fit$estimate,p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2)
          beta_start<-crtb[[1]]
          beta_loc<-crtb[[2]]
          for(i in 1:n_rep){
            id<-sample(1:n1,n1,replace=T)
            signk<-crt_signk(k_hat[id,],p1_sp,p2_sp,model_u2,myEgger0)
            initial_test<-suppressWarnings(est_proc_cont_p3.1_f(beta_start,m_matrix=m_hat[id,],sign_k1=signk$sign_k1,sign_k2=signk$sign_k2,
                                                                f1_matr=data_opt_p3.1$f1_matr[id,],f2_matr=data_opt_p3.1$f2_matr[id,],sigma_prime_list=data_opt_p3.1$sigma_prime_list[id], 
                                                                post_sample_k1=post_sample_k1[id,],post_sample_k2=post_sample_k2[id,],
                                                                p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2,beta_loc=beta_loc))
            if(initial_test%in%c(NA,NaN,Inf,-Inf)){fit1<-list(estimate=NA,code=NA)}else{
              fit1<-tryCatch({nlminb2nlm(nlminb2(f=est_proc_cont_p3.1_f,p=beta_start,m_matrix=m_hat[id,],sign_k1=signk$sign_k1,sign_k2=signk$sign_k2,
                                                 f1_matr=data_opt_p3.1$f1_matr[id,],f2_matr=data_opt_p3.1$f2_matr[id,],sigma_prime_list=data_opt_p3.1$sigma_prime_list[id], 
                                                 post_sample_k1=post_sample_k1[id,],post_sample_k2=post_sample_k2[id,],
                                                 p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2,beta_loc=beta_loc,
                                                 control=nlminb_c_list))},
                             error=function(e){return(list(estimate=NA,code=NA))})
              fit1$estimate<-ext_estimate(fit1$estimate,beta_start0,beta_loc)
              fit1$code<-check_fit(fit1,p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2,lower=control_p3$check_fit_lower,upper=control_p3$check_fit_upper)
            }
            if(fit1$code%in%c(0)){
              out_boot<-rbind(out_boot,trans_par_norm(fit1$estimate,p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2,est_type="c"))
            }
            if(parallel_trace&mc.cores==1){my_moni("Bootstrap repeat",i,n_tar)}
            if(mc.cores!=1&parallel_trace){
              my_proc<-my_moni2(paste0("Child process ",pid,":"),i,n_tar,my_proc,time=T,t0=t0)
            }
            if(nrow(out_boot)>=n_tar){break}
          }
          return(out_boot)
        }
        
        mylist<-list()
        for(j in 1:mc.cores){
          mylist[[j]]<-list(even_allo(n_rep_max*n_boot,mc.cores)[j],j,even_allo(n_boot,mc.cores)[j])
        }
        add_obj_list<-list(var=c("n1","k_hat","p1_sp","p2_sp","model_u2","myEgger0","r_sp",
                                 "m_hat","control_p3","post_sample_k1","post_sample_k2",
                                 "data_opt_p3.1","fit","beta_start0"),
                           env=environment())
        exec_base_func<-function(x){
          suppressWarnings(library(MRprollim,quietly=T))
        }
        myfit<-my_parallel(X=mylist,FUN=my_task,mc.cores=mc.cores,PSOCK=PSOCK,dt=dt,
                           print_message=parallel_trace,add_obj_list=add_obj_list,exec_base_func=exec_base_func,export_parent_func=T,seed=round(runif(1,1,.Machine$integer.max)))
        mybind_list_parallel<-function(list){
          out<-NULL
          for(i in 1:length(list)){
            out<-rbind(out,list[[i]])
          }
          return(out)
        }
        out_boot<-mybind_list_parallel(myfit)
        if(nrow(out_boot)<n_boot){warning("The limit n_rep_max*n_boot was reached, but MR-PROLLIM failed to obtain n_boot repeats.")}
        vcov_boot<-cov(out_boot,use="na.or.complete")
      }
      
      #asymptotic
      if(dt){cat("Calculate asymptotic variance.\r\n")}
      beta_norm<-trans_par_norm(fit$estimate,p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2,est_type="c")
      
      vcov<-tryCatch({est_variance_reMP(f=est_proc_cont_p3.1_f,beta_opt=beta_norm,beta_name=names(beta_norm),
                                        sandwich=control_p3$sandwich,hessian=control_p3$hessian,adj_rs=control_p3$adj_rs,
                                        mc.cores=mc.cores,boot_n=control_p3$se_adj_boot_n,parallel_trace=parallel_trace,
                                        m_matrix=m_hat,sign_k1=signk$sign_k1,sign_k2=signk$sign_k2,
                                        f1_matr=data_opt_p3.1$f1_matr,f2_matr=data_opt_p3.1$f2_matr,sigma_prime_list=data_opt_p3.1$sigma_prime_list, 
                                        post_sample_k1=post_sample_k1,post_sample_k2=post_sample_k2,
                                        p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2,individual=F,vcov_est=T)},error=function(e){return(e)})
      
      se_u1<-p_u1<-NA
      if(!"error"%in%class(vcov)){
        se_u1<-tryCatch({sqrt(vcov_sandw["u1","u1"])},error=function(e){NA})
        p_u1<-pnorm(-abs(beta_norm["u1"]/se_u1))*2
        names(p_u1)<-NULL
      }
      
      #output
      if(identical(as.numeric(p1_sp),1)&identical(as.numeric(p2_sp),0)&identical(myEgger0,T)){
        final_Egger_flag<-T
      }else{
        final_Egger_flag<-F
      }
      out<-list(beta_norm=beta_norm,vcov=vcov,vcov_boot=vcov_boot,se_u1_sandwich=se_u1,p_u1_sandwich=p_u1)
      par<-list(p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2,est_type="c",Egger_info=Egger_info,final_Egger_flag=final_Egger_flag)
      out<-list(estimate=out,parameter=par,maxlik=list(fit_final=fit,fit_gs=fit_gs,fit_nlminb=fit_nlminb),prior_est=fit_k,data=list(m_hat=m_hat,m_sigma=sigma1_matr,k_hat=k_hat,
                                                                                                                                    k_sigma=sigma2_matr,post_sample_k1=post_sample_k1,post_sample_k2=post_sample_k2,
                                                                                                                                    mk_sigma_list=mk_sigma_list,wald_p=wald_p,u_input=u_input,post_messages=post_messages,boot_data=out_boot))
      class(out)<-"MR-PROLLIM output"
      if(dt){cat("Random-effects MR-PROLLIM (Procedure 3) finished.\r\n")}
      return(out)
    }
    
    if(control_p3$nome==T){
      if(!is.null(data_p3_opt)){
        m_hat<-data_p3_opt$m_hat;sigma1_matr<-data_p3_opt$m_sigma
        k_hat<-data_p3_opt$k_hat;sigma2_matr<-data_p3_opt$k_sigma
        mk_sigma_list<-data_p3_opt$mk_sigma_list
        wald_p<-data_p3_opt$wald_p;u_input<-data_p3_opt$u_input
      }
      
      if(control_p3$inspect_data_p3_opt){
        return(list(m_hat=m_hat,m_sigma=sigma1_matr,
                    k_hat=k_hat,k_sigma=sigma2_matr,
                    mk_sigma_list=mk_sigma_list,
                    wald_p=wald_p,u_input=u_input))
      }
      
      myEgger0<-control_p3$Egger
      signk<-crt_signk(k_hat,p1_sp,p2_sp,model_u2,myEgger0)
      if(identical(myEgger0,"auto")){
        Egger_info<-"NULL"
      }else{
        if(myEgger0){Egger_info<-"Egger"}else{
          Egger_info<-"intercept"
        }
      }
      if(control_p3$auto_s1){
        beta_start[5]<-floor(log(var(m_hat[,1])))
        control_gs$lower[5]<-beta_start[5]-control_gs$auto_s1_k
        control_gs$upper[5]<-beta_start[5]+control_gs$auto_s1_k
      }
      
      #check initial value
      if(dt){cat("Check initial value.\r\n")}
      while(T){
        if(identical(p1_sp,1)&identical(p2_sp,0)&is.null(r_sp)&identical(model_u2,T)&identical(myEgger0,"auto")){
          signk<-crt_signk(k_hat,p1_sp,p2_sp,model_u2,F)
          initial_test1<-suppressWarnings(est_proc_cont_p3_f(beta_start,m_matrix=m_hat,sigma1_matr=sigma1_matr,k_matrix=k_hat,
                                                             sign_k1=signk$sign_k1,sign_k2=signk$sign_k2,
                                                             p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2))
          signk<-crt_signk(k_hat,p1_sp,p2_sp,model_u2,T)
          initial_test2<-suppressWarnings(est_proc_cont_p3_f(beta_start,m_matrix=m_hat,sigma1_matr=sigma1_matr,k_matrix=k_hat,
                                                             sign_k1=signk$sign_k1,sign_k2=signk$sign_k2,
                                                             p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2))
          if(initial_test1%in%c(NA,NaN,Inf,-Inf)|initial_test2%in%c(NA,NaN,Inf,-Inf)){initial_test<-NA}else{initial_test<-0}
        }else{
          initial_test<-suppressWarnings(est_proc_cont_p3_f(beta_start,m_matrix=m_hat,sigma1_matr=sigma1_matr,k_matrix=k_hat,
                                                            sign_k1=signk$sign_k1,sign_k2=signk$sign_k2,
                                                            p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2))
        }
        if(ask){
          if(initial_test%in%c(NA,NaN,Inf,-Inf)){
            beta_start<-eval(parse(text=readline("Original beta_start for Procedure 3 is not appropriate.\r\nPlease provide a new one (vector):")))
          }else{break}
        }else{
          if(initial_test%in%c(NA,NaN,Inf,-Inf)){
            warning("beta_start is not appropriate. data_input is returned.")
            return(list(m_hat=m_hat,m_sigma=sigma1_matr,k_hat=k_hat,k_sigma=sigma2_matr,
                        mk_sigma_list=mk_sigma_list,wald_p=wald_p,u_input=u_input))
          }else{break}
        }
      }
      
      est_proc_cont_p3_f0<-function(...){
        out<-MRprollim:::est_proc_cont_p3_f(...)
        if(is.na(out)){return(Inf)}
        out
      }
      if(mc.cores>1&control_gs$global_search&control_gs$gs_type=="genoud"){
        eval(parse(text=paste0(".mrp_p3_cl<-parallel::makePSOCKcluster(",mc.cores,")")),envir=.GlobalEnv)
        add_par_list<-check_genoud_par(...,return_list=T)
        assign(".mrp_p3_f1",withWarnings,envir=.GlobalEnv)
        eval(expression(environment(.mrp_p3_f1)<-.GlobalEnv),envir=.GlobalEnv)
        assign(".mrp_p3_f2",est_proc_cont_p3_f0,envir=.GlobalEnv)
        eval(expression(environment(.mrp_p3_f2)<-.GlobalEnv),envir=.GlobalEnv)
      }else{
        mycl<-F
      }
      if(dt){cat("Start optimizing.\r\n")}
      auto_c<-1
      fit_gs<-fit_nlminb<-list()
      beta_start0<-beta_start
      stopifnot(length(beta_start0)==7)
      while(T){
        crtb<-crt_beta_loc(beta_start0,p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2)
        beta_start<-crtb[[1]]
        beta_loc<-crtb[[2]]
        if(auto_c!=1){
          if(identical(p1_sp,1)&identical(p2_sp,0)&is.null(r_sp)&identical(model_u2,T)&identical(myEgger0,"auto")){
            signk<-crt_signk(k_hat,p1_sp,p2_sp,model_u2,F)
            initial_test1<-suppressWarnings(est_proc_cont_p3_f(beta_start,m_matrix=m_hat,sigma1_matr=sigma1_matr,k_matrix=k_hat,
                                                               sign_k1=signk$sign_k1,sign_k2=signk$sign_k2,
                                                               p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2,beta_loc=beta_loc))
            signk<-crt_signk(k_hat,p1_sp,p2_sp,model_u2,T)
            initial_test2<-suppressWarnings(est_proc_cont_p3_f(beta_start,m_matrix=m_hat,sigma1_matr=sigma1_matr,k_matrix=k_hat,
                                                               sign_k1=signk$sign_k1,sign_k2=signk$sign_k2,
                                                               p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2,beta_loc=beta_loc))
            if(initial_test1%in%c(NA,NaN,Inf,-Inf)|initial_test2%in%c(NA,NaN,Inf,-Inf)){initial_test<-NA}else{initial_test<-0}
          }else{
            initial_test<-suppressWarnings(est_proc_cont_p3_f(beta_start,m_matrix=m_hat,sigma1_matr=sigma1_matr,k_matrix=k_hat,
                                                              sign_k1=signk$sign_k1,sign_k2=signk$sign_k2,
                                                              p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2,beta_loc=beta_loc))
          }
        }
        if(initial_test%in%c(NA,NaN,Inf,-Inf)){
          fit$code<-"inappropriate initial values for the next model"
          warning("Inappropriate initial values for the next model.")
          return(list(fit=fit,fit_nlminb=fit_nlminb,fit_gs=fit_gs,m_hat=m_hat,m_sigma=sigma1_matr,k_hat=k_hat,k_sigma=sigma2_matr,
                      mk_sigma_list=mk_sigma_list,wald_p=wald_p,u_input=u_input))
        }else{
          if(identical(p1_sp,1)&identical(p2_sp,0)&is.null(r_sp)&identical(model_u2,T)&identical(myEgger0,"auto")){
            signk<-crt_signk(k_hat,p1_sp,p2_sp,model_u2,F)
            fit_int<-tryCatch({nlminb2nlm(nlminb2(f=est_proc_cont_p3_f,p=beta_start,m_matrix=m_hat,sigma1_matr=sigma1_matr,k_matrix=k_hat,
                                                  sign_k1=signk$sign_k1,sign_k2=signk$sign_k2,
                                                  p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2,beta_loc=beta_loc,
                                                  control=nlminb_c_list))},
                              error=function(e){return(e)})
            fit_nlminb<-c(fit_nlminb,list(fit_int))
            if(control_gs$global_search&control_gs$global_search_EI){
              if(control_gs$gs_type=="genoud"){
                if(mc.cores>1){
                  .mrp_p3_par<-list(beta_start=beta_start,beta_loc=beta_loc,
                                    print.level=genoud_control$print.level,max.generations=genoud_control$max.generations,pop.size.EI=genoud_control$pop.size.EI,
                                    Domains=cbind(control_gs$lower[beta_loc],control_gs$upper[beta_loc]),
                                    solution.tolerance=genoud_control$solution.tolerance,
                                    gradient.check=genoud_control$gradient.check,
                                    control=update_opt_c(genoud_control$optim_control,beta_loc),
                                    BFGSburnin=genoud_control$BFGSburnin,
                                    balance=genoud_control$balance,
                                    wait.generations=genoud_control$wait.generations,hard.generation.limit=genoud_control$hard.generation.limit,
                                    m_hat=m_hat,sigma1_matr=sigma1_matr,k_hat=k_hat,
                                    signk=signk,
                                    p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2)
                  assign(".mrp_p3_par",.mrp_p3_par,envir=.GlobalEnv)
                  eval(expression(parallel::clusterExport(.mrp_p3_cl,varlist=".mrp_p3_par",envir=.GlobalEnv)),envir=.GlobalEnv)
                  exp0<-join_par(expression(rgenoud::genoud(
                    fn=.mrp_p3_f2,nvars=length(.mrp_p3_par$beta_start),
                    starting.values=.mrp_p3_par$beta_start,
                    print.level=.mrp_p3_par$print.level,max.generations=.mrp_p3_par$max.generations,pop.size=.mrp_p3_par$pop.size.EI,
                    Domains=.mrp_p3_par$Domains,cluster=.mrp_p3_cl,
                    solution.tolerance=.mrp_p3_par$solution.tolerance,
                    gradient.check=.mrp_p3_par$gradient.check,
                    control=.mrp_p3_par$control,
                    BFGSburnin=.mrp_p3_par$BFGSburnin,
                    balance=.mrp_p3_par$balance,
                    wait.generations=.mrp_p3_par$wait.generations,hard.generation.limit=.mrp_p3_par$hard.generation.limit,
                    m_matrix=.mrp_p3_par$m_hat,sigma1_matr=.mrp_p3_par$sigma1_matr,k_matrix=.mrp_p3_par$k_hat,
                    sign_k1=.mrp_p3_par$signk$sign_k1,sign_k2=.mrp_p3_par$signk$sign_k2,
                    p1_sp=.mrp_p3_par$p1_sp,p2_sp=.mrp_p3_par$p2_sp,r_sp=.mrp_p3_par$r_sp,model_u2=.mrp_p3_par$model_u2,beta_loc=.mrp_p3_par$beta_loc
                  )),as.expression(add_par_list))
                  eval(parse(text=paste0(".mrp_p3_fit<-tryCatch({.mrp_p3_f1(",exp0,")},error=function(e){e})")),envir =.GlobalEnv)
                  fit1_int<-get(".mrp_p3_fit",envir=.GlobalEnv)
                }else{
                  fit1_int<-tryCatch({withWarnings(rgenoud::genoud(
                    fn=est_proc_cont_p3_f0,nvars=length(beta_start),
                    starting.values=beta_start,
                    print.level=genoud_control$print.level,max.generations=genoud_control$max.generations,pop.size=genoud_control$pop.size.EI,
                    Domains=cbind(control_gs$lower[beta_loc],control_gs$upper[beta_loc]),cluster=mycl,
                    solution.tolerance=genoud_control$solution.tolerance,
                    gradient.check=genoud_control$gradient.check,
                    control=update_opt_c(genoud_control$optim_control,beta_loc),
                    BFGSburnin=genoud_control$BFGSburnin,
                    balance=genoud_control$balance,
                    wait.generations=genoud_control$wait.generations,hard.generation.limit=genoud_control$hard.generation.limit,
                    m_matrix=m_hat,sigma1_matr=sigma1_matr,k_matrix=k_hat,
                    sign_k1=signk$sign_k1,sign_k2=signk$sign_k2,
                    p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2,beta_loc=beta_loc,...
                  ))},error=function(e){return(e)})
                }
              }else{
                fit1_int<-tryCatch({GenSA::GenSA(fn=est_proc_cont_p3_f0,par=beta_start,
                                                 lower=control_gs$lower[beta_loc],
                                                 upper=control_gs$upper[beta_loc],
                                                 control=update_GenSA(control_gs$GenSA_control,1),
                                                 m_matrix=m_hat,sigma1_matr=sigma1_matr,k_matrix=k_hat,
                                                 sign_k1=signk$sign_k1,sign_k2=signk$sign_k2,
                                                 p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2,beta_loc=beta_loc
                )},error=function(e){return(e)})
              }
              
              fit_gs<-c(fit_gs,list(fit1_int))
              
              if((!"error"%in%class(fit1_int))&(!"error"%in%class(fit_int))){
                if("mywarning"%in%class(fit1_int)){
                  warning_out<-fit1_int$warning[!stringr::str_detect(fit1_int$warning,"(Out of Boundary individual)|(NaNs produced)")]
                  if(length(warning_out)>0){
                    message("Warnings in genoud:\n",paste0(warning_out,collapse="\n"))
                  }
                  fit1_int<-fit1_int$value
                }
                if(fit1_int$value<fit_int$minimum){
                  fit_int<-format_fit1(fit1_int)
                }
              }else{
                message("Errors detected in fit_nlminb or fit_gs.")
              }
            }
            
            signk<-crt_signk(k_hat,p1_sp,p2_sp,model_u2,T)
            fit_egger<-tryCatch({nlminb2nlm(nlminb2(f=est_proc_cont_p3_f,p=beta_start,m_matrix=m_hat,sigma1_matr=sigma1_matr,k_matrix=k_hat,
                                                    sign_k1=signk$sign_k1,sign_k2=signk$sign_k2,
                                                    p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2,beta_loc=beta_loc,
                                                    control=nlminb_c_list))},
                                error=function(e){return(e)})
            fit_nlminb<-c(fit_nlminb,list(fit_egger))
            if(control_gs$global_search&control_gs$global_search_EI){
              if(control_gs$gs_type=="genoud"){
                if(mc.cores>1){
                  .mrp_p3_par<-list(beta_start=beta_start,beta_loc=beta_loc,
                                    print.level=genoud_control$print.level,max.generations=genoud_control$max.generations,pop.size.EI=genoud_control$pop.size.EI,
                                    Domains=cbind(control_gs$lower[beta_loc],control_gs$upper[beta_loc]),
                                    solution.tolerance=genoud_control$solution.tolerance,
                                    gradient.check=genoud_control$gradient.check,
                                    control=update_opt_c(genoud_control$optim_control,beta_loc),
                                    BFGSburnin=genoud_control$BFGSburnin,
                                    balance=genoud_control$balance,
                                    wait.generations=genoud_control$wait.generations,hard.generation.limit=genoud_control$hard.generation.limit,
                                    m_hat=m_hat,sigma1_matr=sigma1_matr,k_hat=k_hat,
                                    signk=signk,
                                    p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2)
                  assign(".mrp_p3_par",.mrp_p3_par,envir=.GlobalEnv)
                  eval(expression(parallel::clusterExport(.mrp_p3_cl,varlist=".mrp_p3_par",envir=.GlobalEnv)),envir=.GlobalEnv)
                  exp0<-join_par(expression(rgenoud::genoud(
                    fn=.mrp_p3_f2,nvars=length(.mrp_p3_par$beta_start),
                    starting.values=.mrp_p3_par$beta_start,
                    print.level=.mrp_p3_par$print.level,max.generations=.mrp_p3_par$max.generations,pop.size=.mrp_p3_par$pop.size.EI,
                    Domains=.mrp_p3_par$Domains,cluster=.mrp_p3_cl,
                    solution.tolerance=.mrp_p3_par$solution.tolerance,
                    gradient.check=.mrp_p3_par$gradient.check,
                    control=.mrp_p3_par$control,
                    BFGSburnin=.mrp_p3_par$BFGSburnin,
                    balance=.mrp_p3_par$balance,
                    wait.generations=.mrp_p3_par$wait.generations,hard.generation.limit=.mrp_p3_par$hard.generation.limit,
                    m_matrix=.mrp_p3_par$m_hat,sigma1_matr=.mrp_p3_par$sigma1_matr,k_matrix=.mrp_p3_par$k_hat,
                    sign_k1=.mrp_p3_par$signk$sign_k1,sign_k2=.mrp_p3_par$signk$sign_k2,
                    p1_sp=.mrp_p3_par$p1_sp,p2_sp=.mrp_p3_par$p2_sp,r_sp=.mrp_p3_par$r_sp,model_u2=.mrp_p3_par$model_u2,beta_loc=.mrp_p3_par$beta_loc
                  )),as.expression(add_par_list))
                  eval(parse(text=paste0(".mrp_p3_fit<-tryCatch({.mrp_p3_f1(",exp0,")},error=function(e){e})")),envir =.GlobalEnv)
                  fit1_egger<-get(".mrp_p3_fit",envir=.GlobalEnv)
                }else{
                  fit1_egger<-tryCatch({withWarnings(rgenoud::genoud(
                    fn=est_proc_cont_p3_f0,nvars=length(beta_start),
                    starting.values=beta_start,
                    print.level=genoud_control$print.level,max.generations=genoud_control$max.generations,pop.size=genoud_control$pop.size.EI,
                    Domains=cbind(control_gs$lower[beta_loc],control_gs$upper[beta_loc]),cluster=mycl,
                    solution.tolerance=genoud_control$solution.tolerance,
                    gradient.check=genoud_control$gradient.check,
                    control=update_opt_c(genoud_control$optim_control,beta_loc),
                    BFGSburnin=genoud_control$BFGSburnin,
                    balance=genoud_control$balance,
                    wait.generations=genoud_control$wait.generations,hard.generation.limit=genoud_control$hard.generation.limit,
                    m_matrix=m_hat,sigma1_matr=sigma1_matr,k_matrix=k_hat,
                    sign_k1=signk$sign_k1,sign_k2=signk$sign_k2,
                    p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2,beta_loc=beta_loc,...
                  ))},error=function(e){return(e)})
                }
              }else{
                fit1_egger<-tryCatch({GenSA::GenSA(fn=est_proc_cont_p3_f0,par=beta_start,
                                                   lower=control_gs$lower[beta_loc],
                                                   upper=control_gs$upper[beta_loc],
                                                   control=update_GenSA(control_gs$GenSA_control,1),
                                                   m_matrix=m_hat,sigma1_matr=sigma1_matr,k_matrix=k_hat,
                                                   sign_k1=signk$sign_k1,sign_k2=signk$sign_k2,
                                                   p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2,beta_loc=beta_loc
                )},error=function(e){return(e)})
              }
              
              fit_gs<-c(fit_gs,list(fit1_egger))
              
              if((!"error"%in%class(fit1_egger))&(!"error"%in%class(fit_egger))){
                if("mywarning"%in%class(fit1_egger)){
                  warning_out<-fit1_egger$warning[!stringr::str_detect(fit1_egger$warning,"(Out of Boundary individual)|(NaNs produced)")]
                  if(length(warning_out)>0){
                    message("Warnings in genoud:\n",paste0(warning_out,collapse="\n"))
                  }
                  fit1_egger<-fit1_egger$value
                }
                if(fit1_egger$value<fit_egger$minimum){
                  fit_egger<-format_fit1(fit1_egger)
                }
              }else{
                message("Errors detected in fit_nlminb or fit_gs.")
              }
            }
            
            if("error"%in%class(fit_int)|"error"%in%class(fit_egger)){
              warning("Errors detected in fit_int or fit_egger. Data_input and nlminb_output are returned.")
              return(list(fit_nlminb=fit_nlminb,fit_gs=fit_gs,m_hat=m_hat,m_sigma=sigma1_matr,k_hat=k_hat,k_sigma=sigma2_matr,
                          mk_sigma_list=mk_sigma_list,wald_p=wald_p,u_input=u_input))
            }
            if(dt){message("Log likelihood for intercept model: ",-fit_int$minimum)}
            if(dt){message("Log likelihood for Egger model: ",-fit_egger$minimum)}
            if(fit_egger$minimum<fit_int$minimum){
              fit<-fit_egger;Egger_info<-"Egger";myEgger0<-T
              signk<-crt_signk(k_hat,p1_sp,p2_sp,model_u2,T)
            }else{
              fit<-fit_int;Egger_info<-"intercept";myEgger0<-F
              signk<-crt_signk(k_hat,p1_sp,p2_sp,model_u2,F)
            }
            if(dt){message("Model selection: to ",Egger_info," model.")}
          }else{
            fit<-tryCatch({nlminb2nlm(nlminb2(f=est_proc_cont_p3_f,p=beta_start,m_matrix=m_hat,sigma1_matr=sigma1_matr,k_matrix=k_hat,
                                              sign_k1=signk$sign_k1,sign_k2=signk$sign_k2,
                                              p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2,beta_loc=beta_loc,
                                              control=nlminb_c_list))},
                          error=function(e){return(e)})
            fit_nlminb<-c(fit_nlminb,list(fit))
            if(control_gs$global_search){
              if(identical(p1_sp,1)&identical(p2_sp,0)){
                pop.size<-genoud_control$pop.size.EI
                maxit.exc<-1
              }else{
                pop.size<-genoud_control$pop.size
                maxit.exc<-2
              }
              
              if(control_gs$gs_type=="genoud"){
                if(mc.cores>1){
                  .mrp_p3_par<-list(beta_start=beta_start,beta_loc=beta_loc,
                                    print.level=genoud_control$print.level,max.generations=genoud_control$max.generations,pop.size=pop.size,
                                    Domains=cbind(control_gs$lower[beta_loc],control_gs$upper[beta_loc]),
                                    solution.tolerance=genoud_control$solution.tolerance,
                                    gradient.check=genoud_control$gradient.check,
                                    control=update_opt_c(genoud_control$optim_control,beta_loc),
                                    BFGSburnin=genoud_control$BFGSburnin,
                                    balance=genoud_control$balance,
                                    wait.generations=genoud_control$wait.generations,hard.generation.limit=genoud_control$hard.generation.limit,
                                    m_hat=m_hat,sigma1_matr=sigma1_matr,k_hat=k_hat,
                                    signk=signk,
                                    p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2)
                  assign(".mrp_p3_par",.mrp_p3_par,envir=.GlobalEnv)
                  eval(expression(parallel::clusterExport(.mrp_p3_cl,varlist=".mrp_p3_par",envir=.GlobalEnv)),envir=.GlobalEnv)
                  exp0<-join_par(expression(rgenoud::genoud(
                    fn=.mrp_p3_f2,nvars=length(.mrp_p3_par$beta_start),
                    starting.values=.mrp_p3_par$beta_start,
                    print.level=.mrp_p3_par$print.level,max.generations=.mrp_p3_par$max.generations,pop.size=.mrp_p3_par$pop.size,
                    Domains=.mrp_p3_par$Domains,cluster=.mrp_p3_cl,
                    solution.tolerance=.mrp_p3_par$solution.tolerance,
                    gradient.check=.mrp_p3_par$gradient.check,
                    control=.mrp_p3_par$control,
                    BFGSburnin=.mrp_p3_par$BFGSburnin,
                    balance=.mrp_p3_par$balance,
                    wait.generations=.mrp_p3_par$wait.generations,hard.generation.limit=.mrp_p3_par$hard.generation.limit,
                    m_matrix=.mrp_p3_par$m_hat,sigma1_matr=.mrp_p3_par$sigma1_matr,k_matrix=.mrp_p3_par$k_hat,
                    sign_k1=.mrp_p3_par$signk$sign_k1,sign_k2=.mrp_p3_par$signk$sign_k2,
                    p1_sp=.mrp_p3_par$p1_sp,p2_sp=.mrp_p3_par$p2_sp,r_sp=.mrp_p3_par$r_sp,model_u2=.mrp_p3_par$model_u2,beta_loc=.mrp_p3_par$beta_loc
                  )),as.expression(add_par_list))
                  eval(parse(text=paste0(".mrp_p3_fit<-tryCatch({.mrp_p3_f1(",exp0,")},error=function(e){e})")),envir =.GlobalEnv)
                  fit1<-get(".mrp_p3_fit",envir=.GlobalEnv)
                }else{
                  fit1<-tryCatch({withWarnings(rgenoud::genoud(
                    fn=est_proc_cont_p3_f0,nvars=length(beta_start),
                    starting.values=beta_start,
                    print.level=genoud_control$print.level,max.generations=genoud_control$max.generations,pop.size=pop.size,
                    Domains=cbind(control_gs$lower[beta_loc],control_gs$upper[beta_loc]),cluster=mycl,
                    solution.tolerance=genoud_control$solution.tolerance,
                    gradient.check=genoud_control$gradient.check,
                    control=update_opt_c(genoud_control$optim_control,beta_loc),
                    BFGSburnin=genoud_control$BFGSburnin,
                    balance=genoud_control$balance,
                    wait.generations=genoud_control$wait.generations,hard.generation.limit=genoud_control$hard.generation.limit,
                    m_matrix=m_hat,sigma1_matr=sigma1_matr,k_matrix=k_hat,
                    sign_k1=signk$sign_k1,sign_k2=signk$sign_k2,
                    p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2,beta_loc=beta_loc,...
                  ))},error=function(e){return(e)})
                }
              }else{
                fit1<-tryCatch({GenSA::GenSA(fn=est_proc_cont_p3_f0,par=beta_start,
                                             lower=control_gs$lower[beta_loc],
                                             upper=control_gs$upper[beta_loc],
                                             control=update_GenSA(control_gs$GenSA_control,maxit.exc),
                                             m_matrix=m_hat,sigma1_matr=sigma1_matr,k_matrix=k_hat,
                                             sign_k1=signk$sign_k1,sign_k2=signk$sign_k2,
                                             p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2,beta_loc=beta_loc
                )},error=function(e){return(e)})
              }
              
              fit_gs<-c(fit_gs,list(fit1))
              
              if((!"error"%in%class(fit1))&(!"error"%in%class(fit))){
                if("mywarning"%in%class(fit1)){
                  warning_out<-fit1$warning[!stringr::str_detect(fit1$warning,"(Out of Boundary individual)|(NaNs produced)")]
                  if(length(warning_out)>0){
                    message("Warnings in genoud:\n",paste0(warning_out,collapse="\n"))
                  }
                  fit1<-fit1$value
                }
                if(fit1$value<fit$minimum){
                  fit<-format_fit1(fit1)
                }
              }else{
                message("Errors detected in fit_nlminb or fit_gs.")
              }
            }
          }
        }
        if("error"%in%class(fit)){
          warning("Errors detected. Already extracted data are returned.")
          return(list(fit=fit,fit_nlminb=fit_nlminb,fit_gs=fit_gs,m_hat=m_hat,m_sigma=sigma1_matr,k_hat=k_hat,k_sigma=sigma2_matr,
                      mk_sigma_list=mk_sigma_list,wald_p=wald_p,u_input=u_input))
        }
        fit$estimate0<-fit$estimate
        fit$estimate<-ext_estimate(fit$estimate,beta_start0,beta_loc)
        
        if(control_p3$model_select==F){break}
        diag_out<-diagnose_fit(fit,p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2,lower=control_p3$check_fit_lower,upper=control_p3$check_fit_upper,myEgger=myEgger0,dt=dt)
        if(identical(diag_out,"pass")){
          if((median(sqrt(sigma1_matr[,1]))*s1_cut_k>sqrt(exp(fit$estimate[5])))&!identical(p1_sp,1)){
            message("Model selection: to Egger or intercept model. (s1 may be too small: ",sqrt(exp(fit$estimate[5])),")")
            diag_out<-list(p1_sp=1,p2_sp=0,r_sp=NULL,model_u2=T,myEgger=myEgger0)
          }else{
            break
          }
        }
        
        if(identical(diag_out,"fail")){fit$code<-"NAs detected by diagnose_fit";break}
        p1_sp<-diag_out$p1_sp
        p2_sp<-diag_out$p2_sp
        r_sp<-diag_out$r_sp
        model_u2<-diag_out$model_u2
        myEgger0<-diag_out$myEgger
        signk<-crt_signk(k_hat,p1_sp,p2_sp,model_u2,myEgger0)
        auto_c<-auto_c+1
      }
      if(mc.cores>1&control_gs$global_search&control_gs$gs_type=="genoud"){
        eval(expression({parallel::stopCluster(.mrp_p3_cl);rm(.mrp_p3_cl,.mrp_p3_f1,.mrp_p3_f2,.mrp_p3_par,.mrp_p3_fit)}),envir=.GlobalEnv)
        gc()
      }
      
      #use nlminb to check the gradients
      if(tryCatch({identical(fit$message$type,"GenSA or genoud")},error=function(e){F})){
        if(dt){cat("Genoud or GenSA finds estimates with a higher likelihood.\r\n")}
        crtb<-crt_beta_loc(fit$estimate,p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2)
        beta_start<-crtb[[1]]
        beta_loc<-crtb[[2]]
        fit0<-tryCatch({nlminb2nlm(nlminb2(f=est_proc_cont_p3_f,p=beta_start,m_matrix=m_hat,sigma1_matr=sigma1_matr,k_matrix=k_hat,
                                           sign_k1=signk$sign_k1,sign_k2=signk$sign_k2,
                                           p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2,beta_loc=beta_loc,
                                           control=nlminb_c_list))},
                       error=function(e){return(e)})
        fit_nlminb<-c(fit_nlminb,list(fit0))
        if("error"%in%class(fit0)){
          fit$code<-"nlminb outputs errors with the solutions given by genoud or GenSA"
        }else{
          fit<-fit0
          fit$estimate0<-fit$estimate
          fit$estimate<-ext_estimate(fit$estimate,beta_start0,beta_loc)
        }
      }
      
      #in case if model_select=F
      fit$code<-check_fit(fit,p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2,lower=control_p3$check_fit_lower,upper=control_p3$check_fit_upper)
      if(fit$code%in%c(0,1)){
        crtb<-crt_beta_loc(fit$estimate,p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2)
        beta_start<-crtb[[1]]
        beta_loc<-crtb[[2]]
        fit0<-tryCatch({suppressWarnings(nlm(f=est_proc_cont_p3_f,p=beta_start,m_matrix=m_hat,sigma1_matr=sigma1_matr,k_matrix=k_hat,
                                             sign_k1=signk$sign_k1,sign_k2=signk$sign_k2,
                                             p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2,beta_loc=beta_loc,
                                             hessian=nlm_c_list$hessian,
                                             fscale=nlm_c_list$fscale,print.level=nlm_c_list$print.level,ndigit=nlm_c_list$ndigit,
                                             gradtol=nlm_c_list$gradtol,stepmax=nlm_c_list$stepmax,steptol=nlm_c_list$steptol[1],
                                             iterlim=nlm_c_list$iterlim,check.analyticals=nlm_c_list$check.analyticals))},
                       error=function(e){return(e)})
        fit_nlminb<-c(fit_nlminb,list(fit0))
        if("error"%in%class(fit0)){
          fit$code<-"nlm outputs errors for the current solution"
        }else{
          if(anyNA(fit0$gradient)){fit0$code<-"NA_gradient"}
          fit0$message<-paste0("nlm code: ",fit0$code)
          if(fit0$code%in%c(1,2,3)){fit0$code<-0}else{fit0$code<-1}
          fit0$estimate0<-fit0$estimate
          fit0$estimate<-ext_estimate(fit0$estimate,beta_start0,beta_loc)
          fit0$code<-check_fit(fit0,p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2,lower=control_p3$check_fit_lower,upper=control_p3$check_fit_upper)
          fit<-fit0
        }
      }
      if(!fit$code%in%c(0)){
        message("Abnormal optimization code; already extracted data are returned.")
        return(list(fit=fit,fit_nlminb=fit_nlminb,fit_gs=fit_gs,m_hat=m_hat,m_sigma=sigma1_matr,k_hat=k_hat,k_sigma=sigma2_matr,
                    mk_sigma_list=mk_sigma_list,wald_p=wald_p,u_input=u_input))
      }
      
      #bootstrap
      vcov_boot<-out_boot<-NULL
      if(boot_se){
        if(dt){cat("Start bootstrapping.\r\n")}
        n1<-nrow(m_hat)
        
        my_task<-function(n_rep){
          t0<-Sys.time()
          pid<-n_rep[[2]]
          n_tar<-n_rep[[3]]
          n_rep<-n_rep[[1]]
          my_proc<-0
          out_boot<-NULL
          
          crtb<-crt_beta_loc(fit$estimate,p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2)
          beta_start<-crtb[[1]]
          beta_loc<-crtb[[2]]
          for(i in 1:n_rep){
            id<-sample(1:n1,n1,replace=T)
            signk<-crt_signk(k_hat[id,],p1_sp,p2_sp,model_u2,myEgger0)
            initial_test<-suppressWarnings(est_proc_cont_p3_f(beta_start,m_matrix=m_hat[id,],sigma1_matr=sigma1_matr[id,],k_matrix=k_hat[id,],
                                                              sign_k1=signk$sign_k1,sign_k2=signk$sign_k2,
                                                              p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2,beta_loc=beta_loc))
            if(initial_test%in%c(NA,NaN,Inf,-Inf)){fit1<-list(estimate=NA,code=NA)}else{
              fit1<-tryCatch({nlminb2nlm(nlminb2(f=est_proc_cont_p3_f,p=beta_start,m_matrix=m_hat[id,],sigma1_matr=sigma1_matr[id,],k_matrix=k_hat[id,],
                                                 sign_k1=signk$sign_k1,sign_k2=signk$sign_k2,
                                                 p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2,beta_loc=beta_loc,
                                                 control=nlminb_c_list))},
                             error=function(e){return(list(estimate=NA,code=NA))})
              fit1$estimate<-ext_estimate(fit1$estimate,beta_start0,beta_loc)
              fit1$code<-check_fit(fit1,p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2,lower=control_p3$check_fit_lower,upper=control_p3$check_fit_upper)
            }
            if(fit1$code%in%c(0)){
              out_boot<-rbind(out_boot,trans_par_norm(fit1$estimate,p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2,est_type="c"))
            }
            if(parallel_trace&mc.cores==1){my_moni("Bootstrap repeat",i,n_tar)}
            if(mc.cores!=1&parallel_trace){
              my_proc<-my_moni2(paste0("Child process ",pid,":"),i,n_tar,my_proc,time=T,t0=t0)
            }
            if(nrow(out_boot)>=n_tar){break}
          }
          return(out_boot)
        }
        
        mylist<-list()
        for(j in 1:mc.cores){
          mylist[[j]]<-list(even_allo(n_rep_max*n_boot,mc.cores)[j],j,even_allo(n_boot,mc.cores)[j])
        }
        add_obj_list<-list(var=c("n1","k_hat","p1_sp","p2_sp","model_u2","myEgger0","r_sp",
                                 "m_hat","sigma1_matr","control_p3","fit","beta_start0"),
                           env=environment())
        exec_base_func<-function(x){
          suppressWarnings(library(MRprollim,quietly=T))
        }
        myfit<-my_parallel(X=mylist,FUN=my_task,mc.cores=mc.cores,PSOCK=PSOCK,dt=dt,
                           print_message=parallel_trace,add_obj_list=add_obj_list,exec_base_func=exec_base_func,export_parent_func=T,seed=round(runif(1,1,.Machine$integer.max)))
        mybind_list_parallel<-function(list){
          out<-NULL
          for(i in 1:length(list)){
            out<-rbind(out,list[[i]])
          }
          return(out)
        }
        out_boot<-mybind_list_parallel(myfit)
        if(nrow(out_boot)<n_boot){warning("The limit n_rep_max*n_boot reached, but MR-PROLLIM failed to obtain n_boot repeats.")}
        vcov_boot<-cov(out_boot,use="na.or.complete")
      }
      
      #asymptotic
      if(dt){cat("Calculate asymptotic variance.\r\n")}
      beta_norm<-trans_par_norm(fit$estimate,p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2,est_type="c")
      
      vcov<-tryCatch({est_variance_ml(f=est_proc_cont_p3_f,beta_opt=beta_norm,n=nrow(m_hat),loc=c(1:7),beta_name=names(beta_norm),name="est_variance_ml",
                                      sandwich=control_p3$sandwich,hessian=control_p3$hessian,
                                      m_matrix=m_hat,sigma1_matr=sigma1_matr,k_matrix=k_hat,
                                      sign_k1=signk$sign_k1,sign_k2=signk$sign_k2,
                                      p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2,
                                      individual=T,vcov_est=T,save_sandw=T)},error=function(e){return(e)})
      
      se_u1<-p_u1<-NA
      if(!"error"%in%class(vcov)){
        se_u1<-tryCatch({sqrt(vcov_sandw["u1","u1"])},error=function(e){NA})
        p_u1<-pnorm(-abs(beta_norm["u1"]/se_u1))*2
        names(p_u1)<-NULL
      }
      
      #output
      if(identical(as.numeric(p1_sp),1)&identical(as.numeric(p2_sp),0)&identical(myEgger0,T)){
        final_Egger_flag<-T
      }else{
        final_Egger_flag<-F
      }
      out<-list(beta_norm=beta_norm,vcov=vcov,vcov_boot=vcov_boot,se_u1_sandwich=se_u1,p_u1_sandwich=p_u1)
      par<-list(p1_sp=p1_sp,p2_sp=p2_sp,r_sp=r_sp,model_u2=model_u2,est_type="c",Egger_info=Egger_info,final_Egger_flag=final_Egger_flag)
      out<-list(estimate=out,parameter=par,maxlik=list(fit_final=fit,fit_gs=fit_gs,fit_nlminb=fit_nlminb),prior_est=NULL,data=list(m_hat=m_hat,m_sigma=sigma1_matr,k_hat=k_hat,k_sigma=sigma2_matr,
                                                                                                                                   mk_sigma_list=mk_sigma_list,wald_p=wald_p,u_input=u_input,boot_data=out_boot))
      class(out)<-"MR-PROLLIM output"
      if(dt){cat("Random-effects MR-PROLLIM (Procedure 3) finished.\r\n")}
      return(out)
    }
  }
}

