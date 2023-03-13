get_sig_snp<-function(x,g,c,c_inherit=T,start=NULL,type=c("c","b"),bi_type=c("log","logit"),
                      p_cut=c(1e-3,1e-5,5e-8),return_dt=T,
                      cd=T,max_unique=10,cd_g_code=T,
                      trinary_only=T,n_min_limit=100,
                      control_limit_c=NULL,scale=T,
                      mc.cores=1,PSOCK=F,parallel_trace=F,dt=T,
                      nlm=T,nlm_control=list(gradtol=1e-8,stepmax=5,steptol=1e-8,iterlim=300),
                      nlminb_control=list(nlminb_control=list(scale=1,eval.max=300,iter.max=300),
                                          nlm_control=list(gradtol=1e-10,stepmax=2,steptol=1e-10))){
  #this function is not designed for fast de novo GWAS
  
  ncol2<-function(x){
    if(is.null(x)){return(0)}else{return(ncol(x))}
  }
  est_indiv_sim<-function(B_solve,M_indiv,loc){
    out<-matrix(NA,nrow=nrow(M_indiv),ncol=length(loc))
    for(i in 1:length(loc)){
      out[,i]<-c(M_indiv%*%(-B_solve[loc[i],]))
    }
    return(out)
  }
  
  type<-match.arg(type)
  bi_type<-match.arg(bi_type)
  
  #check data
  if(cd){
    y<-rep(0,length(x));y[1]<-1
    z<-check_data(x=x,y=y,g=g,c=c,c_inherit=c_inherit,type=type,u_limit=max_unique,twosample_data=NULL,cd_g_code=cd_g_code)
    g<-z[[1]]
    rm(z);gc()
    if(dt){cat("Appropriate input data: true.\r\n")}
  }
  check_start(g,start,"start")
  control_limit_c_org<-list(limit_c=T,dum_loc="auto",quantile=c(0,1),outlier=T)
  control_limit_c<-match.list(control_limit_c,control_limit_c_org)
  check_dum_loc(c,control_limit_c$dum_loc,"control_limit_c$dum_loc")
  update_dum_loc<-function(dum_loc,loc,c_inherit){
    if(is.null(dum_loc)){return(dum_loc)}
    if(identical(dum_loc,"auto")){return(dum_loc)}
    if(c_inherit){
      return(dum_loc)
    }else{
      return(dum_loc[loc])
    }
  }
  
  #trinary g
  if(trinary_only){
    if(dt){"Start trinary g check.\r\n"}
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
      control_limit_c$dum_loc<-update_dum_loc(control_limit_c$dum_loc,which(j),c_inherit)
    }
    stopifnot(ncol(g)!=0)
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
    if(length(control_limit_c$dum_loc)!=length(c)){stop("length(control_limit_c$dum_loc) != length(c).")}
    if(dt){cat("Start modifying the extreme values in c_list.\r\n")}
    for(i in 1:length(c)){
      c[[i]]<-limit_c(c[[i]],dum_loc=control_limit_c$dum_loc[[i]],quantile=control_limit_c$quantile,outlier=control_limit_c$outlier)
    }
  }
  
  #scaling
  if(scale){
    if(dt){cat("Continuous x and continuous c (if appropriate) are automatically scaled.\r\n")}
    if(type=="c"){x<-as.numeric(scale(x))}
    if(!is.null(c)){
      c<-lapply(c,FUN=function(x){
        loc<-apply(x,2,FUN=function(x){length(unique(x))>2})
        x_out<-as.matrix(scale(x))
        x_out[,!loc]<-x[,!loc]
        return(x_out)
      })
    }
  }
  gc()
  
  nlm_c_list_org<-list(hessian=F,fscale=1,print.level=0,ndigit=12,
                       gradtol=1e-8,stepmax=5,steptol=1e-8,iterlim=300,check.analyticals=F)
  nlm_c_list<-match.list(nlm_control,nlm_c_list_org)
  if(dt){cat("Start SNP-Exposure effect calculation.\r\n")}
  my_task<-function(my_loc){
    t0<-Sys.time()
    my_proc<-0
    pid<-my_loc[[2]]
    my_loc<-my_loc[[1]]
    mywarn<-paste("Child process",pid,"done.")
    withCallingHandlers(
      {
        eff_all<-eff<-v<-list()
        p<-rep(NA,length(my_loc))
        
        g<-myselect(g,my_loc)
        start<-start[my_loc]
        n1<-ncol(g)
        gc()
        for(i in 1:n1){
          if(c_inherit){c_m<-c[[1]]}else{
            c_m<-c[[my_loc[i]]]
          }
          loc_m<-which(!is.na(g[,i]))
          x_m<-x[loc_m]
          g_m<-g[loc_m,i]
          c_m<-myselect(c_m,loc_m,type="r")
          if(identical(length(unique(c_m[,1])),1L)){c_m<-NULL}
          
          if(type=="c"){
            g_dum<-crt_dum(g_m)
            d_matr<-cbind(g_dum,c_m,1)
            m<-solve(t(d_matr)%*%d_matr)
            b<-m%*%t(d_matr)%*%x_m
            a1_indiv_c<-est_indiv_sim(B_solve=-m*length(x_m),M_indiv=c(x_m-d_matr%*%b)*d_matr,1:ncol(g_dum))
            vcov<-est_vcov(a1_indiv_c)
            p[i]<-mywald_test(c(b)[1:ncol(g_dum)],vcov)[2]
            if(return_dt){
              eff[[i]]<-c(b)[1:ncol(g_dum)]
              v[[i]]<-vcov
              eff_all[[i]]<-c(b)
            }
          }else{
            if(bi_type=="log"){
              g_dum<-crt_dum(g_m)
              if(is.null(start)){start_m<-c(rep(0,ncol(g_dum)),rep(0,ncol2(c_m)),-1)}else{
                start_m<-start[[i]][1:length(c(rep(0,ncol(g_dum)),rep(0,ncol2(c_m)),-1))]
              }
              
              if(!nlm){
                fit<-mynlminb(f=f_log_linear,p=start_m,x=x_m,g=g_dum,c=c_m,name=colnames(g)[i],control=nlminb_control)
              }else{
                fit<-tryCatch({suppressWarnings(nlm(f=f_log_linear,p=start_m,x=x_m,g=g_dum,c=c_m,
                                                    hessian=nlm_c_list$hessian,
                                                    fscale=nlm_c_list$fscale,print.level=nlm_c_list$print.level,ndigit=nlm_c_list$ndigit,
                                                    gradtol=nlm_c_list$gradtol,stepmax=nlm_c_list$stepmax,steptol=nlm_c_list$steptol[1],
                                                    iterlim=nlm_c_list$iterlim,check.analyticals=nlm_c_list$check.analyticals))},error=function(e){e})
                if("error"%in%class(fit)){
                  warning(paste0(colnames(g)[i]," failed in the effect calculation (nlm gives an error)"))
                  fit<-list(estimate=rep(NA,length(start_m)))
                }else{
                  if(anyNA(fit$gradient)|fit$code%in%c(4,5)){
                    warning(paste0(colnames(g)[i]," failed in the effect calculation (nlm code: ",fit$code,")"))
                    fit$estimate<-rep(NA,length(fit$estimate))
                  }
                }
              }
              vcov<-sandwich_log_linear(1:ncol(g_dum),fit$estimate,x=x_m,g=g_dum,c=c_m,name=colnames(g)[i])
              if(anyNA(vcov)){fit$estimate<-rep(NA,length(fit$estimate))}
              p[i]<-mywald_test(fit$estimate[1:ncol(g_dum)],vcov)[2]
              if(return_dt){
                eff[[i]]<-fit$estimate[1:ncol(g_dum)]
                v[[i]]<-vcov
                eff_all[[i]]<-fit$estimate
              }
            }else{
              fit<-glm(x_m~cbind(g_dum,c_m),family=binomial())
              eff0<-c(coef(fit)[-1],coef(fit)[1])
              vcov<-vcov(fit)[(1:ncol(g_dum))+1,(1:ncol(g_dum))+1]
              names(eff0)<-NULL
              rownames(vcov)<-colnames(vcov)<-NULL
              p[i]<-mywald_test(eff0[1:ncol(g_dum)],vcov)[2]
              if(return_dt){
                eff[[i]]<-eff0[1:ncol(g_dum)]
                v[[i]]<-vcov
                eff_all[[i]]<-eff0
              }
            }
          }
          if(parallel_trace&mc.cores==1){my_moni("SNP",i,n1)}
          if(mc.cores!=1&parallel_trace){
            my_proc<-my_moni2(paste0("Child process ",pid,":"),i,n1,my_proc,time=T,t0=t0)
          }
        }
        names(p)<-colnames(g)
        if(return_dt){names(eff)<-names(v)<-names(eff_all)<-colnames(g)}
      }
      ,warning=function(w){mywarn<<-c(mywarn,w$message);invokeRestart("muffleWarning")}
    )
    if(T){
      if(length(mywarn)>1|parallel_trace){message_parallel(mywarn)}
    }
    return(list(p=p,eff=eff,vcov=v,eff_all=eff_all))
  }
  
  mynum<-cut_num(1:ncol(g),mc.cores)
  mycheck<-"pass"
  exec_base_func<-function(x){
    suppressWarnings(library(MRprollim,quietly=T))
  }
  myfit<-withCallingHandlers({my_parallel(X=mynum,FUN=my_task,mc.cores=mc.cores,PSOCK=PSOCK,dt=dt,
                                          print_message=parallel_trace,exec_base_func=exec_base_func,export_parent=T,seed=NULL)},warning=function(w){mycheck<<-w})
  if((!identical(mycheck,"pass"))&mc.cores!=1){
    warning("An error occurred. Output of my_parallel with errors is returned.")
    message(mycheck)
    class(myfit)<-"myerror"
    return(myfit)
  }
  
  mybind_list_parallel<-function(list_parallel){
    out1<-out2<-out3<-out4<-NULL
    for(i in 1:length(list_parallel)){
      out1<-c(out1,list_parallel[[i]][[1]])
      out2<-c(out2,list_parallel[[i]][[2]])
      out3<-c(out3,list_parallel[[i]][[3]])
      out4<-c(out4,list_parallel[[i]][[4]])
    }
    return(list(out1,out2,out3,out4))
  }
  myfit<-mybind_list_parallel(myfit)
  p<-myfit[[1]]
  sig_snp<-list()
  for(i in 1:length(p_cut)){
    sig_snp[[i]]<-names(p)[which(p<p_cut[i])]
  }
  
  return(list(sig_snp=sig_snp,p=p,eff=myfit[[2]],vcov=myfit[[3]],eff_all=myfit[[4]]))
}
