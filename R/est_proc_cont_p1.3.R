est_proc_cont_p1.3<-function(x,y,g_nohp,c,c_inherit,start,mc.cores,PSOCK,parallel_trace,dt,twosample_data=NULL,nlminb_control,return_b1_indiv=F){
  ncol2<-function(x){
    if(is.null(x)){return(0)}else{return(ncol(x))}
  }
  crt_mymatrix<-function(matr,loc,n){
    out<-matrix(NA,ncol=ncol(matr),nrow=n)
    out[loc,]<-matr
    return(out)
  }
  
  my_task<-function(my_loc){
    t0<-Sys.time()
    my_proc<-0
    pid<-my_loc[[2]]
    my_loc<-my_loc[[1]]
    mywarn<-paste("Child process",pid,"done.")
    withCallingHandlers(
      {
        eff<-se<-NULL
        w1<-w2<-z<-NULL
        g_nohp<-myselect(g_nohp,my_loc)
        if(!is.null(twosample_data)){g_x<-myselect(twosample_data$g,my_loc)}
        start<-start[my_loc]
        n1<-ncol(g_nohp)
        for(i in 1:n1){
          if(c_inherit){c_m<-c[[1]]}else{
            c_m<-c[[my_loc[i]]]
          }
          loc_m<-which(!is.na(g_nohp[,i]))
          y_m<-y[loc_m]
          if(is.null(twosample_data)){x_m<-x[loc_m]}
          g_m<-g_nohp[loc_m,i]
          c_m<-myselect(c_m,loc_m,type="r")
          
          if(is.null(twosample_data)){
            if(identical(as.numeric(length(unique(c_m[,1]))),1)){
              c_m<-NULL
            }
            g_dum<-crt_dum(g_m)
            d_matr<-cbind(g_dum,c_m,1)
            a1<-as.numeric(solve(t(d_matr)%*%d_matr)%*%t(d_matr)%*%x_m)[1:ncol(g_dum)]
            if(is.null(start)){start_m<-c(0,rep(0,ncol2(c_m)),-1)}else{
              start_m<-start[[i]][1:length(c(0,rep(0,ncol2(c_m)),-1))]
            }
            fit1<-mynlminb(f=f_cont_c_mr,p=start_m,y=y_m,g=g_m,g_dum=g_dum,c=c_m,a1=a1,b4_prior=c(0),name=colnames(g_nohp)[i],control=nlminb_control)
            z_indiv<-sandwich_cont_c_mr(ncol(g_dum)+1,beta_opt=fit1$estimate,y=y_m,g=g_m,g_dum=g_dum,c=c_m,a1=a1,x=x_m,b4_prior=c(0),return_indiv=T,name=colnames(g_nohp)[i])
            vcov<-est_vcov(z_indiv)
            if(anyNA(vcov)){fit1$estimate<-rep(NA,length(fit1$estimate))}
          }else{
            if(c_inherit){c_m2<-twosample_data$c[[1]]}else{
              c_m2<-twosample_data$c[[my_loc[i]]]
            }
            loc_m2<-which(!is.na(g_x[,i]))
            x_m2<-x[loc_m2]
            g_m2<-g_x[loc_m2,i]
            c_m2<-myselect(c_m2,loc_m2,type="r")
            if(identical(as.numeric(length(unique(c_m2[,1]))),1)){
              c_m2<-NULL
            }
            if(identical(as.numeric(length(unique(c_m[,1]))),1)){
              c_m<-NULL
            }
            g_dum2<-crt_dum(g_m2)
            g_dum<-crt_dum(g_m)
            d_matr2<-cbind(g_dum2,c_m2,1)
            b<-as.numeric(solve(t(d_matr2)%*%d_matr2)%*%t(d_matr2)%*%x_m2)
            a1_2<-b[1:ncol(g_dum2)]
            if(is.null(start)){start_m<-c(0,rep(0,ncol2(c_m)),-1)}else{
              start_m<-start[[i]][1:length(c(0,rep(0,ncol2(c_m)),-1))]
            }
            fit1<-mynlminb(f=f_cont_c_mr,p=start_m,y=y_m,g=g_m,g_dum=g_dum,c=c_m,a1=a1_2,b4_prior=c(0),name=colnames(g_nohp)[i],control=nlminb_control)
            fit_indiv<-sandwich_cont_c_mr(ncol(g_dum)+1,beta_opt=fit1$estimate,y=y_m,g=g_m,g_dum=g_dum,c=c_m,a1=a1_2,a1_indiv=matrix(a1_2,ncol=length(a1_2),nrow=length(y_m),byrow=T),
                                          b4_prior=c(0),return_indiv=T,return_w=T,name=colnames(g_nohp)[i])
            
            a1_indiv<-est_indiv(B=-t(d_matr2)%*%d_matr2/length(x_m2),M_indiv=c(x_m2-d_matr2%*%b)*d_matr2,1:length(a1_2),a1_2)
            vcov<-est_vcov(fit_indiv[[1]])+fit_indiv[[2]]%*%est_vcov(a1_indiv)%*%t(fit_indiv[[2]])
            if(anyNA(vcov)){fit1$estimate<-rep(NA,length(fit1$estimate))}
            z_indiv<-fit_indiv[[1]]
          }
          
          se<-cbind(se,sqrt(vcov))
          eff<-cbind(eff,fit1$estimate[1])
          
          if(return_b1_indiv){
            z_m<-matrix(NA,ncol=ncol(z_indiv),nrow=length(y))
            if(!anyNA(z_indiv)){z_m[loc_m,]<-z_indiv}
            z<-cbind(z,z_m[,1])
            
            if(!is.null(twosample_data)){
              w1<-c(w1,list(crt_mymatrix(a1_indiv,loc_m2,length(x))))
              w2<-c(w2,list(fit_indiv[[2]][1,]))
            }
          }
          
          if(parallel_trace&mc.cores==1){my_moni("SNP",i,n1)}
          if(mc.cores!=1&parallel_trace){
            my_proc<-my_moni2(paste0("Child process ",pid,":"),i,n1,my_proc,time=T,t0=t0)
          }
        }
        colnames(eff)<-colnames(se)<-colnames(g_nohp)
        if(return_b1_indiv){
          colnames(z)<-colnames(g_nohp)
        }
      }
      ,warning=function(w){mywarn<<-c(mywarn,w$message);invokeRestart("muffleWarning")}
    )
    if(T){
      if(length(mywarn)>1|parallel_trace){message_parallel(mywarn)}
    }
    return(list(eff,se,z,w1,w2))
  }
  
  mynum<-cut_num(1:ncol(g_nohp),mc.cores)
  exec_base_func<-function(x){
    suppressWarnings(library(MRprollim,quietly=T))
  }
  mycheck<-"pass"
  myfit<-withCallingHandlers({my_parallel(X=mynum,FUN=my_task,mc.cores=mc.cores,PSOCK=PSOCK,dt=dt,
                                          print_message=parallel_trace,export_parent=T,exec_base_func=exec_base_func,seed=NULL)},warning=function(w){mycheck<<-w})
  if((!identical(mycheck,"pass"))&mc.cores!=1){
    warning("An error occurred. Output of my_parallel with errors is returned.")
    message(mycheck)
    class(myfit)<-"myerror"
    return(myfit)
  }
  mybind_list_parallel<-function(list_parallel){
    out1<-out2<-out3<-out4<-out5<-NULL
    for(i in 1:length(list_parallel)){
      out1<-cbind(out1,list_parallel[[i]][[1]])
      out2<-cbind(out2,list_parallel[[i]][[2]])
      out3<-cbind(out3,list_parallel[[i]][[3]])
      out4<-c(out4,list_parallel[[i]][[4]])
      out5<-c(out5,list_parallel[[i]][[5]])
    }
    return(list(out1,out2,out3,out4,out5))
  }
  myfit<-mybind_list_parallel(myfit)
  return(list(eff=myfit[[1]],se=myfit[[2]],b1_indiv=myfit[[3]],w1=myfit[[4]],w2=myfit[[5]]))
}
