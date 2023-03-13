est_proc_cont_p1.2<-function(x,y,g_unsuit,c,c_inherit,start,
                             b1_a1_list=NULL,w_a1_list=NULL,
                             b1_sp,b1_indiv_matr,se_indiv,se_m,
                             mc.cores,PSOCK,parallel_trace,dt,twosample_data=NULL,nlminb_control){
  ncol2<-function(x){
    if(is.null(x)){return(0)}else{return(ncol(x))}
  }
  crt_mymatrix<-function(matr,loc,n){
    out<-matrix(NA,ncol=ncol(matr),nrow=n)
    out[loc,]<-matr
    return(out)
  }
  
  stopifnot(is.matrix(g_unsuit))
  
  my_task<-function(my_loc){
    t0<-Sys.time()
    my_proc<-0
    pid<-my_loc[[2]]
    my_loc<-my_loc[[1]]
    mywarn<-paste("Child process",pid,"done.")
    
    withCallingHandlers(
      {
        wald<-NULL
        out<-list()
        se<-eff<-NULL
        if(!is.null(twosample_data)){g_x<-myselect(twosample_data$g,my_loc)}
        g_unsuit<-myselect(g_unsuit,my_loc)
        start<-start[my_loc]
        n1<-ncol(g_unsuit)
        for(i in 1:n1){
          if(c_inherit){c_m<-c[[1]]}else{
            c_m<-c[[my_loc[i]]]
          }
          loc_m<-which(!is.na(g_unsuit[,i]))
          y_m<-y[loc_m]
          if(is.null(twosample_data)){x_m<-x[loc_m]}
          g_m<-g_unsuit[loc_m,i]
          c_m<-myselect(c_m,loc_m,type="r")
          
          if(is.null(twosample_data)){
            if(identical(as.numeric(length(unique(c_m[,1]))),1)){
              c_m<-NULL
            }
            g_dum<-crt_dum(g_m)
            d_matr<-cbind(g_dum,c_m,1)
            a1<-as.numeric(solve(t(d_matr)%*%d_matr)%*%t(d_matr)%*%x_m)[1:ncol(g_dum)]
            if(is.null(start)){start_m<-c(rep(0,ncol(g_dum)),rep(0,ncol2(c_m)),-1)}else{
              start_m<-start[[i]][1:length( c(rep(0,ncol(g_dum)),rep(0,ncol2(c_m)),-1) )]
            }
            
            fit1<-mynlminb(f=f_cont_c_mr,p=start_m,y=y_m,g=g_dum,g_dum=g_dum,c=c_m,a1=a1,b4_prior=rep(1,ncol(g_dum)),b1_sp=b1_sp,name=colnames(g_unsuit)[i],control=nlminb_control)
            vcov<-sandwich_cont_c_mr((ncol(g_dum)+2):(2*ncol(g_dum)+1),beta_opt=fit1$estimate,
                                     y=y_m,g=g_dum,g_dum=g_dum,c=c_m,a1=a1,x=x_m,
                                     b4_prior=rep(1,ncol(g_dum)),
                                     b1_sp=b1_sp,b1_sp_indiv=rep(b1_sp,length(y_m)),
                                     loc_indiv=loc_m,se=se_indiv,se_m=se_m,b1_indiv_matr=b1_indiv_matr,
                                     name=colnames(g_unsuit)[i])
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
            if(is.null(start)){start_m<-c(rep(0,ncol(g_dum)),rep(0,ncol2(c_m)),-1)}else{
              start_m<-start[[i]][1:length( c(rep(0,ncol(g_dum)),rep(0,ncol2(c_m)),-1) )]
            }
            fit1<-mynlminb(f=f_cont_c_mr,p=start_m,y=y_m,g=g_dum,g_dum=g_dum,c=c_m,a1=a1_2,b4_prior=rep(1,ncol(g_dum)),b1_sp=b1_sp,name=colnames(g_unsuit)[i],control=nlminb_control)
            
            a1_indiv<-est_indiv(B=-t(d_matr2)%*%d_matr2/length(x_m2),M_indiv=c(x_m2-d_matr2%*%b)*d_matr2,1:length(a1_2),a1_2)
            
            vcov<-sandwich_cont_c_mr((ncol(g_dum)+2):(2*ncol(g_dum)+1),beta_opt=fit1$estimate,
                                     y=y_m,g=g_dum,g_dum=g_dum,c=c_m,a1=a1_2,a1_indiv=matrix(a1_2,ncol=length(a1_2),nrow=length(y_m),byrow=T),
                                     b4_prior=rep(1,ncol(g_dum)),
                                     b1_sp=b1_sp,b1_sp_indiv=rep(b1_sp,length(y_m)),
                                     a1_matrix=crt_mymatrix(a1_indiv,loc_m2,length(x)),b1_a1_list=b1_a1_list,w_a1_list=w_a1_list,
                                     loc_indiv=loc_m,se=se_indiv,se_m=se_m,b1_indiv_matr=b1_indiv_matr,
                                     name=colnames(g_unsuit)[i])
            if(anyNA(vcov)){fit1$estimate<-rep(NA,length(fit1$estimate))}
          }
          
          eff<-fit1$estimate[1:ncol(g_dum)]
          wald<-cbind(wald,mywald_test(eff,vcov))
          out[[i]]<-list(eff=eff,vcov=vcov)
          if(parallel_trace&mc.cores==1){my_moni("SNP",i,n1)}
          if(mc.cores!=1&parallel_trace){
            my_proc<-my_moni2(paste0("Child process ",pid,":"),i,n1,my_proc,time=T,t0=t0)
          }
        }
        names(out)<-colnames(g_unsuit)
        colnames(wald)<-colnames(g_unsuit)
        rownames(wald)<-c("chisq","p")
      }
      ,warning=function(w){mywarn<<-c(mywarn,w$message);invokeRestart("muffleWarning")}
    )
    if(T){
      if(length(mywarn)>1|parallel_trace){message_parallel(mywarn)}
    }
    return(list(wald,out))
  }
  
  mynum<-cut_num(1:ncol(g_unsuit),mc.cores)
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
    out1<-out2<-NULL
    for(i in 1:length(list_parallel)){
      out1<-cbind(out1,list_parallel[[i]][[1]])
      out2<-c(out2,list_parallel[[i]][[2]])
    }
    return(list(out1,out2))
  }
  myfit<-mybind_list_parallel(myfit)
  wald<-myfit[[1]]
  out<-myfit[[2]]
  return(list(wald_test=wald,hp_vcov=out))
}
