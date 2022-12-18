select_suit<-function(x,g,c,c_inherit,twosample_data=NULL,mc.cores,PSOCK,dt,parallel_trace,p_cut=0.1,adj_m="bonferroni"){
  #c is NULL or list(matrix)
  select_suit2<-function(x,g,c){
    #c is NULL or matrix
    crt_dum<-function(x,ref=0){
      x_u=sort(unique(x),decreasing=F)
      x_u=x_u[-which(x_u==ref)]
      out=NULL
      for(i in 1:(length(x_u))){
        out=cbind(out,as.numeric(x==(x_u[i])))
      }
      colnames(out)=x_u
      return(out)
    }
    get_p<-function(eff,vcov){
      eff_2<-eff[1]*2-eff[2]
      se_2<-sqrt(vcov[1,1]*4+vcov[2,2]-4*vcov[1,2])
      return(pnorm(-abs(eff_2/se_2))*2)
    }
    if(identical(length(unique(c[,1])),1L)){
      c<-NULL
    }
    g_dum<-crt_dum(g)
    fit<-lm(x~cbind(g_dum,c))
    p<-get_p(coef(fit)[2:(ncol(g_dum)+1)],vcov(fit)[2:(ncol(g_dum)+1),2:(ncol(g_dum)+1)])
    return(p)
  }

  if(!is.null(twosample_data)){
    g<-twosample_data$g
    c<-twosample_data$c
  }
  
  my_task<-function(my_loc){
    pid<-my_loc[[2]]
    my_loc<-my_loc[[1]]
    mywarn<-paste("Child process",pid,"done.")
    withCallingHandlers(
      {
        g<-myselect(g,my_loc)
        n1<-ncol(g)
        out<-c()
        for(i in 1:n1){
          if(c_inherit){c_m<-c[[1]]}else{
            c_m<-c[[my_loc[i]]]
          }
          loc_m<-which(!is.na(g[,i]))
          x_m<-x[loc_m]
          g_m<-g[loc_m,i]
          c_m<-myselect(c_m,loc_m,type="r")
          out[i]<-select_suit2(x=x_m,g=g_m,c=c_m)
        }
      }
      ,warning=function(w){mywarn<<-c(mywarn,w$message);invokeRestart("muffleWarning")}
    )
    if(T){
      if(length(mywarn)>1){message_parallel(mywarn)}
    }
    return(out)
  }
  mynum<-cut_num(1:ncol(g),mc.cores)
  add_obj_list<-list(var=c("g","c","x","c_inherit"),
                     env=environment())
  exec_base_func<-function(x){
    Sys.sleep(x/10)
    library(MRprollim,quietly=T)
  }
  mycheck<-"pass"
  myfit<-withCallingHandlers({my_parallel(X=mynum,FUN=my_task,mc.cores=mc.cores,PSOCK=PSOCK,dt=dt,
                                          print_message=parallel_trace,export_parent_func=T,add_obj_list=add_obj_list,exec_base_func=exec_base_func)},warning=function(w){mycheck<<-w})
  if((!identical(mycheck,"pass"))&mc.cores!=1){
    warning("An error occurred. Output of my_parallel with errors is returned.")
    message(mycheck)
    class(myfit)<-"myerror"
    return(myfit)
  }
  p<-unlist(myfit)
  j<-p.adjust(p,adj_m)<=p_cut
  if(anyNA(j)){message("NAs detected in select_suit.");j[is.na(j)]<-F}
  return(j)
}
