est_proc_bi_p1_dl<-function(x,y,g_suit,c,c_inherit,dum_loc_list,start,mc.cores,PSOCK,parallel_trace,dt,
                            k1_td_p=0.1,k1_td_adj_m="bonferroni",delta_p_cut=0.01,
                            p2_cut=0.1,p2_cut_adj_m="none",p2_cut2=0.05,p2_cut2_adj_m="none",p3_cut=0.05,p3_cut_adj_m="bonferroni",boot_n=10000,
                            select_list=NULL,
                            d=.Machine$double.eps^(1/3),n_random=NULL,n_max=2^16,twosnp_p_cut=1e-5,
                            p_limit,est_type,data_k_all_list=NULL,crt_data3_list=NULL,nlminb_control=list()){
  stopifnot(p3_cut_adj_m%in%c("bonferroni","bonf","none"))
  #functions
  calcu_b1_b4<-function(est_mkp,vcov_mkp,est_indiv,boot_n=10000,delta_p_cut,length_all,nonNA_loc,d){
    if(anyNA(est_mkp)|anyNA(vcov_mkp)){
      return(
        list(fb11=c(est=NA,se=NA),
             fb12=c(est=NA,se=NA),
             expb41=c(est=NA,se=NA),
             expb42=c(est=NA,se=NA),
             fb101_indiv=rep(NA,length_all),
             fb102_indiv=rep(NA,length_all),
             delta_p=NA,j=F)
      )
    }
    
    f<-function(est_mkp){
      m1<-est_mkp[1]
      m2<-est_mkp[2]
      k1<-est_mkp[3]
      k2<-est_mkp[4]
      p_td<-est_mkp[5]
      
      k1_td<-(exp(k1)-1)*p_td
      k2_td<-(exp(k2)-1)*p_td
      m_td<-exp(2*m1-m2)
      
      delta<-m_td^2*k2_td^2-4*m_td*k1_td*k2_td+4*m_td*k1_td^2
      if(delta<0){
        return(list(s1=c(NA,NA),s2=c(NA,NA)))
      }
      
      fb11<-(m_td*k2_td-2*k1_td-sqrt(delta))/(2*k1_td^2)
      fb12<-(m_td*k2_td-2*k1_td+sqrt(delta))/(2*k1_td^2)
      expb41<-exp(m1)/(k1_td*fb11+1)
      expb42<-exp(m1)/(k1_td*fb12+1)
      return(
        list(s1=c(fb11,expb41),s2=c(fb12,expb42))
      )
    }
    s0<-f(est_mkp)
    fb101<-s0$s1[1]
    fb102<-s0$s2[1]
    expb401<-s0$s1[2]
    expb402<-s0$s2[2]
    
    if(anyNA(fb101)){
      return(
        list(fb11=c(est=NA,se=NA),
             fb12=c(est=NA,se=NA),
             expb41=c(est=NA,se=NA),
             expb42=c(est=NA,se=NA),
             fb101_indiv=rep(NA,length_all),
             fb102_indiv=rep(NA,length_all),
             delta_p=NA,j=F)
      )
    }
    
    boot_data<-mvrnorm(boot_n,est_mkp,vcov_mkp)
    
    m1<-boot_data[,1]
    m2<-boot_data[,2]
    k1<-boot_data[,3]
    k2<-boot_data[,4]
    p_td<-boot_data[,5]
    
    k1_td<-(exp(k1)-1)*p_td
    k2_td<-(exp(k2)-1)*p_td
    m_td<-exp(2*m1-m2)
    
    #k1_td0<-(exp(est_mkp[3])-1)*est_mkp[5]
    #k1_td_p<-pnorm(-abs(k1_td0)/sd(k1_td))*2
    
    delta<-m_td^2*k2_td^2-4*m_td*k1_td*k2_td+4*m_td*k1_td^2
    delta_p<-mean(delta<0)
    
    fb101_indiv<-fb102_indiv<-rep(NA,length_all)
    v1<-v2<-v3<-v4<-NA
    if(delta_p<delta_p_cut){
      j<-T
      der1<-der2<-der3<-der4<-rep(NA,5)
      for(i in 1:5){
        est_mkp2<-est_mkp1<-est_mkp
        est_mkp1[i]<-est_mkp[i]-d
        est_mkp2[i]<-est_mkp[i]+d
        s1<-f(est_mkp1)
        s2<-f(est_mkp2)
        der1[i]<-(s2$s1[1]-s1$s1[1])/(2*d)
        der2[i]<-(s2$s2[1]-s1$s2[1])/(2*d)
        der3[i]<-(s2$s1[2]-s1$s1[2])/(2*d)
        der4[i]<-(s2$s2[2]-s1$s2[2])/(2*d)
      }
      est_indiv<-scale(est_indiv,scale=F,center=T)
      n<-nrow(est_indiv)
      v_ind1<-est_indiv%*%der1+fb101
      v_ind2<-est_indiv%*%der2+fb102
      v1<-var(v_ind1)/n
      v2<-var(v_ind2)/n
      fb101_indiv[nonNA_loc]<-v_ind1
      fb102_indiv[nonNA_loc]<-v_ind2
      v3<-var(est_indiv%*%der3)/n
      v4<-var(est_indiv%*%der4)/n
      if(anyNA(c(v1,v2,v3,v4))){j<-F}
    }else{
      j<-F
    }
    
    #fb11<-suppressWarnings((m_td*k2_td-2*k1_td-sqrt(delta))/(2*k1_td^2))
    #fb12<-suppressWarnings((m_td*k2_td-2*k1_td+sqrt(delta))/(2*k1_td^2))
    #expb41<-suppressWarnings(exp(m1)/(k1_td*fb11+1))
    #expb42<-suppressWarnings(exp(m1)/(k1_td*fb12+1))

    return(
      list(fb11=c(est=fb101,se=sqrt(v1)),
           fb12=c(est=fb102,se=sqrt(v2)),
           expb41=c(est=expb401,se=sqrt(v3)),
           expb42=c(est=expb402,se=sqrt(v4)),
           fb101_indiv=fb101_indiv,
           fb102_indiv=fb102_indiv,
           delta_p=delta_p,j=j)
    )
  }
  get_lowest_q<-function(est_list,se_list,est_list2,se_list2,fb1_indiv_list,n_random=100,n_max=2^16,p_cut=1e-5){
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
    cochran_q<-function(eff,se){
      loc<-which(!is.na(eff))
      eff<-eff[loc]
      se<-se[loc]
      w<-1/se^2/sum(1/se^2)
      eff_m<-sum(w*eff)
      q=sum((eff-eff_m)^2/se^2)
      return(q)
    }
    get_choice_matr<-function(list,n_max=2^16){
      for(i in 1:length(list)){
        if(i==1){out<-cbind(list[[i]]);next}
        out<-reshape2(out,length(list[[i]]))
        out<-cbind(out,list[[i]])
        if(nrow(out)>n_max){stop("Too many combinations.")}
      }
      return(out)
    }
    
    n1<-length(est_list)
    n2<-sum(unlist(lapply(est_list,FUN=function(x){length(x)==2})))
    m<-ceiling(n2/log(n_max,2))
    if(m==0){m<-1}
    if(is.null(n_random)){
      m0<-floor(n1/m)
      if(m0==0){m0<-1}
      n_random<-ceiling(log(p_cut,(n1-m0)/(n1-1)))
      if(n_random==0){n_random<-1}
    }
    loc_list<-list()
    for(i in 1:n1){
      loc_list[[i]]<-1:length(est_list[[i]])
    }
    
    des_random<-loc_random<-list()
    for(i in 1:n_random){
      des<-cut_num(sample(1:n1,n1),m)
      loc_matr_m<-list()
      for(j in 1:m){
        est_matr<-get_choice_matr(est_list[des[[j]]],n_max)
        se_matr<-get_choice_matr(se_list[des[[j]]],n_max)
        loc_matr<-get_choice_matr(loc_list[des[[j]]],n_max)
        q<-rep(NA,nrow(est_matr))
        for(k in 1:nrow(est_matr)){
          q[k]<-cochran_q(est_matr[k,],se_matr[k,])
        }
        loc_matr_m[[j]]<-loc_matr[which.min(q),]
      }
      des_random[[i]]<-unlist(des)
      loc_random[[i]]<-unlist(loc_matr_m)
    }
    
    q<-rep(NA,n_random)
    for(i in 1:n_random){
      o1<-est_list[des_random[[i]]]
      o2<-loc_random[[i]]
      o3<-se_list[des_random[[i]]]
      eff<-se<-rep(NA,n1)
      for(j in 1:n1){
        eff[j]<-o1[[j]][o2[j]]
        se[j]<-o3[[j]][o2[j]]
      }
      q[i]<-cochran_q(eff,se)
    }
    final<-which.min(q)
    
    o1<-est_list[des_random[[final]]]
    o12<-est_list2[des_random[[final]]]
    o2<-loc_random[[final]]
    o3<-se_list[des_random[[final]]]
    o32<-se_list2[des_random[[final]]]
    o4<-fb1_indiv_list[des_random[[final]]]
    name<-names(est_list)[des_random[[final]]]
    eff<-se<-eff2<-se2<-rep(NA,n1)
    fb1_indiv<-NULL
    for(j in 1:n1){
      eff[j]<-o1[[j]][o2[j]]
      se[j]<-o3[[j]][o2[j]]
      eff2[j]<-o12[[j]][o2[j]]
      se2[j]<-o32[[j]][o2[j]]
      fb1_indiv<-cbind(fb1_indiv,o4[[j]][,o2[j]])
    }
    ivw<-meta_eff(eff,se)$est
    loc1<-des_random[[final]]
    loc2<-o2
    
    names(eff)<-names(se)<-names(eff2)<-names(se2)<-colnames(fb1_indiv)<-name
    return(list(ivw=ivw,q=q[final],eff=eff,se=se,eff2=eff2,se2=se2,fb1_indiv=fb1_indiv,loc1=loc1,loc2=loc2))
  }
  cover_p<-function(est,se,lower=0,upper=1,p=1e-5){
    if(anyNA(est)|anyNA(se)){return(F)}
    u<-est+se*abs(qnorm(p))
    l<-est-se*abs(qnorm(p))
    return((lower<l)&(upper>u))
  }
  crt_data0<-function(y,x,g,c,dum_loc,start=NULL,p_limit=1e-5,name,k_vcov_r=T,data_k_all=NULL,nlminb_control=list()){
    loc_m<-which(!is.na(g))
    y<-y[loc_m]
    x<-x[loc_m]
    g<-g[loc_m]
    
    j<-F
    if(is.null(c)){
      o1<-x[g==1]
      o2<-x[g==2]
      o3<-x[g==0]
      o11<-y[g==1]
      o22<-y[g==2]
      o33<-y[g==0]
      o<-x[g==0]
      o[o33==0]<-NA
      
      j<-unlist(lapply(list(o1,o2,o3,o,o11,o22,o33),FUN=function(x){length(unique2(x))%in%c(1,0)}))
      
      if(!any(j)){
        out1<-log(mean(o1))-log(mean(o3))#k1
        out2<-log(mean(o2))-log(mean(o3))#k2
        out3<-mean(o,na.rm=T)#p
        out4<-log(mean(o11))-log(mean(o33))#m1
        out5<-log(mean(o22))-log(mean(o33))#m2
        
        m1<-1/(mean(o3)^2)*var(o3)/length(o3)
        m2<-(-1/mean(o3)*cov(o3,o,use="complete.obs")/length(o3))
        vcov<-rbind(
          c( 1/(mean(o1)^2)*var(o1)/length(o1)+m1, m1, m2),
          c( m1, 1/(mean(o2)^2)*var(o2)/length(o2)+m1, m2),
          c( m2, m2, var(o,na.rm=T)/sum(!is.na(o)))
        )
        
        k_hat<-rbind(c(out1,out2,out3))
        m_hat<-rbind(c(out4,out5))
        m3<-1/(mean(o33)^2)*var(o33)/length(o33)
        m_sigma<-rbind(c(1/(mean(o11)^2)*var(o11)/length(o11)+m3,m3,1/(mean(o22)^2)*var(o22)/length(o22)+m3))
        
        m4<-1/mean(o3)/mean(o33)*cov(o3,o33)/length(o3)
        m5<-(-1/mean(o33)*cov(o,o33,use="na.or.complete")/length(o33))
        vcov_mkp<-rbind(
          c(1/mean(o1)/mean(o11)*cov(o1,o11)/length(o1)+m4, m4 , m5),
          c(m4, 1/mean(o2)/mean(o22)*cov(o2,o22)/length(o2)+m4 , m5)
        )
        
        mkp_vcov<-matrix(NA,5,5)
        mkp_vcov[1,1:2]<-m_sigma[1,1:2]
        mkp_vcov[2,1:2]<-m_sigma[1,2:3]
        mkp_vcov[3:5,3:5]<-vcov
        mkp_vcov[1:2,3:5]<-vcov_mkp
        mkp_vcov[3:5,1:2]<-t(vcov_mkp)
        
        f<-function(x,loc,type){
          p0<-mean(loc)
          x0<-mean(x[loc])
          XI<-as.numeric(loc)
          XE<-x*XI
          if(type=="log"){
            return(XE/p0/x0-1/p0*XI+log(x0))
          }
          if(type=="x"){
            return(1/p0*XE-x0/p0*XI+x0)
          }
          return(NA)
        }
        
        mkp_ind<-cbind(f(y,g==1,"log")-f(y,g==0,"log"),
                       f(y,g==2,"log")-f(y,g==0,"log"),
                       f(x,g==1,"log")-f(x,g==0,"log"),
                       f(x,g==2,"log")-f(x,g==0,"log"),
                       f(x,g==0&y==1,"x"))
        
        j<-anyNA(m_hat)|anyNA(m_sigma)|anyNA(k_hat)|anyNA(vcov)|anyNA(mkp_vcov)|anyNA(mkp_ind)|anyNA(loc_m)|(!cover_p(out3,sqrt(diag(vcov)[3]),p=p_limit))
      }
    }
    if(!is.null(c)){
      c<-myselect(c,loc_m,"r")
      g_dum<-crt_dum(g)
      
      if(is.null(data_k_all)){
        if(is.null(start)){start_k<-c(rep(0,ncol(g_dum)),rep(0,ncol(c)),-1)}else{
          start_k<-start[1:length(c(rep(0,ncol(g_dum)),rep(0,ncol(c)),-1))]
        }
        fit<-mynlminb(f=f_log_linear,p=start_k,x=x,g=g_dum,c=c,name=name,control=nlminb_control)
      }else{
        fit<-list(estimate=data_k_all)
      }
      
      k_indiv<-sandwich_log_linear(1:2,fit$estimate,x=x,g=g_dum,c=c,return_indiv=T,name=name)
      if(anyNA(k_indiv)|is.null(k_indiv)){
        m_hat<-rbind(c(NA,NA))
        m_sigma<-rbind(c(NA,NA,NA))
        k_hat<-rbind(rep(NA,3))
        vcov<-matrix(NA,ncol=3,nrow=3)
        mkp_vcov<-matrix(NA,ncol=5,nrow=5)
        mkp_ind<-NULL
        return(list(m_hat=m_hat,m_sigma=m_sigma,k_hat=k_hat,k_sigma=vcov,mkp_sigma=mkp_vcov,mkp_ind=mkp_ind,nonNA_loc=loc_m))
      }
      
      x0<-rep(NA,length(x))
      loc<-which(g==0&y==1)
      x0[loc]<-x[loc]
      out0<-mean(x0,na.rm=T)
      out0_ind<-x*y*as.numeric(g==0)*length(g)/length(loc)-as.numeric(g==0&y==1)*out0/length(loc)*length(g)+out0
      
      vcov<-matrix(NA,3,3)
      vcov[1:2,1:2]<-est_vcov(k_indiv)
      cov1<-cov(k_indiv[,1],x0,use="na.or.complete")/length(x)
      cov2<-cov(k_indiv[,2],x0,use="na.or.complete")/length(x)
      var1<-var(x0,na.rm=T)/length(loc)
      vcov[,3]<-c(cov1,cov2,var1)
      vcov[3,1:2]<-c(cov1,cov2)
      
      loc1<-which(g==1)
      y1<-y[loc1];n1<-length(y1)
      loc2<-which(g==2)
      y2<-y[loc2];n2<-length(y2)
      loc3<-which(g==0)
      o3<-myselect(c,loc3,"r")
      y0<-y[loc3];n3<-length(y0)
      out3<-mean(y0)
      n<-length(y)
      
      o11<-my_allocate_bi(myselect(c,loc1,"r"),o3,dum_loc=dum_loc,cpp=cpp)
      out1<-mean(o11*y1)
      out1_indiv<-(y1-out1)*o11+out1
      
      o22<-my_allocate_bi(myselect(c,loc2,"r"),o3,dum_loc=dum_loc,cpp=cpp)
      out2<-mean(o22*y2)
      out2_indiv<-(y2-out2)*o22+out2
      
      out1_indiv2<-rep(0,length(y))
      out1_indiv2[g==1]<-out1_indiv
      m1_indiv<-log(out1/out3)+g_dum[,1]*(out1_indiv2/out1-1)*n/n1-as.numeric(g==0)*(y/out3-1)*n/n3
      
      out2_indiv2<-rep(0,length(y))
      out2_indiv2[g==2]<-out2_indiv
      m2_indiv<-log(out2/out3)+g_dum[,2]*(out2_indiv2/out2-1)*n/n2-as.numeric(g==0)*(y/out3-1)*n/n3
      
      mkp_ind<-cbind(m1_indiv,m2_indiv,k_indiv[,1],k_indiv[,2],out0_ind)
      mkp_vcov<-est_vcov(mkp_ind)
      if(k_vcov_r){mkp_vcov[3:5,3:5]<-vcov}
      
      m_hat<-(c(log(out1/out3),log(out2/out3)))
      m_hat[m_hat%in%c(Inf,-Inf)]<-NA
      m_hat<-rbind(m_hat)
      m_sigma<-rbind(c(mkp_vcov[1,1],mkp_vcov[1,2],mkp_vcov[2,2]))
      k_hat<-rbind(c(fit$estimate[1:2],out0))
      
      j<-anyNA(m_hat)|anyNA(m_sigma)|anyNA(k_hat)|anyNA(vcov)|anyNA(mkp_vcov)|anyNA(mkp_ind)|anyNA(loc_m)|(!cover_p(out0,sqrt(var1),p=p_limit))
    }
    
    
    if(sum(j)!=0){
      m_hat<-rbind(c(NA,NA))
      m_sigma<-rbind(c(NA,NA,NA))
      k_hat<-rbind(rep(NA,3))
      vcov<-matrix(NA,ncol=3,nrow=3)
      mkp_vcov<-matrix(NA,ncol=5,nrow=5)
      mkp_ind<-NULL
    }
    return(list(m_hat=m_hat,m_sigma=m_sigma,k_hat=k_hat,k_sigma=vcov,mkp_sigma=mkp_vcov,mkp_ind=mkp_ind,nonNA_loc=loc_m))
  }
  my_allocate_bi<-function(c_matr1,c_matr2,dum_loc=NULL,...){
    if(!is.null(dum_loc)){
      cont_loc<-(1:ncol(c_matr1))[-dum_loc]
      if(length(cont_loc)>0){
        c_matr1_cont<-myselect(c_matr1,cont_loc)
        c_matr2_cont<-myselect(c_matr2,cont_loc)
      }
      gc1<-mypaste(cbind(c_matr1[,dum_loc]))
      gc2<-mypaste(cbind(c_matr2[,dum_loc]))
      gc1_u<-unique(gc1)
      if(length(gc1_u)>=100){
        warning("Too many different combinations of categorical C. Please check dum_loc_list.")
      }
      k_adj<-length(gc2)/sum(gc2%in%gc1_u)
      allo_out<-rep(0,nrow(c_matr1))
      for(i in 1:length(gc1_u)){
        loc1<-which(gc1==gc1_u[i])
        loc2<-which(gc2==gc1_u[i])
        if(length(cont_loc)>0){
          allo_out[loc1]<-my_allocate_cpp2(t(myselect(c_matr1_cont,loc1,"r")),myselect(c_matr2_cont,loc2,"r"))
        }else{
          allo_out[loc1]<-length(loc2)/length(loc1)
        }
      }
      allo_out<-allo_out*k_adj
    }else{
      c_matr1_t<-t(c_matr1)
      nr<-nrow(c_matr1_t)
      nc<-ncol(c_matr1_t)
      if(T){
        allo_out<-my_allocate_cpp2(c_matr1_t,c_matr2)
      }else{
        for(i in 1:nrow(c_matr2)){
          d<-.colSums((c_matr1_t-c_matr2[i,])^2,nr,nc)
          loc<-which(d==min(d))
          allo_out[loc]<-allo_out[loc]+1/length(loc)
        }
      }
    }
    return(allo_out/nrow(c_matr2)*nrow(c_matr1))
  }
  check_select_list<-function(fit_b1_b4,select_list){
    if(is.null(select_list)){return("pass")}
    if(!is.list(select_list)){return("select_list should be a named list or NULL.")}
    name<-names(select_list)
    if(is.null(name)){return("select_list should be a named list or NULL.")}
    if(length(unique(name))!=length(select_list)){return("names for select_list should be unique.")}
    if(any(!name%in%names(fit_b1_b4))){return("names of select_list not detected in names(dl_data2)")}
    if( any( !unlist(lapply(select_list,FUN=function(x){is.vector(x)&length(x)>=1})) ) ){return("select_list should be a list containing vectors c(1), c(2), or c(1,2)")}
    "pass"
  }
  
  stopifnot(is.matrix(g_suit))
  if(!is.null(crt_data3_list)){
    if(!identical(names(crt_data3_list),colnames(g_suit))){
      loc<-which(colnames(g_suit)%in%names(crt_data3_list))
      g_suit<-myselect(g_suit,loc)
      stopifnot(identical(names(crt_data3_list),colnames(g_suit)))
      data_k_all_list<-data_k_all_list[loc]
      start<-start[loc]
      if(!c_inherit){
        c<-c[loc]
        dum_loc_list<-dum_loc_list[loc]
      }
      stopifnot(ncol(g_suit)!=0)
      gc()
    }
  }
  
  if(is.null(crt_data3_list)){
    my_task<-function(my_loc){
      my_proc<-0
      pid<-my_loc[[2]]
      my_loc<-my_loc[[1]]
      mywarn<-paste("Child process",pid,"done.")
      t0<-Sys.time()
      withCallingHandlers(
        {
          out<-list()
          g_suit<-myselect(g_suit,my_loc)
          start<-start[my_loc]
          data_k_all_list<-data_k_all_list[my_loc]
          n1<-ncol(g_suit)
          for(i in 1:n1){
            if(c_inherit){c_m<-c[[1]];dum_loc_m<-dum_loc_list[[1]]}else{
              c_m<-c[[my_loc[i]]]
              dum_loc_m<-dum_loc_list[[my_loc[i]]]
            }
            
            fit<-crt_data0(y=y,x=x,g=g_suit[,i],c=c_m,dum_loc=dum_loc_m,start=start[[i]],
                           p_limit=p_limit,name=colnames(g_suit)[i],k_vcov_r=T,data_k_all=data_k_all_list[[i]],
                           nlminb_control=nlminb_control)
            out<-c(out,list(fit))
            if(parallel_trace&mc.cores==1){my_moni("SNP",i,n1)}
            if(mc.cores!=1&parallel_trace){
              my_proc<-my_moni2(paste0("Child process ",pid,":"),i,n1,my_proc,time=T,t0=t0)
            }
          }
          names(out)<-colnames(g_suit)
        }
        ,warning=function(w){mywarn<<-c(mywarn,w$message);invokeRestart("muffleWarning")}
      )
      if(T){
        if(length(mywarn)>1|parallel_trace){message_parallel(mywarn)}
      }
      return(out)
    }
    mynum<-cut_num(1:ncol(g_suit),mc.cores)
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
      out<-list()
      for(i in 1:length(list_parallel)){
        out<-c(out,list_parallel[[i]])
      }
      return(out)
    }
    myfit<-mybind_list_parallel(myfit)
  }else{
    myfit<-crt_data3_list
  }
  
  calcu_k1_td_p<-function(myfit){
    p<-NA
    for(i in 1:length(myfit)){
      eff<-myfit[[i]]$k_hat[1,]
      vcov<-myfit[[i]]$k_sigma
      k1_td<-(exp(eff[1])-1)*eff[3]
      se<-sqrt(vcov[1,1]*(eff[3]*exp(eff[1]))^2+(exp(eff[1])-1)^2*vcov[3,3]+2*vcov[1,3]*(exp(eff[1])-1)*eff[3]*exp(eff[1]))
      p[i]<-pnorm(-abs(k1_td/se))*2
    }
    return(p)
  }
  loc<-p.adjust(calcu_k1_td_p(myfit),k1_td_adj_m)<k1_td_p
  loc[is.na(loc)]<-F
  if(sum(loc)<=1){
    warning("The number of SNPs that passed the k1_td_p filtration is <= 1, unsuitable for the subsequent analysis.")
    out<-list(dl_data1=myfit)
    class(out)<-"myerror"
    return(out)
  }
  
  fit_b1_b4<-list()
  loc1<-which(loc)
  for(i in 1:length(loc1)){
    eff<-c(myfit[[loc1[i]]]$m_hat[1,],myfit[[loc1[i]]]$k_hat[1,])
    vcov<-myfit[[loc1[i]]]$mkp_sigma
    fit_b1_b4[[i]]<-calcu_b1_b4(est_mkp=eff,vcov_mkp=vcov,est_indiv=myfit[[loc1[i]]]$mkp_ind,boot_n=boot_n,
                                delta_p_cut=delta_p_cut,length_all=nrow(g_suit),nonNA_loc=myfit[[loc1[i]]]$nonNA_loc,d=d)
  }
  names(fit_b1_b4)<-names(myfit)[loc1]
  loc2<-unlist(lapply(fit_b1_b4,FUN=function(x){x$j}))
  name_suit<-names(fit_b1_b4)[loc2]
  
  if(sum(loc2)<=1){
    warning("The number of SNPs that passed the delta_p_cut filtration is <= 1, unsuitable for the subsequent analysis.")
    out<-list(dl_data1=myfit,dl_data2=fit_b1_b4)
    class(out)<-"myerror"
    return(out)
  }
  
  if(est_type=="p2"){
    f<-function(fit_b1_b4){
      f1<-function(x){
        pnorm(-abs((x[1]-1)/x[2]))*2
      }
      out<-NULL
      for(i in 1:length(fit_b1_b4)){
        out<-c(out,c(f1(fit_b1_b4[[i]]$expb41),f1(fit_b1_b4[[i]]$expb42)))
      }
      return(out)
    }
    p2<-f(fit_b1_b4)
    p2<-matrix(p.adjust(p2,p2_cut_adj_m)<p2_cut,ncol=2,byrow=T)
    loc3<-apply(p2,1,FUN=function(x){!identical(x,c(T,T))})&loc2
  }else{
    loc3<-loc2
  }
  if(sum(loc3)<=1){
    warning("The number of SNPs that passed the p2_cut filtration (stage 1) is <= 1, unsuitable for the subsequent analysis.")
    out<-list(dl_data1=myfit,dl_data2=fit_b1_b4)
    class(out)<-"myerror"
    return(out)
  }
  
  j<-check_select_list(fit_b1_b4[loc3],select_list)
  if(!identical(j,"pass")){
    warning(j)
    out<-list(dl_data1=myfit,dl_data2=fit_b1_b4[loc3])
    class(out)<-"myerror"
    return(out)
  }

  get_data<-function(fit_b1_b4,loc,p_cut,adj_m,select_list){
    est_list<-se_list<-fb1_indiv<-est_list2<-se_list2<-list()
    if(!is.null(select_list)){
      name<-names(select_list)
      for(i in 1:length(name)){
        s0<-select_list[[i]]
        s<-fit_b1_b4[[which(names(fit_b1_b4)==name[i])]]
        est_list[[i]]<-c(s$fb11[1],s$fb12[1])[s0]
        se_list[[i]]<-c(s$fb11[2],s$fb12[2])[s0]
        est_list2[[i]]<-c(s$expb41[1],s$expb42[1])[s0]
        se_list2[[i]]<-c(s$expb41[2],s$expb42[2])[s0]
        fb1_indiv[[i]]<-cbind(cbind(s$fb101_indiv,s$fb102_indiv)[,s0])
      }
      return(list(est_list=est_list,se_list=se_list,fb1_indiv=fb1_indiv,est_list2=est_list2,se_list2=se_list2))
    }
    
    f<-function(fb1,expb4,se1,se2){
      p1<-pnorm((fb1-1)/se1,lower.tail=F)
      p2<-pnorm(expb4/se2,lower.tail=T)
      c(p1,p2)
    }
    loc<-which(loc)
    p<-NULL
    for(i in 1:length(loc)){
      s<-fit_b1_b4[[loc[i]]]
      est_list[[i]]<-c(s$fb11[1],s$fb12[1])
      se_list[[i]]<-c(s$fb11[2],s$fb12[2])
      fb1_indiv[[i]]<-cbind(s$fb101_indiv,s$fb102_indiv)
      est_list2[[i]]<-c(s$expb41[1],s$expb42[1])
      se_list2[[i]]<-c(s$expb41[2],s$expb42[2])
      p<-c(p,c(f(s$fb11[1],s$expb41[1],s$fb11[2],s$expb41[2]),f(s$fb12[1],s$expb42[1],s$fb12[2],s$expb42[2])))
    }
    
    if(adj_m%in%c("bonferroni","bonf")){
      p<-matrix(p*length(loc)*2<p_cut,ncol=4,byrow=T)
    }else{
      p<-matrix(p<p_cut,ncol=4,byrow=T)
    }
    for(i in 1:length(loc)){
      j<-c(!any(p[i,1:2]),!any(p[i,3:4]))
      est_list[[i]]<-est_list[[i]][j]
      se_list[[i]]<-se_list[[i]][j]
      fb1_indiv[[i]]<-cbind(fb1_indiv[[i]][,j])
      est_list2[[i]]<-est_list2[[i]][j]
      se_list2[[i]]<-se_list2[[i]][j]
    }
    names(est_list)<-names(se_list)<-names(fb1_indiv)<-names(est_list2)<-names(se_list2)<-names(fit_b1_b4)[loc]
    
    loc1<-unlist(lapply(est_list,FUN=function(x){length(x)!=0}))
    return(list(est_list=est_list[loc1],se_list=se_list[loc1],est_list2=est_list2[loc1],se_list2=se_list2[loc1],fb1_indiv=fb1_indiv[loc1]))
  }
  fit1<-get_data(fit_b1_b4,loc3,p3_cut,p3_cut_adj_m,select_list)
  #return(fit1)
  
  if(length(fit1$est_list)<=1){
    warning("The number of SNPs that passed the p3_cut (or select_list) filtration is <= 1, unsuitable for the subsequent analysis.")
    out<-list(dl_data1=myfit,dl_data2=fit_b1_b4[loc3])
    class(out)<-"myerror"
    return(out)
  }
  fit2<-get_lowest_q(fit1$est_list,fit1$se_list,fit1$est_list2,fit1$se_list2,fit1$fb1_indiv,n_random,n_max,twosnp_p_cut)
  if(est_type=="p2"){
    p22<-pnorm(-abs((fit2$eff2-1)/fit2$se2))*2
    p22_loc<-!p.adjust(p22,p2_cut2_adj_m)<p2_cut2
    if(sum(p22_loc)==0){
      warning("The number of SNPs that passed the p2_cut filtration (stage 2) is = 0, unsuitable for the subsequent analysis.")
      out<-list(dl_data1=myfit,dl_data2=fit_b1_b4[loc3])
      class(out)<-"myerror"
      return(out)
    }
    p22_loc<-which(p22_loc)
    fit2$eff<-fit2$eff[p22_loc]
    fit2$se<-fit2$se[p22_loc]
    fit2$eff2<-fit2$eff2[p22_loc]
    fit2$se2<-fit2$se2[p22_loc]
    fit2$fb1_indiv<-myselect(fit2$fb1_indiv,p22_loc)
  }
  
  if(length(fit2$eff)>1){
    q<-cochran_q(fit2$eff,fit2$se)
    p<-pchisq(q,length(fit2$eff)-1,lower.tail=F)
    if(p<0.05){message(paste("Heterogeneity detected.\r\nQ_stat:",round(q,2)," P_value:",formatC(p,digits=2,format="e")))}
  }else{
    q<-NA
    p<-NA
  }
  return(list(eff=fit2$eff,se=fit2$se,fb1_indiv=fit2$fb1_indiv,exph_est=fit2$eff2,exph_se=fit2$se2,name_suit=name_suit,dl_data1=myfit,dl_data2=fit_b1_b4[loc3],Q_stage1=c(Q=q,p.value=p)))
}
