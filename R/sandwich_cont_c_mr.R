sandwich_cont_c_mr<-function(loc,beta_opt,y,g,g_dum,c,a1,a1_indiv="auto",x=NULL,b1_sp=NULL,b1_sp_indiv=NULL,b4_prior,
                             return_indiv=F,name=NULL,loc_indiv=NULL,se=NULL,se_m=NULL,b1_indiv_matr=NULL,
                             a1_matrix=NULL,b1_a1_list=NULL,w_a1_list=NULL,
                             direct_M=F,return_w=F,d=.Machine$double.eps^(1/3)){
  mymeat_cont_mr2<-function(mymeat_cont_mr_indiv,loc_indiv,loc_b1,se,se_m,b1_indiv_matr){
    M<-mymeat_cont_mr_indiv
    n<-nrow(M)
    n1<-ncol(M)
    M5<-matrix(NA,nrow=n1,ncol=n1)
    se_loc<-which(!is.na(se))
    
    if(length(se_loc)==0|anyNA(se_m)|identical(unique(c(M)),NA)){print("s");return(M5)}
    se<-se[se_loc]
    b1_indiv_matr<-myselect(b1_indiv_matr,se_loc)
    
    for(i in 1:n1){
      for(j in 1:n1){
        M6<-cbind(M[,i],M[,j])
        M6<-na.omit(M6)
        if(nrow(M6)==0){M5[i,j]<-0}else{
          M5[i,j]<-cov(M6[,1],M6[,2])*nrow(M6)/n/n
        }
      }
    }
    
    sum.na<-function(x){sum(!is.na(x))}
    N<-nrow(b1_indiv_matr)
    w<-1/se^2/sum(1/se^2)
    
    myout<-NA
    for(i in 1:n1){
      o<-rep(NA,N)
      o[loc_indiv]<-M[,i]
      out<-0
      for(j in 1:ncol(b1_indiv_matr)){
        M6<-cbind(o,b1_indiv_matr[,j])
        M6<-na.omit(M6)
        if(nrow(M6)==0){o1<-0}else{
          o1<-cov(M6[,1],M6[,2])
        }
        out<-out+o1*w[j]/n/sum.na(b1_indiv_matr[,j])*nrow(M6)
      }
      myout[i]<-out
    }
    myout[loc_b1]<-se_m^2
    M5[loc_b1,]<-myout
    M5[,loc_b1]<-myout
    return(M5)
  }
  mymeat_cont_mr3<-function(mymeat_cont_mr_indiv,loc_indiv,a1_matrix,loc_b1,se,se_m,b1_indiv_matr,b1_a1_list,w_a1_list){
    M<-mymeat_cont_mr_indiv
    n<-nrow(M)
    n1<-ncol(M)
    M5<-matrix(NA,nrow=n1,ncol=n1)
    se_loc<-which(!is.na(se))
    
    if(length(se_loc)==0|anyNA(se_m)|identical(unique(c(M)),NA)){return(M5)}
    se<-se[se_loc]
    b1_indiv_matr<-myselect(b1_indiv_matr,se_loc)
    b1_a1_list<-b1_a1_list[se_loc]
    w_a1_list<-w_a1_list[se_loc]
    
    for(i in 1:n1){
      for(j in 1:n1){
        M6<-cbind(M[,i],M[,j])
        M6<-na.omit(M6)
        if(nrow(M6)==0){M5[i,j]<-0}else{
          M5[i,j]<-cov(M6[,1],M6[,2])*nrow(M6)/n/n
        }
      }
    }
    
    sum.na<-function(x){sum(!is.na(x))}
    calcu_cor<-function(a_matr1,a_matr2,w){
      out<-rep(0,ncol(a_matr1))
      for(i in 1:ncol(a_matr1)){
        o<-cbind(a_matr1[,i],a_matr2%*%w)
        o<-na.omit(o)
        if(nrow(o)==0){out[i]<-0}else{
          out[i]<-cov(o[,1],o[,2])/sum.na(a_matr1[,1])/sum.na(a_matr2[,1])*nrow(o)
        }
      }
      return(out)
    }
    
    N<-nrow(b1_indiv_matr)
    w<-1/se^2/sum(1/se^2)
    
    myout<-NA
    for(i in 1:n1){
      o<-rep(NA,N)
      o[loc_indiv]<-M[,i]
      out<-0
      for(j in 1:ncol(b1_indiv_matr)){
        M6<-cbind(o,b1_indiv_matr[,j])
        M6<-na.omit(M6)
        if(nrow(M6)==0){o1<-0}else{
          o1<-cov(M6[,1],M6[,2])
        }
        out<-out+o1*w[j]/n/sum.na(b1_indiv_matr[,j])*nrow(M6)
      }
      myout[i]<-out
    }
    myout[loc_b1]<-se_m^2
    M5[loc_b1,]<-myout
    M5[,loc_b1]<-myout
    M5[1:ncol(a1_matrix),1:ncol(a1_matrix)]<-est_vcov(na.omit(a1_matrix))
    b1_cor<-rep(0,ncol(a1_matrix))
    for(j in 1:ncol(b1_indiv_matr)){
      b1_cor<-b1_cor+calcu_cor(a1_matrix,b1_a1_list[[j]],w_a1_list[[j]])*w[j]
    }
    M5[loc_b1,1:ncol(a1_matrix)]<-b1_cor
    M5[1:ncol(a1_matrix),loc_b1]<-b1_cor
    return(M5)
  }
  
  if(anyNA(beta_opt)){
    if(return_w){
      return(list(matrix(NA,nrow=length(loc),ncol=length(loc)),matrix(NA,nrow=length(loc),ncol=length(a1))))
    }else{
      return(matrix(NA,ncol=length(loc),nrow=length(loc)))
    }
  }
  if(identical(a1_indiv,"auto")){
    d_matr<-cbind(g_dum,c,1)
    b<-solve(t(d_matr)%*%d_matr)%*%t(d_matr)%*%x
    a1_indiv_c<-est_indiv(B=-t(d_matr)%*%d_matr/length(x),M_indiv=c(x-d_matr%*%b)*d_matr,1:length(a1),rep(0,length(a1)))
  }else{
    a1_indiv_c<-t(t(a1_indiv)-a1)
  }
  
  n<-length(y)
  n1<-length(c(a1,b1_sp))
  n2<-length(beta_opt)
  B1<-cbind(diag(rep(-1,n1)),matrix(0,ncol=n2,nrow=n1))
  
  beta0<-c(a1,b1_sp,beta_opt)
  B2<-matrix(NA,nrow=n2,ncol=n1+n2)
  for(i in 1:(n1+n2)){
    beta1<-beta2<-beta0
    beta1[i]<-beta1[i]+d
    beta2[i]<-beta2[i]-d
    
    b1_sp1<-b1_sp2<-NULL
    a11<-beta1[1:(n1)]
    a12<-beta2[1:(n1)]
    if(!is.null(b1_sp)){
      b1_sp1<-beta1[n1]
      b1_sp2<-beta2[n1]
      a11<-beta1[1:(n1-1)]
      a12<-beta2[1:(n1-1)]
    }
    f1<-f_cont_c_mr(beta=beta1[-(1:n1)],y=y,g=g,g_dum=g_dum,c=c,a1=a11,b1_sp=b1_sp1,b4_prior=b4_prior)
    f1<-attr(f1,"gradient")
    f2<-f_cont_c_mr(beta=beta2[-(1:n1)],y=y,g=g,g_dum=g_dum,c=c,a1=a12,b1_sp=b1_sp2,b4_prior=b4_prior)
    f2<-attr(f2,"gradient")
    B2[,i]<-(f1-f2)/2/d/n
  }
  B<-rbind(B1,B2)
  
  M_indiv<-cbind(a1_indiv_c,b1_sp_indiv-b1_sp,f_cont_c_mr(beta_opt,y=y,g=g,g_dum=g_dum,c=c,a1=a1,b1_sp=b1_sp,b4_prior=b4_prior,der_indiv=T,grad=T))
  if(!is.null(b1_sp)&!direct_M){
    if(is.null(b1_a1_list)){
      M<-mymeat_cont_mr2(mymeat_cont_mr_indiv=M_indiv,loc_indiv=loc_indiv,loc_b1=length(a1)+1,
                         se=se,se_m=se_m,b1_indiv_matr=b1_indiv_matr)
    }else{
      M<-mymeat_cont_mr3(mymeat_cont_mr_indiv=M_indiv,loc_indiv=loc_indiv,a1_matrix=a1_matrix,loc_b1=length(a1)+1,
                         se=se,se_m=se_m,b1_indiv_matr=b1_indiv_matr,b1_a1_list=b1_a1_list,w_a1_list=w_a1_list)
    }
    vcov<-tryCatch({solve(B)%*%M%*%t(solve(B))},error=function(e){"e"})
    if(identical(vcov,"e")){
      warning(paste(name,"failed in the calculation of sandwich variance estimator."))
      vcov<-matrix(NA,ncol=ncol(B),nrow=nrow(B))
    }
    return(vcov[loc,loc])
  }
  fit<-est_indiv(B=B,M_indiv=M_indiv,loc=loc,eff=beta0[loc],name=name)
  
  if(return_indiv){
    if(return_w){
      inv_B<-tryCatch({solve(B)},error=function(e){return(matrix(NA,nrow=nrow(B),ncol=ncol(B)))})
      return(list(fit,myselect(myselect(inv_B,loc,"r"),1:length(a1))))
    }else{return(fit)}
  }
  out<-est_vcov(fit)
  colnames(out)<-rownames(out)<-NULL
  return(out)
}
