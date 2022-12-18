est_variance_reMP<-function(f,beta_opt,...,beta_name=NULL,name="est_variance_reMP",sandwich=T,hessian=T,
                            adj_rs=F,boot_n=10000,mc.cores=1,parallel_trace=F,cpp=T,
                            d1=.Machine$double.eps^(1/4),d2=.Machine$double.eps^(1/4),d3=.Machine$double.eps^(1/3)){
  
  n1<-length(beta_opt)
  names(beta_opt)<-NULL
  
  der1_list<-list()
  for(i in 1:n1){
    beta1<-beta2<-beta_opt
    beta1[i]<-beta_opt[i]+d3
    beta2[i]<-beta_opt[i]-d3
    der1_list[[i]]<-(f(beta1,...,indiv_p=T)-f(beta2,...,indiv_p=T))/2/d3
  }
  
  #is a matrix because of posterior samples
  p_matr<-f(beta_opt,...,indiv_p=T)
  
  der2_list<-list()
  for(i in 1:n1){
    for(j in 1:n1){
      beta1<-beta_opt
      beta1[j]<-beta_opt[j]+d2
      beta12<-beta11<-beta1
      beta11[i]<-beta1[i]+d1
      beta12[i]<-beta1[i]-d1
      der1<-(f(beta11,...,indiv_p=T)-f(beta12,...,indiv_p=T))/2/d1
      
      beta2<-beta_opt
      beta2[j]<-beta_opt[j]-d2
      beta22<-beta21<-beta2
      beta21[i]<-beta2[i]+d1
      beta22[i]<-beta2[i]-d1
      der2<-(f(beta21,...,indiv_p=T)-f(beta22,...,indiv_p=T))/2/d1
      
      der2_list[[n1*(i-1)+j]]<-(der1-der2)/2/d2
    }
  }
  
  loc<-unlist(lapply(der1_list,FUN=function(x){
    sum(x==0)!=length(x)
  }))
  loc2<-unlist(lapply(der2_list,FUN=function(x){
    sum(x==0)!=length(x)
  }))
  
  n1<-sum(loc)
  der1_list<-der1_list[loc]
  der2_list<-der2_list[loc2]
  stopifnot(length(der1_list)==n1&length(der2_list)==(n1*n1))
  
  f<-function(){
    nr<-nrow(p_matr)
    nc<-ncol(p_matr)
    p_b<-rep(NA,nr)
    der1_matr<-matrix(NA,nrow=nr,ncol=n1)
    der2_array<-array(NA,dim=c(nr,n1,n1))
    for(j in 1:nr){
      p_b[j]<-mean(p_matr[j,])
      for(k in 1:n1){
        der1_matr[j,k]<-mean(der1_list[[k]][j,])
        for(l in 1:n1){
          der2_array[j,k,l]<-mean(der2_list[[n1*(k-1)+l]][j,])
        }
      }
    }
    der1<-NULL
    der2<-matrix(NA,n1,n1)
    der2_indiv<-matrix(NA,nrow=nr,ncol=n1*n1)
    for(j in 1:n1){
      der1<-cbind(der1,-1/p_b*der1_matr[,j])
      for(k in 1:n1){
        der2_indiv[,n1*(j-1)+k]<- -(-1/(p_b^2)*der1_matr[,j]*der1_matr[,k]+1/p_b*der2_array[,j,k])
        der2[j,k]<-mean(der2_indiv[,n1*(j-1)+k])
      }
    }
    invh<-tryCatch({solve(der2)/nr},error=function(e){matrix(NA,n1,n1)})
    sandw<-tryCatch({solve(der2)%*%(cov(der1)/nr)%*%t(solve(der2))},error=function(e){matrix(NA,n1,n1)})
    return(list(der1,der2_indiv,invh,sandw))
  }
  
  if(!is.null(name)){
    name<-paste0(name,":")
  }
  
  out_cov<-f()
  der1<-out_cov[[1]]
  der2_indiv<-out_cov[[2]]
  invh<-out_cov[[3]]
  out_cov<-sandw<-out_cov[[4]]
  
  if(anyNA(invh)|anyNA(sandw)){
    warning(name,": NAs (or non-invertible hessian) detected in computing asymptotic variance estimates.")
    out_f<-matrix(NA,n1,n1)
    colnames(out_f)<-rownames(out_f)<-beta_name[loc]
    return(out_f)
  }
  
  nr<-nrow(p_matr)
  
  if(hessian&sandwich){
    if(out_cov[1,1]<invh[1,1]&!any(diag(invh)<0)){
      out_cov<-invh
    }
  }
  
  if(hessian&!sandwich){
    out_cov<-invh
  }
  
  colnames(out_cov)<-rownames(out_cov)<-beta_name[loc]
  if(!adj_rs){return(out_cov)}
  
  if(cpp){
    myboot<-function(n_rep){
      out<-my_boot_adj(n_rep[[1]],n1,p_matr,der1_list,der2_list)
      myout<-matrix(NA,nrow=n1,ncol=n_rep[[1]])
      for(i in 1:(n_rep[[1]])){
        myout[,i]<-tryCatch({-solve(out[[2]][[i]])%*%(out[[1]][[i]])},error=function(e){cbind(rep(NA,n1))})
      }
      return(myout)
    }
  }else{
    myboot<-function(n_rep){
      t0<-Sys.time()
      my_proc<-0
      pid<-n_rep[[2]]
      n_rep<-n_rep[[1]]
      
      nr<-nrow(p_matr)
      nc<-ncol(p_matr)
      out<-matrix(NA,nrow=n1,ncol=n_rep)
      for(i in 1:n_rep){
        p_b<-rep(NA,nr)
        der1_matr<-matrix(NA,nrow=nr,ncol=n1)
        der2_array<-array(NA,dim=c(nr,n1,n1))
        for(j in 1:nr){
          id<-sample(1:nc,nc,replace=T)
          p_b[j]<-mean(p_matr[j,id])
          for(k in 1:n1){
            der1_matr[j,k]<-mean(der1_list[[k]][j,id])
            for(l in 1:n1){
              der2_array[j,k,l]<-mean(der2_list[[n1*(k-1)+l]][j,id])
            }
          }
        }
        der1<-rep(NA,n1)
        der2<-matrix(NA,n1,n1)
        for(j in 1:n1){
          der1[j]<- -mean(1/p_b*der1_matr[,j])
          for(k in 1:n1){
            der2[j,k]<- -mean(-1/(p_b^2)*der1_matr[,j]*der1_matr[,k]+1/p_b*der2_array[,j,k])
          }
        }
        out[,i]<- tryCatch({-solve(der2)%*%der1},error=function(e){rep(NA,n1)})
        if(mc.cores==1&parallel_trace){my_moni("iter",i,n_rep)}
        if(mc.cores!=1&parallel_trace){
          my_proc<-my_moni2(paste0("Child process ",pid,":"),i,n_rep,my_proc,time=T,t0=t0)
        }
      }
      return(out)
    }
  }
  
  stopifnot(boot_n>=mc.cores)
  out<-mclapply(even_allo2(boot_n,mc.cores),FUN=myboot,mc.cores=mc.cores)
  boot_out<-NULL
  for(i in 1:length(out)){
    boot_out<-cbind(boot_out,out[[i]])
  }
  boot_out<-t(boot_out)
  boot_out<-na.omit(boot_out)
  add_cov<-cov(boot_out,use="na.or.complete")
  if(nrow(boot_out)/boot_n<0.9){warning(paste(name,"Only",nrow(boot_out),"repeats in adj_rs procedure are valid."))}
  return(out_cov+add_cov)
}
