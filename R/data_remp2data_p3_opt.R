data_remp2data_p3_opt<-function(data_remp){
  stopifnot(identical(class(data_remp),"MR-PROLLIM output"))
  if(is.null(data_remp$parameter$est_type)){stop("data_remp should be a random-effects MR-PROLLIM output")}
  if(!is.null(data_remp$prior_est)){
    if(data_remp$parameter$est_type=="c"){
      m_hat<-data_remp$data$m_hat
      sigma1_matr<-data_remp$data$m_sigma
      k_hat<-data_remp$data$k_hat
      sigma2_matr<-data_remp$data$k_sigma
      mk_sigma_list<-data_remp$data$mk_sigma_list
      post_sample_k1<-data_remp$data$post_sample_k1;post_sample_k2<-data_remp$data$post_sample_k2
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
      wald_p<-data_remp$data$wald_p;u_input<-data_remp$data$u_input;post_messages<-data_remp$data$post_messages
      fit_k<-data_remp$prior_est
      
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
    
    if(data_remp$parameter$est_type=="b"){
      m_hat<-data_remp$data$m_hat
      sigma1_matr<-data_remp$data$m_sigma
      k_hat<-data_remp$data$k_hat
      sigma2_list<-data_remp$data$k_sigma
      mkp_sigma_list<-data_remp$data$mkp_sigma_list
      post_sample_k1<-data_remp$data$post_sample_k1;post_sample_k2<-data_remp$data$post_sample_k2;post_sample_p<-data_remp$data$post_sample_p
      post_sample_k1_exp1<-exp(post_sample_k1)-1
      post_sample_k2_exp1<-exp(post_sample_k2)-1
      prepare_p3.1_bi<-function(k_hat,post_sample_k1,post_sample_k2,post_sample_p,mkp_sigma_list){
        out1<-out2<-matrix(NA,nrow=nrow(post_sample_k1),ncol=ncol(post_sample_k1))
        out3<-list()
        for(i in 1:nrow(k_hat)){
          s<-mkp_sigma_list[[i]][1:2,3:5]%*%solve(mkp_sigma_list[[i]][3:5,3:5])
          x<-cbind(k_hat[i,1]-post_sample_k1[i,],k_hat[i,2]-post_sample_k2[i,],k_hat[i,3]-post_sample_p[i,])%*%t(s)
          out1[i,]<-x[,1]
          out2[i,]<-x[,2]
          out3[[i]]<-mkp_sigma_list[[i]][1:2,1:2]-s%*%mkp_sigma_list[[i]][3:5,1:2]
        }
        return(list(g1_matr=out1,g2_matr=out2,sigma_prime_list=out3))
      }
      data_opt_p3.1<-prepare_p3.1_bi(k_hat,post_sample_k1,post_sample_k2,post_sample_p,mkp_sigma_list)
      wald_p<-data_remp$data$wald_p;u_input<-data_remp$data$u_input;post_messages<-data_remp$data$post_messages
      fit_k<-data_remp$prior_est
      
      return(list(m_hat=m_hat,
                  m_sigma=sigma1_matr,
                  k_hat=k_hat,
                  k_sigma=sigma2_list,
                  mkp_sigma_list=mkp_sigma_list,
                  post_sample_k1=post_sample_k1,post_sample_k2=post_sample_k2,post_sample_p=post_sample_p,
                  post_sample_k1_exp1=post_sample_k1_exp1,post_sample_k2_exp1=post_sample_k2_exp1,
                  data_opt_p3.1=data_opt_p3.1,
                  wald_p=wald_p,u_input=u_input,post_messages=post_messages,
                  fit_k=fit_k))
    }
  }else{
    if(data_remp$parameter$est_type=="c"){
      m_hat<-data_remp$data$m_hat;sigma1_matr<-data_remp$data$m_sigma
      k_hat<-data_remp$data$k_hat;sigma2_matr<-data_remp$data$k_sigma
      mk_sigma_list<-data_remp$data$mk_sigma_list
      wald_p<-data_remp$data$wald_p;u_input<-data_remp$data$u_input
      
      return(list(m_hat=m_hat,m_sigma=sigma1_matr,
                  k_hat=k_hat,k_sigma=sigma2_matr,
                  mk_sigma_list=mk_sigma_list,
                  wald_p=wald_p,u_input=u_input))
    }
    
    if(data_remp$parameter$est_type=="b"){
      m_hat<-data_remp$data$m_hat;sigma1_matr<-data_remp$data$m_sigma
      k_hat<-data_remp$data$k_hat;sigma2_list<-data_remp$data$k_sigma
      mkp_sigma_list<-data_remp$data$mkp_sigma_list
      wald_p<-data_remp$data$wald_p;u_input<-data_remp$data$u_input
      
      return(list(m_hat=m_hat,m_sigma=sigma1_matr,
                  k_hat=k_hat,k_sigma=sigma2_list,
                  mkp_sigma_list=mkp_sigma_list,
                  wald_p=wald_p,u_input=u_input))
    }
  }
  NULL
}
