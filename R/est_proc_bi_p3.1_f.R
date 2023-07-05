est_proc_bi_p3.1_f<-function(beta,m_matrix,post_sample_k1,post_sample_k1_exp1,post_sample_k2,post_sample_k2_exp1,post_sample_p,
                             g1_matr,g2_matr,sigma_prime_list,sign_k1,sign_k2,
                             p1_sp=NULL,p2_sp=NULL,r_sp=NULL,model_u2=T,log_appr=2,individual=F,vcov_est=F,indiv_p=F,cpp=T,t_b1=F,beta_loc=NULL){
  if(!is.null(beta_loc)){
    beta0<-rep(NA,7)
    beta0[beta_loc]<-beta
    beta<-beta0
  }
  if(vcov_est){
    b1_t<-beta[1]
    if(is.null(p1_sp)){p1<-beta[2]}else{p1<-p1_sp}
    if(is.null(p2_sp)){p2<-beta[3]}else{p2<-p2_sp}
    u1<-beta[4]
    s1_2<-beta[5]^2
    if(is.null(r_sp)){r<-tanh(beta[6])}else{r<-r_sp}
    if(model_u2){u2<-beta[7]}else{u2<-0}
  }else{
    if(t_b1){
      b1_t<-beta[1]
    }else{
      b1<-beta[1]
      b1_t<-(1-exp(-b1))
    }
    if(is.null(p1_sp)){p1<-1/(1+1/exp(beta[2]))}else{p1<-p1_sp}
    if(is.null(p2_sp)){p2<-1/(1+1/exp(beta[3]))}else{p2<-p2_sp}
    u1<-beta[4]
    s1_2<-exp(beta[5])
    if(is.null(r_sp)){r<-1/(1+1/exp(beta[6]))}else{r<-r_sp}
    if(model_u2){u2<-beta[7]}else{u2<-0}
  }
  s_h<-matrix(c(s1_2,2*s1_2*r,2*s1_2*r,4*s1_2),2)
  
  if(indiv_p){
    if(cpp){
      if(is.null(p1_sp)){p1_sp<- -1}
      if(is.null(p2_sp)){p2_sp<- -1}
      p<-loop_p3_1_bi_indiv_p(p1_sp, p2_sp, b1_t, u1, u2, p1, p2, log_appr, m_matrix, 
                              s_h, post_sample_k1, post_sample_k1_exp1, post_sample_k2, 
                              post_sample_k2_exp1, post_sample_p, g1_matr, g2_matr, sigma_prime_list, sign_k1, sign_k2)
    }else{
      pd<-pun<-pno<-matrix(NA,nrow=nrow(post_sample_k1),ncol=ncol(post_sample_k1))
      for(i in 1:nrow(m_matrix)){
        if(identical(log_appr,2)){
          z1<-(exp(post_sample_k1[i,])-1)*b1_t*post_sample_p[i,]+1
          z2<-(exp(post_sample_k2[i,])-1)*b1_t*post_sample_p[i,]+1
          z1_m<-mean(z1);z2_m<-mean(z2)
          z1_z<-z1/z1_m-1
          z2_z<-z2/z2_m-1
          z1_log<-log(z1_m)+z1_z-(z1_z)^2
          z2_log<-log(z2_m)+z2_z-(z2_z)^2
        }else if(identical(log_appr,1)){
          z1<-(exp(post_sample_k1[i,])-1)*b1_t*post_sample_p[i,]+1
          z2<-(exp(post_sample_k2[i,])-1)*b1_t*post_sample_p[i,]+1
          z1_m<-mean(z1);z2_m<-mean(z2)
          z1_log<-log(z1_m)+z1/z1_m-1
          z2_log<-log(z2_m)+z2/z2_m-1
        }else{
          z1_log<-log((exp(post_sample_k1[i,])-1)*b1_t*post_sample_p[i,]+1)
          z2_log<-log((exp(post_sample_k2[i,])-1)*b1_t*post_sample_p[i,]+1)
        }
        
        #double pleiotropy
        if(identical(p2_sp,0))(pd[i,]<-0)else{
          pd[i,]<-(mvdnorm3(m_matrix[i,],
                              cbind(u1*post_sample_k1[i,]+u2*sign_k1[i]+g1_matr[i,]+z1_log,u1*post_sample_k2[i,]+2*u2*sign_k2[i]+g2_matr[i,]+z2_log),
                              s_h+sigma_prime_list[[i]]))
        }
        
        #uncorrelated
        if(identical(p2_sp,1)){pun[i]<-0}else{
          pun[i,]<-(mvdnorm3(m_matrix[i,],
                               cbind(u2*sign_k1[i]+g1_matr[i,]+z1_log,2*u2*sign_k2[i]+g2_matr[i,]+z2_log),
                               s_h+sigma_prime_list[[i]]))
        }
        
        #no pleiotropy
        if(identical(p1_sp,1)){pno[i]<-0}else{
          pno[i,]<-(mvdnorm3(m_matrix[i,],
                               cbind(z1_log+g1_matr[i,],z2_log+g2_matr[i,]),
                               sigma_prime_list[[i]]))
        }
      }
      p<-p1*p2*pd+p1*(1-p2)*pun+(1-p1)*pno
    }
    return(p)
  }
  
  if(cpp){
    if(is.null(p1_sp)){p1_sp<- -1}
    if(is.null(p2_sp)){p2_sp<- -1}
    p<-loop_p3_1_bi(p1_sp, p2_sp, b1_t, u1, u2, p1, p2, log_appr, m_matrix, 
                    s_h, post_sample_k1, post_sample_k1_exp1, post_sample_k2, 
                    post_sample_k2_exp1, post_sample_p, g1_matr, g2_matr, sigma_prime_list, sign_k1, sign_k2)
  }else{
    pd<-pun<-pno<-rep(0,nrow(m_matrix))
    for(i in 1:nrow(m_matrix)){
      if(identical(log_appr,2)){
        z1<-(exp(post_sample_k1[i,])-1)*b1_t*post_sample_p[i,]+1
        z2<-(exp(post_sample_k2[i,])-1)*b1_t*post_sample_p[i,]+1
        z1_m<-mean(z1);z2_m<-mean(z2)
        z1_z<-z1/z1_m-1
        z2_z<-z2/z2_m-1
        z1_log<-log(z1_m)+z1_z-(z1_z)^2
        z2_log<-log(z2_m)+z2_z-(z2_z)^2
      }else if(identical(log_appr,1)){
        z1<-(exp(post_sample_k1[i,])-1)*b1_t*post_sample_p[i,]+1
        z2<-(exp(post_sample_k2[i,])-1)*b1_t*post_sample_p[i,]+1
        z1_m<-mean(z1);z2_m<-mean(z2)
        z1_log<-log(z1_m)+z1/z1_m-1
        z2_log<-log(z2_m)+z2/z2_m-1
      }else{
        z1_log<-log((exp(post_sample_k1[i,])-1)*b1_t*post_sample_p[i,]+1)
        z2_log<-log((exp(post_sample_k2[i,])-1)*b1_t*post_sample_p[i,]+1)
      }
      
      #double pleiotropy
      if(identical(p2_sp,0))(pd[i]<-0)else{
        pd[i]<-mean(mvdnorm3(m_matrix[i,],
                               cbind(u1*post_sample_k1[i,]+u2*sign_k1[i]+g1_matr[i,]+z1_log,u1*post_sample_k2[i,]+2*u2*sign_k2[i]+g2_matr[i,]+z2_log),
                               s_h+sigma_prime_list[[i]]))
      }
      
      #uncorrelated
      if(identical(p2_sp,1)){pun[i]<-0}else{
        pun[i]<-mean(mvdnorm3(m_matrix[i,],
                                cbind(u2*sign_k1[i]+g1_matr[i,]+z1_log,2*u2*sign_k2[i]+g2_matr[i,]+z2_log),
                                s_h+sigma_prime_list[[i]]))
      }
      
      #no pleiotropy
      if(identical(p1_sp,1)){pno[i]<-0}else{
        pno[i]<-mean(mvdnorm3(m_matrix[i,],
                                cbind(z1_log+g1_matr[i,],z2_log+g2_matr[i,]),
                                sigma_prime_list[[i]]))
      }
    }
    p<-p1*p2*pd+p1*(1-p2)*pun+(1-p1)*pno
  }
  
  if(individual){return(-log(p))}
  return(-sum(log(p)))
}
