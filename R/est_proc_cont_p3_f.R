est_proc_cont_p3_f<-function(beta,m_matrix,sigma1_matr,k_matrix,sign_k1,sign_k2,
                             p1_sp=NULL,p2_sp=NULL,r_sp=NULL,model_u2=T,individual=F,vcov_est=F,beta_loc=NULL){
  if(!is.null(beta_loc)){
    beta0<-rep(NA,7)
    beta0[4]<-0
    beta0[beta_loc]<-beta
    beta<-beta0
  }
  
  if(vcov_est){
    b1<-beta[1]
    if(is.null(p1_sp)){p1<-beta[2]}else{p1<-p1_sp}
    if(is.null(p2_sp)){p2<-beta[3]}else{p2<-p2_sp}
    u1<-beta[4]
    s1_2<-beta[5]^2
    if(is.null(r_sp)){r<-tanh(beta[6])}else{r<-r_sp}
    if(model_u2){u2<-beta[7]}else{u2<-0}
  }else{
    b1<-beta[1]
    if(is.null(p1_sp)){p1<-1/(1+1/exp(beta[2]))}else{p1<-p1_sp}
    if(is.null(p2_sp)){p2<-1/(1+1/exp(beta[3]))}else{p2<-p2_sp}
    u1<-beta[4]
    s1_2<-exp(beta[5])
    if(is.null(r_sp)){r<-1/(1+1/exp(beta[6]))}else{r<-r_sp}
    if(model_u2){u2<-beta[7]}else{u2<-0}
  }
  
  u21<-u2*sign_k1
  u22<-2*u2*sign_k2
  
  u1_m<-(b1+u1)*k_matrix[,1]+u21
  u2_m<-(b1+u1)*k_matrix[,2]+u22
  s1<-sqrt(s1_2+sigma1_matr[,1])
  s2<-sqrt(4*s1_2+sigma1_matr[,3])
  rho<-(2*r*s1_2+sigma1_matr[,2])/s1/s2
  
  u1_m2<-(b1)*k_matrix[,1]+u21
  u2_m2<-(b1)*k_matrix[,2]+u22
  s1.2<-s1
  s2.2<-s2
  rho.2<-rho
  
  u1_m3<-u1_m2-u21
  u2_m3<-u2_m2-u22
  s1.3<-sqrt(sigma1_matr[,1])
  s2.3<-sqrt(sigma1_matr[,3])
  rho.3<-(sigma1_matr[,2])/s1.3/s2.3
  
  p<-p1*p2*mvdnorm(x1=m_matrix[,1],x2=m_matrix[,2],u1=u1_m,u2=u2_m,
                   s1=s1,s2=s2,r=rho)+
    p1*(1-p2)*mvdnorm(x1=m_matrix[,1],x2=m_matrix[,2],u1=u1_m2,u2=u2_m2,
                      s1=s1.2,s2=s2.2,r=rho.2)+
    (1-p1)*mvdnorm(x1=m_matrix[,1],x2=m_matrix[,2],u1=u1_m3,u2=u2_m3,
                   s1=s1.3,s2=s2.3,r=rho.3)
  
  if(individual){return(-log(p))}
  return(-sum(log(p)))
}
