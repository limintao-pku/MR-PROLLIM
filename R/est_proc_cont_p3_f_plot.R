est_proc_cont_p3_f_plot<-function(beta,m_matrix,sigma1_matr,k_matrix,sign_k1,sign_k2,
                                  p1_sp=NULL,p2_sp=NULL,r_sp=NULL,model_u2=T){
  b1<-beta[1]
  if(is.null(p1_sp)){p1<-1/(1+1/exp(beta[2]))}else{p1<-p1_sp}
  if(is.null(p2_sp)){p2<-1/(1+1/exp(beta[3]))}else{p2<-p2_sp}
  u1<-beta[4]
  s1_2<-exp(beta[5])
  if(is.null(r_sp)){r<-1/(1+1/exp(beta[6]))}else{r<-r_sp}
  if(model_u2){u2<-beta[7]}else{u2<-0}
  
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
  
  pd<-mvdnorm(x1=m_matrix[,1],x2=m_matrix[,2],u1=u1_m,u2=u2_m,
              s1=s1,s2=s2,r=rho)
  pun<-mvdnorm(x1=m_matrix[,1],x2=m_matrix[,2],u1=u1_m2,u2=u2_m2,
               s1=s1.2,s2=s2.2,r=rho.2)
  pno<-mvdnorm(x1=m_matrix[,1],x2=m_matrix[,2],u1=u1_m3,u2=u2_m3,
               s1=s1.3,s2=s2.3,r=rho.3)
  p<-p1*p2*pd+p1*(1-p2)*pun+(1-p1)*pno
  return(cbind((1-p1)*pno/p,p1*(1-p2)*pun/p,p1*p2*pd/p))
}
