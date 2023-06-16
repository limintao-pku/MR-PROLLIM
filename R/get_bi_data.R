get_bi_data<-function(n,PAHP=F){
  #this function generates genotype and phenotype data for MR-PROLLIM. The exposure is binary.
  n_snp<-125
  
  c<-rbinom(n,1,0.5)
  u10<-0.9*rnorm(n)+sqrt(1-0.9^2)*c*2
  u20<-0.9*rnorm(n)+sqrt(1-0.9^2)*c*2
  
  k0<-runif(n_snp,-0.05,0.05);pg<-runif(n_snp,0.1,0.4)
  g<-crt_MR_g(n_snp,n,matrix(pg,nrow=n,ncol=n_snp,byrow=T)+
                t(t(matrix(c,nrow=n,ncol=n_snp))*k0))
  
  g3<-g2<-g
  colnames(g)<-paste0("rs",1:n_snp)
  
  if(T){
    sk_2<-0.15^2
    out_k<-MASS::mvrnorm(n_snp,c(0,0),matrix(c(sk_2*1,sk_2*2*0.9,sk_2*2*0.9,sk_2*4),2))
    k1<-out_k[,1]
    k2<-out_k[,2]
  }
  
  b1<-(-0.3)
  p1<-0.5
  p2<-0.5
  u1<-0.75
  s1_2<-0.12^2
  if(PAHP){
    rho<-1
    k3<-2*k1
  }else{
    rho<-0.9
    k3<-k2
  }
  
  u2<-0
  
  for(j in 1:n_snp){
    q1<-rbinom(1,1,p1)
    q2<-rbinom(1,1,p2)
    out_h<-q1*(MASS::mvrnorm(1,c(u2,2*u2),matrix(c(s1_2,2*rho*s1_2,2*rho*s1_2,4*s1_2),2))+
                 c(k1[j],k3[j])*u1*q2)
    
    g2[which(g[,j]==1),j]<-k1[j]
    g2[which(g[,j]==2),j]<-k2[j]
    
    g3[which(g[,j]==1),j]<-out_h[1]
    g3[which(g[,j]==2),j]<-out_h[2]
  }
  
  if(T){
    px<-g2%*%rep(1,n_snp)+0.2*c+0.3*u10
    px<-exp(px-mean(px)-2.2)
    #print(quantile(px,c(0.50,0.99,0.999,0.9999,1)))
    px[px>=1]<-1
    x<-rbinom(n,1,px)
    #print(mean(x))
  }
  
  if(T){
    py<-b1*x+g3%*%rep(1,n_snp)+0.3*u10+0.3*u20
    py<-exp(py-mean(py)-1.8)
    #print(quantile(py,c(0.50,0.99,0.999,0.9999,1)))
    py[py>=1]<-1
    y<-rbinom(n,1,py)
    #print(mean(y))
  }
  
  g[sample(1:length(g),round(length(g)*0.05))]<-NA
  
  return(list(x=x,y=y,g=g,c=c))
}