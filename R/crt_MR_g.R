crt_MR_g<-function(n_snp,n_sample,p_minor){
  #this function generates genotype data
  p_minor<-rbind(p_minor)
  out<-matrix(NA,nrow=n_sample,ncol=n_snp)
  for(i in 1:n_snp){
    x<-rbinom(n_sample,1,1-(1-p_minor[,i])^2)
    y<-rbinom(n_sample,1,p_minor[,i]^2/(1-(1-p_minor[,i])^2))+1L
    x[x==1]<-y[x==1]
    out[,i]<-x
  }
  return(out)
}