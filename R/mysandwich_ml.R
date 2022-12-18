mysandwich_ml<-function(f,beta_opt,n,loc=c(1),name,d1=.Machine$double.eps^(1/4),d2=.Machine$double.eps^(1/4),d3=.Machine$double.eps^(1/3),hessian=T,...){
  #f is an individual function
  n1<-length(beta_opt)
  mybread<-matrix(NA,ncol=n1,nrow=n1)
  mymeat_indiv<-matrix(NA,ncol=n1,nrow=n)
  for(i in 1:n1){
    for(j in 1:n1){
      beta1<-beta_opt
      beta1[j]<-beta_opt[j]+d2
      beta12<-beta11<-beta1
      beta11[i]<-beta1[i]+d1
      beta12[i]<-beta1[i]-d1
      l1<-sum(f(beta11,...))
      l2<-sum(f(beta12,...))
      der1<-(l1-l2)/(2*d1)/n
      
      beta2<-beta_opt
      beta2[j]<-beta_opt[j]-d2
      beta22<-beta21<-beta2
      beta21[i]<-beta2[i]+d1
      beta22[i]<-beta2[i]-d1
      l1<-sum(f(beta21,...))
      l2<-sum(f(beta22,...))
      der2<-(l1-l2)/(2*d1)/n
      
      mybread[i,j]<-(der1-der2)/(2*d2)
    }
  }
  loc1<-which(!apply(mybread,1,FUN=function(x){sum(x==0)==length(x)}))
  mybread<-cbind(mybread[(loc1),(loc1)])
  colnames(mybread)<-rownames(mybread)<-name[loc1]
  
  loc2<-rep(F,n1)
  loc2[loc]<-T
  loc2<-loc2[loc1]
  loc<-which(loc2)
  
  if(hessian){return(solve(mybread*n))}
  for(i in 1:n1){
    beta2<-beta1<-beta_opt
    beta1[i]<-beta_opt[i]+d3
    beta2[i]<-beta_opt[i]-d3
    mymeat_indiv[,i]<-(f(beta1,...)-f(beta2,...))/(2*d3)
  }
  mymeat_indiv<-cbind(mymeat_indiv[,loc1])
  colnames(mymeat_indiv)<-name[loc1]
  #return((solve(mybread)%*%est_vcov(mymeat_indiv)%*%t(solve(mybread))))
  #return(est_vcov(mymeat_indiv))
  fit<-est_indiv(B=mybread,M_indiv=mymeat_indiv,loc=loc,eff=beta_opt[loc1][loc])
  out<-est_vcov(fit)
  rownames(out)<-colnames(out)<-colnames(mymeat_indiv)[loc]
  return(out)
}

