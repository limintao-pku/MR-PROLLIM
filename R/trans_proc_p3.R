trans_proc_p3<-function(beta,p1_sp=NULL,p2_sp=NULL,r_sp=NULL,model_u2=F){
  b1<-beta[1]
  if(is.null(p1_sp)){p1<-1/(1+1/exp(beta[2]))}else{p1<-p1_sp}
  if(is.null(p2_sp)){p2<-1/(1+1/exp(beta[3]))}else{p2<-p2_sp}
  u1<-beta[4]
  s1_2<-exp(beta[5])
  if(is.null(r_sp)){r<-1/(1+1/exp(beta[6]))}else{r<-r_sp}
  if(model_u2){u2<-beta[7]}else{u2<-0}
  
  return(c(b1=b1,p1=p1,p2=p2,u1=u1,s1_2=s1_2,r=r,u2=u2))
}

