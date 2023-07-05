est_indiv<-function(B,M_indiv,loc,eff,name=NULL){
  if(anyNA(eff)){return(matrix(NA,ncol=length(loc),nrow=length(loc)))}
  for(i in 1:ncol(M_indiv)){
    M_indiv[which(is.na(M_indiv[,i])),i]<-0
  }
  B1<-tryCatch({solve(B)},error=function(e){"e"})
  if(identical(B1,"e")|anyNA(B)){
    warning(paste(name,"(one SNP) failed in the calculation of sandwich variance estimator."))
    return(matrix(NA,ncol=length(loc),nrow=length(loc)))
  }
  
  out<-NULL
  for(i in 1:length(loc)){
    z<-c(M_indiv%*%(-B1[loc[i],]))
    z<-z+eff[i]
    out<-cbind(out,z)
  }
  return(out)
}

