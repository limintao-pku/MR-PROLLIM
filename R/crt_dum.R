crt_dum<-function(x,ref=0){
  x_u=sort(unique(x),decreasing=F)
  x_u=x_u[-which(x_u==ref)]
  out=NULL
  for(i in 1:(length(x_u))){
    out=cbind(out,as.numeric(x==(x_u[i])))
  }
  colnames(out)=x_u
  return(out)
}

