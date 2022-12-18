cut_num<-function(x,n){
  if(length(x)<n){n<-length(x)}
  out<-list()
  n1<-length(x)
  x1<-cut(1:n1,breaks=seq(0,n1,length.out=n+1),labels=1:n)
  for(i in 1:n){
    out[[i]]<-list(x[which(x1==i)],i)
  }
  return(out)
}

