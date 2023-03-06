check_start<-function(g,start,name){
  if(is.null(start)){return(NULL)}
  if(!is.list(start)){stop(name," should be a list(vectors) or NULL.")}
  if(!is.vector(start[[1]])){stop(name," should be a list(vectors) or NULL.")}
  if(ncol(g)!=length(start)){stop("length(",name,") should be equal to ncol(g).")}
  NULL
}
