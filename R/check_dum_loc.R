check_dum_loc<-function(c,dum_loc,name){
  if(is.null(dum_loc)){return(NULL)}
  if(identical(dum_loc,"auto")){return(NULL)}
  if(!is.list(dum_loc)){stop(name," should be a list(vectors/NULLs), NULL, or 'auto'.")}
  if(length(c)!=length(dum_loc)){stop("length(",name,") should be equal to length(c).")}
  NULL
}
