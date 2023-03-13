meta_eff<-function(eff,se,name=NULL){
  if(length(eff)==0){
    return(list(est=c(NA,NA),name=NULL))
  }
  id<-which(!is.na(eff))
  name<-name[id]
  eff<-eff[id];se<-se[id]
  eff_m<-tryCatch({sum(1/(se^2)/sum(1/(se^2))*eff)},warning=function(w){return("w")})
  if(identical(eff_m,"w")){stop("An error occurred in the IVW")}
  se_m<-sqrt(sum(( 1/(se^2)/sum(1/(se^2)) )^2*(se^2)))
  return(list(est=c(eff_m,se_m),name=name))
}
