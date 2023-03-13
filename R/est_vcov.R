est_vcov<-function(est_indiv_out){
  indiv<-scale(est_indiv_out,scale=F)
  out<-t(indiv)%*%indiv/(nrow(indiv)-1)/nrow(indiv)
  colnames(out)<-rownames(out)<-NULL
  return(out)
}

