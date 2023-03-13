est_proc_bi_p1.3_dl<-function(dl_data1_nohp,cor_correct,length_all,parallel_trace){
  eff<-se<-eff_indiv<-NULL
  for(i in 1:length(dl_data1_nohp)){
    s<-dl_data1_nohp[[i]]
    fit<-dl_nohp(est_mkp=c(s$m_hat[1,],s$k_hat[1,]),vcov_mkp=NULL,mkp_indiv=s$mkp_ind,name=names(dl_data1_nohp)[i],nonNA_loc=s$nonNA_loc,length_all=length_all)
    eff<-cbind(eff,fit$est)
    se<-cbind(se,fit$se)
    if(cor_correct){
      eff_indiv<-cbind(eff_indiv,fit$est_indiv)
    }
    if(parallel_trace){my_moni("SNP",i,length(dl_data1_nohp))}
  }
  colnames(eff)<-colnames(se)<-colnames(eff_indiv)<-names(dl_data1_nohp)
  return(list(eff=eff,se=se,eff_indiv=eff_indiv))
}
