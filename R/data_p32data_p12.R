data_p32data_p12<-function(data_p3,p_cut=5e-8,adj_m="none"){
  loc<-which(p.adjust(data_p3$wald_p,adj_m)<p_cut)
  out<-list()
  for(i in 1:length(loc)){
    out<-c(out,list(list(m_hat=rbind(data_p3$m_hat[loc[i],]),
                         m_sigma=rbind(data_p3$m_sigma[loc[i],]),
                         k_hat=rbind(data_p3$k_hat[loc[i],]),
                         k_sigma=data_p3$k_sigma[[loc[i]]],
                         mkp_sigma=data_p3$mkp_sigma_list[[loc[i]]],
                         mkp_ind=data_p3$mkp_ind_list[[loc[i]]],
                         nonNA_loc=data_p3$nonNA_loc_list[[loc[i]]])))
  }
  names(out)<-rownames(data_p3$m_hat)[loc]
  return(out)
}
