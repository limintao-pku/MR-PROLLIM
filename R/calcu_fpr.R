calcu_fpr<-function(p1,p2,summary_data,x,g,c,c_inherit=T,start=NULL,type=c("c","b"),
                    cd=T,cd_g_code=F,max_unique=10,n_min_limit=100,
                    control_limit_c=NULL,scale=T,
                    mc.cores=1,PSOCK=F,parallel_trace=F,dt=T,
                    nlm=T,nlm_control=list(gradtol=1e-8,steptol=1e-6,stepmax=5,iterlim=100),
                    nlminb_control=list(rel.tol=1e-12,sing.tol=1e-12,step.min=0.8,eval.max=300,iter.max=300),
                    rm_sig_p_cut=0.05,rm_sig_adj_m="bonferroni",cover=0.9999,n_boot=10000,get_sig_snp_data=NULL){
  if(!is.data.frame(summary_data)){stop("summary_data should be a data.frame.")}
  if(any(!c("snp","eff","se")%in%colnames(summary_data))){stop("summary_data should contain columns 'snp', 'eff', and 'se'.")}
  if(!identical(summary_data$snp,colnames(g))){stop("summary_data$snp should be identical to colnames(g).")}
  if(any(duplicated(summary_data$snp))){stop("Duplicated SNP names detected.")}
  summary_data$p<-pnorm(-abs(summary_data$eff/summary_data$se))*2
  
  #function
  calcu_c_prob<-function(p1,p2,cov_m,cover=0.9999,n=10000){
    get_rotate_sample<-function(s_s,mychisq,n=10000,center=c(0,0)){
      if(mychisq==0){return(list(sample=matrix(center,nrow=n,ncol=2,byrow=T),area=0))}
      A=s_s[1,1]/mychisq
      B=s_s[2,2]/mychisq
      C=s_s[1,2]/mychisq
      
      if(abs(A-B)/abs(A)<1e-3){
        r<-pi/4
        c2<-C/(sin(r)*cos(r))
        c3<-A+B
        a<-sqrt(2/(c3+c2))
        b<-sqrt(2/(c3-c2))
      }else{
        f<-function(x,v){-sin(x)*cos(x)/cos(2*x)-v}
        r<-uniroot(f=f,interval=c(-pi/4,pi/4),v=C/(B-A),extendInt ="no")$root
        b<-sqrt((-cos(2*r))/(A*sin(r)^2-B*cos(r)^2))
        a<-sqrt((-cos(2*r))/(-A*cos(r)^2+B*sin(r)^2))
      }
      
      mysample<-NULL
      while(T){
        mysample1<-cbind(runif(round(n*1.5),-a,a),runif(round(n*1.5),-b,b))
        id<-mysample1[,1]^2/a^2+mysample1[,2]^2/b^2<=1
        mysample<-rbind(mysample,mysample1[id,])
        if(nrow(mysample)>=n){break}
      }
      mysample<-mysample[1:n,]
      myrotate<-function(x,y,r,center){
        cbind(x*cos(r)-y*sin(r)+center[1],y*cos(r)+x*sin(r)+center[2])
      }
      return(list(sample=myrotate(mysample[,1],mysample[,2],r,center),area=pi*a*b))
    }

    v1<-cov_m[1,1]
    cov2<-cov_m[2:3,2:3]
    cov2_s<-solve(cov2)
    
    p3<-(1-cover)*p2
    mychisq3<-qchisq(1-p3,2)
    mychisq2<-qchisq(1-p2,2)
    if(is.infinite(mychisq3)|is.infinite(mychisq2)){stop("p2 or (1-cover) may be too small.")}
    
    s3<-NULL
    while(T){
      s1<-get_rotate_sample(cov2_s,mychisq3,n,c(0,0))
      j<-apply(s1$sample,1,FUN=function(x){t(x)%*%cov2_s%*%cbind(x)>mychisq2})
      s3<-rbind(s3,s1$sample[j,])
      if(nrow(s3)>=n){break}
    }
    p_max<-max(mvdnorm3(s3[1:n,],c(0,0),cov2))*1.5
    
    s3<-NULL
    while(T){
      s1<-get_rotate_sample(cov2_s,mychisq3,n,c(0,0))
      j<-apply(s1$sample,1,FUN=function(x){t(x)%*%cov2_s%*%cbind(x)>mychisq2})
      if(sum(j)==0){next}
      s2<-rbind(s1$sample[j,])
      p<-mvdnorm3(s2,c(0,0),cov2)/p_max
      p[p>=1]<-1
      j2<-rbinom(sum(j),1,p)==1
      s3<-rbind(s3,s2[j2,])
      if(nrow(s3)>=n){break}
    }
    s3<-s3[1:n,]
    
    q1<-abs(qnorm(p1/2)*sqrt(v1))
    out<-mean(apply(s3,1,FUN=function(x){
      u<-rbind(cov_m[1,2:3])%*%cov2_s%*%cbind(x)
      s<-sqrt(v1-rbind(cov_m[1,2:3])%*%cov2_s%*%cbind(cov_m[2:3,1]))
      pnorm(-q1,u,s)+1-pnorm(q1,u,s)
    }))
    out2<-out*cover+1*(1-cover)
    if(out2/out>1.1){warning("cover may be too small.")}
    out2
  }

  if(is.null(get_sig_snp_data)){
    mydata1<-get_sig_snp(x,g,c,c_inherit,start,type,bi_type="log",
                         p_cut=5e-8,return_dt=T,
                         cd,cd_g_code,max_unique,
                         trinary_only=T,n_min_limit,
                         control_limit_c,scale,
                         mc.cores,PSOCK,parallel_trace,dt,
                         nlm,nlm_control,nlminb_control)
  }else{
    mydata1<-get_sig_snp_data
  }
  
  name_re<-names(mydata1$eff)[str_detect(names(mydata1$eff),"_recoded")]
  name_re<-str_remove_all(name_re,"_recoded")
  if(length(name_re)>0){
    for(i in 1:length(name_re)){
      loc<-which(summary_data$snp==name_re[i])
      summary_data$eff[loc]<- -summary_data$eff[loc]
      summary_data$snp[loc]<-paste0(summary_data$snp[loc],"_recoded")
    }
  }
  
  summary_data<-myselect(summary_data,which(summary_data$snp%in%names(mydata1$eff)),"r")
  stopifnot(identical(summary_data$snp,names(mydata1$eff)))
  
  loc_sig<-which(p.adjust(mydata1$p,rm_sig_adj_m)<rm_sig_p_cut|p.adjust(summary_data$p,rm_sig_adj_m)<rm_sig_p_cut)
  if(length(loc_sig)>0){
    if(dt){cat(length(loc_sig)," significant SNPs are removed.\r\n")}
    for(i in 1:length(loc_sig)){
      mydata1$eff[[loc_sig[i]]]<-rep(NA,2)
      mydata1$vcov[[loc_sig[i]]]<-matrix(NA,2,2)
    }
  }
  loc_nonNA<-unlist(lapply(mydata1$eff,FUN=function(x){!anyNA(x)}))&(!is.na(summary_data$p))
  if(dt){cat(sum(loc_nonNA)," SNPs will be used.\r\n")}
  stopifnot(sum(loc_nonNA)>=2)

  modify_mydata2<-function(mydata2){
    for(i in 1:length(mydata2$eff)){
      mydata2$eff[[i]]<-mydata2$eff[[i]]/sqrt(diag(mydata2$vcov[[i]]))
      mydata2$vcov[[i]]<-cov2cor(mydata2$vcov[[i]])
    }
    mydata2
  }
  m_mean<-function(x){
    x11<-mean(unlist(lapply(x,FUN=function(x){x[1,1]})))
    x22<-mean(unlist(lapply(x,FUN=function(x){x[2,2]})))
    x12<-mean(unlist(lapply(x,FUN=function(x){x[1,2]})))
    matrix(c(x11,x12,x12,x22),2)
  }
  mydata2<-modify_mydata2(mydata1)
  sigma23<-m_mean(mydata2$vcov)
  
  modify_sum_data<-function(sum_data){
    sum_data$eff<-sum_data$eff/sum_data$se
    sum_data$se<-rep(1,length(sum_data$se))
    sum_data
  }
  summary_data<-modify_sum_data(summary_data)
  v1<-1
  
  mycor<-cor(cbind(summary_data$eff,unlist(lapply(mydata2$eff,FUN=function(x){x[1]})),unlist(lapply(mydata2$eff,FUN=function(x){x[2]}))),use="na.or.complete")
  
  vcov_f<-matrix(NA,3,3)
  vcov_f[1,1]<-v1
  vcov_f[2:3,2:3]<-sigma23
  vcov_f[1,2:3]<-vcov_f[2:3,1]<-sqrt(v1*diag(sigma23))*mycor[1,2:3]
  
  if(dt){cat("Start bootstrapping.\r\n")}
  out1<-out2<-matrix(NA,nrow=length(p1),ncol=length(p2))
  for(i in 1:length(p1)){
    for(j in 1:length(p2)){
      out1[i,j]<-calcu_c_prob(p1[i],p2[j],vcov_f,cover,n_boot)
      out2[i,j]<-out1[i,j]*p2[j]
    }
  }
  
  out<-list(p_fpr=out2,p_conditional=out1,par=list(vcov_f=vcov_f,cover=cover,n_boot=n_boot,n_snp=sum(loc_nonNA)),get_sig_snp_data=mydata1)
  return(out)
}
