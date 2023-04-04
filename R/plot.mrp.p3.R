plot.mrp.p3<-function(est_out,ci_cover=0.95,use_trans=T,trans_lower=0,reorder=NULL,
                      emph_loc=NULL,emph_col="purple",
                      cluster_type=3,cluster_color=c("black","blue","red"),bg_color="white",
                      control_plot1=list(),control_plot2=list(),expr1=NULL,expr2=NULL,
                      control_lines=list(),control_points=list(),control_abline=list(),
                      sub_fig_label=NULL,
                      sub_fig_line=3,sub_fig_k=0.08,sub_fig_cex=1,
                      control_legend=NULL,
                      interactive=F,obj_name=NULL){
  if(!identical(class(est_out),"MR-PROLLIM output")){
    stop("est_out should be an MR-PROLLIM output.")
  }
  if(!cluster_type%in%c(2,3)){stop("cluster_type should be 2 or 3.")}
  
  trans<-function(x,lower,use_trans=T){
    if(!use_trans){
      return(x)
    }
    if(length(unique(x))==1){return(rep(1,length(x)))}
    return((x-min(x))/(max(x)-min(x))*(1-lower)+lower)
  }
  get_color<-function(color1="white",color2,x){
    myrgb<-function(matr){
      rgb(matr[,1],matr[,2],matr[,3],maxColorValue=255)
    }
    myrgb(colorRamp(c(color1,color2),bias=1)(x))
  }
  match.list2<-function(list1,list2){
    list3<-c(list1,list2)
    if(""%in%names(list3)){stop("Control_list should contain names for all arguments.")}
    list3[!duplicated(names(list3))]
  }
  my_eval<-function(code_start,code_list,code_end=")",env){
    exp2char<-function(exp){
      stringr::str_remove_all(paste0(capture.output(exp),collapse=""),"(^expression\\()|(\\)$)")
    }
    eval(parse(text=paste0(code_start,exp2char(as.expression(code_list)),code_end)),envir=env)
  }
  myround<-function(x,d=3){
    format(round(x,d),digits=d,nsmall=d)
  }
  get_legend<-function(cluster,color,color2,min,mid,max,draw_flag=rep(1,length(cluster)),
                       mar=c(2,2,2,2),x_extend=0.05,width=0.4,
                       border_col="gray",border_lwd=1,
                       scale_cex=0.7,scale_offset=0.2,
                       cluster_cex=0.8,c_x_adj=0,c_y_adj=-2,cluster_center=T,
                       point_cex=2,point_lwd=2,p_y_adj=0,
                       title=NULL,title_cex=1,t_x_adj=0,t_y_adj=-2){
    par(mar=mar)
    plot(y=c(0,100),x=c(0-x_extend,1+x_extend),xaxs='i',yaxs='i',type='n',ann=F,axes=F)
    n<-length(cluster)
    xl<-(0:(n-1))/n
    xr<-xl+1/n*width
    
    for(i in 1:n){
      if(draw_flag[i]==0){next}
      col<-color[[i]]
      rect(xleft=rep(xl[i],100),ybottom=seq(25,60-0.35,length.out=100),xright=rep(xr[i],100),ytop=seq(25,60-0.35,length.out=100)+0.35,col=col,border=col)
      
      rect(xl[i],25,xr[i],60,col="transparent",border=border_col,lwd=border_lwd)
      
      text(xr[i],25,labels=min[i],pos=4,cex=scale_cex,offset=scale_offset)
      text(xr[i],25/2+60/2,labels=mid[i],pos=4,cex=scale_cex,offset=scale_offset)
      text(xr[i],60,labels=max[i],pos=4,cex=scale_cex,offset=scale_offset)
      
      if(cluster_center){
        text(xl[i]/2+xr[i]/2+c_x_adj,71.25+c_y_adj,labels=cluster[i],cex=cluster_cex,font=2)
      }else{
        text(xl[i]+c_x_adj,71.25+c_y_adj,labels=cluster[i],pos=4,cex=cluster_cex,offset=0,font=2)
      }
      
      points(xl[i]/2+xr[i]/2,63.75+p_y_adj,pch=21,col=color2[i],cex=point_cex,lwd=point_lwd)
    }
    
    if(is.null(title)){
      if(sum(draw_flag==1)>1){
        title<-paste0("Estimated posterior\nprobabilities ","by\nclusters")
      }else{
        title<-"Estimated posterior\nprobabilities"
      }
    }
    text(0+t_x_adj,78.75+t_y_adj,labels=title,pos=4,cex=title_cex,offset=0,font=2)
  }
  my_abline<-function(b0,b1,range,...){
    corner<-par("usr")
    l<-max(corner[1],range[1])
    r<-min(corner[2],range[2])
    lines(c(l,r),c(l*b1+b0,r*b1+b0),...)
  }
  my_abline_ggplot<-function(b0,b1,range,obj_ggplot,...){
    p<-ggplot2::ggplot_build(obj_ggplot)
    corner<-p$layout$panel_params[[1]]$x.range
    l<-max(corner[1],range[1])
    r<-min(corner[2],range[2])
    #ggplot2::geom_line(aes(x=a,y=b),data=data.frame(a=c(l,r),b=c(l*b1+b0,r*b1+b0)),...)
    #ggplot2::annotate("segment",x=l,xend=r,y=l*b1+b0,yend=r*b1+b0,...)
    ggplot2::geom_segment(x=l,xend=r,y=l*b1+b0,yend=r*b1+b0,...)
  }
  
  fit<-est_out$maxlik$fit_final
  if(is.null(fit)){stop("est_out should be a final output of random-effects MR-PROLLIM.")}
  fit_data<-est_out$data
  par<-est_out$parameter
  
  crt_signk<-function(k_hat,p1_sp,p2_sp,model_u2,Egger){
    stopifnot(!is.null(Egger))
    stopifnot(!anyNA(k_hat))
    if(identical(Egger,"auto")){Egger<-T}
    j<-p1_sp==1&p2_sp==0&model_u2&Egger
    if(length(j)==0){j<-F}
    if(j){
      egger<-T
    }else{
      egger<-F
    }
    if(egger){
      sign_k1<-as.integer(2*(as.numeric(k_hat[,1]>=0)-0.5))
      sign_k2<-as.integer(2*(as.numeric(k_hat[,2]>=0)-0.5))
    }else{
      sign_k1<-sign_k2<-rep(1L,nrow(k_hat))
    }
    return(list(sign_k1=sign_k1,sign_k2=sign_k2))
  }
  signk<-crt_signk(fit_data$k_hat,par$p1_sp,par$p2_sp,par$model_u2,par$final_Egger_flag)
  
  if(par$est_type=="b"){
    if(is.null(fit_data$post_sample_k1)){
      #bi_p3_nome
      post_prob<-est_proc_bi_p3_f_plot(fit$estimate,m_matrix=fit_data$m_hat,sigma1_matr=fit_data$m_sigma,k_matrix=fit_data$k_hat,
                                       sign_k1=signk$sign_k1,sign_k2=signk$sign_k2,
                                       p1_sp=par$p1_sp,p2_sp=par$p2_sp,r_sp=par$r_sp,model_u2=par$model_u2,t_b1=par$t_b1)
    }else{
      #bi_p3
      post_sample_k1_exp1<-exp(fit_data$post_sample_k1)-1
      post_sample_k2_exp1<-exp(fit_data$post_sample_k2)-1
      
      prepare_p3.1_bi<-function(k_hat,post_sample_k1,post_sample_k2,post_sample_p,mkp_sigma_list){
        out1<-out2<-matrix(NA,nrow=nrow(post_sample_k1),ncol=ncol(post_sample_k1))
        out3<-list()
        for(i in 1:nrow(k_hat)){
          s<-mkp_sigma_list[[i]][1:2,3:5]%*%solve(mkp_sigma_list[[i]][3:5,3:5])
          x<-cbind(k_hat[i,1]-post_sample_k1[i,],k_hat[i,2]-post_sample_k2[i,],k_hat[i,3]-post_sample_p[i,])%*%t(s)
          out1[i,]<-x[,1]
          out2[i,]<-x[,2]
          out3[[i]]<-mkp_sigma_list[[i]][1:2,1:2]-s%*%mkp_sigma_list[[i]][3:5,1:2]
        }
        return(list(g1_matr=out1,g2_matr=out2,sigma_prime_list=out3))
      }
      data_opt_p3.1<-prepare_p3.1_bi(fit_data$k_hat,fit_data$post_sample_k1,fit_data$post_sample_k2,fit_data$post_sample_p,fit_data$mkp_sigma_list)
      
      post_prob<-est_proc_bi_p3.1_f_plot(fit$estimate,m_matrix=fit_data$m_hat,
                                         log_appr=par$log_appr,sign_k1=signk$sign_k1,sign_k2=signk$sign_k2,
                                         post_sample_k1=fit_data$post_sample_k1,post_sample_k2=fit_data$post_sample_k2,post_sample_p=fit_data$post_sample_p,
                                         post_sample_k1_exp1=post_sample_k1_exp1,post_sample_k2_exp1=post_sample_k2_exp1,
                                         g1_matr=data_opt_p3.1$g1_matr,g2_matr=data_opt_p3.1$g2_matr,sigma_prime_list=data_opt_p3.1$sigma_prime_list,
                                         p1_sp=par$p1_sp,p2_sp=par$p2_sp,r_sp=par$r_sp,model_u2=par$model_u2,t_b1=par$t_b1)
    }
    b1_td<-est_out$estimate$beta_norm[1]
    k1<-fit_data$k_hat[,1]
    k2<-fit_data$k_hat[,2]
    o1<-(exp(fit_data$k_hat[,1])-1)*fit_data$k_hat[,3]*b1_td+1
    o2<-(exp(fit_data$k_hat[,2])-1)*fit_data$k_hat[,3]*b1_td+1
    h1<-fit_data$m_hat[,1]-log(o1)
    h2<-fit_data$m_hat[,2]-log(o2)
    
    if(!is.null(fit_data$post_sample_k1)){
      vcov1<-lapply(fit_data$mkp_sigma_list,FUN=function(x){x[c(1,3,5),c(1,3,5)]})
      a1<-exp(fit_data$k_hat[,1])*fit_data$k_hat[,3]*b1_td/o1
      a2<-(exp(fit_data$k_hat[,1])-1)*b1_td/o1
      se1<-NA
      for(i in 1:length(a1)){
        o<-c(1,-a1[i],-a2[i])
        se1[i]<-sqrt(t(o)%*%vcov1[[i]]%*%o)
      }
      
      vcov2<-lapply(fit_data$mkp_sigma_list,FUN=function(x){x[c(2,4,5),c(2,4,5)]})
      a1<-exp(fit_data$k_hat[,2])*fit_data$k_hat[,3]*b1_td/o2
      a2<-(exp(fit_data$k_hat[,2])-1)*b1_td/o2
      se2<-NA
      for(i in 1:length(a1)){
        o<-c(1,-a1[i],-a2[i])
        se2[i]<-sqrt(t(o)%*%vcov2[[i]]%*%o)
      }
    }else{
      se1<-sqrt(fit_data$m_sigma[,1])
      se2<-sqrt(fit_data$m_sigma[,3])
    }
    se_h<-cbind(se1,se2)
    se_k<-matrix(unlist(lapply(fit_data$k_sigma,FUN=function(x){sqrt(c(x[1,1],x[2,2]))})),ncol=2,byrow=T)
  }
  
  if(par$est_type=="c"){
    if(is.null(fit_data$post_sample_k1)){
      #cont_p3_nome
      post_prob<-est_proc_cont_p3_f_plot(fit$estimate,m_matrix=fit_data$m_hat,sigma1_matr=fit_data$m_sigma,k_matrix=fit_data$k_hat,
                                         sign_k1=signk$sign_k1,sign_k2=signk$sign_k2,
                                         p1_sp=par$p1_sp,p2_sp=par$p2_sp,r_sp=par$r_sp,model_u2=par$model_u2)
    }else{
      #cont_p3
      prepare_p3.1<-function(k_hat,post_sample_k1,post_sample_k2,mk_sigma_list){
        out1<-out2<-matrix(NA,nrow=nrow(post_sample_k1),ncol=ncol(post_sample_k1))
        out3<-list()
        for(i in 1:nrow(k_hat)){
          s<-mk_sigma_list[[i]][1:2,3:4]%*%solve(mk_sigma_list[[i]][3:4,3:4])
          x<-cbind(k_hat[i,1]-post_sample_k1[i,],k_hat[i,2]-post_sample_k2[i,])%*%t(s)
          out1[i,]<-x[,1]
          out2[i,]<-x[,2]
          out3[[i]]<-mk_sigma_list[[i]][1:2,1:2]-s%*%mk_sigma_list[[i]][3:4,1:2]
        }
        return(list(f1_matr=out1,f2_matr=out2,sigma_prime_list=out3))
      }
      data_opt_p3.1<-prepare_p3.1(fit_data$k_hat,fit_data$post_sample_k1,fit_data$post_sample_k2,fit_data$mk_sigma_list)
      
      post_prob<-est_proc_cont_p3.1_f_plot(fit$estimate,m_matrix=fit_data$m_hat,sign_k1=signk$sign_k1,sign_k2=signk$sign_k2,
                                           f1_matr=data_opt_p3.1$f1_matr,f2_matr=data_opt_p3.1$f2_matr,sigma_prime_list=data_opt_p3.1$sigma_prime_list, 
                                           post_sample_k1=fit_data$post_sample_k1,post_sample_k2=fit_data$post_sample_k2,
                                           p1_sp=par$p1_sp,p2_sp=par$p2_sp,r_sp=par$r_sp,model_u2=par$model_u2)
    }
    b1<-est_out$estimate$beta_norm[1]
    k1<-fit_data$k_hat[,1]
    k2<-fit_data$k_hat[,2]
    o1<-fit_data$k_hat[,1]*b1
    o2<-fit_data$k_hat[,2]*b1
    h1<-fit_data$m_hat[,1]-o1
    h2<-fit_data$m_hat[,2]-o2
    
    if(!is.null(fit_data$post_sample_k1)){
      vcov1<-lapply(fit_data$mk_sigma_list,FUN=function(x){x[c(1,3),c(1,3)]})
      se1<-NA
      for(i in 1:length(vcov1)){
        o<-c(1,-b1)
        se1[i]<-sqrt(t(o)%*%vcov1[[i]]%*%o)
      }
      
      vcov2<-lapply(fit_data$mk_sigma_list,FUN=function(x){x[c(2,4),c(2,4)]})
      se2<-NA
      for(i in 1:length(vcov2)){
        o<-c(1,-b1)
        se2[i]<-sqrt(t(o)%*%vcov2[[i]]%*%o)
      }
    }else{
      se1<-sqrt(fit_data$m_sigma[,1])
      se2<-sqrt(fit_data$m_sigma[,3])
    }
    se_h<-cbind(se1,se2)
    se_k<-sqrt(fit_data$k_sigma[,c(1,3)])
  }
  
  par_p2<-est_out$estimate$beta_norm[3]
  par_u1<-est_out$estimate$beta_norm[4]
  par_u2<-est_out$estimate$beta_norm[7]
  
  if(cluster_type==2){
    post_prob<-cbind(post_prob[,1],post_prob[,2]+post_prob[,3])
  }
  
  rownames(post_prob)<-rownames(fit_data$m_hat)
  
  if(cluster_type==2){
    cluster<-c("No HP","HP")[apply(post_prob,1,FUN=which.max)]
  }else{
    cluster<-c("No HP","Uncor HP","Cor HP")[apply(post_prob,1,FUN=which.max)]
  }
  
  post_prob2<-apply(post_prob,1,FUN=max)
  
  col1<-col2<-rep(NA,length(cluster))
  
  if(cluster_type==2){
    col1[cluster=="No HP"]<-cluster_color[1]
    col1[cluster=="HP"]<-cluster_color[2]
    c_u<-unique(cluster)
    for(i in 1:length(c_u)){
      col2[cluster==c_u[i]]<-get_color(bg_color,
                                       cluster_color[which(c("No HP","HP")==c_u[i])],
                                       trans(post_prob2[cluster==c_u[i]],trans_lower,use_trans))
    }
  }else{
    col1[cluster=="No HP"]<-cluster_color[1]
    col1[cluster=="Uncor HP"]<-cluster_color[2]
    col1[cluster=="Cor HP"]<-cluster_color[3]
    c_u<-unique(cluster)
    for(i in 1:length(c_u)){
      col2[cluster==c_u[i]]<-get_color(bg_color,
                                       cluster_color[which(c("No HP","Uncor HP","Cor HP")==c_u[i])],
                                       trans(post_prob2[cluster==c_u[i]],trans_lower,use_trans))
    }
  }
  if(!is.null(emph_loc)){
    stopifnot(is.numeric(emph_loc))
    col2[emph_loc]<-col1[emph_loc]<-emph_col
  }
  
  old<-par(no.readonly=T)
  layout(matrix(c(1,1,2,2,3),ncol=5))
  
  control_plot1_org<-list(xlab=expression("Estimate of"~italic(k)[1]),ylab=expression("Estimate of"~italic(h)[1]),main="")
  control_plot1<-match.list2(control_plot1,control_plot1_org)
  control_lines_org<-list(lty=1,lwd=1,col="gray")
  control_lines<-match.list2(control_lines,control_lines_org)
  control_points_org<-list(pch=21,cex=1.2,lwd=1.6)
  control_points<-match.list2(control_points,control_points_org)
  control_abline_org<-list(lty=2,lwd=1)
  control_abline<-match.list2(control_abline,control_abline_org)
  
  if(is.null(reorder)){
    rd<-1:length(k1)
  }else{
    stopifnot(length(reorder)==length(k1))
    stopifnot(is.numeric(reorder))
    rd<-reorder
  }
  
  my_eval("plot(x=k1[rd],y=h1[rd],type=\"n\",",control_plot1,")",environment())
  corner<-par("usr")
  mtext(sub_fig_label[1],side=2,outer=F,line=sub_fig_line,las=1,at=(corner[4]-corner[3])*sub_fig_k+corner[4],cex=sub_fig_cex)
  if(!is.null(expr1)){
    eval(expr1,envir=environment())
  }
  if(!is.null(ci_cover)){
    q<-abs(qnorm((1-ci_cover)/2))
    upper<-h1+q*se_h[,1]
    lower<-h1-q*se_h[,1]
    for(i in 1:length(h1)){
      my_eval("lines(x=c(k1[rd[i]],k1[rd[i]]),y=c(lower[rd[i]],upper[rd[i]]),",control_lines,")",environment())
    }
    
    upper<-k1+q*se_k[,1]
    lower<-k1-q*se_k[,1]
    for(i in 1:length(k1)){
      my_eval("lines(y=c(h1[rd[i]],h1[rd[i]]),x=c(lower[rd[i]],upper[rd[i]]),",control_lines,")",environment())
    }
  }
  
  if(par_p2!=1&cluster_type==3){
    if(par$final_Egger_flag){
      my_eval("my_abline(b1=0,b0=par_u2,",c(list(range=c(0,Inf)),control_abline,list(col=cluster_color[2])),")",environment())
      my_eval("my_abline(b1=0,b0=-par_u2,",c(list(range=c(-Inf,0)),control_abline,list(col=cluster_color[2])),")",environment())
    }else{
      my_eval("my_abline(b1=0,b0=par_u2,",c(list(range=c(-Inf,Inf)),control_abline,list(col=cluster_color[2])),")",environment())
    }
  }
  if(par_u1!=0&par_p2!=0&cluster_type==3){
    my_eval("abline(b=par_u1,a=par_u2,",c(control_abline,list(col=cluster_color[3])),")",environment())
  }
  my_eval("points(x=k1[rd],y=h1[rd],bg=col2[rd],col=col1[rd],",control_points,")",environment())
  
  control_plot2_org<-list(xlab=expression("Estimate of"~italic(k)[2]),ylab=expression("Estimate of"~italic(h)[2]),main="")
  control_plot2<-match.list2(control_plot2,control_plot2_org)
  
  my_eval("plot(x=k2[rd],y=h2[rd],type=\"n\",",control_plot2,")",environment())
  corner<-par("usr")
  mtext(sub_fig_label[2],side=2,outer=F,line=sub_fig_line,las=1,at=(corner[4]-corner[3])*sub_fig_k+corner[4],cex=sub_fig_cex)
  if(!is.null(expr2)){
    eval(expr2,envir=environment())
  }
  if(!is.null(ci_cover)){
    q<-abs(qnorm((1-ci_cover)/2))
    upper<-h2+q*se_h[,2]
    lower<-h2-q*se_h[,2]
    for(i in 1:length(h2)){
      my_eval("lines(x=c(k2[rd[i]],k2[rd[i]]),y=c(lower[rd[i]],upper[rd[i]]),",control_lines,")",environment())
    }
    
    upper<-k2+q*se_k[,2]
    lower<-k2-q*se_k[,2]
    for(i in 1:length(k2)){
      my_eval("lines(y=c(h2[rd[i]],h2[rd[i]]),x=c(lower[rd[i]],upper[rd[i]]),",control_lines,")",environment())
    }
  }
  
  if(par_p2!=1&cluster_type==3){
    if(par$final_Egger_flag){
      my_eval("my_abline(b1=0,b0=2*par_u2,",c(list(range=c(0,Inf)),control_abline,list(col=cluster_color[2])),")",environment())
      my_eval("my_abline(b1=0,b0=-2*par_u2,",c(list(range=c(-Inf,0)),control_abline,list(col=cluster_color[2])),")",environment())
    }else{
      my_eval("my_abline(b1=0,b0=2*par_u2,",c(list(range=c(-Inf,Inf)),control_abline,list(col=cluster_color[2])),")",environment())
    }
  }
  if(par_u1!=0&par_p2!=0&cluster_type==3){
    my_eval("abline(b=par_u1,a=2*par_u2,",c(control_abline,list(col=cluster_color[3])),")",environment())
  }
  my_eval("points(x=k2[rd],y=h2[rd],bg=col2[rd],col=col1[rd],",control_points,")",environment())
  
  control_legend_org<-list(draw_legend=T,mar=c(2,0,2,2),x_extend=0.05,width=0.4,
                           border_col="gray",border_lwd=1,
                           scale_cex=0.8,scale_offset=0.2,
                           cluster_cex=0.9,c_x_adj=0,c_y_adj=-2,cluster_center=T,
                           point_cex=2,point_lwd=2,p_y_adj=0,
                           title=NULL,title_cex=1,t_x_adj=0,t_y_adj=-1)
  c_l<-match.list(control_legend,control_legend_org)
  
  if(c_l$draw_legend){
    get_legend_data1<-function(x,cluster_type){
      rr<-function(x){
        out<-rep("null",3)
        loc<-which(x!="null")
        if(length(loc)==1){
          out[2]<-x[loc]
        }
        if(length(loc)==2){
          out[1]<-x[loc[1]]
          out[3]<-x[loc[2]]
        }
        if(length(loc)==3){
          out<-x
        }
        return(out)
      }
      out<-rep("null",3)
      if(cluster_type==2){
        if("No HP"%in%x){out[1]<-"No HP"}
        if("HP"%in%x){out[2]<-"HP"}
      }else{
        if("No HP"%in%x){out[1]<-"No HP"}
        if("Uncor HP"%in%x){out[2]<-"Uncor HP"}
        if("Cor HP"%in%x){out[3]<-"Cor HP"}
      }
      rr(out)
    }
    get_legend_data2<-function(x,cluster,post_prob2,cluster_color,bg_color,cluster_type,trans_lower,use_trans){
      max_data<-min_data<-mid_data<-rep(0,3)
      color_l<-list()
      color2_l<-rep("black",3)
      x_u<-x[x!="null"]
      draw_flag<-as.numeric(x!="null")
      for(i in 1:length(x_u)){
        loc<-which(x==x_u[i])
        val<-post_prob2[cluster==x_u[i]]
        max_data[loc]<-max(val)
        min_data[loc]<-min(val)
        mid_data[loc]<-(max(val)+min(val))/2
        
        val2<-trans(val,trans_lower,use_trans)
        if(cluster_type==2){
          color_l[[loc]]<-get_color(bg_color,
                                    cluster_color[which(c("No HP","HP")==x_u[i])],
                                    seq(min(val2),max(val2),length.out=100))
          color2_l[loc]<-cluster_color[which(c("No HP","HP")==x_u[i])]
        }else{
          color_l[[loc]]<-get_color(bg_color,
                                    cluster_color[which(c("No HP","Uncor HP","Cor HP")==x_u[i])],
                                    seq(min(val2),max(val2),length.out=100))
          color2_l[loc]<-cluster_color[which(c("No HP","Uncor HP","Cor HP")==x_u[i])]
        }
      }
      return(list(color=color_l,color2=color2_l,min=min_data,mid=mid_data,max=max_data,draw_flag=draw_flag))
    }
    
    cluster_l<-get_legend_data1(cluster,cluster_type)
    data_l<-get_legend_data2(cluster_l,cluster,post_prob2,cluster_color,bg_color,cluster_type,trans_lower,use_trans)
    
    get_legend(cluster_l,color=data_l$color,color2=data_l$color2,
               min=myround(data_l$min,2),mid=myround(data_l$mid,2),max=myround(data_l$max,2),
               draw_flag=data_l$draw_flag,c_l$mar,c_l$x_extend,c_l$width,
               c_l$border_col,c_l$border_lwd,
               c_l$scale_cex,c_l$scale_offset,
               c_l$cluster_cex,c_l$c_x_adj,c_l$c_y_adj,c_l$cluster_center,
               c_l$point_cex,c_l$point_lwd,c_l$p_y_adj,
               c_l$title,c_l$title_cex,c_l$t_x_adj,c_l$t_y_adj)
  }
  
  par(old)
  
  if(interactive){
    if(!"ggplot2"%in%(.packages())){
      library(ggplot2)
    }
    if(!"plotly"%in%(.packages())){
      library(plotly)
    }
    label<-apply(post_prob,1,FUN=function(x){paste0("(",paste0(myround(x,2),collapse=", "),")")})
    label<-paste(rownames(post_prob),label)
    
    p01<-ggplot(data=data.frame(x=k1[rd],y=h1[rd],label=label[rd]),mapping=aes(x=x,y=y,key=label))
    
    if(!is.null(ci_cover)){
      q<-abs(qnorm((1-ci_cover)/2))
      upper<-h1+q*se_h[,1]
      lower<-h1-q*se_h[,1]
      p1_line<-geom_errorbar(mapping=aes(ymin=lower,ymax=upper),
                             data=data.frame(x=k1[rd],y=h1[rd],upper=upper[rd],lower=lower[rd]),
                             colour="gray")
      
      upper<-k1+q*se_k[,1]
      lower<-k1-q*se_k[,1]
      p12_line<-geom_errorbar(mapping=aes(xmin=lower,xmax=upper),
                              data=data.frame(x=k1[rd],y=h1[rd],upper=upper[rd],lower=lower[rd]),
                              colour="gray")
    }else{
      p1_line<-NULL
      p12_line<-NULL
    }
    
    p01<-p01+p1_line+p12_line
    
    if(par_p2!=1&cluster_type==3){
      if(par$final_Egger_flag){
        p01<-p01+my_abline_ggplot(b0=par_u2,b1=0,range=c(0,Inf),obj_ggplot=p01,lty=2,color=cluster_color[2])+
          my_abline_ggplot(b0=-par_u2,b1=0,range=c(-Inf,0),obj_ggplot=p01,lty=2,color=cluster_color[2])
      }else{
        p01<-p01+my_abline_ggplot(b0=par_u2,b1=0,range=c(-Inf,Inf),obj_ggplot=p01,lty=2,color=cluster_color[2])
      }
    }
    
    if(par_u1!=0&par_p2!=0&cluster_type==3){
      p01<-p01+geom_abline(slope=par_u1,intercept=par_u2,lty=2,color=cluster_color[3])
    }
    
    p1_point<-geom_point(size=2,stroke=0.6,shape=21,colour=col1[rd],fill=col2[rd])
    p1<-p01+p1_point+
      xlab("Estimate of k_1")+
      ylab("Estimate of h_1")+
      theme_bw()+theme(legend.position='none')
    
    p02<-ggplot(data=data.frame(x=k2[rd],y=h2[rd],label=label[rd]),mapping=aes(x=x,y=y,key=label))
    
    if(!is.null(ci_cover)){
      q<-abs(qnorm((1-ci_cover)/2))
      upper<-h2+q*se_h[,2]
      lower<-h2-q*se_h[,2]
      p2_line<-geom_errorbar(mapping=aes(ymin=lower,ymax=upper),
                             data=data.frame(x=k2[rd],y=h2[rd],upper=upper[rd],lower=lower[rd]),
                             colour="gray")
      
      upper<-k2+q*se_k[,2]
      lower<-k2-q*se_k[,2]
      p22_line<-geom_errorbar(mapping=aes(xmin=lower,xmax=upper),
                              data=data.frame(x=k2[rd],y=h2[rd],upper=upper[rd],lower=lower[rd]),
                              colour="gray")
      
    }else{
      p2_line<-NULL
      p22_line<-NULL
    }
    
    p02<-p02+p2_line+p22_line
    
    if(par_p2!=1&cluster_type==3){
      if(par$final_Egger_flag){
        p02<-p02+my_abline_ggplot(b0=2*par_u2,b1=0,range=c(0,Inf),obj_ggplot=p02,lty=2,color=cluster_color[2])+
          my_abline_ggplot(b0=-2*par_u2,b1=0,range=c(-Inf,0),obj_ggplot=p02,lty=2,color=cluster_color[2])
      }else{
        p02<-p02+my_abline_ggplot(b0=2*par_u2,b1=0,range=c(-Inf,Inf),obj_ggplot=p02,lty=2,color=cluster_color[2])
      }
    }
    
    if(par_u1!=0&par_p2!=0&cluster_type==3){
      p02<-p02+geom_abline(slope=par_u1,intercept=2*par_u2,lty=2,color=cluster_color[3])
    }
    
    p2_point<-geom_point(size=2,stroke=0.6,shape=21,colour=col1[rd],fill=col2[rd])
    p2<-p02+p2_point+
      xlab("Estimate of k_2")+
      ylab("Estimate of h_2")+
      theme_bw()+theme(legend.position='none')
    p1<-ggplotly(p1)
    p2<-ggplotly(p2)
    print(subplot(p1,p2,titleX=T,titleY=T,margin=0.05))
  }
  
  if(!is.null(obj_name)){
    names(cluster)<-rownames(post_prob)
    assign(obj_name,list(post_prob=post_prob,cluster=cluster,
                         data1=list(k1=k1,h1=h1,k1_se=se_k[,1],h1_se=se_h[,1]),
                         data2=list(k2=k2,h2=h2,k2_se=se_k[,2],h2_se=se_h[,2])),envir=.GlobalEnv)
  }
  invisible()
}