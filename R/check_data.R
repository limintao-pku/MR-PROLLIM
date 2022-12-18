check_data<-function(x,y,g,c=NULL,c_inherit,twosample_data=NULL,type=c("b","c"),u_limit=10,cd_g_code=T){
  type<-match.arg(type)
  unique2<-function(x){
    x<-na.omit(x)
    unique(x)
  }
  table2<-function(x){
    y<-table(x)
    return(y/sum(y))
  }
  recode_g_f<-function(g,force=F){
    if(force){
      g1<-g
      g1[g==2]<-0L
      g1[g==0]<-2L
      return(g1)
    }
    m<-mean(g,na.rm=T)/2
    if(m>0.5){
      g1<-g
      g1[g==2]<-0L
      g1[g==0]<-2L
      return(list(g1,T))
    }else{
      return(list(g,F))
    }
  }
  
  if(!is.null(twosample_data)){
    if((!is.list(twosample_data))|sum(c("g","c")%in%names(twosample_data))!=2){
      stop("twosample_data should be a named list containing 'g' matrix and 'c' list(matrices). 'c' can be NULL.")
    }
  }
  
  #check c: NAs, unique value, format, length, etc.
  check_c<-function(g,c,c_inherit,u_limit){
    if(!is.null(c)){
      if((!is.list(c))|(!is.matrix(c[[1]]))){stop("c should be a list containing c matirx")}
      stopifnot(nrow(g)==nrow(c[[1]]))
      
      mycheck2<-function(m,u_limit){
        if(anyNA(m)){stop("NAs in c_list detected")}
        j<-F
        for(i in 1:ncol(m)){
          x<-length(unique(m[,i]))
          if(x<=1){stop("constant c detected in c_list")}
          if(x>=3&x<=u_limit){j<-T}
        }
        if(j){warning("There seems to be a categorical c in c_list that has not been transformed to dummy variables.")}
        NULL
      }
      mywarn<-F
      withCallingHandlers({lapply(c,FUN=mycheck2,u_limit=10)},
                          warning=function(w){mywarn<<-T;invokeRestart("muffleWarning")})
      if(mywarn){message("There seems to be a categorical c in c_list that has not been transformed to dummy variables.")}
      
      if(c_inherit){
        if(length(c)!=1){stop("c_list should contain only 1 matirx if c_inherit = T")}
      }else{
        if(!identical(length(c),ncol(g))){stop("length(c) and ncol(g) should be equal if c_inherit = F")}
      }
    }
    return(NULL)
  }
  check_c(g=g,c=c,c_inherit=c_inherit,u_limit)
  check_c(g=twosample_data$g,c=twosample_data$c,c_inherit=c_inherit,u_limit)
  
  #check NAs
  if(anyNA(x)|anyNA(y)){stop("NAs detected in x or y")}
  
  #check y
  if(!identical(as.numeric(sort(unique(y))),c(0,1))|is.factor(y)){stop("y should be binary and coded with numeric/integer 0 and 1.")}
  
  #check x
  if(type=="b"){
    if(!identical(as.numeric(sort(unique(x))),c(0,1))|is.factor(x)){stop("x should be binary and coded with numeric/integer 0 and 1.")}
  }
  
  #check g
  check_g<-function(g,g2,u_limit,cd_g_code){
    if(is.null(g2)){
      stopifnot(is.matrix(g))
      if(is.null(colnames(g))){stop("Column names for g_matrix are required")}
      if(length(unique(colnames(g)))!=ncol(g)){stop("Column names for g_matrix are not unique")}
      
      m<-NULL
      for(i in 1:ncol(g)){
        x<-sort(unique2(g[,i]))
        if(!any(x==0)){stop("One SNP is not properly coded")}
        if(length(x)>=u_limit){stop("One SNP has too many unique values")}
        if(cd_g_code){
          if(identical(as.numeric(x),c(0,1,2))){
            y<-recode_g_f(g[,i])
            if(y[[2]]){m<-c(m,i);g[,i]<-y[[1]]}
          }
        }
      }
      if(length(m)>0){
        message(length(m)," SNPs are recoded. Corresponding SNP names are attached with '_recoded'")
        colnames(g)[m]<-paste0(colnames(g)[m],"_recoded")
      }
      return(list(g))
    }else{
      stopifnot(is.matrix(g));stopifnot(is.matrix(g2))
      if(is.null(colnames(g))){stop("Column names for g_matrix are required")}
      if(length(unique(colnames(g)))!=ncol(g)){stop("Column names for g_matrix are not unique")}
      if(!identical(colnames(g),colnames(g2))){stop("Column names for g_matrix and tsd_g_matrix are not identical")}
      
      m<-NULL
      for(i in 1:ncol(g)){
        x<-sort(unique2(g[,i]))
        x2<-sort(unique2(g2[,i]))
        if(!identical(as.numeric(x),as.numeric(x2))){stop("SNPs in g_matrix and tsd_g_matrix are not identically coded.")}
        if(!any(x==0)){stop("One SNP is not properly coded")}
        if(length(x)>=u_limit){stop("One SNP has too many unique values")}
        if(any(abs(table2(g[,i])-table2(g2[,i]))>0.1)){message(colnames(g)[i]," exhibits substantially different genotype frequencies in g_matrix and tsd_g_matrix.")}
        if(cd_g_code){
          if(identical(as.numeric(x),c(0,1,2))){
            y<-recode_g_f(g[,i])
            if(y[[2]]){m<-c(m,i);g[,i]<-y[[1]];g2[,i]<-recode_g_f(g2[,i],T)}
          }
        }
      }
      if(length(m)>0){
        message(length(m),"SNPs are recoded. Corresponding SNP names are attached with '_recoded'")
        colnames(g)[m]<-colnames(g2)[m]<-paste0(colnames(g)[m],"_recoded")
      }
      return(list(g,g2))
    }
  }
  z<-check_g(g,twosample_data$g,u_limit,cd_g_code)
  g<-z[[1]]
  if(!is.null(twosample_data)){
    twosample_data$g<-z[[2]]
  }  
  
  #check length
  stopifnot(nrow(g)==length(y))
  if(is.null(twosample_data)){stopifnot(nrow(g)==length(x))}else{stopifnot(nrow(twosample_data$g)==length(x))}
  
  return(list(g,twosample_data))
}

