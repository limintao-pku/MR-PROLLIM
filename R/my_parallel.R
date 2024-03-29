my_parallel<-function(X,FUN,...,mc.cores=1,PSOCK=F,dt=T,print_message=T,
                      cl=NULL,stopcl=T,exec_base_func=NULL,export_parent=F,export_parent_func=F,
                      add_obj_list=NULL,outfile=NULL,seed=NULL){
  cores<-mc.cores
  stopifnot(cores>=1)
  stopifnot(length(X)<=cores)
  if(length(X)<cores){cores<-length(X)}
  if(cores>parallel::detectCores()){message("cores > parallel::detectCores()")}
  
  if(cores==1){
    out<-parallel::mclapply(X=X,FUN=FUN,...,mc.cores=cores)
    return(out)
  }
  
  get_function<-function(env){
    name<-ls(envir=env)
    j<-NA
    for(i in 1:length(name)){
      j[i]<-"function"%in%eval(parse(text=paste0("class(",name[i],")")),envir=env)
    }
    list(var=name[j],env=env)
  }
  if(Sys.info()[1]=="Windows"|PSOCK){
    if(Sys.info()[1]=="Windows"){
      tmpdir<-tempdir()
      if(!is.null(outfile)){
        temp1<-outfile
      }else{
        temp1<-tempfile(tmpdir=tmpdir,fileext=".txt")
      }
      if(print_message){
        temp2<-tempfile(tmpdir=tmpdir,fileext=".ps1")
        temp3<-tempfile(tmpdir=tmpdir,fileext=".bat")
        writeLines(paste0("Get-Content ",temp1," -wait"),con=temp2)
        writeLines(paste0("powershell -ExecutionPolicy RemoteSigned -File ",temp2,"\npause"),con=temp3)
      }
    }else{
      if(!is.null(outfile)){
        if(print_message){
          temp1<-""
          message("Print_message = T, and outfile will not be used.")
        }else{
          temp1<-outfile
        }
      }else{
        if(print_message){
          temp1<-""
        }else{
          temp1<-tempfile()
        }
      }
    }
    
    if(is.null(cl)){
      cl<-parallel::makePSOCKcluster(cores,outfile=temp1)
      if(!is.null(seed)){
        #to be consistent with mclapply
        myclusterSetRNGStream<-function(cl=NULL,iseed=NULL){
          cl <- parallel:::defaultCluster(cl)
          oldseed <- if (exists(".Random.seed", envir = .GlobalEnv, 
                                inherits = FALSE)) 
            get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
          else NULL
          RNGkind("L'Ecuyer-CMRG")
          if (!is.null(iseed)) 
            set.seed(iseed)
          nc <- length(cl)
          seeds <- vector("list", nc)
          seeds[[1L]] <- parallel::nextRNGStream(.Random.seed)
          for (i in seq_len(nc - 1L)) seeds[[i + 1L]] <- parallel::nextRNGStream(seeds[[i]])
          if (!is.null(oldseed)) 
            assign(".Random.seed", oldseed, envir = .GlobalEnv)
          else rm(.Random.seed, envir = .GlobalEnv)
          for (i in seq_along(cl)) {
            expr <- substitute(assign(".Random.seed", seed, envir = .GlobalEnv), 
                               list(seed = seeds[[i]]))
            parallel:::sendCall(cl[[i]], eval, list(expr))
          }
          parallel:::checkForRemoteErrors(lapply(cl, parallel:::recvResult))
          invisible()
        }
        myclusterSetRNGStream(cl,seed)
      }
      j<-(!is.null(exec_base_func))|export_parent|export_parent_func|(!is.null(add_obj_list))
      if(j){
        if(F){cat("Start exporting objects to clusters.\r\n")}
      }
      if(!is.null(exec_base_func)){
        #if(dt){cat("Start exec_base_func\r\n")}
        invisible(parallel::parLapply(cl=cl,X=1:cores,fun=exec_base_func))
        #if(dt){cat("Exec_base_func done.\r\n")}
      }
      if(export_parent){
        #if(dt){cat("Start export_parent\r\n")}
        parallel::clusterExport(cl,ls(envir=parent.frame()),envir=parent.frame())
        #if(dt){cat("Export_parent done.\r\n")}
      }
      if(export_parent_func){
        par_func<-get_function(parent.frame())
        #if(dt){cat("Start export_parent_func\r\n")}
        parallel::clusterExport(cl,par_func$var,envir=par_func$env)
        #if(dt){cat("Export_parent_func done.\r\n")}
      }
      if(!is.null(add_obj_list)){
        #if(dt){cat("Start add_obj_list\r\n")}
        parallel::clusterExport(cl,varlist=add_obj_list$var,envir=add_obj_list$env)
        #if(dt){cat("Add_obj_list done.\r\n")}
      }
      if(j){
        if(F){cat("Exporting done.\r\n")}
      }
    }
    
    if(print_message&Sys.info()[1]=="Windows"){
      shell.exec(temp3)
    }
    out<-tryCatch({parallel::parLapply(cl=cl,X=X,fun=FUN,...)},error=function(e){warning("Errors occured in parLapply");e})
    if(stopcl){
      parallel::stopCluster(cl)
      Sys.sleep(0.3)
    }
    if(!print_message){
      parlog<-suppressWarnings(readLines(con=temp1))
      if(length(parlog)>cores){
        message("There seems to be printed information (e.g., warnings) in the logfile.")
        message("Logfile: ",temp1)
      }
    }
  }else{
    if(!is.null(seed)){
      myseed0<-tryCatch({get(".Random.seed",envir=.GlobalEnv,inherits=FALSE)},error=function(e){"e"})
      RNGkind("L'Ecuyer-CMRG")
      set.seed(seed)
    }
    out<-parallel::mclapply(X=X,FUN=FUN,...,mc.cores=cores)
    if(!is.null(seed)){
      if(!identical(myseed0,"e")){
        assign(".Random.seed",myseed0,envir=.GlobalEnv)
      }else{
        rm(.Random.seed,envir=.GlobalEnv)
      }
    }
  }
  return(out)
}
