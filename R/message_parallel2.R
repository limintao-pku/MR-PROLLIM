message_parallel2<-function(begin,x,time=F,t0=NULL){
  if(Sys.info()[1]=="Windows"){
    if(time){
      cat(paste(begin,x,"|",capture.output(Sys.time()-t0),"\n"))
    }else{
      cat(paste(begin,x,"\n"))
    }
  }else{
    if(time){
      system( sprintf('echo "%s"',paste(begin,x,"|",capture.output(Sys.time()-t0))) )
    }else{
      system( sprintf('echo "%s"',paste(begin,x)) )
    }
  }
}
