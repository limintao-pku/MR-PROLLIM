message_parallel<-function(x){
  if(Sys.info()[1]=="Windows"){
    cat(paste0(paste0(paste("Messages/Warnings:",x),collapse="\n"),"\n"))
  }else{
    system(sprintf('echo "%s"',paste0(paste("Messages/Warnings:",x),collapse="\n")))
  }
}
