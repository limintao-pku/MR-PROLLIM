str_remove_all<-function(string,pattern,fixed=F){
  gsub(pattern=pattern,replacement="",x=string,ignore.case=F,perl=F,
       fixed=fixed,useBytes=F)
}