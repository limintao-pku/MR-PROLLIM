str_detect<-function(string,pattern,fixed=F){
  grepl(pattern=pattern,x=string,ignore.case=F,perl=F,
        fixed=fixed,useBytes=F)
}