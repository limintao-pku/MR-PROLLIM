match.list<-function(list1,list2,vector=F){
  if(is.null(list1)){return(list2)}
  if(!vector){
    if(!is.list(list1)){stop("control_list should be a list.")}
    if(is.null(names(list1))){stop("control_list should be a named list.")}
  }
  na2<-names(list2)
  na1<-names(list1)
  for(i in 1:length(list1)){
    loc<-which(na2==na1[i])
    if(length(loc)!=1){stop(paste(na1[i],"is an unknown parameter."))}
    if(is.list(list1)){
      if(is.null(list1[[i]])){list2[loc]<-list(NULL)}else{list2[[loc]]<-list1[[i]]}
    }else{list2[loc]<-list1[i]}
  }
  return(list2)
}

