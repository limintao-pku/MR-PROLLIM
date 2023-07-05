nlminb2nlm<-function(nlminb_out){
  if(!is.list(nlminb_out)){return(nlminb_out)}
  list(minimum=nlminb_out$objective,estimate=nlminb_out$par,gradient="nlminb_out",
       code=nlminb_out$convergence,iterations=nlminb_out$evaluations,message=nlminb_out$message)
}

