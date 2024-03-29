nlminb2<-function(f,p,...,gradient=F,scale=1,control=list(),lower=-Inf,upper=Inf) 
{
  objective<-f
  start<-p
  par <- setNames(as.double(start), names(start))
  n <- length(par)
  iv <- integer(78 + 3 * n)
  v <- double(130 + (n * (n + 27))/2)
  .Call(stats:::C_port_ivset, 2, iv, v)
  if (length(control)) {
    nms <- names(control)
    if (!is.list(control) || is.null(nms)) 
      stop("'control' argument must be a named list")
    pos <- pmatch(nms, names(stats:::port_cpos))
    if (any(nap <- is.na(pos))) {
      warning(sprintf(ngettext(length(nap), "unrecognized control element named %s ignored", 
                               "unrecognized control elements named %s ignored"), 
                      paste(sQuote(nms[nap]), collapse = ", ")), domain = NA)
      pos <- pos[!nap]
      control <- control[!nap]
    }
    ivpars <- pos <= 4
    vpars <- !ivpars
    if (any(ivpars)) 
      iv[stats:::port_cpos[pos[ivpars]]] <- as.integer(unlist(control[ivpars]))
    if (any(vpars)) 
      v[stats:::port_cpos[pos[vpars]]] <- as.double(unlist(control[vpars]))
  }
  if(!gradient){
    obj <- quote(objective(.par, ...))
  }else{
    obj <- quote(objective(.par, ..., grad=F))
  }
  rho <- new.env(parent = environment())
  assign(".par", par, envir = rho)
  grad <- hess <- low <- upp <- NULL
  if (gradient) {
    grad <- quote(attr(objective(.par, ..., grad=T),"gradient"))
  }
  if (any(lower != -Inf) || any(upper != Inf)) {
    low <- rep_len(as.double(lower), length(par))
    upp <- rep_len(as.double(upper), length(par))
  }
  else low <- upp <- numeric()
  .Call(stats:::C_port_nlminb, obj, grad, hess, rho, low, upp, d = rep_len(as.double(scale), 
                                                                   length(par)), iv, v)
  iv1 <- iv[1L]
  list(par = get(".par", envir = rho), objective = v[10L], 
       convergence = (if (iv1 %in% 3L:6L) 0L else 1L), iterations = iv[31L], 
       evaluations = c(`function` = iv[6L], gradient = iv[30L]), 
       message = if (19 <= iv1 && iv1 <= 43) {
         if (any(B <- iv1 == stats:::port_cpos)) sprintf("'control' component '%s' = %g, is out of range", 
                                                 names(stats:::port_cpos)[B], v[iv1]) else sprintf("V[IV[1]] = V[%d] = %g is out of range (see PORT docu.)", 
                                                                                           iv1, v[iv1])
       } else stats:::port_msg(iv1))
}

