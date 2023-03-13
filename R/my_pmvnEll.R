my_pmvnEll<-function(chisq,sigma,mu,e,x0,lower.tail=TRUE){
  L <- .Internal(La_chol(e, T, -1))
  Esqrt <- L[, order(attr(L, "pivot"))]
  xmu1 <- Esqrt %*% (x0 - mu)
  sigma1 <- Esqrt %*% sigma %*% t(Esqrt)

  z <- .Internal(La_rs(sigma1,F))
  ord <- rev(seq_along(z$values))
  S1eig <- list(values=z$values[ord],vectors=z$vectors[,ord])
  
  xmu2 <- t(S1eig$vectors) %*% xmu1
  ncp <- xmu2^2/S1eig$values
  r<-length(S1eig$values)
  cqf <- if (r == 1L) {
    pchisq(chisq/S1eig$values, df = 1, ncp = ncp, lower.tail = FALSE)
  } else {
    1-.C("ruben", lambda = as.double(S1eig$values), h = rep(1L,r), 
       delta = as.double(ncp), r = r, q = as.double(chisq), 
       mode = as.double(1), maxit = 100000L, eps = as.double(10^(-10)), 
       dnsty = as.double(0), ifault = 0L, 
       res = as.double(0), PACKAGE = "CompQuadForm")$res
  }
  if (lower.tail) {return(1 - cqf)} else {return(cqf)}
}

