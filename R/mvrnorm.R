mvrnorm<-function(n,mu,Sigma,tol=1e-6){ 
  p <- length(mu)
  eS <- my_eigen(Sigma)
  ev <- eS$values
  if (!all(ev >= -tol * abs(ev[1L]))) 
    stop("'Sigma' is not positive definite")
  X <- matrix(rnorm(p * n), n)
  X <- drop(mu) + eS$vectors %*% diag(sqrt(pmax(ev, 0)), p) %*% t(X)
  if (n == 1) 
    drop(X)
  else t(X)
}

