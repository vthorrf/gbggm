dmnorm    <- function(x, mu=rep(0, ncol(x)), sigma, log=FALSE){
  k <- ncol(sigma)
  x <- t(x)
  density <- {{-.5} *diag(t(x-mu)%*%solve(sigma)%*%(x-mu))} - log(sqrt(((2*pi)^k)*det(sigma)))
  return(if(log) density else exp(density))
}