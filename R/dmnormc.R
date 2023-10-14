dmnormc <- function(x, mu=NULL, L, log=TRUE) {
  if(is.null(mu)) mu <- colMeans(x)
  d <- ncol(x)
  Mu <- matrix(mu, byrow=T, nrow=nrow(x), ncol=ncol(x))
  z <- solve(L) %*% t(x - Mu)
  logDensity <- {{-d/2}*log(2*pi)} +
    {-.5 * 2 * log(diag(L))} +
    {-.5 * colSums(z * z)}
  if(log) {
    return(logDensity)
  } else {
    return(exp(logDensity))
  }
}