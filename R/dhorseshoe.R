dhorseshoe <- function (x, lambda, tau, log = FALSE) {
  dens <- dnorm(x, 0, lambda * tau, log = log)
  return(dens)
}