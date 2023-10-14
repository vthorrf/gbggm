dtruncnorm <- function(x, mean=0, sd=1, log=FALSE) {
  dtrunc(x, "norm", 0, Inf, log=log, mean=mean, sd=sd)
}