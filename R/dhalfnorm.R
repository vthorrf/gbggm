dhalfnorm <- function(x, sd = 1, log = FALSE) {
  dtrunc(x, "norm", 0, Inf, log=log, mean=0, sd=sd)
}