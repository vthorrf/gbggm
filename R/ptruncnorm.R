ptruncnorm <- function(x, mean=0, sd=1) {
  ptrunc(x, "norm", 0, Inf, mean=mean, sd=sd)
}