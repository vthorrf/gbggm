ptruncnorm <- function(x, mean=0, sd=1) {
  ptrunc(x, dist = "norm", a = 0, b = Inf, mean=mean, sd=sd)
}