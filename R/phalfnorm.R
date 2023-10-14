phalfnorm <- function(x, sd=1) {
  ptrunc(x, "norm", 0, Inf, mean=0, sd=sd)
  #{erf(x/{sd*sqrt(2)}) + 1}/2
}