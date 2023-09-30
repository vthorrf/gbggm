phalfnorm <- function(x, mean=0, sd=1) {
  {erf({x-mean}/{sd*sqrt(2)}) + 1}/2
}