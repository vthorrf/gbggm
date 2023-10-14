dhalfcauchy <- function(x, scale = 1, log = FALSE) {
  dtrunc(x, "cauchy", 0, Inf, log=log, location=0, scale=scale)
}