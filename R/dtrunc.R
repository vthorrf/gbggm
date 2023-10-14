dtrunc <- function (x, dist, a = -Inf, b = Inf, log = FALSE, ...) {
  if (a >= b) {
    stop("Lower bound a is not less than upper bound b.")
  }
  if (any(x < a) | any(x > b)) {
    stop("At least one value of x is outside of the bounds.")
  }
  pdf <- get(paste("d", dist, sep = ""), mode = "function")
  cdf <- get(paste("p", dist, sep = ""), mode = "function")
  if (log) {
    dens <- pdf(x, log=TRUE, ...) - log(cdf(b, ...) - cdf(a, ...))
  }
  else {
    dens <- pdf(x, ...)/(cdf(b, ...) - cdf(a, ...))
  }
  return(dens)
}