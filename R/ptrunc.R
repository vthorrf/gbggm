ptrunc <- function (x, dist, a = -Inf, b = Inf, ...) {
  if (a >= b) {
    stop("Lower bound a is not less than upper bound b.")
  }
  if (any(x < a) | any(x > b)) {
    stop("At least one value of x is outside of the bounds.")
  }
  lb <- rep(a, length(x))
  ub <- rep(b, length(x))
  cdf <- get(paste("p", dist, sep = ""), mode = "function")
  p <- cdf(apply(cbind(apply(cbind(x, ub), 1, min), lb), 1, max), ...)
  p <- p - cdf(lb, ...)
  p <- p/{cdf(ub, ...) - cdf(lb, ...)}
  return(p)
}