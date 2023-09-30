ptrunc <- function (x, spec, a = -Inf, b = Inf, ...) {
  if (a >= b) 
    stop("Lower bound a is not less than upper bound b.")
  if (any(x < a) | any(x > b)) 
    stop("At least one instance of (x < a) or (x > b) found.")
  p <- x
  aa <- rep(a, length(x))
  bb <- rep(b, length(x))
  G <- get(paste("p", spec, sep = ""), mode = "function")
  p <- G(apply(cbind(apply(cbind(x, bb), 1, min), aa), 1, max), 
         ...)
  p <- p - G(aa, ...)
  p <- p/{
    G(bb, ...) - G(aa, ...)
  }
  return(p)
}