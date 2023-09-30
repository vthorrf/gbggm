dtrunc <- function (x, spec, a = -Inf, b = Inf, log = FALSE, ...) {
  if (a >= b) 
    stop("Lower bound a is not less than upper bound b.")
  if (any(x < a) | any(x > b)) 
    stop("At least one instance of (x < a) or (x > b) found.")
  dens <- rep(0, length(x))
  g <- get(paste("d", spec, sep = ""), mode = "function")
  G <- get(paste("p", spec, sep = ""), mode = "function")
  if (log == TRUE) {
    dens <- g(x, log = TRUE, ...) - log(G(b, ...) - G(a, 
                                                      ...))
  }
  else {
    dens <- g(x, ...)/(G(b, ...) - G(a, ...))
  }
  return(dens)
}