dhalfcauchy <- function(x, location = 0, scale = 1, log = FALSE) {
  y <- x; y[y < 0] <- -Inf
  y <- as.vector(y)
  location <- as.vector(location)
  scale <- as.vector(scale)
  density <- if(log) {
    log(2*dcauchy(x=y, location=location, scale=scale, log=F))
  } else {
    2*dcauchy(x=y, location=location, scale=scale, log=F)
  }
  return(density)
}