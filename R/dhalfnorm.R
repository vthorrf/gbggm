dhalfnorm <- function(x, mean = 0, sd = 1, log = FALSE) {
  y <- x; y[y < 0] <- -Inf
  y <- as.vector(y)
  location <- as.vector(mean)
  scale <- as.vector(sd)
  density <- if(log) {
    log(2*dnorm(x=y, mean=location, sd=scale, log=F))
  } else {
    2*dnorm(x=y, mean=location, sd=scale, log=F)
  }
  return(density)
}