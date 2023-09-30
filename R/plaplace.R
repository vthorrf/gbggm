plaplace <- function (q, location = 0, scale = 1) {
  q <- as.vector(q)
  location <- as.vector(location)
  scale <- as.vector(scale)
  if (any(scale <= 0)) 
    stop("The scale parameter must be positive.")
  z <- {
    q - location
  }/scale
  NN <- max(length(q), length(location), length(scale))
  q <- rep(q, len = NN)
  location <- rep(location, len = NN)
  p <- q
  temp <- which(q < location)
  p[temp] <- 0.5 * exp(z[temp])
  temp <- which(q >= location)
  p[temp] <- 1 - 0.5 * exp(-z[temp])
  return(p)
}