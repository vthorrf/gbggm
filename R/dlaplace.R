dlaplace  <- function(x, location=0, scale=1, log=FALSE) {
  x <- as.vector(x)
  location <- as.vector(location)
  scale <- as.vector(scale)
  density <- {-abs(x - location)/scale} - log(2 * scale)
  return(if(log) density else exp(density))
}