dlomax    <- function(x, location=0, scale=1, shape=1, log=FALSE) {
  x <- as.vector(x)
  location <- as.vector(location)
  scale <- as.vector(scale)
  shape <- as.vector(shape)
  density <- log(.5) + log(shape) - log(scale) - {{shape+1}*log(1+{abs(x-location)/scale})}
  return(if(log) density else exp(density))
}