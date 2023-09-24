d3t       <- function(x, location=0, scale=1, df=length(x), log=FALSE) {
  x <- as.vector(x)
  location <- rep(as.vector(location), length=length(x))
  scale <- rep(as.vector(scale), length=length(x))
  df <- rep(as.vector(df), length=length(x))
  inter <- lgamma({df + 1}/2) - lgamma(df/2) - log(sqrt(pi * df) * scale)
  density <- inter + log({1 + {1/df} * {{{x - location}/scale}^2}} ^ {-{df + 1}/2})
  return(if(log) density else exp(density))
}