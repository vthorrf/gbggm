dhypersec <- function(x, mu=0, sigma=1, log=TRUE) {
  density <- -log(2) - log(sigma) - log( cosh( 0.5*pi*(x-mu)/sigma ) ) 
  if(log) {
    return(density)
  } else {
    return(exp(density))
  }
}