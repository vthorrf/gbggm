dNEG <- function(x, mu=0, k=1, theta=1, log=T) {
  kappa <- {{2^k} * k}/{{theta*sqrt(pi)}*gamma(k+.5)}
  density <- kappa * exp({{x-mu}^2}/{4*{theta^2}}) * paracyl(x, mu, k, theta)
  if(log) {
    return(log(density))
  } else {
    density
  }
}
