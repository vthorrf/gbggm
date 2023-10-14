paracyl <- function(x, mu=0, k=1, theta=1) {
  vals <- x
  for(i in seq_along(vals)) {
    vals[i] <- integrate(paracyl_weight, 0, Inf, x=x[i], mu=mu, k=k, theta=theta)$value
  }
  result <- {1/gamma({2*k}+1)} * exp(-{{x-mu}^2}/{4*{theta^2}}) * vals
  return(result)
}