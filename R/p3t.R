p3t <- function (q, mu = 0, sigma = 1, nu = 10, lower.tail = TRUE, log.p = FALSE) {
  q <- as.vector(q)
  mu <- as.vector(mu)
  sigma <- as.vector(sigma)
  nu <- as.vector(nu)
  if (any(sigma <= 0)) 
    stop("The sigma parameter must be positive.")
  if (any(nu <= 0)) 
    stop("The nu parameter must be positive.")
  NN <- max(length(q), length(mu), length(sigma), length(nu))
  q <- rep(q, len = NN)
  mu <- rep(mu, len = NN)
  sigma <- rep(sigma, len = NN)
  nu <- rep(nu, len = NN)
  p <- pt({
    q - mu
  }/sigma, df = nu, lower.tail = lower.tail, log.p = log.p)
  temp <- which(nu > 1e+06)
  p[temp] <- pnorm(q[temp], mu[temp], sigma[temp], lower.tail = lower.tail, 
                   log.p = log.p)
  return(p)
}