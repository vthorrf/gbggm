dkaniadakis <- function(x, mu=0, beta=1, kappa=.5, log=F) {
  # Checks
  if({kappa >= 1} | {kappa <= 0}) stop("kappa should be in the interval (0, 1)")
  if(beta <= 0) stop("beta should be positive")
  # PDF
  beta <- 1/beta
  y    <- -beta * (x - mu) * (x - mu)
  M    <- {1/{2*kappa}}+.25
  m    <- {1/{2*kappa}}-.25
  Zkappa   <- log(sqrt({2*beta*kappa}/pi)) +
              log(1+{.5*kappa}) +
              lgamma(M) - 
              lgamma(m)
  expKappa <- log(expKappa(y, kappa))
  # Result
  if(log) {
    Zkappa + expKappa
  } else {
    exp(Zkappa + expKappa)
  }
}

expKappa <- function(x, kappa=1) {
  if(kappa == 0) {
    exp(x)
  } else {
    {sqrt(1+{kappa*kappa}*{x*x}) + kappa*x}^{1/kappa}
  }
}