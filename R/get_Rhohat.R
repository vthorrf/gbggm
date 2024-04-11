get_Rhohat <- function(fit, Data) {
  n_par <- {Data$V * {Data$V-1}}/2
  Rho_hat <- diag(Data$V)
  Rho_hat[lower.tri(Rho_hat)] <- fit$Monitor[which.max(fit$LP),{n_par+1}:{n_par*2}]
  Rho_hat[upper.tri(Rho_hat)] <- t(Rho_hat)[upper.tri(Rho_hat)]
  return(Rho_hat)
}