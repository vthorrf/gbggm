get_Rhat <- function(fit, Data) {
  n_par <- {Data$V * {Data$V-1}}/2
  R_hat <- diag(Data$V)
  R_hat[lower.tri(R_hat)] <- fit$Monitor[which.max(fit$LP),1:n_par]
  R_hat[upper.tri(R_hat)] <- t(R_hat)[upper.tri(R_hat)]
  return(R_hat)
}