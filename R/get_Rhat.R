get_Rhat <- function(fit, Data) {
  alpha <- colMeans(fit$posterior[,grep("alpha",colnames(fit$posterior))])
  C_hat <- diag(Data$V)
  C_hat[lower.tri(C_hat)] <- alpha
  norms <- apply(C_hat, 1, Data$euclidean)
  L_hat <- t(t(C_hat) %*% diag(1/norms))
  R_hat <- L_hat %*% t(L_hat)
  return(R_hat)
}