gbggm <- function(data, reg, cor, sparse, ...) {
  ## Check what regularization method to use
  if(is.null(reg)) reg <- "normal"
  if(!{reg %in% c("normal", "laplace", "logistic", "cauchy", "t", "lomax")}) stop("Unknown regularization method!")
  
  ## Check what correlation method to estimate
  if(is.null(cor)) cor <- "pearson"
  if(!{cor %in% c("pearson","spearman","poly")}) stop("Unknown correlation method!")
  
  ## Check if the model is sparse of not
  if(is.null(sparse)) sparse <- FALSE
  if(!is.logical(sparse)) stop("The `sparse` argument should be logical.")
  
  ## Set information
  V <- ncol(data) # Number of variables
  N <- nrow(data) # Sample size
  
  ## Generate data list, random initial values and select the model
  if(sparse) {
    Data <- dataSparseFUN(X=data, V=V, N=N, reg=reg, cor=cor)
    Model <- sparseFUN(reg=reg, cor=cor)
  } else {
    Data <- dataListFUN(X=data, V=V, N=N, reg=reg, cor=cor)
    Model <- modelFUN(reg=reg, cor=cor)
  }
  Initial.Values <- Data$PGF(Data)
  
  ## Fit the model and generate summaries
  fit <- MCMC(Model=Model, Data=Data, Initial.Values=Initial.Values, ...)
  R_hat   <- get_Rhat(fit, Data)
  Rho_hat <- pcor(R_hat)
  sample_R_hat <- fit$Monitor[,1:{ncol(fit$Monitor)/3}]
  sample_Rho_hat <- fit$Monitor[,{{ncol(fit$Monitor)/3}+1}:{{ncol(fit$Monitor)/3}*2}]
  # Get matrix of probability of inclusion if the estimated model is sparse
  if(sparse) {
    sample_Delta_hat <- fit$Monitor[,{{{ncol(fit$Monitor)/3}*2}+1}:ncol(fit$Monitor)]
    Delta_hat <- diag(V)
    Delta_hat[lower.tri(Delta_hat)] <- colMeans(sample_Delta_hat)
    Delta_hat[upper.tri(Delta_hat)] <- t(Delta_hat)[upper.tri(Delta_hat)]
  }
  
  ## Final results
  if(sparse) {
    Results <- list("R_hat"=R_hat, "Rho_hat"=Rho_hat, "Delta_hat"=Delta_hat,
                    "sample_R_hat"=sample_R_hat,
                    "sample_Rho_hat"=sample_Rho_hat,
                    "sample_Delta_hat"=sample_Delta_hat,
                    "full_output"=fit)
  } else {
    Results <- list("R_hat"=R_hat, "Rho_hat"=Rho_hat,
                    "sample_R_hat"=sample_R_hat,
                    "sample_Rho_hat"=sample_Rho_hat,
                    "full_output"=fit)
  }
  class(Results) <- "bggm"
  return(Results)
}