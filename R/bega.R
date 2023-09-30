bega <- function(data, reg, cor,...) {
  ## Check what regularization method to use
  if(is.null(reg)) reg <- "normal"
  if(!{reg %in% c("normal", "laplace", "logistic", "cauchy", "t", "lomax")}) stop("Unknown regularization method!")
  
  ## Check what correlation method to estimate
  if(is.null(cor)) cor <- "pearson"
  if(!{cor %in% c("pearson","spearman","poly")}) stop("Unknown correlation method!")

  ## Set information
  V <- ncol(data) # Number of variables
  N <- nrow(data) # Sample size
  
  ## Generate data list, random initial values and select the model
  Data <- dataSparseFUN(data=data, V=V, N=N, reg=reg, cor=cor)
  Model <- sparseFUN(reg=reg, cor=cor)
  Initial.Values <- Data$PGF(Data)
  
  ## Fit the model and generate summaries
  fit <- MCMC(Model=Model, Data=Data, Initial.Values=Initial.Values, ...)
  #R_hat   <- get_Rhat(fit, Data)
  #n_cor <- {ncol(R_hat)*{ncol(R_hat)-1}}/2
  #Rho_hat <- pcor(R_hat)
  #sample_R_hat <- fit$Monitor[,1:n_cor]
  #sample_Rho_hat <- fit$Monitor[,{n_cor+1}:{n_cor*2}]
  ## Get matrix of probability of inclusion if the estimated model is sparse
  #sample_Delta_hat <- fit$Monitor[,{{n_cor*2}+1}:{n_cor*3}]
  #Delta_hat <- diag(V)
  #Delta_hat[lower.tri(Delta_hat)] <- colMeans(sample_Delta_hat)
  #Delta_hat[upper.tri(Delta_hat)] <- t(Delta_hat)[upper.tri(Delta_hat)]
  ## Get the communities
  #sample_communities <- fit$Monitor[,{{n_cor*3}+1}:ncol(fit$Monitor)]
  #communities <- apply(sample_communities, 2, median)
  
  ## Final results
  Results <- list(#"R_hat"=R_hat, "Rho_hat"=Rho_hat, "Delta_hat"=Delta_hat,
                  #"wc"=communities,
                  #"sample_R_hat"=sample_R_hat,
                  #"sample_Rho_hat"=sample_Rho_hat,
                  #"sample_Delta_hat"=sample_Delta_hat,
                  #"sample_wc"=sample_communities,
                  "full_output"=fit)
  class(Results) <- "bega"
  return(Results)
}