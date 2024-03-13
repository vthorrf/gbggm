bggm <- function(data, reg=NULL, cor=NULL, sparse=NULL, method=NULL, full=FALSE, ...) {
  ## Check what regularization method to use
  if(is.null(reg)) reg <- "normal"
  if(!{reg %in% c("normal", "laplace", "logistic", "cauchy",
                  "hypersec", "t", "lomax", "NEG")}) stop("Unknown regularization method!")
  
  ## Check what correlation method to estimate
  if(is.null(cor)) cor <- "pearson"
  if(!{cor %in% c("pearson","spearman","poly")}) stop("Unknown correlation method!")
  
  ## Check if the model is sparse of not
  if(is.null(sparse)) sparse <- FALSE
  if(!is.logical(sparse)) stop("The `sparse` argument should be logical.")
  
  ## Check what estimationg method to use
  if(is.null(method)) method <- "LA"
  if(!{method %in% c("LA","MCMC")}) stop("Unknown estimation method!")
  
  ## Set information
  V <- ncol(data) # Number of variables
  N <- nrow(data) # Sample size
  
  ## Generate data list, random initial values and select the model
  if(sparse) {
    Data <- dataSparseFUN(data=data, V=V, N=N, reg=reg, cor=cor)
    Model <- sparseFUN(reg=reg, cor=cor)
  } else {
    Data <- dataListFUN(data=data, V=V, N=N, reg=reg, cor=cor)
    Model <- modelFUN(reg=reg, cor=cor)
  }
  Initial.Values <- Data$PGF(Data)
  
  ## Fit the model and generate summaries
  if(method == "LA") {
    fit <- LA(Model=Model, Data=Data, Initial.Values=Initial.Values, ...)
  } else {
    fit <- MCMC(Model=Model, Data=Data, Initial.Values=Initial.Values, ...)
  }
  # Get matrix of probability of inclusion if the estimated model is sparse
  if(sparse) {
    Rho_hat   <- get_Rhat(fit, Data)
    Rho_hat[!diag(ncol(Rho_hat))] <- -Rho_hat[!diag(ncol(Rho_hat))]
    R_hat <- tryCatch(corp(Rho_hat), error=function(e) NA)
    if(any(is.na(R_hat))) {
      R_hat <- get_Rhat(fit, Data)
      Rho_hat <- pcor(R_hat)
      sample_R_hat <- fit$Monitor[,1:ncol(R_hat)]
      sample_Rho_hat <- fit$Monitor[,{ncol(R_hat)+1}:{ncol(R_hat)*2}]
    } else {
      sample_R_hat <- -fit$Monitor[,{ncol(R_hat)+1}:{ncol(R_hat)*2}]
      sample_Rho_hat <- -fit$Monitor[,1:ncol(R_hat)]
    }
    sample_Delta_hat <- fit$Monitor[,{{ncol(R_hat)*2}+1}:ncol(fit$Monitor)]
    Delta_hat <- diag(V)
    Delta_hat[lower.tri(Delta_hat)] <- colMeans(sample_Delta_hat)
    Delta_hat[upper.tri(Delta_hat)] <- t(Delta_hat)[upper.tri(Delta_hat)]
  } else {
    R_hat   <- get_Rhat(fit, Data)
    Rho_hat <- pcor(R_hat)
    sample_R_hat <- fit$Monitor[,1:ncol(R_hat)]
    sample_Rho_hat <- fit$Monitor[,{ncol(R_hat)+1}:{ncol(R_hat)*2}]
  }
  
  ## Final results
  if(sparse) {
    Results <- list("R_hat"=R_hat, "Rho_hat"=Rho_hat, "Delta_hat"=Delta_hat,
                    "sample_R_hat"=sample_R_hat,
                    "sample_Rho_hat"=sample_Rho_hat,
                    "sample_Delta_hat"=sample_Delta_hat)
    if(full) {
      Results$full_output <- fit
    }
  } else {
    Results <- list("R_hat"=R_hat, "Rho_hat"=Rho_hat,
                    "sample_R_hat"=sample_R_hat,
                    "sample_Rho_hat"=sample_Rho_hat)
    if(full) {
      Results$full_output <- fit
    }
  }
  class(Results) <- "bggm"
  return(Results)
}