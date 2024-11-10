bggm <- function(data, reg=NULL, cor=NULL, sparse=NULL, est=NULL, full=FALSE, ...) {
  ## Check what regularization method to use
  if(is.null(reg)) reg <- "normal"
  if(!{reg %in% c("normal", "laplace", "logistic", "cauchy", "kaniadakis",
                  "hypersec", "t", "lomax", "NEG")}) stop("Unknown regularization method!")
  
  ## Check what correlation method to estimate
  if(is.null(cor)) cor <- "pearson"
  if(!{cor %in% c("pearson","spearman","poly")}) stop("Unknown correlation method!")
  
  ## Check if the model is sparse of not
  if(is.null(sparse)) sparse <- FALSE
  if(!is.logical(sparse)) stop("The `sparse` argument should be logical.")
  
  ## Check what estimationg method to use
  if(is.null(est)) est <- "LA"
  if(!{est %in% c("LA","MCMC")}) stop("Unknown estimation method!")
  
  ## Set information
  V <- ncol(data) # Number of variables
  N <- nrow(data) # Sample size
  
  ## Generate data list, random initial values and select the model
  if(sparse) {
    Data <- dataSparseFUN(data=data, V=V, N=N, reg=reg, cor=cor)
  } else {
    Data <- dataListFUN(data=data, V=V, N=N, reg=reg, cor=cor)
  }
  Model <- modelFUN(reg=reg, cor=cor, sparse=sparse)
  Initial.Values <- Data$PGF(Data)
  
  ## Fit the model and generate summaries
  if(est == "LA") {
    fit <- LA(Model=Model, Data=Data, Initial.Values=Initial.Values, ...)
  } else {
    fit <- MCMC(Model=Model, Data=Data, Initial.Values=Initial.Values, ...)
  }
  # Statistical fit
  K      <- length(Data$parm.names)
  Sigma  <- cov(fit$posterior)
  s      <- rankMatrix(Sigma)
  #lambda <- eigen(Sigma)$values
    # AIC
  AIC    <- min(fit$AIC)
  AICc   <- min(fit$Dev) + {2*K*{N/{N-K-1}}}
  CAIC   <- min(fit$Dev) + K * {log(N) + 1}
    # BIC
  BIC    <- min(fit$Dev) + K * log(N)
  SABIC  <- min(fit$Dev) + K * log({N+2}/24)
    # DIC
  DIC    <- fit$DIC$DIC
  ICOMP  <- min(fit$Dev) + 2 * {{{s/2}*log(sum(diag(Sigma))/s)} - {.5*log(det(Sigma))}}
  #IFIM   <- min(fit$Dev) + s * log(mean(lambda)/exp(mean(log(lambda))))
  indices <- c(AIC, AICc, CAIC, BIC, SABIC, DIC, ICOMP)#, IFIM)
  names(indices) <- c("AIC", "AICc", "CAIC", "BIC", "SABIC", "DIC", "ICOMP")#, "IFIM")
  # Get matrix of probability of inclusion if the estimated model is sparse
  n_par <- {Data$V * {Data$V-1}}/2
  R_hat   <- get_Rhat(fit, Data)
  Rho_hat <- get_Rhohat(fit, Data)
  sample_R_hat <- fit$Monitor[,1:n_par]
  sample_Rho_hat <- fit$Monitor[,{n_par+1}:{n_par*2}]
  if(sparse) {
    sample_Delta_hat <- fit$Monitor[,{{n_par*2}+1}:ncol(fit$Monitor)]
    Delta_hat <- diag(V)
    Delta_hat[lower.tri(Delta_hat)] <- sample_Delta_hat[which.max(fit$LP),]
    Delta_hat[upper.tri(Delta_hat)] <- t(Delta_hat)[upper.tri(Delta_hat)]
  }
  if(cor == "poly") {
    sample_Tau_hat <- fit$posterior[,grep("delta", Data$parm.names)]
    Tau_hat <- fit$posterior[which.max(fit$LP),grep("delta", Data$parm.names)]
  }
  
  ## Final results
  if(sparse) {
    Results <- list("R_hat"=R_hat, "Rho_hat"=Rho_hat, "Delta_hat"=Delta_hat,
                    "sample_R_hat"=sample_R_hat,
                    "sample_Rho_hat"=sample_Rho_hat,
                    "sample_Delta_hat"=sample_Delta_hat,
                    "fit"=indices)
  } else {
    Results <- list("R_hat"=R_hat, "Rho_hat"=Rho_hat,
                    "sample_R_hat"=sample_R_hat,
                    "sample_Rho_hat"=sample_Rho_hat,
                    "fit"=indices)
  }
  if(full) {
    Results$full_output <- fit
  }
  if(cor == "poly") {
    Results$sample_Tau_hat <- sample_Tau_hat
    Results$Tau_hat <- Tau_hat
  }
  class(Results) <- "bggm"
  return(Results)
}