modelNEG <- function(cor) {
  if(cor %in% c("spearman", "pearson", "s", "p")) {
    Model <- function(parm, Data){
      
      ## Prior parameters
      lambda <- exp(parm[Data$pos.lambda])
      tau    <- exp(parm[Data$pos.tau])
      alpha  <- parm[Data$pos.alpha]
      
      ### Log-Likelihood
      C_hat <- diag(Data$V)
      C_hat[lower.tri(C_hat)] <- alpha
      norms <- apply(C_hat, 1, euclidean)
      L_hat <- t(t(C_hat) %*% diag(1/norms))
      R_hat <- L_hat %*% t(L_hat)
      LL <- sum( dmnorm(Data$X, sigma=R_hat, log=T) )
      
      ### Clustering of (dis)connected nodes
      Rho_hat <- pcor(R_hat)[lower.tri(R_hat)]
      
      ### Log-Priors
      lambda.prior <- sum( dhalfcauchy(lambda, 1, log=T) )
      tau.prior    <- sum( dhalfcauchy(tau, 1, log=T) )
      alpha.prior  <- sum( dNEG(solve(L_hat)[lower.tri(L_hat)], 0, 1/lambda, tau, log=T) )
      Lpp <- lambda.prior + alpha.prior + tau.prior
      
      ### Estimates
      yhat <- c(MASS::mvrnorm(Data$N, rep(0, Data$V), R_hat))
      
      ### Log-Posterior
      LP <- Lpp + LL
      
      ### Output
      Monitor=c(R_hat[lower.tri(R_hat)], Rho_hat)
      Modelout <- list(LP=LP, Dev=-2*LL, Monitor=Monitor, parm=parm, yhat=yhat)
      return(Modelout)
    }
  } else if(cor == "poly") {
    Model <- function(parm, Data){
      
      ## Prior parameters
      lambda <- exp(parm[Data$pos.lambda])
      tau    <- exp(parm[Data$pos.tau])
      alpha  <- parm[Data$pos.alpha]
      delta  <- parm[Data$pos.delta]
      
      ### Log-Likelihood
      C_hat <- diag(Data$V)
      C_hat[lower.tri(C_hat)] <- alpha
      norms <- apply(C_hat, 1, euclidean)
      L_hat <- t(t(C_hat) %*% diag(1/norms))
      R_hat <- L_hat %*% t(L_hat)
      taus <- lapply(unique(Data$id.delta), function(g) delta[which(Data$id.delta == g)])
      LL <- sum(dpoly(Data$X, R=R_hat, taus=taus) * {Data$N + Data$n_par + Data$n_thr} * .5)
      
      ### Clustering of (dis)connected nodes
      Rho_hat   <- pcor(R_hat)[lower.tri(R_hat)]
      
      ### Log-Priors
      lambda.prior <- sum( dhalfcauchy(lambda, 1, log=T) )
      tau.prior    <- sum( dhalfcauchy(tau, 1, log=T) )
      alpha.prior  <- sum( dNEG(solve(L_hat)[lower.tri(L_hat)], 0, 1/lambda, tau, log=T) )
      delta.prior  <- sum( dnorm(delta, 0, 1, log=T) )
      Lpp <- lambda.prior + alpha.prior + tau.prior + delta.prior
      
      ### Estimates
      yhat <- c(MASS::mvrnorm(Data$N, rep(0, Data$V), R_hat))
      
      ### Log-Posterior
      LP <- Lpp + LL
      
      ### Output
      Monitor=c(R_hat[lower.tri(R_hat)], Rho_hat)
      Modelout <- list(LP=LP, Dev=-2*LL, Monitor=Monitor, parm=parm, yhat=yhat)
      return(Modelout)
    }
  } else{
    stop("Unkown correlation method!")
  }
  return(Model)
}