sparseFUN <- function(reg, cor) {
  if(cor %in% c("pearson","spearman")) {
    if(reg %in% c("normal", "laplace", "logistic", "cauchy", "hypersec")) {
      Model <- function(parm, Data){
        
        ## Prior parameters
        lambda <- exp(parm[Data$pos.lambda])
        gamma  <- exp(parm[Data$pos.gamma])
        sigma  <- exp(parm[Data$pos.sigma])
        alpha  <- parm[Data$pos.alpha]
        
        ### Log-Priors
        lambda.prior <- sum( dhalfcauchy(lambda, 1, log=T) )
        gamma.prior  <- sum( dhalfcauchy(gamma, 1, log=T) )
        sigma.prior  <- sum( dhalfcauchy(sigma, 1, log=T) )
        alpha.prior  <- sum( Data$density(alpha, 0, 1/lambda, log=T) )
        Lpp <- lambda.prior + gamma.prior + sigma.prior + alpha.prior
        
        ### Log-Likelihood
        kappa   <- 1 - {ptruncnorm(abs(alpha), 0, sigma) - ptruncnorm(abs(alpha), lambda, sigma)}
        keep <- ptruncnorm(abs(alpha), lambda, sigma)/{1-ptruncnorm(abs(alpha), 0, sigma)}
        C_hat   <- Imat <- diag(Data$V)
        C_hat[lower.tri(C_hat)] <- alpha
        norms <- apply(C_hat, 1, euclidean)
        L_hat   <- t(t(C_hat) %*% diag(1/norms))
        Rho_hat <- {{-{L_hat %*% t(L_hat)}} + 2*Imat}
        R_hat   <- corp(Rho_hat)
        LL <- sum( dmnorm(Data$X, sigma=R_hat, log=T) ) + sum(log(dhorseshoe(alpha, 1-kappa, gamma)))
        
        ### Estimates
        yhat <- c(MASS::mvrnorm(Data$N, rep(0, Data$V), R_hat))
        
        ### Log-Posterior
        LP <- Lpp + LL
        
        ### Output
        Monitor=c(R_hat[lower.tri(R_hat)], Rho_hat[lower.tri(R_hat)], {keep > {1/30}}*1)
        Modelout <- list(LP=LP, Dev=-2*LL, Monitor=Monitor, parm=parm, yhat=yhat)
        return(Modelout)
      }
    } else if (reg %in% c("t", "lomax", "NEG")) {
      Model <- function(parm, Data){
        
        ## Prior parameters
        lambda <- exp(parm[Data$pos.lambda])
        gamma  <- exp(parm[Data$pos.gamma])
        tau    <- exp(parm[Data$pos.tau])
        sigma  <- exp(parm[Data$pos.sigma])
        alpha  <- parm[Data$pos.alpha]
        
        ### Log-Priors
        lambda.prior <- sum( dhalfcauchy(lambda, 1, log=T) )
        gamma.prior  <- sum( dhalfcauchy(gamma, 1, log=T) )
        tau.prior    <- sum( dhalfcauchy(tau, 1, log=T) )
        sigma.prior  <- sum( dhalfcauchy(sigma, 1, log=T) )
        alpha.prior  <- sum( Data$density(alpha, 0, 1/lambda, tau, log=T) )
        Lpp <- lambda.prior + tau.prior + gamma.prior + sigma.prior + alpha.prior
        
        ### Log-Likelihood
        kappa   <- 1 - {ptruncnorm(abs(alpha), 0, sigma) - ptruncnorm(abs(alpha), lambda, sigma)}
        keep <- ptruncnorm(abs(alpha), lambda, sigma)/{1-ptruncnorm(abs(alpha), 0, sigma)}
        C_hat   <- Imat <- diag(Data$V)
        C_hat[lower.tri(C_hat)] <- alpha
        norms <- apply(C_hat, 1, euclidean)
        L_hat   <- t(t(C_hat) %*% diag(1/norms))
        Rho_hat <- {{-{L_hat %*% t(L_hat)}} + 2*Imat}
        R_hat   <- corp(Rho_hat)
        LL <- sum( dmnorm(Data$X, sigma=R_hat, log=T) ) + sum(log(dhorseshoe(alpha, 1-kappa, gamma)))
        
        ### Estimates
        yhat <- c(MASS::mvrnorm(Data$N, rep(0, Data$V), R_hat))
        
        ### Log-Posterior
        LP <- Lpp + LL
        
        ### Output
        Monitor=c(R_hat[lower.tri(R_hat)], Rho_hat[lower.tri(R_hat)], {keep > {1/30}}*1)
        Modelout <- list(LP=LP, Dev=-2*LL, Monitor=Monitor, parm=parm, yhat=yhat)
        return(Modelout)
      }
    } else {
      stop("Unknown regularization prior!")
    }
  } else if(cor == "poly") {
    if(reg %in% c("normal", "laplace", "logistic", "cauchy", "hypersec")) {
      Model <- function(parm, Data){
        
        ## Prior parameters
        lambda <- exp(parm[Data$pos.lambda])
        gamma  <- exp(parm[Data$pos.gamma])
        sigma  <- exp(parm[Data$pos.sigma])
        alpha  <- parm[Data$pos.alpha]
        delta  <- parm[Data$pos.delta]
        
        ### Log-Priors
        lambda.prior <- sum( dhalfcauchy(lambda, 1, log=T) )
        gamma.prior  <- sum( dhalfcauchy(gamma, 1, log=T) )
        sigma.prior  <- sum( dhalfcauchy(sigma, 1, log=T) )
        alpha.prior  <- sum( Data$density(alpha, 0, 1/lambda, log=T) )
        delta.prior  <- sum( dnorm(delta, 0, 1, log=T) )
        Lpp <- lambda.prior + delta.prior + gamma.prior + sigma.prior + alpha.prior
        
        ### Log-Likelihood
        kappa   <- 1 - {ptruncnorm(abs(alpha), 0, sigma) - ptruncnorm(abs(alpha), lambda, sigma)}
        keep <- ptruncnorm(abs(alpha), lambda, sigma)/{1-ptruncnorm(abs(alpha), 0, sigma)}
        C_hat   <- Imat <- diag(Data$V)
        C_hat[lower.tri(C_hat)] <- alpha
        norms <- apply(C_hat, 1, euclidean)
        L_hat   <- t(t(C_hat) %*% diag(1/norms))
        Rho_hat <- {{-{L_hat %*% t(L_hat)}} + 2*Imat}
        R_hat   <- corp(Rho_hat)
        taus <- lapply(unique(Data$id.delta), function(g) delta[which(Data$id.delta == g)])
        LL <- sum(dpoly(Data$X, R=R_hat, taus=taus)$loglik[lower.tri(R_hat)] *
                    {Data$N + Data$n_par + Data$n_thr}) + sum(log(dhorseshoe(alpha, 1-kappa, gamma)))
        
        ### Estimates
        yhat <- c(MASS::mvrnorm(Data$N, rep(0, Data$V), R_hat))
        
        ### Log-Posterior
        LP <- Lpp + LL
        
        ### Output
        Monitor=c(R_hat[lower.tri(R_hat)], Rho_hat[lower.tri(R_hat)], {keep > {1/30}}*1)
        Modelout <- list(LP=LP, Dev=-2*LL, Monitor=Monitor, parm=parm, yhat=yhat)
        return(Modelout)
      }
    } else if (reg %in% c("t", "lomax", "NEG")) {
      Model <- function(parm, Data){
        
        ## Prior parameters
        lambda <- exp(parm[Data$pos.lambda])
        gamma  <- exp(parm[Data$pos.gamma])
        tau    <- exp(parm[Data$pos.tau])
        sigma  <- exp(parm[Data$pos.sigma])
        alpha  <- parm[Data$pos.alpha]
        delta  <- parm[Data$pos.delta]
        
        ### Log-Priors
        lambda.prior <- sum( dhalfcauchy(lambda, 1, log=T) )
        gamma.prior  <- sum( dhalfcauchy(gamma, 1, log=T) )
        tau.prior    <- sum( dhalfcauchy(tau, 1, log=T) )
        sigma.prior  <- sum( dhalfcauchy(sigma, 1, log=T) )
        alpha.prior  <- sum( Data$density(alpha, 0, 1/lambda, tau, log=T) )
        delta.prior  <- sum( dnorm(delta, 0, 1, log=T) )
        Lpp <- lambda.prior + tau.prior + delta.prior + gamma.prior + sigma.prior + alpha.prior
        
        ### Log-Likelihood
        kappa   <- 1 - {ptruncnorm(abs(alpha), 0, sigma) - ptruncnorm(abs(alpha), lambda, sigma)}
        keep <- ptruncnorm(abs(alpha), lambda, sigma)/{1-ptruncnorm(abs(alpha), 0, sigma)}
        C_hat   <- Imat <- diag(Data$V)
        C_hat[lower.tri(C_hat)] <- alpha
        norms <- apply(C_hat, 1, euclidean)
        L_hat   <- t(t(C_hat) %*% diag(1/norms))
        Rho_hat <- {{-{L_hat %*% t(L_hat)}} + 2*Imat}
        R_hat   <- corp(Rho_hat)
        taus <- lapply(unique(Data$id.delta), function(g) delta[which(Data$id.delta == g)])
        LL <- sum(dpoly(Data$X, R=R_hat, taus=taus)$loglik[lower.tri(R_hat)] *
                    {Data$N + Data$n_par + Data$n_thr}) + sum(log(dhorseshoe(alpha, 1-kappa, gamma)))
        
        ### Estimates
        yhat <- c(MASS::mvrnorm(Data$N, rep(0, Data$V), R_hat))
        
        ### Log-Posterior
        LP <- Lpp + LL
        
        ### Output
        Monitor=c(R_hat[lower.tri(R_hat)], Rho_hat[lower.tri(R_hat)], {keep > {1/30}}*1)
        Modelout <- list(LP=LP, Dev=-2*LL, Monitor=Monitor, parm=parm, yhat=yhat)
        return(Modelout)
      }
    } else {
      stop("Unknown regularization prior!")
    }
  } else {
    stop("Unkown correlation method")
  }
  return(Model)
}