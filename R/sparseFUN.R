sparseFUN <- function(reg, cor) {
  if(cor %in% c("pearson","spearman")) {
    if(reg %in% c("normal", "laplace", "logistic", "cauchy", "hypersec")) {
      Model <- function(parm, Data){
        
        ## Prior parameters
        lambda <- exp(parm[Data$pos.lambda])
        gamma  <- exp(parm[Data$pos.gamma])
        sigma  <- exp(parm[Data$pos.sigma])
        theta  <- exp(parm[Data$pos.theta])
        alpha  <- parm[Data$pos.alpha]
        
        ### Log-Likelihood
        ## Data
        C_hat <- diag(Data$V)
        C_hat[lower.tri(C_hat)] <- alpha
        norms <- apply(C_hat, 1, euclidean)
        L_hat <- t(t(C_hat) %*% diag(1/norms))
        R_hat <- L_hat %*% t(L_hat)
        d_LL <- sum( dmnorm(Data$X, sigma=R_hat, log=T) )
        ## Sparsity
        Rho_hat <- pcor(R_hat)[lower.tri(R_hat)]
        Rho_tra <- tanh(theta*atanh(Rho_hat))
        kappa   <- ptrunc(atanh(abs(Rho_tra)), "norm", 0, Inf, mean=0, sd=sigma)
        pi <- 1 - kappa
        s_LL <- sum(log( {pi * dnorm(atanh(Rho_hat), mean=0, sd=gamma)} + {kappa*kappa}))
        ## Total Likelihood
        LL <- d_LL + s_LL
        
        ### Log-Priors
        lambda.prior <- sum( dhalfcauchy(lambda, 1, log=T) )
        gamma.prior  <- sum( dhalfcauchy(gamma, 1, log=T) )
        sigma.prior  <- sum( dhalfcauchy(sigma, 1, log=T) )
        alpha.prior  <- sum( Data$density(solve(L_hat)[lower.tri(L_hat)], 0, 1/lambda, log=T) )
        theta.prior  <- sum( dhalfcauchy(theta, 1e2, log=T) )
        Lpp <- lambda.prior + gamma.prior + sigma.prior + alpha.prior + theta.prior
        
        ### Estimates
        yhat <- c(MASS::mvrnorm(Data$N, rep(0, Data$V), R_hat))
        
        ### Log-Posterior
        LP <- Lpp + LL
        
        ### Output
        Monitor=c(R_hat[lower.tri(R_hat)], Rho_hat, kappa)
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
        theta  <- exp(parm[Data$pos.theta])
        alpha  <- parm[Data$pos.alpha]
        
        ### Log-Likelihood
        ## Data
        C_hat <- diag(Data$V)
        C_hat[lower.tri(C_hat)] <- alpha
        norms <- apply(C_hat, 1, euclidean)
        L_hat <- t(t(C_hat) %*% diag(1/norms))
        R_hat <- L_hat %*% t(L_hat)
        d_LL <- sum( dmnorm(Data$X, sigma=R_hat, log=T) )
        ## Sparsity
        Rho_hat <- pcor(R_hat)[lower.tri(R_hat)]
        Rho_tra <- tanh(theta*atanh(Rho_hat))
        kappa   <- ptrunc(atanh(abs(Rho_tra)), "norm", 0, Inf, mean=0, sd=sigma)
        pi <- 1 - kappa
        s_LL <- sum(log( {pi * dnorm(atanh(Rho_hat), mean=0, sd=gamma)} + {kappa*kappa}))
        ## Total Likelihood
        LL <- d_LL + s_LL
        
        ### Log-Priors
        lambda.prior <- sum( dhalfcauchy(lambda, 1, log=T) )
        gamma.prior  <- sum( dhalfcauchy(gamma, 1, log=T) )
        tau.prior    <- sum( dhalfcauchy(tau, 1, log=T) )
        sigma.prior  <- sum( dhalfcauchy(sigma, 1, log=T) )
        alpha.prior  <- sum( Data$density(solve(L_hat)[lower.tri(L_hat)], 0, 1/lambda, tau, log=T) )
        theta.prior  <- sum( dhalfcauchy(theta, 1e2, log=T) )
        Lpp <- lambda.prior + gamma.prior + tau.prior + sigma.prior + alpha.prior + theta.prior
        
        ### Estimates
        yhat <- c(MASS::mvrnorm(Data$N, rep(0, Data$V), R_hat))
        
        ### Log-Posterior
        LP <- Lpp + LL
        
        ### Output
        Monitor=c(R_hat[lower.tri(R_hat)], Rho_hat, kappa)
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
        theta  <- exp(parm[Data$pos.theta])
        alpha  <- parm[Data$pos.alpha]
        delta  <- parm[Data$pos.delta]
        
        ### Log-Likelihood
        ## Data
        C_hat <- diag(Data$V)
        C_hat[lower.tri(C_hat)] <- alpha
        norms <- apply(C_hat, 1, euclidean)
        L_hat <- t(t(C_hat) %*% diag(1/norms))
        R_hat <- L_hat %*% t(L_hat)
        taus <- lapply(unique(Data$id.delta), function(g) delta[which(Data$id.delta == g)])
        d_LL <- sum(dpoly(Data$X, R=R_hat, taus=taus)$loglik[lower.tri(R_hat)] *
                      {Data$N + Data$n_par + Data$n_thr})
        ## Sparsity
        Rho_hat <- pcor(R_hat)[lower.tri(R_hat)]
        Rho_tra <- tanh(theta*atanh(Rho_hat))
        kappa   <- ptrunc(atanh(abs(Rho_tra)), "norm", 0, Inf, mean=0, sd=sigma)
        pi <- 1 - kappa
        s_LL <- sum(log( {pi * dnorm(atanh(Rho_hat), mean=0, sd=gamma)} + {kappa*kappa}))
        ## Total Likelihood
        LL <- d_LL + s_LL

        ### Log-Priors
        lambda.prior <- sum( dhalfcauchy(lambda, 1, log=T) )
        gamma.prior  <- sum( dhalfcauchy(gamma, 1, log=T) )
        sigma.prior  <- sum( dhalfcauchy(sigma, 1, log=T) )
        alpha.prior  <- sum( Data$density(solve(L_hat)[lower.tri(L_hat)], 0, 1/lambda, log=T) )
        theta.prior  <- sum( dhalfcauchy(theta, 1e2, log=T) )
        delta.prior  <- sum( dnorm(delta, 0, 1, log=T) )
        Lpp <- lambda.prior + gamma.prior + delta.prior + sigma.prior + alpha.prior + theta.prior
        
        ### Estimates
        yhat <- c(MASS::mvrnorm(Data$N, rep(0, Data$V), R_hat))
        
        ### Log-Posterior
        LP <- Lpp + LL
        
        ### Output
        Monitor=c(R_hat[lower.tri(R_hat)], Rho_hat, kappa)
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
        theta  <- exp(parm[Data$pos.theta])
        alpha  <- parm[Data$pos.alpha]
        delta  <- parm[Data$pos.delta]
        
        ### Log-Likelihood
        ## Data
        C_hat <- diag(Data$V)
        C_hat[lower.tri(C_hat)] <- alpha
        norms <- apply(C_hat, 1, euclidean)
        L_hat <- t(t(C_hat) %*% diag(1/norms))
        R_hat <- L_hat %*% t(L_hat)
        taus <- lapply(unique(Data$id.delta), function(g) delta[which(Data$id.delta == g)])
        d_LL <- sum(dpoly(Data$X, R=R_hat, taus=taus)$loglik[lower.tri(R_hat)] *
                      {Data$N + Data$n_par + Data$n_thr})
        ## Sparsity
        Rho_hat <- pcor(R_hat)[lower.tri(R_hat)]
        Rho_tra <- tanh(theta*atanh(Rho_hat))
        kappa   <- ptrunc(atanh(abs(Rho_tra)), "norm", 0, Inf, mean=0, sd=sigma)
        pi <- 1 - kappa
        s_LL <- sum(log( {pi * dnorm(atanh(Rho_hat), mean=0, sd=gamma)} + {kappa*kappa}))
        ## Total Likelihood
        LL <- d_LL + s_LL
        
        ### Log-Priors
        lambda.prior <- sum( dhalfcauchy(lambda, 1, log=T) )
        gamma.prior  <- sum( dhalfcauchy(gamma, 1, log=T) )
        tau.prior    <- sum( dhalfcauchy(tau, 1, log=T) )
        sigma.prior  <- sum( dhalfcauchy(sigma, 1, log=T) )
        alpha.prior  <- sum( Data$density(solve(L_hat)[lower.tri(L_hat)], 0, 1/lambda, log=T) )
        theta.prior  <- sum( dhalfcauchy(theta, 1e2, log=T) )
        delta.prior  <- sum( dnorm(delta, 0, 1, log=T) )
        Lpp <- lambda.prior + gamma.prior + tau.prior + delta.prior + sigma.prior + alpha.prior + theta.prior
        
        ### Estimates
        yhat <- c(MASS::mvrnorm(Data$N, rep(0, Data$V), R_hat))
        
        ### Log-Posterior
        LP <- Lpp + LL
        
        ### Output
        Monitor=c(R_hat[lower.tri(R_hat)], Rho_hat, kappa)
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