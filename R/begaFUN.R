begaFUN <- function(reg, cor) {
  if(cor %in% c("pearson","spearman")) {
    if(reg %in% c("normal", "laplace", "logistic", "cauchy")) {
      Model <- function(parm, Data){
        
        ## Prior parameters
        lambda <- exp(parm[Data$pos.lambda])
        gamma  <- exp(parm[Data$pos.gamma])
        sigma  <- exp(parm[Data$pos.sigma])
        alpha  <- parm[Data$pos.alpha]
        
        ### Log-Priors
        lambda.prior <- sum( Data$DHC(lambda, 0, 1, log=T) )
        gamma.prior  <- sum( Data$DHN(gamma, 0, 1, log=T) )
        sigma.prior  <- sum( Data$DHC(sigma, 0, 1, log=T) )
        alpha.prior  <- sum( Data$density(alpha, 0, 1/lambda, log=T) )
        Lpp <- lambda.prior + alpha.prior + gamma.prior + sigma.prior
        
        ### Log-Likelihood
        C_hat <- diag(Data$V)
        C_hat[lower.tri(C_hat)] <- alpha
        norms <- apply(C_hat, 1, Data$euclidean)
        L_hat <- t(t(C_hat) %*% diag(1/norms))
        R_hat <- L_hat %*% t(L_hat)
        LL <- sum( dmnorm(Data$X, sigma=R_hat, log=T) )
        
        ### Clustering of (dis)connected nodes
        net       <- pcor(R_hat)
        Rho_hat   <- atanh(net[lower.tri(R_hat)])
        weight    <- pnorm(abs(Rho_hat), mean=-gamma, sd=sigma) - pnorm(abs(Rho_hat), mean=gamma, sd=sigma)
        clusterLL <- sum( dnorm(abs(Rho_hat), {1-weight}*gamma, sigma, log=T) )
        
        ### Louvain likelihood
        mu        <- {abs(Rho_hat) > {gamma/2}}*1
        net[lower.tri(R_hat)] <- mu
        net[upper.tri(net)] <- t(net)[upper.tri(net)]
        graph <- suppressWarnings(suppressMessages(graph_from_adjacency_matrix(abs(net), mode="undirected",
                                                   weighted=T, add.colnames = FALSE)))
        communities <- cluster_louvain(graph, weights=NULL, resolution=1)
        wc <- communities$membership
        mod <- communities$modularity[which.max(communities$modularity)]
        comLL <- dbeta({mod + 1} / 2, 2, 2, log=T)
        
        ### Estimates
        yhat <- c(MASS::mvrnorm(Data$N, rep(0, Data$V), R_hat))
        
        ### Log-Posterior
        LP <- Lpp + LL + clusterLL + comLL
        
        ### Output
        Monitor=c(R_hat[lower.tri(R_hat)], tanh(Rho_hat), 1-weight, wc)
        Modelout <- list(LP=LP, Dev=-2*LL, Monitor=Monitor, parm=parm, yhat=yhat)
        return(Modelout)
      }
    } else if (reg %in% c("t", "lomax")) {
      Model <- function(parm, Data){
        
        ## Prior parameters
        lambda <- exp(parm[Data$pos.lambda])
        gamma  <- exp(parm[Data$pos.gamma])
        tau    <- exp(parm[Data$pos.tau])
        sigma  <- exp(parm[Data$pos.sigma])
        alpha  <- parm[Data$pos.alpha]
        
        ### Log-Priors
        lambda.prior <- sum( Data$DHC(lambda, 0, 1, log=T) )
        gamma.prior  <- sum( Data$DHN(gamma, 0, 1, log=T) )
        tau.prior    <- sum( Data$DHC(tau, 0, 1, log=T) )
        sigma.prior  <- sum( Data$DHC(sigma, 0, 1, log=T) )
        alpha.prior  <- sum( Data$density(alpha, 0, 1/lambda, tau, log=T) )
        Lpp <- lambda.prior + alpha.prior + gamma.prior + tau.prior + sigma.prior
        
        ### Log-Likelihood
        C_hat <- diag(Data$V)
        C_hat[lower.tri(C_hat)] <- alpha
        norms <- apply(C_hat, 1, Data$euclidean)
        L_hat <- t(t(C_hat) %*% diag(1/norms))
        R_hat <- L_hat %*% t(L_hat)
        LL <- sum( dmnorm(Data$X, sigma=R_hat, log=T) )
        
        ### Clustering of (dis)connected nodes
        net       <- pcor(R_hat)
        Rho_hat   <- atanh(net[lower.tri(R_hat)])
        weight    <- pnorm(abs(Rho_hat), mean=-gamma, sd=sigma) - pnorm(abs(Rho_hat), mean=gamma, sd=sigma)
        clusterLL <- sum( dnorm(abs(Rho_hat), {1-weight}*gamma, sigma, log=T) )
        
        ### Louvain likelihood
        mu        <- {abs(Rho_hat) > {gamma/2}}*1
        net[lower.tri(R_hat)] <- mu
        net[upper.tri(net)] <- t(net)[upper.tri(net)]
        graph <- suppressWarnings(suppressMessages(graph_from_adjacency_matrix(abs(net), mode="undirected",
                                                   weighted=T, add.colnames = FALSE)))
        communities <- cluster_louvain(graph, weights=NULL, resolution=1)
        wc <- communities$membership
        mod <- communities$modularity[which.max(communities$modularity)]
        comLL <- dbeta({mod + 1} / 2, 2, 2, log=T)
        
        ### Estimates
        yhat <- c(MASS::mvrnorm(Data$N, rep(0, Data$V), R_hat))
        
        ### Log-Posterior
        LP <- Lpp + LL + clusterLL + comLL
        
        ### Output
        Monitor=c(R_hat[lower.tri(R_hat)], tanh(Rho_hat), 1-weight)
        Modelout <- list(LP=LP, Dev=-2*LL, Monitor=Monitor, parm=parm, yhat=yhat)
        return(Modelout)
      }
    } else {
      stop("Unknown regularization prior!")
    }
  } else if(cor == "poly") {
    if(reg %in% c("normal", "laplace", "logistic", "cauchy")) {
      Model <- function(parm, Data){
        
        ## Prior parameters
        lambda <- exp(parm[Data$pos.lambda])
        gamma  <- exp(parm[Data$pos.gamma])
        sigma  <- exp(parm[Data$pos.sigma])
        alpha  <- parm[Data$pos.alpha]
        delta  <- parm[Data$pos.delta]
        
        ### Log-Priors
        lambda.prior <- sum( Data$DHC(lambda, 0, 1, log=T) )
        gamma.prior  <- sum( Data$DHN(gamma, 0, 1, log=T) )
        sigma.prior  <- sum( Data$DHC(sigma, 0, 1, log=T) )
        alpha.prior  <- sum( Data$density(alpha, 0, 1/lambda, log=T) )
        delta.prior  <- sum( dnorm(delta, 0, 1, log=T) )
        Lpp <- lambda.prior + alpha.prior + gamma.prior + delta.prior + sigma.prior
        
        ### Log-Likelihood
        C_hat <- diag(Data$V)
        C_hat[lower.tri(C_hat)] <- alpha
        norms <- apply(C_hat, 1, Data$euclidean)
        L_hat <- t(t(C_hat) %*% diag(1/norms))
        R_hat <- L_hat %*% t(L_hat)
        taus <- lapply(unique(Data$id.delta), function(g) delta[which(Data$id.delta == g)])
        LL <- sum(dpoly(Data$X, R=R_hat, taus=taus)$loglik[lower.tri(R_hat)] * {Data$N + Data$n_par + Data$n_thr})
        
        ### Clustering of (dis)connected nodes
        net       <- pcor(R_hat)
        Rho_hat   <- atanh(net[lower.tri(R_hat)])
        weight    <- pnorm(abs(Rho_hat), mean=-gamma, sd=sigma) - pnorm(abs(Rho_hat), mean=gamma, sd=sigma)
        clusterLL <- sum( dnorm(abs(Rho_hat), {1-weight}*gamma, sigma, log=T) )
        
        ### Louvain likelihood
        mu        <- {abs(Rho_hat) > {gamma/2}}*1
        net[lower.tri(R_hat)] <- mu
        net[upper.tri(net)] <- t(net)[upper.tri(net)]
        graph <- suppressWarnings(suppressMessages(graph_from_adjacency_matrix(abs(net), mode="undirected",
                                                   weighted=T, add.colnames = FALSE)))
        communities <- cluster_louvain(graph, weights=NULL, resolution=1)
        wc <- communities$membership
        mod <- communities$modularity[which.max(communities$modularity)]
        comLL <- dbeta({mod + 1} / 2, 2, 2, log=T)
        
        ### Estimates
        yhat <- c(MASS::mvrnorm(Data$N, rep(0, Data$V), R_hat))
        
        ### Log-Posterior
        LP <- Lpp + LL + clusterLL + comLL
        
        ### Output
        Monitor=c(R_hat[lower.tri(R_hat)], tanh(Rho_hat), 1-weight)
        Modelout <- list(LP=LP, Dev=-2*LL, Monitor=Monitor, parm=parm, yhat=yhat)
        return(Modelout)
      }
    } else if (reg %in% c("t", "lomax")) {
      Model <- function(parm, Data){
        
        ## Prior parameters
        lambda <- exp(parm[Data$pos.lambda])
        gamma  <- exp(parm[Data$pos.gamma])
        tau    <- exp(parm[Data$pos.tau])
        sigma  <- exp(parm[Data$pos.sigma])
        alpha  <- parm[Data$pos.alpha]
        delta  <- parm[Data$pos.delta]
        
        ### Log-Priors
        lambda.prior <- sum( Data$DHC(lambda, 0, 1, log=T) )
        gamma.prior  <- sum( Data$DHN(gamma, 0, 1, log=T) )
        tau.prior    <- sum( Data$DHC(tau, 0, 1, log=T) )
        sigma.prior  <- sum( Data$DHC(sigma, 0, 1, log=T) )
        alpha.prior  <- sum( Data$density(alpha, 0, 1/lambda, tau, log=T) )
        delta.prior  <- sum( dnorm(delta, 0, 1, log=T) )
        Lpp <- lambda.prior + alpha.prior + gamma.prior + tau.prior + delta.prior + sigma.prior
        
        ### Log-Likelihood
        C_hat <- diag(Data$V)
        C_hat[lower.tri(C_hat)] <- alpha
        norms <- apply(C_hat, 1, Data$euclidean)
        L_hat <- t(t(C_hat) %*% diag(1/norms))
        R_hat <- L_hat %*% t(L_hat)
        taus <- lapply(unique(Data$id.delta), function(g) delta[which(Data$id.delta == g)])
        LL <- sum(dpoly(Data$X, R=R_hat, taus=taus)$loglik[lower.tri(R_hat)] * {Data$N + Data$n_par + Data$n_thr})
        
        ### Clustering of (dis)connected nodes
        net       <- pcor(R_hat)
        Rho_hat   <- atanh(net[lower.tri(R_hat)])
        weight    <- pnorm(abs(Rho_hat), mean=-gamma, sd=sigma) - pnorm(abs(Rho_hat), mean=gamma, sd=sigma)
        clusterLL <- sum( dnorm(abs(Rho_hat), {1-weight}*gamma, sigma, log=T) )
        
        ### Louvain likelihood
        mu        <- {abs(Rho_hat) > {gamma/2}}*1
        net[lower.tri(R_hat)] <- mu
        net[upper.tri(net)] <- t(net)[upper.tri(net)]
        graph <- suppressWarnings(suppressMessages(graph_from_adjacency_matrix(abs(net), mode="undirected",
                                                   weighted=T, add.colnames = FALSE)))
        communities <- cluster_louvain(graph, weights=NULL, resolution=1)
        wc <- communities$membership
        mod <- communities$modularity[which.max(communities$modularity)]
        comLL <- dbeta({mod + 1} / 2, 2, 2, log=T)
        
        ### Estimates
        yhat <- c(MASS::mvrnorm(Data$N, rep(0, Data$V), R_hat))
        
        ### Log-Posterior
        LP <- Lpp + LL + clusterLL + comLL
        
        ### Output
        Monitor=c(R_hat[lower.tri(R_hat)], tanh(Rho_hat), 1-weight)
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