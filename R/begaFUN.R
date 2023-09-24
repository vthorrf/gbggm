begaFUN <- function(reg, cor) {
  if(cor %in% c("pearson","spearman")) {
    if(reg %in% c("normal", "laplace", "logistic", "cauchy")) {
      Model <- function(parm, Data){
        
        ## Prior parameters
        lambda <- exp(parm[Data$pos.lambda])
        alpha  <- parm[Data$pos.alpha]
        
        ### Log-Priors
        lambda.prior <- sum( Data$DHC(lambda, 0, 1, log=T) )
        alpha.prior  <- sum( Data$density(alpha, 0, 1/lambda, log=T) )
        Lpp <- lambda.prior + alpha.prior
        
        ### Log-Likelihood
        C_hat <- diag(Data$V)
        C_hat[lower.tri(C_hat)] <- alpha
        norms <- apply(C_hat, 1, Data$euclidean)
        L_hat <- t(t(C_hat) %*% diag(1/norms))
        R_hat <- L_hat %*% t(L_hat)
        LL <- sum( dmnorm(Data$X, sigma=R_hat, log=T) )
        
        ### Clustering of (dis)connected nodes
        net <- pcor(R_hat)
        Rho_hat <- net[lower.tri(net)]
        
        ### Louvain likelihood
        graph <- suppressWarnings(suppressMessages(graph_from_adjacency_matrix(abs(net), mode="undirected",
                                  weighted=T, add.colnames = FALSE)))
        communities <- cluster_louvain(graph, weights=NULL, resolution=.5)
        wc <- communities$memberships[which.max(communities$modularity),]
        mod <- communities$modularity[which.max(communities$modularity)]
        comLL <- dbeta({mod + 1} / 2, 2, 1, log=T)
        
        ### Estimates
        yhat <- c(MASS::mvrnorm(Data$N, rep(0, Data$V), R_hat))
        
        ### Log-Posterior
        LP <- Lpp + LL + comLL
        
        ### Output
        Monitor=c(R_hat[lower.tri(R_hat)], Rho_hat, wc)
        Modelout <- list(LP=LP, Dev=-2*LL, Monitor=Monitor, parm=parm, yhat=yhat)
        return(Modelout)
      }
    } else if (reg %in% c("t", "lomax")) {
      Model <- function(parm, Data){
        
        ## Prior parameters
        lambda <- exp(parm[Data$pos.lambda])
        tau    <- exp(parm[Data$pos.tau])
        alpha  <- parm[Data$pos.alpha]
        
        ### Log-Priors
        lambda.prior <- sum( Data$DHC(lambda, 0, 1, log=T) )
        tau.prior    <- sum( Data$DHC(tau, 0, 1, log=T) )
        alpha.prior  <- sum( Data$density(alpha, 0, 1/lambda, tau, log=T) )
        Lpp <- lambda.prior + alpha.prior + tau.prior
        
        ### Log-Likelihood
        C_hat <- diag(Data$V)
        C_hat[lower.tri(C_hat)] <- alpha
        norms <- apply(C_hat, 1, Data$euclidean)
        L_hat <- t(t(C_hat) %*% diag(1/norms))
        R_hat <- L_hat %*% t(L_hat)
        LL <- sum( dmnorm(Data$X, sigma=R_hat, log=T) )
        
        ### Clustering of (dis)connected nodes
        net <- pcor(R_hat)
        Rho_hat <- net[lower.tri(net)]
        
        ### Louvain likelihood
        graph <- suppressWarnings(suppressMessages(graph_from_adjacency_matrix(abs(net), mode="undirected",
                                                                               weighted=T, add.colnames = FALSE)))
        communities <- cluster_louvain(graph, weights=NULL, resolution=.5)
        wc <- communities$memberships[which.max(communities$modularity),]
        mod <- communities$modularity[which.max(communities$modularity)]
        comLL <- dbeta({mod + 1} / 2, 2, 1, log=T)
        
        ### Estimates
        yhat <- c(MASS::mvrnorm(Data$N, rep(0, Data$V), R_hat))
        
        ### Log-Posterior
        LP <- Lpp + LL + comLL
        
        ### Output
        Monitor=c(R_hat[lower.tri(R_hat)], Rho_hat, wc)
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
        alpha  <- parm[Data$pos.alpha]
        delta  <- parm[Data$pos.delta]
        
        ### Log-Priors
        lambda.prior <- sum( Data$DHC(lambda, 0, 1, log=T) )
        alpha.prior  <- sum( Data$density(alpha, 0, 1/lambda, log=T) )
        delta.prior  <- sum( dnorm(delta, 0, 1, log=T) )
        Lpp <- lambda.prior + alpha.prior + delta.prior
        
        ### Log-Likelihood
        C_hat <- diag(Data$V)
        C_hat[lower.tri(C_hat)] <- alpha
        norms <- apply(C_hat, 1, Data$euclidean)
        L_hat <- t(t(C_hat) %*% diag(1/norms))
        R_hat <- L_hat %*% t(L_hat)
        taus <- lapply(unique(Data$id.delta), function(g) delta[which(Data$id.delta == g)])
        LL <- sum(dpoly(Data$X, R=R_hat, taus=taus)$loglik[lower.tri(R_hat)] * {Data$N + Data$n_par + Data$n_thr})
        
        ### Clustering of (dis)connected nodes
        net <- pcor(R_hat)
        Rho_hat <- net[lower.tri(net)]
        
        ### Louvain likelihood
        graph <- suppressWarnings(suppressMessages(graph_from_adjacency_matrix(abs(net), mode="undirected",
                                                                               weighted=T, add.colnames = FALSE)))
        communities <- cluster_louvain(graph, weights=NULL, resolution=.5)
        wc <- communities$memberships[which.max(communities$modularity),]
        mod <- communities$modularity[which.max(communities$modularity)]
        comLL <- dbeta({mod + 1} / 2, 2, 1, log=T)
        
        ### Estimates
        yhat <- c(MASS::mvrnorm(Data$N, rep(0, Data$V), R_hat))
        
        ### Log-Posterior
        LP <- Lpp + LL + comLL
        
        ### Output
        Monitor=c(R_hat[lower.tri(R_hat)], Rho_hat, wc)
        Modelout <- list(LP=LP, Dev=-2*LL, Monitor=Monitor, parm=parm, yhat=yhat)
        return(Modelout)
      }
    } else if (reg %in% c("t", "lomax")) {
      Model <- function(parm, Data){
        
        ## Prior parameters
        lambda <- exp(parm[Data$pos.lambda])
        tau    <- exp(parm[Data$pos.tau])
        alpha  <- parm[Data$pos.alpha]
        delta  <- parm[Data$pos.delta]
        
        ### Log-Priors
        lambda.prior <- sum( Data$DHC(lambda, 0, 1, log=T) )
        tau.prior    <- sum( Data$DHC(tau, 0, 1, log=T) )
        alpha.prior  <- sum( Data$density(alpha, 0, 1/lambda, tau, log=T) )
        delta.prior  <- sum( dnorm(delta, 0, 1, log=T) )
        Lpp <- lambda.prior + alpha.prior + tau.prior + delta.prior
        
        ### Log-Likelihood
        C_hat <- diag(Data$V)
        C_hat[lower.tri(C_hat)] <- alpha
        norms <- apply(C_hat, 1, Data$euclidean)
        L_hat <- t(t(C_hat) %*% diag(1/norms))
        R_hat <- L_hat %*% t(L_hat)
        taus <- lapply(unique(Data$id.delta), function(g) delta[which(Data$id.delta == g)])
        LL <- sum(dpoly(Data$X, R=R_hat, taus=taus)$loglik[lower.tri(R_hat)] * {Data$N + Data$n_par + Data$n_thr})
        
        ### Clustering of (dis)connected nodes
        net <- pcor(R_hat)
        Rho_hat <- net[lower.tri(net)]
        
        ### Louvain likelihood
        graph <- suppressWarnings(suppressMessages(graph_from_adjacency_matrix(abs(net), mode="undirected",
                                                                               weighted=T, add.colnames = FALSE)))
        communities <- cluster_louvain(graph, weights=NULL, resolution=.5)
        wc <- communities$memberships[which.max(communities$modularity),]
        mod <- communities$modularity[which.max(communities$modularity)]
        comLL <- dbeta({mod + 1} / 2, 2, 1, log=T)
        
        ### Estimates
        yhat <- c(MASS::mvrnorm(Data$N, rep(0, Data$V), R_hat))
        
        ### Log-Posterior
        LP <- Lpp + LL + comLL
        
        ### Output
        Monitor=c(R_hat[lower.tri(R_hat)], Rho_hat, wc)
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