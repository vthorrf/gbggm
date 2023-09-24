dataSparseFUN <- function(X, V, N, reg, cor) {
  if(cor == "pearson") {
    if(reg %in% c("normal", "laplace", "logistic", "cauchy")) {
      # Parameter names
      n_par <- {{V*{V-1}}/2}+2
      parm.names <- c("lambda", "gamma",
                      paste("alpha",sprintf(paste0("%0",nchar(V)+1,"d"),1:{n_par-2}),sep="_"))
      # Parameters positions
      pos.lambda <- grep("lambda", parm.names)
      pos.gamma  <- grep("gamma" , parm.names)
      pos.alpha  <- grep("alpha" , parm.names)
      # Probability Generating Function
      PGF <- function(Data) {
        lambda <- qnorm(pgamma(1,1e-2,1e-2))
        gamma  <- rnorm(1)
        alpha  <- rnorm(Data$n_par - 2)
        return(c(lambda, gamma, alpha))
      }
      # Density method
      if(reg=="normal") {
        density <- dnorm
      } else if(reg=="laplace") {
        density <- dlaplace
      } else if(reg=="logistic") {
        density <- dlogis
      } else if(reg=="cauchy") {
        density <- dcauchy
      } else stop("Unknown regularization prior!")
      # Datalist
      Data <- list( X=scale(X), V=V, N=N, n_par=n_par, parm.names=parm.names, pos.gamma=pos.gamma,
                    pos.lambda=pos.lambda, pos.alpha=pos.alpha, density=density, PGF=PGF,
                    DHN=dhalfnorm, DHC=dhalfcauchy, euclidean=euclidean )
    } else if(reg %in% c("t", "lomax")) {
      # Parameter names
      n_par <- {{V*{V-1}}/2}+3
      parm.names <- c("lambda", "gamma", "tau",
                      paste("alpha",sprintf(paste0("%0",nchar(V)+1,"d"),1:{n_par-3}),sep="_"))
      # Parameters positions
      pos.lambda <- grep("lambda", parm.names)
      pos.gamma  <- grep("gamma" , parm.names)
      pos.tau    <- grep("tau"   , parm.names)
      pos.alpha  <- grep("alpha" , parm.names)
      # Probability Generating Function
      PGF <- function(Data) {
        lambda <- qnorm(pgamma(1,1e-2,1e-2))
        gamma  <- rnorm(1)
        tau    <- rnorm(1)
        alpha  <- rnorm(Data$n_par - 3)
        return(c(lambda, gamma, tau, alpha))
      }
      # Density method
      if(reg=="t") {
        density <- d3t
      } else if(reg=="lomax") {
        density <- dlomax
      } else stop("Unknown regularization prior!")
      # Datalist
      Data <- list( X=scale(X), V=V, N=N, n_par=n_par, parm.names=parm.names, pos.gamma=pos.gamma,
                    pos.lambda=pos.lambda, pos.tau=pos.tau, pos.alpha=pos.alpha, density=density,
                    PGF=PGF, DHN=dhalfnorm, DHC=dhalfcauchy, euclidean=euclidean )
      
    } else {
      stop("Unknown regularization prior!")
    }
  } else if(cor == "spearman") {
    if(reg %in% c("normal", "laplace", "logistic", "cauchy")) {
      # Parameter names
      n_par <- {{V*{V-1}}/2}+2
      parm.names <- c("lambda", "gamma",
                      paste("alpha",sprintf(paste0("%0",nchar(V)+1,"d"),1:{n_par-2}),sep="_"))
      # Parameters positions
      pos.lambda <- grep("lambda", parm.names)
      pos.gamma  <- grep("gamma" , parm.names)
      pos.alpha  <- grep("alpha" , parm.names)
      # Probability Generating Function
      PGF <- function(Data) {
        lambda <- qnorm(pgamma(1,1e-2,1e-2))
        gamma  <- rnorm(1)
        alpha  <- rnorm(Data$n_par - 2)
        return(c(lambda, gamma, alpha))
      }
      # Density method
      if(reg=="normal") {
        density <- dnorm
      } else if(reg=="laplace") {
        density <- dlaplace
      } else if(reg=="logistic") {
        density <- dlogis
      } else if(reg=="cauchy") {
        density <- dcauchy
      } else stop("Unknown regularization prior!")
      # Datalist
      Data <- list( X=scale(apply(X,2,rank)), V=V, N=N, n_par=n_par, parm.names=parm.names,
                    pos.gamma=pos.gamma, pos.lambda=pos.lambda, pos.alpha=pos.alpha,
                    density=density, PGF=PGF, DHN=dhalfnorm, DHC=dhalfcauchy, euclidean=euclidean )
    } else if(reg %in% c("t", "lomax")) {
      # Parameter names
      n_par <- {{V*{V-1}}/2}+3
      parm.names <- c("lambda", "gamma", "tau",
                      paste("alpha",sprintf(paste0("%0",nchar(V)+1,"d"),1:{n_par-3}),sep="_"))
      # Parameters positions
      pos.lambda <- grep("lambda", parm.names)
      pos.gamma  <- grep("gamma" , parm.names)
      pos.tau    <- grep("tau"   , parm.names)
      pos.alpha  <- grep("alpha" , parm.names)
      # Probability Generating Function
      PGF <- function(Data) {
        lambda <- qnorm(pgamma(1,1e-2,1e-2))
        gamma  <- rnorm(1)
        tau    <- rnorm(1)
        alpha  <- rnorm(Data$n_par - 3)
        return(c(lambda, gamma, tau, alpha))
      }
      # Density method
      if(reg=="t") {
        density <- d3t
      } else if(reg=="lomax") {
        density <- dlomax
      } else stop("Unknown regularization prior!")
      # Datalist
      Data <- list( X=scale(apply(X,2,rank)), V=V, N=N, n_par=n_par, parm.names=parm.names,
                    pos.gamma=pos.gamma, pos.lambda=pos.lambda, pos.tau=pos.tau,
                    pos.alpha=pos.alpha, density=density, PGF=PGF, DHN=dhalfnorm,
                    DHC=dhalfcauchy, euclidean=euclidean )
    } else {
      stop("Unknown regularization prior!")
    }
  } else if(cor == "poly") {
    # Threshold parameters
    thresholds <- sapply(1:ncol(X), function(g) length(unique(X[,g]))) - 1
    if(any(thresholds > 9)) stop("At least one variable has 10 categories or more.")
    varPt <- cbind(1:ncol(X), thresholds)
    n_thr <- sum(thresholds)
    # Choose regularization method
    if(reg %in% c("normal", "laplace", "logistic", "cauchy")) {
      # Parameter names
      n_par <- {{V*{V-1}}/2}+2
      parm.names <- c("lambda", "gamma",
                      paste("alpha",sprintf(paste0("%0",nchar(V)+1,"d"),1:{n_par-2}),sep="_"),
                      unlist(sapply(1:nrow(varPt), function(g) {
                        sapply(1:varPt[g,2], function(r) {
                          paste0("delta_var",sprintf(paste0("%0",nchar(V)+1,"d"),g),"_",r)
                          
                        })
                      })))
      id.delta <- unlist(lapply(strsplit(parm.names[grep("delta",parm.names)],"_"), function(g) g[2]))
      # Parameters positions
      pos.lambda <- grep("lambda", parm.names)
      pos.gamma  <- grep("gamma" , parm.names)
      pos.alpha  <- grep("alpha" , parm.names)
      pos.delta  <- grep("delta" , parm.names)
      # Probability Generating Function
      PGF <- function(Data) {
        lambda <- qnorm(pgamma(1,1e-2,1e-2))
        gamma  <- rnorm(1)
        alpha  <- rnorm(Data$n_par - 2)
        delta  <- rnorm(Data$n_thr) 
        return(c(lambda, gamma, alpha, delta))
      }
      # Density method
      if(reg=="normal") {
        density <- dnorm
      } else if(reg=="laplace") {
        density <- dlaplace
      } else if(reg=="logistic") {
        density <- dlogis
      } else if(reg=="cauchy") {
        density <- dcauchy
      } else stop("Unknown regularization prior!")
      # Datalist
      if(min(X) > 0) X <- X - min(X)
      if(min(X) < 0) stop("Minimum value should be 0, but a there is at least value below 0.")
      Data <- list( X=X, V=V, N=N, n_par=n_par, n_thr=n_thr, parm.names=parm.names,
                    pos.gamma=pos.gamma, pos.lambda=pos.lambda, pos.alpha=pos.alpha,
                    pos.delta=pos.delta, density=density, PGF=PGF, DHN=dhalfnorm,
                    DHC=dhalfcauchy, euclidean=euclidean, id.delta=id.delta )
    } else if(reg %in% c("t", "lomax")) {
      # Parameter names
      n_par <- {{V*{V-1}}/2}+3
      parm.names <- c("lambda", "gamma", "tau",
                      paste("alpha",sprintf(paste0("%0",nchar(V)+1,"d"),1:{n_par-3}),sep="_"),
                      unlist(sapply(1:nrow(varPt), function(g) {
                        sapply(1:varPt[g,2], function(r) {
                          paste0("delta_var",sprintf(paste0("%0",nchar(V)+1,"d"),g),"_",r)
                        })
                      })))
      id.delta <- unlist(lapply(strsplit(parm.names[grep("delta",parm.names)],"_"), function(g) g[2]))
      # Parameters positions
      pos.lambda <- grep("lambda", parm.names)
      pos.gamma  <- grep("gamma" , parm.names)
      pos.tau    <- grep("tau"   , parm.names)
      pos.alpha  <- grep("alpha" , parm.names)
      pos.delta  <- grep("delta" , parm.names)
      # Probability Generating Function
      PGF <- function(Data) {
        lambda <- qnorm(pgamma(1,1e-2,1e-2))
        gamma  <- rnorm(1)
        tau    <- rnorm(1)
        alpha  <- rnorm(Data$n_par - 3)
        delta  <- rnorm(Data$n_thr)
        return(c(lambda, gamma, tau, alpha, delta))
      }
      # Density method
      if(reg=="t") {
        density <- d3t
      } else if(reg=="lomax") {
        density <- dlomax
      } else stop("Unknown regularization prior!")
      # Datalist
      if(min(X) > 0) X <- X - min(X)
      if(min(X) < 0) stop("Minimum value should be 0, but a there is at least value below 0.")
      Data <- list( X=X, V=V, N=N, n_par=n_par, n_thr=n_thr, parm.names=parm.names,
                    pos.gamma=pos.gamma, pos.lambda=pos.lambda, pos.tau=pos.tau,
                    pos.alpha=pos.alpha, pos.delta=pos.delta, density=density,
                    PGF=PGF, DHN=dhalfnorm, DHC=dhalfcauchy, euclidean=euclidean,
                    id.delta=id.delta )
      
    } else {
      stop("Unknown regularization prior!")
    }
  } else {
    stop("Unkown correlation method")
  }
  return(Data)
}