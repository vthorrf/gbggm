dataListFUN <- function(data, V, N, reg, cor) {
  if(cor == "pearson") {
    X <- scale(as.matrix(data))
  } else if(cor == "spearman") {
    X <- scale(apply(as.matrix(data),2,rank))
  } else if(cor == "poly") {
    X <- as.matrix(data)
    if(min(X) > 0) X <- X - min(X)
    if(min(X) < 0) stop("Minimum value in 'data' should be 0, but there is at least one value below 0.")
  }
  if(cor %in% c("pearson", "spearman")) {
    if(reg %in% c("normal", "laplace", "logistic", "cauchy", "hypersec")) {
      # Parameter names
      n_par <- {{V*{V-1}}/2}+1
      parm.names <- c("lambda",
                      paste("alpha",sprintf(paste0("%0",nchar(V)+1,"d"),1:{n_par-1}),sep="_"))
      # Parameters positions
      pos.lambda <- grep("lambda", parm.names)
      pos.alpha  <- grep("alpha" , parm.names)
      # Probability Generating Function
      PGF <- function(Data) {
        lambda <- qnorm(pgamma(1,1e-2,1e-2))
        alpha  <- rnorm(Data$n_par - 1)
        return(c(lambda, alpha))
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
      } else if(reg=="hypersec") {
        density <- dhypersec
      } else stop("Unknown regularization prior!")
      # Datalist
      Data <- list( X=X, V=V, N=N, n_par=n_par,
                    parm.names=parm.names, pos.lambda=pos.lambda,
                    pos.alpha=pos.alpha, density=density, PGF=PGF )
    } else if(reg %in% c("t", "lomax", "NEG")) {
      # Parameter names
      n_par <- {{V*{V-1}}/2}+2
      parm.names <- c("lambda", "tau",
                      paste("alpha",sprintf(paste0("%0",nchar(V)+1,"d"),1:{n_par-2}),sep="_"))
      # Parameters positions
      pos.lambda <- grep("lambda", parm.names)
      pos.tau    <- grep("tau"   , parm.names)
      pos.alpha  <- grep("alpha" , parm.names)
      # Probability Generating Function
      PGF <- function(Data) {
        lambda <- qnorm(pgamma(1,1e-2,1e-2))
        tau    <- rnorm(1)
        alpha  <- rnorm(Data$n_par - 2)
        return(c(lambda, tau, alpha))
      }
      # Density method
      if(reg=="t") {
        density <- d3t
      } else if(reg=="lomax") {
        density <- dlomax
      } else if(reg=="NEG") {
        density <- dNEG
      } else stop("Unknown regularization prior!")
      # Datalist
      Data <- list( X=X, V=V, N=N, n_par=n_par, parm.names=parm.names,
                    pos.lambda=pos.lambda, pos.tau=pos.tau, pos.alpha=pos.alpha,
                    density=density, PGF=PGF )
      
    } else {
      stop("Unknown regularization prior!")
    }
  } else if(cor == "poly") {
    # Threshold parameters
    thresholds <- sapply(1:ncol(X), function(g) length(unique(X[,g]))) - 1
    if(any(thresholds > 9)) warning("At least one variable has 10 categories or more.")
    varPt <- cbind(1:ncol(X), thresholds)
    n_thr <- sum(thresholds)
    # Choose regularization method
    if(reg %in% c("normal", "laplace", "logistic", "cauchy", "hypersec")) {
      # Parameter names
      n_par <- {{V*{V-1}}/2}+1
      parm.names <- c("lambda",
                      paste("alpha",sprintf(paste0("%0",nchar(V)+1,"d"),1:{n_par-1}),sep="_"),
                      unlist(sapply(1:nrow(varPt), function(g) {
                        sapply(1:varPt[g,2], function(r) {
                          paste0("delta_var",sprintf(paste0("%0",nchar(V)+1,"d"),g),"_",r)
                          
                        })
                      })))
      id.delta <- unlist(lapply(strsplit(parm.names[grep("delta",parm.names)],"_"), function(g) g[2]))
      # Parameters positions
      pos.lambda <- grep("lambda", parm.names)
      pos.alpha  <- grep("alpha" , parm.names)
      pos.delta  <- grep("delta" , parm.names)
      # Probability Generating Function
      PGF <- function(Data) {
        lambda <- qnorm(pgamma(1,1e-2,1e-2))
        alpha  <- rnorm(Data$n_par - 1)
        delta  <- rnorm(Data$n_thr) 
        return(c(lambda, alpha, delta))
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
      } else if(reg=="hypersec") {
        density <- dhypersec
      } else stop("Unknown regularization prior!")
      # Datalist
      Data <- list( X=X, V=V, N=N, n_par=n_par, n_thr=n_thr, parm.names=parm.names,
                    pos.lambda=pos.lambda, pos.alpha=pos.alpha,
                    pos.delta=pos.delta, density=density, PGF=PGF, id.delta=id.delta )
    } else if(reg %in% c("t", "lomax", "NEG")) {
      # Parameter names
      n_par <- {{V*{V-1}}/2}+2
      parm.names <- c("lambda", "tau",
                      paste("alpha",sprintf(paste0("%0",nchar(V)+1,"d"),1:{n_par-2}),sep="_"),
                      unlist(sapply(1:nrow(varPt), function(g) {
                        sapply(1:varPt[g,2], function(r) {
                          paste0("delta_var",sprintf(paste0("%0",nchar(V)+1,"d"),g),"_",r)
                        })
                      })))
      id.delta <- unlist(lapply(strsplit(parm.names[grep("delta",parm.names)],"_"), function(g) g[2]))
      # Parameters positions
      pos.lambda <- grep("lambda", parm.names)
      pos.tau    <- grep("tau"   , parm.names)
      pos.alpha  <- grep("alpha" , parm.names)
      pos.delta  <- grep("delta" , parm.names)
      # Probability Generating Function
      PGF <- function(Data) {
        lambda <- qnorm(pgamma(1,1e-2,1e-2))
        tau    <- rnorm(1)
        alpha  <- rnorm(Data$n_par - 2)
        delta  <- rnorm(Data$n_thr)
        return(c(lambda, tau, alpha, delta))
      }
      # Density method
      if(reg=="t") {
        density <- d3t
      } else if(reg=="lomax") {
        density <- dlomax
      } else if(reg=="NEG") {
        density <- dNEG
      } else stop("Unknown regularization prior!")
      # Datalist
      Data <- list( X=X, V=V, N=N, n_par=n_par, n_thr=n_thr, parm.names=parm.names,
                    pos.lambda=pos.lambda, pos.tau=pos.tau,
                    pos.alpha=pos.alpha, pos.delta=pos.delta,
                    density=density, PGF=PGF, id.delta=id.delta )
      
    } else {
      stop("Unknown regularization prior!")
    }
  } else {
    stop("Unkown correlation method")
  }
  return(Data)
}