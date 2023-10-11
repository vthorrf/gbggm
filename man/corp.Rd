\name{corp}
\alias{corp}
\title{Correlation matrix from partial correlation matrix}
\description{
Takes a partial correlation matrix as input and computes the corresponding correlation matrix.
}
\usage{
corp(Rho)
}
\arguments{
   \item{Rho}{A partial correlation matrix.}
}
\value{
A correlation matrix.
}
\examples{##### Test with a toy simulation
## Load the package
require(gbggm)

### Data====
# Seed for reproducibility
set.seed(3)
# Sample size
N <- 200
# Fixed correlations
b1 <- .3; b2 <- .7
# Correlation matrix implied by the chain DAG: x1 -> x2 -> x3
Sigma <- matrix(c(    1, b1, b1*b2,
                     b1,  1,    b2,
                  b1*b2, b2,     1), byrow=T, nrow=3, ncol=3)
# The implied partial correlation matrix
pTrue <- pcor(Sigma)
# Recover the correlation matrix
round(gbggm::corp(pTrue), .Machine$double.eps) == round(Sigma, .Machine$double.eps)
}
