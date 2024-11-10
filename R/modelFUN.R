modelFUN <- function(reg, cor, sparse) {
  if(reg == "normal") {
    Model <- modelNorm(cor=cor, sparse=sparse)
  } else if(reg == "laplace") {
    Model <- modelLaplace(cor=cor, sparse=sparse)
  } else if(reg == "logistic") {
    Model <- modelLogistic(cor=cor, sparse=sparse)
  } else if(reg == "cauchy") {
    Model <- modelCauchy(cor=cor, sparse=sparse)
  } else if(reg == "hypersec") {
    Model <- modelHypersec(cor=cor, sparse=sparse)
  } else if(reg == "t") {
    Model <- modelT(cor=cor, sparse=sparse)
  } else if(reg == "lomax") {
    Model <- modelLomax(cor=cor, sparse=sparse)
  } else if(reg == "kaniadakis") {
    Model <- modelKaniadakis(cor=cor, sparse=sparse)
  } else if(reg == "NEG") {
    Model <- modelNEG(cor=cor, sparse=sparse)
  } else {
    stop("Unkown regularization method!")
  }
  return(Model)
}