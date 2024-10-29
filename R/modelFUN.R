modelFUN <- function(reg, cor) {
  if(reg == "normal") {
    Model <- modelNorm(cor=cor)
  } else if(reg == "laplace") {
    Model <- modelLaplace(cor=cor)
  } else if(reg == "logistic") {
    Model <- modelLogistic(cor=cor)
  } else if(reg == "cauchy") {
    Model <- modelCauchy(cor=cor)
  } else if(reg == "hypersec") {
    Model <- modelHypersec(cor=cor)
  } else if(reg == "t") {
    Model <- modelT(cor=cor)
  } else if(reg == "lomax") {
    Model <- modelLomax(cor=cor)
  } else if(reg == "kaniadakis") {
    Model <- modelKaniadakis(cor=cor)
  } else if(reg == "NEG") {
    Model <- modelNEG(cor=cor)
  } else {
    stop("Unkown regularization method!")
  }
  return(Model)
}