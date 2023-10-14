paracyl_weight <- function(w, x, mu, k, theta) {
  {w^{2*k}} * exp({{-.5}*{w^2}}-{{abs(x-mu)/theta}*w})
}