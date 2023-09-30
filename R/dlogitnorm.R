# https://stats.stackexchange.com/questions/494924/what-is-the-pdf-of-the-ratio-of-two-random-log-normal-variables
dlogitnorm <- function(x, mean=0, scale=1, log=F) {
  if(log) {
    log({1/{scale * sqrt(2 * pi) * x * (1 - x)}} * 
        exp(-{log(x/{1 - x}) - mean}^2/{2 * scale^2}))
  } else {
    {1/{scale * sqrt(2 * pi) * x * (1 - x)}} * 
      exp(-{log(x/{1 - x}) - mean}^2/{2 * scale^2})
  }
}
