group_inclusion <- function(x, gamma, sigma) {
  1 - ptruncnorm(abs(x), 0, sigma) + ptruncnorm(abs(x), gamma, sigma)
}