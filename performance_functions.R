# This file contains functions to estimate the simulation study performance measures

# bias estimate
bias <- function(theta.hat, theta) {
  # theta.hat: point estimates; theta: truth
  nsim <- length(theta.hat)
  est <- sum(theta.hat)/nsim - theta
  return(est)
}

# Monte Carlo standard error of bias estimate
bias.mcse <- function(theta.hat) {
  # theta.hat: point estimates
  nsim <- length(theta.hat)
  tmp <- sum((theta.hat - mean(theta.hat))^2)
  mcse <- sqrt(1/(nsim*(nsim-1))*tmp)
  return(mcse)
}

# Empirical standard error
empse <- function(theta.hat) {
  # theta.hat: point estimates
  nsim <- length(theta.hat)
  tmp <- sum((theta.hat - mean(theta.hat))^2)
  est <- sqrt(tmp/(nsim-1))
  return(est)
}

# Monte Carlo standard error of empirical standard error estimate
empse.mcse <- function(empse, nsim) {
  # empse: empirical standard error estimate; nsim: number of simulations
  mcse <- empse/(sqrt(2*(nsim-1)))
  return(mcse)
}

# mean square error estimate
mse <- function(theta.hat, theta) {
  # theta.hat: point estimates; theta: truth
  nsim <- length(theta.hat)
  est <- sum((theta.hat-theta)^2)/nsim
  return(est)
}

# Monte Carlo standard error of mean square error estimate
mse.mcse <- function(theta.hat, theta) {
  # theta.hat: point estimates; theta: truth
  nsim <- length(theta.hat)
  tmp <- (theta.hat-theta)^2
  mse.est <- sum(tmp)/nsim
  mcse <- sqrt(sum((tmp - mse.est)^2)/(nsim*(nsim-1)))
  return(mcse)
}

# coverage estimate
coverage <- function(theta.hat.low, theta.hat.upp, theta) {
  # theta.hat.low: lower bound of interval estimate
  # theta.hat.upp; upper bound of interval estimate
  # theta: truth
  nsim <- length(theta.hat.low)
  est <- sum(ifelse(theta>=theta.hat.low & theta<=theta.hat.upp,1,0))/nsim
  return(est)
}

# Monte Carlo standard error of coverage estimate
coverage.mcse <- function(coverage, nsim) {
  # coverage: coverage estimate; nsim: number of simulations 
  mcse <- sqrt((coverage*(1-coverage))/nsim)
  return(mcse)
}
