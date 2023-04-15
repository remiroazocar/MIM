# Example R code implementing the standardization methods on a simulated example

# for non-parametric bootstrap in standard parametric model-based standardization (G-computation)
library("boot") 
# for the synthesis stage in multiple imputation marginalization
library("rstanarm")

# set seed for reproducibility
set.seed(555) 

# setwd("C:/Users/gndtl/OneDrive - Bayer/Personal Data/New Papers/MIM/Example")

### Standard parametric model-based standardization (G-computation) ###

rm(list=ls())
ipd.index <- read.csv("index_IPD.csv") # load patient-level dataset of comparative index study
ipd.target <- read.csv("target_IPD.csv") # load target covariate individual-level dataset

# settings
resamples <- 1000 # number of resamples in non-parametric bootstrap

# function to be bootstrapped
gcomp.ml <- function(data, indices) {
  index.dat = data[indices,]
  # logistic outcome model fitted to index RCT subject-level data using maximum likelihood estimation
  outcome.model <- glm(y~(trt*x1)+(trt*x2), data=index.dat, family=binomial)
  # "counterfactual" target datasets ("target" covariates assumed fixed)
  target.t1 <- target.t0 <- ipd.target
  # intervene on assigned treatment
  target.t1$trt <- 1 # target dataset where everyone is assigned active treatment
  target.t0$trt <- 0 # target dataset where everyone is assigned placebo
  # predict potential individual-level event probabilities, conditional on treatment and covariates
  hat.mu.1.i <- predict(outcome.model, type="response", newdata=target.t1)
  hat.mu.0.i <- predict(outcome.model, type="response", newdata=target.t0)
  # marginal mean probability predictions under each treatment
  hat.mu.1 <- mean(hat.mu.1.i)
  hat.mu.0 <- mean(hat.mu.0.i)
  # transform from probability to linear predictor scale to estimate marginal log odds ratio
  hat.Delta <- log(hat.mu.1/(1-hat.mu.1)) - log(hat.mu.0/(1-hat.mu.0))  
  # hat.Delta <- qlogis(hat.mu.1) - qlogis(hat.mu.0)
  return(hat.Delta)
} 
# non-parametric bootstrap
boot.object <- boot::boot(data=ipd.index, statistic=gcomp.ml, R=resamples)
# bootstrap mean of marginal log odds ratio estimate
hat.Delta <- mean(boot.object$t)
# bootstrap variance of marginal log odds ratio estimate   
hat.var.Delta <- var(boot.object$t)
# Interval estimates derived from the relevant percentiles across the bootstrap resamples
conf.ints <- quantile(boot.object$t, probs = c(.025, .975), type = 6)
lci.Delta <- conf.ints[1] # lower interval bound
uci.Delta <- conf.ints[2] # upper interval bound

### Multiple imputation marginalization (MIM) ###

rm(list=ls())
ipd.index <- read.csv("index_IPD.csv") # load patient-level dataset of comparative index study
ipd.target <- read.csv("target_IPD.csv") # load target covariate individual-level dataset

# settings
M <- 1000 # number of syntheses/imputations used in analysis stage (high for low Monte Carlo error)
n.chains <- 2 # number of Markov chains for MCMC sampler in synthesis stage
warmup <- 2000 # number of discarded warmup/burn-in iterations per chain for MCMC sampler
iters <- 4000 # number of total iterations per chain for MCMC sampler (including warmup)

## SYNTHESIS STAGE ##
# first-stage logistic regression model fitted to index RCT using MCMC (Stan)
outcome.model <- stan_glm(y~(trt*x1)+(trt*x2), data=ipd.index, family=binomial, algorithm="sampling",
                          iter=iters, warmup=warmup, chains=n.chains, 
                          # thin to use M independent samples in analysis stage 
                          thin=(n.chains*(iters-warmup))/M) 
# create augmented target dataset
target.t1 <- target.t0 <- ipd.target
target.t1$trt <- 1 # assign active treatment in synthesis
target.t0$trt <- 0 # assign control in synthesis
aug.target <- rbind(target.t0, target.t1)
# complete syntheses by drawing binary outcomes from their posterior predictive distribution
y.star <- posterior_predict(outcome.model, newdata=aug.target)

## ANALYSIS STAGE ##
# fit second-stage regression to each synthesis using maximum-likelihood estimation  
reg2.fits <- lapply(1:M, function(m) glm(y.star[m,]~trt, data=aug.target, family=binomial))
# treatment effect point estimates given by treatment coefficient for each synthesis
hats.delta <- unlist(lapply(reg2.fits, function(fit) coef(fit)["trt"][[1]]))  
# point estimates for the variance in each synthesis
hats.v <- unlist(lapply(reg2.fits, function(fit) vcov(fit)["trt", "trt"]))
# quantities originally defined by Rubin (1987) for multiple imputation
bar.delta <- mean(hats.delta) # average of treatment effect point estimates
bar.v <- mean(hats.v) # "within" variance (average of variance point estimates) 
b <- var(hats.delta) # "between" variance (sample variance of point estimates)
# pooling: average of point estimates is marginal log odds ratio
hat.Delta <- bar.delta
# pooling: use combining rules to estimate the variance
hat.var.Delta <- (1+(1/M))*b-bar.v 
# Wald-type interval estimates constructed using t-distribution with nu degrees of freedom
nu <- (M-1)*(1+bar.v/((1+1/M)*b))^2
lci.Delta <- hat.Delta + qt(0.025, df=nu)*sqrt(hat.var.Delta)
uci.Delta <- hat.Delta + qt(0.975, df=nu)*sqrt(hat.var.Delta)
