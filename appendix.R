# This file implements the simulation study in Additional file 1 (Supplementary Appendix)
# MIM standardizes over the index study, for which there are some missing outcome values

rm(list=ls())

# setwd("C:/Users/gndtl/OneDrive - Bayer/Personal Data/New Papers/MIM") 

## load the function by Phillippo et al. (2020) used to simulate the data
source("./simulation_functions.R")
## load the functions used to compute the performance measures
source("./performance_functions.R")

## load packages
# for simulating covariates from a multivariate Gaussian copula
if(!require("copula")) {install.packages("copula"); library(copula)}
# for detect cores
if(!require("parallel")) {install.packages("parallel"); library(parallel)}
# for parallel cluster
if(!require("doSNOW")) {install.packages("doSNOW"); library(doSNOW)}
# for the first-stage outcome regression and for drawing outcomes from the PPD in MIM
if(!require("rstanarm")) {install.packages("rstanarm"); library(rstanarm)}
# to plot the results of the simulation study
if(!require("ggplot2")) {install.packages("ggplot2"); library(ggplot2)}
# to plot the results of the simulation study
if(!require("ggridges")) {install.packages("ggridges"); library(ggridges)}
# to plot the results of the simulation study
if(!require("gridExtra")) {install.packages("gridExtra"); library(gridExtra)}

#### SPECIFY SIMULATION STUDY SETUP #####

set.seed(444) # set random seed for reproducibility

N.sim <- 1000 # number of simulated Monte Carlo replicates
N <- 750 # number of participants in index RCT
prog_strength <- 2 # main prognostic coefficients are prog_strength x covariate sd
inter_strength <- 1 # interaction coefficients are inter_strength x covariate sd
rho <- 0.15 # pairwise correlation coefficients for the covariates
prop.missing <- c(0.1, 0.2, 0.3, 0.4) # proportion of missing subject-level outcomes
scenarios <- length(prop.missing) # number of simulation scenarios

#### GENERATE THE SIMULATED DATA #####

for (i in 1:scenarios) {
  print(i)
  # simulate patient-level data for the index RCT and the target covariate dataset
  datasets <- replicate(n=N.sim, expr=simulateData(N_index=N, N_target=100, # target will not be used in analysis 
                                                   prog_strength=prog_strength, inter_strength,
                                                   overlap=0.5, # will not be used in analysis 
                                                   X_corr=rho, 
                                                   missing=TRUE, prop.missing=prop.missing[i],
                                                   truth=FALSE),
                        simplify=FALSE)
  file.id <- paste0("prop_missing_", prop.missing[i])                          
  save(datasets, file=paste0("Data/", file.id, ".RData")) 
}

## load simulated datasets for all scenarios
datasets.all <- vector(mode="list", scenarios)
for (i in 1:scenarios) {
  file.id <- paste0("prop_missing_", prop.missing[i])  
  load(paste0("Data/", file.id, ".RData"))
  datasets.all[[i]] <- datasets
}

#### PERFORM MIM ON THE SIMULATED DATASETS ####

## Settings
M <- 1000 # number of syntheses/imputations used in analysis stage (high for low Monte Carlo error)
n.chains <- 2 # number of Markov chains for MCMC sampler in synthesis stage
warmup <- 2000 # number of discarded warmup/burn-in iterations per chain for MCMC sampler
iters <- 4000 # number of total iterations per chain for MCMC sampler (including warmup)

## Multiple imputation marginalization (MIM)
mim.wrapper <- function(datasets, M, n.chains, warmup, iters) {
  # Inputs: datasets - index RCT and target covariate datasets (target not used)
  # M - number of syntheses used in analysis stage (high for low Monte Carlo error)  
  # n.chains, warmup, iters - info for MCMC sampler
  ipd.index <- datasets[[1]]
  ipd.target <- datasets[[1]] # target assumed to be the index RCT
  ## SYNTHESIS STAGE ##
  # first-stage logistic regression model fitted to index RCT using MCMC (Stan)
  outcome.model <- stan_glm(y~(trt*x1)+(trt*x2), data=ipd.index, family=binomial, algorithm="sampling",
                            iter=iters, warmup=warmup, chains=n.chains, 
                            # thin to use M independent samples in analysis stage 
                            thin=(n.chains*(iters-warmup))/M) 
  # create augmented target dataset
  target.t1 <- target.t0 <- ipd.target[,c("x1","x2")]
  target.t1$trt <- 1
  target.t0$trt <- 0
  aug.target <- rbind(target.t0, target.t1)
  # complete syntheses by drawing binary outcomes from their posterior predictive distribution
  y.star <- posterior_predict(outcome.model, newdata=aug.target)
  ## ANALYSIS STAGE ##
  # fit second-stage regression to each synthesis using maximum-likelihood estimation  
  reg2.fits <- lapply(1:M, function(m) glm(y.star[m,]~trt, data=aug.target, family=binomial))
  # treatment effect point estimates in each synthesis
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
  list(hat.Delta, hat.var.Delta, lci.Delta, uci.Delta)
}

#### ANALYSIS OF SIMULATED DATASETS ####

# set up cluster for parallel computing
num.cores <- detectCores()-1
cluster <- makeCluster(num.cores, type="SOCK", outfile="")
registerDoSNOW(cluster)
pb <- txtProgressBar(max=N.sim, style=3) # progress bar
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)

# combine lists in parallelisation
comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}

## for each missingness scenario, run MIM for all replicates in parallel
for(i in 1:scenarios) {
  datasets <- datasets.all[[i]]
  file.id <- paste0("prop_missing_", prop.missing[i])  
  ### Parametric model-based standardization using MIM
  mim.results <- foreach(j=1:N.sim, .combine='comb', .multicombine=TRUE,
                         .init=list(list(),list(),list(),list()), .options.snow=opts,
                         .packages=c("rstanarm")) %dopar% {
                           results <- mim.wrapper(datasets[[j]], M, n.chains, warmup, iters)
                           return(results)
                         }
  means <- unlist(mim.results[[1]])
  variances <- unlist(mim.results[[2]])
  lcis <- unlist(mim.results[[3]])
  ucis <- unlist(mim.results[[4]])
  save(means, file=paste0("Results/MIM/means_", file.id, ".RData"))
  save(variances, file=paste0("Results/MIM/variances_", file.id, ".RData"))
  save(lcis, file=paste0("Results/MIM/lcis_", file.id, ".RData"))
  save(ucis, file=paste0("Results/MIM/ucis_", file.id, ".RData"))
}

stopCluster(cluster)

#### COMPUTE TRUE MARGINAL ESTIMANDS #### 

true.log.OR <- rep(NA, scenarios)
for (i in 1:scenarios) {
  N_cohort <- 2000000 # size of simulated "potential" cohort
  # compute the true index study marginal estimand via simulation
  truth.datasets <- simulateData(N_index=N_cohort, N_target=N_cohort, prog_strength=prog_strength,
                                 inter_strength=inter_strength, overlap=0.5, X_corr=rho, truth=TRUE)
  truth.target <- truth.datasets[[1]]
  p1.target <- mean(truth.target[which(truth.target$trt==1),]$y)
  p0.target <- mean(truth.target[which(truth.target$trt==0),]$y)
  # true index study marginal log OR for active treatment vs. control
  true.log.OR[i] <- log((p1.target/(1-p1.target))/(p0.target/(1-p0.target)))  
} 

#### PROCESS RESULTS OF THE SIMULATION STUDY ####

# function that computes the performance measures
process.metrics <- function(means, lcis, ucis, truth) {
  # means: point estimates for a given scenario; 
  # lcis: lower interval estimates; ucis: upper interval estimates; truth: true estimand value
  N.sim <- length(means) # number of simulations
  # bias with mcse
  bias.metric <- bias(means, truth)
  bias.metric.mcse <- bias.mcse(means)
  # empirical standard error with mcse
  ese <- empse(means)
  ese.mcse <- empse.mcse(ese, N.sim) 
  # mean square error with mcse
  mse.metric <- mse(means, truth) 
  mse.metric.mcse <- mse.mcse(means, truth)
  # coverage with mcse
  cov <- coverage(lcis, ucis, truth)
  cov.mcse <- coverage.mcse(cov, length(means))
  list(bias.metric, bias.metric.mcse, ese, ese.mcse, mse.metric, mse.metric.mcse,
       cov, cov.mcse)
} 

# generate table to display performance metrics 
display.table <- as.data.frame(matrix(nrow=scenarios, ncol=5)) 
display.table[,1] <- prop.missing 
colnames(display.table) <- c(expression(paste("Missingness proportion (", pi, ")")), 
                             "Bias", "ESE", "MSE", "Coverage")

# fills in results table with performance metrics for each missingness scenario
for (i in 1:scenarios) {
  file.id <- paste0("prop_missing_", prop.missing[i]) 
  load(file=paste0("Results/MIM/means_", file.id, ".RData"))
  load(file=paste0("Results/MIM/variances_", file.id, ".RData"))
  load(file=paste0("Results/MIM/lcis_", file.id, ".RData"))
  load(file=paste0("Results/MIM/ucis_", file.id, ".RData")) 
  mim.metrics <- process.metrics(means, lcis, ucis, truth=true.log.OR[i])
  display.table[i,2:5] <- cbind(Bias=paste0(format(round(mim.metrics[[1]],digits=3),nsmall=3)," (",
                                            format(round(mim.metrics[[2]],digits=3),nsmall=3),")"),
                                ESE=paste0(format(round(mim.metrics[[3]],digits=3),nsmall=3)," (",
                                           format(round(mim.metrics[[4]],digits=3),nsmall=3),")"),
                                MSE=paste0(format(round(mim.metrics[[5]],digits=3),nsmall=3)," (",
                                           format(round(mim.metrics[[6]],digits=3),nsmall=3),")"),
                                Coverage=paste0(format(round(mim.metrics[[7]],digits=3),nsmall=3)," (",
                                                format(round(mim.metrics[[8]],digits=3),nsmall=3),")"))
}  
