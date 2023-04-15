# This file specifies the setup of the simulation study, generates the simulated data 
# and performs the standardization methods on the simulated data

rm(list=ls())

# setwd("C:/Users/gndtl/OneDrive - Bayer/Personal Data/New Papers/MIM") 

## load the function by Phillippo et al. (2020) used to simulate the data
source("./simulation_functions.R")
## load the functions used to compute the performance measures in the simulation study
source("./performance_functions.R")

## load packages
# for simulating covariates from a multivariate Gaussian copula when simulating the data
if(!require("copula")) {install.packages("copula"); library(copula)}
# for detect cores
if(!require("parallel")) {install.packages("parallel"); library(parallel)}
# for parallel cluster
if(!require("doSNOW")) {install.packages("doSNOW"); library(doSNOW)}
# for non-parametric bootstrap in standard parametric model-based standardization (G-computation)
if(!require("boot")) {install.packages("boot"); library(boot)}
# for the first-stage outcome regression and for drawing outcomes from the PPD in MIM
if(!require("rstanarm")) {install.packages("rstanarm"); library(rstanarm)}
# to plot the results of the simulation study
if(!require("ggplot2")) {install.packages("ggplot2"); library(ggplot2)}
# to plot the results of the simulation study
if(!require("ggridges")) {install.packages("ggridges"); library(ggridges)}
# to plot the results of the simulation study
if(!require("gridExtra")) {install.packages("gridExtra"); library(gridExtra)}

set.seed(444) # set random seed for reproducibility

## Define simulation study parameters
N.sim <- 1000 # number of simulated Monte Carlo replicates per scenario/method
N <- c(500,1000,2000) # number of participants in index RCT
N_tar <- 2000 # number of subjects in target covariate dataset
prog_strength <- 2 # main prognostic coefficients are prog_strength x covariate sd
inter_strength <- 1 # interaction coefficients are inter_strength x covariate sd
kappa <- c(0.5,1) # proxy overlap parameter
rho <- 0.15 # pairwise correlation coefficients for the covariates
# parameter combinations for each simulation study scenario
pc <- expand.grid(N=N, overlap=kappa)
scenarios <- nrow(pc) # number of simulation scenarios

#### SIMULATE DATA #####

for (i in 1:scenarios) {
  print(i)
  # simulate patient-level data for the index RCT and the target covariate dataset
  datasets <- replicate(n=N.sim, expr=simulateData(N_index=pc$N[i], N_target=N_tar, 
                                                   prog_strength=prog_strength, inter_strength,
                                                   overlap=pc$overlap[i], X_corr=rho, truth=FALSE),
                        simplify=FALSE)
  file.id <- paste0("N_", pc$N[i], "_overlap_", pc$overlap[i])                          
  save(datasets, file=paste0("Data/", file.id, ".RData")) 
}

## load simulated datasets for all scenarios
datasets.all <- vector(mode="list", scenarios)
for (i in 1:scenarios) {
  file.id <- paste0("N_", pc$N[i], "_overlap_", pc$overlap[i])
  load(paste0("Data/", file.id, ".RData"))
  load(paste0("Data/", file.id, ".RData"))
  datasets.all[[i]] <- datasets
}

#### STANDARDIZATION METHODOLOGIES ####

## settings for the standardization methodologies 
# standard version of model-based standardization
resamples <- 1000 # number of resamples in non-parametric bootstrap
# multiple imputation marginalization (MIM)
M <- 1000 # number of syntheses/imputations used in analysis stage (high for low Monte Carlo error)
n.chains <- 2 # number of Markov chains for MCMC sampler in synthesis stage
warmup <- 2000 # number of discarded warmup/burn-in iterations per chain for MCMC sampler
iters <- 4000 # number of total iterations per chain for MCMC sampler (including warmup)

## Standard parametric model-based standardization (G-computation) with maximum-likelihood estimation and bootstrapping
gcomp.ml.wrapper <- function(datasets, resamples) {
  # Inputs: datasets - index RCT and target covariate datasets
  # resamples - number of resamples for non-parametric bootstrap
  ipd.index <- datasets[[1]]
  ipd.target <- datasets[[2]]
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
  lci.Delta <- unname(conf.ints[1]) # lower interval bound
  uci.Delta <- unname(conf.ints[2]) # upper interval bound
  list(hat.Delta, hat.var.Delta, lci.Delta, uci.Delta)
}  

## Multiple imputation marginalization (MIM)
mim.wrapper <- function(datasets, M, n.chains, warmup, iters) {
  # Inputs: datasets - index RCT and target covariate datasets
  # M - number of syntheses used in analysis stage (high for low Monte Carlo error)  
  # n.chains, warmup, iters - info for MCMC sampler
  ipd.index <- datasets[[1]]
  ipd.target <- datasets[[2]]
  ## SYNTHESIS STAGE ##
  # first-stage logistic regression model fitted to index RCT using MCMC (Stan)
  outcome.model <- stan_glm(y~(trt*x1)+(trt*x2), data=ipd.index, family=binomial, algorithm="sampling",
                            iter=iters, warmup=warmup, chains=n.chains, 
                            # thin to use M independent samples in analysis stage 
                            thin=(n.chains*(iters-warmup))/M) 
  # create augmented target dataset
  target.t1 <- target.t0 <- ipd.target
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

## for each scenario, run model-based standardization methods for all replicates in parallel
for(i in 1:scenarios) {
  datasets <- datasets.all[[i]]
  file.id <- paste0("N_", pc$N[i], "_overlap_", pc$overlap[i])  
  ### Standard parametric model-based standardization (G-computation) with maximum-likelihood estimation and bootstrapping
  gcomp.ml.results <- foreach(j=1:N.sim, .combine='comb', .multicombine=TRUE,
                              .init=list(list(),list(),list(),list()), .options.snow=opts,
                              .packages=c("boot")) %dopar% {
                                results <- gcomp.ml.wrapper(datasets[[j]], resamples)
                                return(results)
                              }
  means <- unlist(gcomp.ml.results[[1]])
  variances <- unlist(gcomp.ml.results[[2]])
  lcis <- unlist(gcomp.ml.results[[3]])
  ucis <- unlist(gcomp.ml.results[[4]]) 
  save(means, file=paste0("Results/GcompML/means_", file.id, ".RData"))
  save(variances, file=paste0("Results/GcompML/variances_", file.id, ".RData"))
  save(lcis, file=paste0("Results/GcompML/lcis_", file.id, ".RData"))
  save(ucis, file=paste0("Results/GcompML/ucis_", file.id, ".RData"))
  ### Parametric model-based standardization using multiple imputation marginalization (MIM)
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
  # compute true marginal estimand in the target covariate distribution via simulation
  truth.datasets <- simulateData(N_index=N_cohort, N_target=N_cohort, prog_strength=prog_strength,
                                 inter_strength=inter_strength, overlap=pc$overlap[i], X_corr=rho, truth=TRUE)
  truth.target <- truth.datasets[[2]]
  p1.target <- mean(truth.target[which(truth.target$trt==1),]$y)
  p0.target <- mean(truth.target[which(truth.target$trt==0),]$y)
  # true marginal log OR for active treatment vs. control in the target covariate distribution
  true.log.OR[i] <- log((p1.target/(1-p1.target))/(p0.target/(1-p0.target)))  
} 

#### PROCESS RESULTS OF THE SIM. STUDY, COMPUTE AND PLOT SIMULATION STUDY METRICS ####

# function that computes the performance measures for a given method
process.metrics <- function(means, lcis, ucis, truth) {
  # means: point estimates for a given method/scenario; 
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

simulation.metrics <- as.data.frame(matrix(nrow=2*scenarios, ncol=11)) # store performance metrics here
colnames(simulation.metrics) <- c("N", "Overlap", "Method", "Bias", "Bias.MCSE", "ESE", "ESE.MCSE", 
                                  "MSE", "MSE.MCSE", "Cov", "Cov.MCSE")
                                  
ate.table <- as.data.frame(matrix(nrow=2*scenarios*N.sim, ncol=4)) # store all marginal log OR point estimates here
colnames(ate.table) <- c("N", "Overlap", "Method", "ATE")

# insert performance metrics and point estimates for each method/scenario into the data frames
for (i in 1:scenarios) {
  # row indices for the methods in the data frame of performance metrics
  mbs.id <- (2*i)-1
  mim.id <- 2*i 
  # row indices for the methods in the data frame of marginal log OR point estimates
  mbs.id.ate <- ((mbs.id-1)*N.sim+1):(mbs.id*N.sim) 
  mim.id.ate <- ((mim.id-1)*N.sim+1):(mim.id*N.sim)
  # column 1 is the number of subjects and column 2 is the overlap level in the simulation scenario
  simulation.metrics[mbs.id:mim.id, 1] <- pc$N[i]
  simulation.metrics[mbs.id:mim.id, 2] <- pc$overlap[i]
  ate.table[mbs.id.ate, 1] <- pc$N[i]
  ate.table[mbs.id.ate, 2] <- pc$overlap[i]
  ate.table[mim.id.ate, 1] <- pc$N[i]
  ate.table[mim.id.ate, 2] <- pc$overlap[i]  
  file.id <- paste0("N_", pc$N[i], "_overlap_", pc$overlap[i]) 
  # standard parametric model-based standardization (G-computation) with maximum-likelihood estimation and bootstrapping
  simulation.metrics[mbs.id, 3] <- "Standard"
  load(file=paste0("Results/GcompML/means_", file.id, ".RData"))
  load(file=paste0("Results/GcompML/variances_", file.id, ".RData"))
  load(file=paste0("Results/GcompML/lcis_", file.id, ".RData"))
  load(file=paste0("Results/GcompML/ucis_", file.id, ".RData"))
  mbs.metrics <- process.metrics(means, lcis, ucis, truth=true.log.OR[i])
  simulation.metrics[mbs.id, 4:11] <- unlist(mbs.metrics)
  ate.table[mbs.id.ate, 3] <- "Standard"
  ate.table[mbs.id.ate, 4] <- means
  # multiple imputation marginalization (MIM)
  simulation.metrics[mim.id, 3] <- "MIM"
  load(file=paste0("Results/MIM/means_", file.id, ".RData"))
  load(file=paste0("Results/MIM/variances_", file.id, ".RData"))
  load(file=paste0("Results/MIM/lcis_", file.id, ".RData"))
  load(file=paste0("Results/MIM/ucis_", file.id, ".RData")) 
  mim.metrics <- process.metrics(means, lcis, ucis, truth=true.log.OR[i])
  simulation.metrics[mim.id, 4:11] <- unlist(mim.metrics)
  ate.table[mim.id.ate, 3] <- "MIM"
  ate.table[mim.id.ate, 4] <- means
}  

# function generates ridgeline plot for point estimates for a specific scenario
plot.results <- function(scenario) {
  # scenario: index for the scenario
  scenario.ates <- subset(ate.table, N==pc$N[scenario] & Overlap==pc$overlap[scenario]) # point estimates for a scenario
  ridge.plot <- ggplot(scenario.ates, aes(x=ATE, y=Method, fill=Method)) +
    geom_density_ridges(alpha=0.65) +
    geom_vline(xintercept=true.log.OR[scenario], linetype="dashed", color ="red") +
    scale_x_continuous(limits=c(-2.8, -0.6)) +
    scale_y_discrete(limits=c("MIM", "Standard")) +
    theme_classic() + 
    theme(legend.position = "none", 
          axis.text.y = element_text(color="grey20", size=9, face ="plain"),
          axis.text.x = element_text(color="grey20", size=7, face="plain"),
          plot.title = element_text(color="grey20", size=10, face ="plain"),
          axis.title.y = element_blank(),
          axis.title.x = element_text(color="grey20", size=9, face="plain")) +
    scale_fill_brewer(palette="Dark2") + xlab("Point estimates")
  # a ridgeline plot with the spread of point estimates is returned
  return(ridge.plot)
}  

# Function generates results table with performance metrics for a specific scenario
table.results <- function(scenario) {
  # scenario: index for the scenario
  scenario.metrics <- subset(simulation.metrics, N==pc$N[scenario] & Overlap==pc$overlap[scenario]) # metrics for a scenario
  display.table <- cbind(Method=scenario.metrics$Method,
                         Bias=paste0(format(round(scenario.metrics$Bias,digits=3),nsmall=3)," (",
                                     format(round(scenario.metrics$Bias.MCSE,digits=3),nsmall=3),")"),
                         ESE=paste0(format(round(scenario.metrics$ESE,digits=3),nsmall=3)," (",
                                    format(round(scenario.metrics$ESE.MCSE,digits=3),nsmall=3),")"),
                         MSE=paste0(format(round(scenario.metrics$MSE,digits=3),nsmall=3)," (",
                                    format(round(scenario.metrics$MSE.MCSE,digits=3),nsmall=3),")"),
                         Coverage=paste0(format(round(scenario.metrics$Cov,digits=3),nsmall=3)," (",
                                    format(round(scenario.metrics$Cov.MCSE,digits=3),nsmall=3),")"))
  table.grob <- tableGrob(display.table, theme=ttheme_minimal(base_size=7))
  # a table with the performance measures and corresponding MCSEs is returned
  return(table.grob)
}

# ridgeline plot for each scenario
ridge.plot.s1 <- plot.results(scenario=1) + ggtitle(expression(paste(italic(N), "=500, limited overlap (", kappa, "=0.5)")))  
ridge.plot.s2 <- plot.results(scenario=2) + ggtitle(expression(paste(italic(N), "=1000, limited overlap (", kappa, "=0.5)")))
ridge.plot.s3 <- plot.results(scenario=3) + ggtitle(expression(paste(italic(N), "=2000, limited overlap (", kappa, "=0.5)")))
ridge.plot.s4 <- plot.results(scenario=4) + ggtitle(expression(paste(italic(N), "=500, full overlap (", kappa, "=1)")))
ridge.plot.s5 <- plot.results(scenario=5) + ggtitle(expression(paste(italic(N), "=1000, full overlap (", kappa, "=1)")))
ridge.plot.s6 <- plot.results(scenario=6) + ggtitle(expression(paste(italic(N), "=2000, full overlap (", kappa, "=1)")))

# table of results for each scenario
table.grob.s1 <- table.results(scenario=1)  
table.grob.s2 <- table.results(scenario=2)  
table.grob.s3 <- table.results(scenario=3)  
table.grob.s4 <- table.results(scenario=4)  
table.grob.s5 <- table.results(scenario=5)  
table.grob.s6 <- table.results(scenario=6)  

# figure with ridgeline plot and table of results for each scenario
ridge.grid <- arrangeGrob(ridge.plot.s1, table.grob.s1, ridge.plot.s2, table.grob.s2,
                          ridge.plot.s3, table.grob.s3, ridge.plot.s4, table.grob.s4,
                          ridge.plot.s5, table.grob.s5, ridge.plot.s6, table.grob.s6, 
                          ncol=2, widths=c(0.8,1.2)) 
ggsave(file="Figure2.pdf", plot=ridge.grid, width=170, height=225, units="mm", dpi = 300)
