#########################################################################################
# Simulation functions to generate the simulation study data
# 
# Based on "Assessing the performance of population adjustment methods for anchored 
# indirect comparisons: A simulation study."
# by Phillippo DM, Dias S, Ades AE, Welton NJ
#########################################################################################

# Defines the function used to simulate the data for the simulation study
simulateData <- function(N_index, # number of subjects in the index study
                         N_target, # number of subjects in the target dataset
                         prog_strength, # determines the strength of the main "prognostic" coefficients
                         inter_strength, # determines the strength of the interaction coefficients
                         overlap, # determines the level of covariate overlap
                         X_corr, # pairwise correlation coefficient for covariates
                         # truth=TRUE means that the returned data frame will be used to compute
                         # the true marginal estimands in the index and target (not for the sim. study)
                         truth=FALSE) { 
  n_X <- 2 # number of baseline covariates
  # Mean and standard deviation of normally-distributed marginal covariates in index RCT
  x1_index_mean <- 1
  x1_index_sd <- 0.5
  x2_index_mean <- 0.5
  x2_index_sd <- 0.2
  # Set mean and sd of the target normally-distributed marginal covariates based on overlap level
  x1_target_mean <- x1_index_mean * (1.1 + (1 - overlap)^n_X)
  x1_target_sd <- x1_index_sd * 0.75
  x2_target_mean <- x2_index_mean * (1.1 + (1 - overlap)^n_X)
  x2_target_sd <- x2_index_sd * 0.75  
  # Set intercept of outcome-generating logistic model
  b_0 <- -0.5
  # Set conditional treatment effect of active treatment vs. placebo at baseline (x=0)
  b_t <- -1.5
  # Set main "prognostic" effects (in log OR scale)
  b_11 <- prog_strength * x1_index_sd
  b_12 <- prog_strength * x2_index_sd
  # Set treatment-covariate interactions (conditional effect measure modification)
  b_21 <- inter_strength * x1_index_sd
  b_22 <- inter_strength * x2_index_sd  
  simulateIPD <- function(N, x1_mean, x1_sd, x2_mean, x2_sd, truth) {
    # Simulate covariates -------------------------------------------------------
    # Use the inverse CDF approach to incorporate correlations with a Gaussian copula
    cop <- normalCopula(X_corr, dim = n_X, dispstr = "un")  
    u <- cCopula(matrix(runif(n_X*N), ncol=n_X), cop, inverse=TRUE)
    colnames(u) <- paste0("u", 1:n_X)
    x1 <- qnorm(u[,"u1"], x1_mean, x1_sd)
    x2 <- qnorm(u[,"u2"], x2_mean, x2_sd)
    if (truth==FALSE) {
      # Treatment -----------------------------------------------------------------
      trt <- rep(c(1,0), each=N/2) # 1:1 treatment allocation ratio
      # Simulate outcomes ---------------------------------------------------------
      # Generate binary outcomes with logit link
      LP <- ifelse(trt==0, b_0 + x1*b_11 + x2*b_12, b_0 + x1*b_11 + x2*b_12 + b_t + x1*b_21 + x2*b_22) 
      y  <- ifelse(runif(N) > 1-plogis(LP), 1, 0)
      ipd <- data.frame(x1=x1, x2=x2, trt=trt, y=y)
    } else {
      # Simulate outcomes --------------------------------------------------------- 
      LP1 <- b_0 + x1*b_11 + x2*b_12 + b_t + x1*b_21 + x2*b_22
      LP0 <- b_0 + x1*b_11 + x2*b_12
      # two potential "counterfactual" datasets under each treatment 
      ipd.t0 <- ipd.t1 <- data.frame(x1=x1, x2=x2)
      ipd <- rbind(ipd.t1, ipd.t0)
      LP <- c(LP1, LP0)
      ipd$trt <- rep(c(1,0), each=N) # N is simulated cohort size (two copies of the cohort are concatenated)
      # Potential subject-level binary outcomes under active treatment and control
      ipd$y <- ifelse(runif(2*N) > 1-plogis(LP), 1, 0) # simulate binary outcomes with logit link
    }
    return(ipd)
  }
  ipd_index <- simulateIPD(N=N_index, x1_mean=x1_index_mean, x1_sd=x1_index_sd, 
                           x2_mean=x2_index_mean, x2_sd=x2_index_sd, truth=truth)
  ipd_target <- simulateIPD(N=N_target, x1_mean=x1_target_mean, x1_sd=x1_target_sd, 
                            x2_mean=x2_target_mean, x2_sd=x2_target_sd, truth=truth)
  if(truth==FALSE) {
    # individual-level treatment and outcomes in target assumed unavailable or not relevant
    ipd_target <- subset(ipd_target, select = -c(trt,y))
  }
  return(list(ipd_index=ipd_index, ipd_target=ipd_target))
}
