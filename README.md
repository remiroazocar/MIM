# Model-based standardization using multiple imputation: Code

### Antonio Remiro-Azócar, Anna Heath, Gianluca Baio
### *remiroantonio@gmail.com*
### *2023*

This repository contains the R code used for my paper [Model-based standardization using multiple imputation][1], co-authored with [Prof. Gianluca Baio][2] and [Prof. Anna Heath][3]. 

## Utilizing the Scripts

In order to use this repository, the user must first download a copy to their local machine. The user must set the working directory to the location 
where the download was made, and open and run the `main.R` script. This script specifies the settings of the simulation study, generates the data, 
performs the standardization methods (saving the point estimates, variances and interval estimates to the `"./Results/"` subdirectory), 
computes the relevant performance metrics, and graphs the results of the simulation study. The `simulation_functions.R` script contains user-defined
functions, adapted from [Phillippo et al. (2020)][4], to generate the simulation study data. The `performance_functions.R` script contains user-defined functions to evaluate the performance measures of interest. Run the `appendix.R` script to reproduce the simulation study in Additional file 1 (Supplementary Appendix). 

The `doSNOW` package is used to parallelize the performance of the methods in the simulation study, distributing the tasks to different cores of the computer. 

The code was prepared in `RStudio` using `R` version `4.1.1` in a Windows architecture, with a 64-bit operating system. The following packages and versions were used:

* `boot 1.3.28` for the non-parametric bootstrap in standard parametric model-based standardization (G-computation) 
* `copula 1.0.1` to simulate covariates from a multivariate Gaussian copula when simulating the data
* `doSNOW 1.0.19` used in combination with `foreach()` to start up local clusters that distribute parallel tasks to different cores
* `ggplot2 3.3.5` to plot the simulation study results (Figure 2 in the article)
* `ggridges 0.5.3` to plot the simulation study results (Figure 2 in the article)
* `gridExtra 2.3` to plot the simulation study results (Figure 2 in the article)
* `parallel 4.1.1` to detect the number of CPU cores
* `rstanarm 2.21.1` for the synthesis stage (fitting the first-stage outcome regression and drawing outcomes from the posterior predictive distribution) in multiple imputation marginalization

[1]: https://doi.org/10.48550/arXiv.2210.01757
[2]: https://gianluca.statistica.it/
[3]: https://sites.google.com/site/annaheathstats/
[4]: https://doi.org/10.1002/sim.8759
