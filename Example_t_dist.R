###------------------------------------------------------------------###
###------------------------------------------------------------------###
###             Code for one-shot simulation study                   ###
###                    under t-distribution                          ###
###------------------------------------------------------------------###
###------------------------------------------------------------------###
rm(list=ls())


## load R files and packages
library(coda)
source("Sampler_t_dist.R")


## simulation settings 
n <- 50     # sample size
true_alpha <- 1   # true degrees of freedom parameter
true_theta <- 3    # true location parameter
true_tau <- 1   # true scale parameter


## data generation 
x <- rt(n, df=2*true_alpha) * sqrt(true_tau) + true_theta


## posterior computation (using standard gamma prior) 
start_time <- proc.time()[3]
DA_pos <- DA_sampler(x)
time_DA <- proc.time()[3] - start_time

start_time <- proc.time()[3]
AMH_pos <- AMH_sampler(x)
time_AMH <- proc.time()[3] - start_time

# effective sample size
apply(DA_pos, 2, effectiveSize)
apply(AMH_pos, 2, effectiveSize)

# effective sample size (per unit second )
apply(DA_pos, 2, effectiveSize) / time_DA
apply(AMH_pos, 2, effectiveSize) / time_AMH


# posterior means
apply(DA_pos, 2, mean)
apply(AMH_pos, 2, mean)

c(true_theta, true_tau^2, true_alpha)   # true value




## posterior computation (using truncated gamma prior ) 
start_time <- proc.time()[3]
DA_pos_trunc <- DA_sampler(x, alpha_min=0.5)
time_DA_trunc <- proc.time()[3] - start_time

start_time <- proc.time()[3]
AMH_pos_trunc <- AMH_sampler(x, alpha_min=0.5)
time_AMH_trunc <- proc.time()[3] - start_time


# effective sample size (per unit second )
apply(DA_pos_trunc, 2, effectiveSize) / time_DA_trunc
apply(AMH_pos_trunc, 2, effectiveSize) / time_AMH_trunc


# posterior means
apply(DA_pos_trunc, 2, mean)
apply(AMH_pos_trunc, 2, mean)

c(true_theta, true_tau^2, true_alpha)   # true value

