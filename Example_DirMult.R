###------------------------------------------------------------------###
###------------------------------------------------------------------###
###             Code for one-shot simulation study                   ###
###           under Dirichlet-Multinomial distribution               ###
###------------------------------------------------------------------###
###------------------------------------------------------------------###
rm(list=ls())

## load R files and packages
library(coda)
source("Sampler_DirMult.R")


## simulation settings 
L <- 10    # number of items
n <- 100     # sample size
N <- 500    # size in multinomial distribution 
true_alpha <- rep(1/10, L)   # true shape parameters of Dirichlet distribution  


## data generation 
ep <- rdirichlet(n, true_alpha)
x <- matrix(NA, n, L)
for(i in 1:n){
  x[i,] <- rmultinom(1, size=N, prob=ep[i,])
}


## posterior computation 
start_time <- proc.time()[3]
DAN_pos <- DAN_sampler(x=x)
time_DAN <- proc.time()[3] - start_time

start_time <- proc.time()[3]
DAP_pos <- DAP_sampler(x=x)
time_DAP <- proc.time()[3] - start_time

start_time <- proc.time()[3]
DAPT_pos <- DAPT_sampler(x=x)
time_DAPT <- proc.time()[3] - start_time

start_time <- proc.time()[3]
ERG_pos <- ERG_sampler(x=x)
time_ERG <- proc.time()[3] - start_time



## effective sample size
apply(DAN_pos, 2, effectiveSize)
apply(DAP_pos, 2, effectiveSize)
apply(DAPT_pos, 2, effectiveSize)
apply(ERG_pos, 2, effectiveSize)


## effective sample size (per unit second )
apply(DAN_pos, 2, effectiveSize) / time_DAN
apply(DAP_pos, 2, effectiveSize) / time_DAP
apply(DAPT_pos, 2, effectiveSize) / time_DAPT
apply(ERG_pos, 2, effectiveSize) / time_ERG    


## posterior means
apply(DAN_pos, 2, mean)
apply(DAP_pos, 2, mean)
apply(DAPT_pos, 2, mean)
apply(ERG_pos, 2, mean)

