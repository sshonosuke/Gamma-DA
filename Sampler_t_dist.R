###------------------------------------------------------------------###
###------------------------------------------------------------------###
###                Code implementing MCMC algorithms                 ###
###                      for t-distribution                          ###
###------------------------------------------------------------------###
###------------------------------------------------------------------###
library(truncdist)
library(zipfR)

###------------------------------------------------------------------###
###         MCMC algorithm for fitting t-distribution                ###
###               (The new data augmentation)                        ###
###------------------------------------------------------------------###
## INPUT 
# x: data (vector) 
# prior_theta: mean and standard deviation of normal prior for 'theta' 
# prior_tau: shape and rate of inverse-gamma prior for 'tau' 
# prior_alpha: shape and rate of gamma prior for 'alpha' 
# alpha_min: lower bound of alpha (optional)
# mc: MCMC length 
# burn: length of burn-in period 
## OUTPUT
# matrix of posterior samples of ('theta', 'tau', 'alpha')

DA_sampler <- function(x, prior_theta=c(0, sqrt(0.1)), prior_tau=c(0.1, 0.1), 
                       prior_alpha=c(0.1, 0.1), alpha_min=NULL, mc=5000, burn=1000){
  ## preparation 
  n <- length(x)
  a <- prior_theta[2]
  b <- prior_theta[1]
  c <- prior_tau[1]
  d <- prior_tau[2]
  a0 <- prior_alpha[1]
  b0 <- prior_alpha[2]
  
  ## function (Metropolis-Hastings ratio)
  MHratio <- function(new, old){
    (new - 1/2) * log(new) - lgamma(new) - new - (old - 1/2) * log(old) + lgamma(old) + old
  }
  
  ## initial values 
  if(is.null(alpha_min)){ 
    al <- 1 
  }else{
    al <- alpha_min + 1 
  }
  w <- rep(1, n)
  th <- ta <- NA
  rho <- rep(NA, n - 1)
  
  ## matrices to store posterior samples
  Pos_sample <- matrix(NA, mc, 3)
  dimnames(Pos_sample)[[2]] <- c("theta", "tau", "alpha")
  
  ## MCMC iterations
  for(r in 1:mc){
    # tau (scale) 
    tau_shape <- n/2 + c
    tau_rate <- 0.5*(a^2*b^2 + sum(w*x^2) - (a^2*b + sum(w*x))^2/(a^2 + sum(w))) + d
    ta <- 1 / rgamma(1, shape=tau_shape, rate=tau_rate)
    Pos_sample[r, 2] <- ta
    
    # theta (location)
    th_mean <- (a^2*b + sum(w*x)) / (a^2 + sum(w))
    th_sd <- sqrt(ta/(a^2 + sum(w)))
    th <- rnorm(1, mean=th_mean, sd=th_sd)
    Pos_sample[r, 1] <- th
    
    # w (latent variable)
    w <- rgamma(n, shape=al+1/2, rate=al+(x-th)^2/(2*ta))
    
    # rho (latent variable for full conditional of 'alpha')
    rho <- rbeta(n-1, shape1 = al+((2:n)-1)/n, shape2=(n-(2:n)+1)/n)
    
    # alpha (degrees of freedom)
    alpha_shape <- n + 1/2 - 1 + a0
    alpha_rate <- b0 - n + sum(w-log(w)) + sum(log(1/rho))
    if(is.null(alpha_min)){
      al_proposal <- rgamma(1, shape=alpha_shape, rate=alpha_rate)   # proposal
    }else{
      al_proposal <- al
      try( al_proposal <- rtrunc(1, spec="gamma", a=alpha_min, b=Inf, shape=alpha_shape, rate=alpha_rate)  )  # proposal
    }
    log_uniform <- log(runif(1))
    if(log_uniform <= (-1)/(12*n*al_proposal)){   # MH step 
      al <- al_proposal
    }
    else{
      A_proposal <- n * al_proposal
      A_old <- n * al
      log_ratio <- MHratio(new=A_proposal, old=A_old)
      if(log_uniform <= log_ratio){   # accept & reject
        al <- al_proposal
      }
    }
    Pos_sample[r, 3] <- al
  }
  
  ## Summary
  Pos_sample <- Pos_sample[-(1:burn),]
  return(Pos_sample)
}
  




###------------------------------------------------------------------###
###         MCMC algorithm for fitting t-distribution                ###
###         (approximated MH algorithm by Miller (2019))             ###
###------------------------------------------------------------------###
## INPUT 
# x: data (vector) 
# prior_theta: mean and precision of normal prior for 'theta' 
# prior_tau: shape and rate of inverse-gamma prior for 'tau' 
# prior_alpha: shape and rate of gamma prior for 'alpha' 
# alpha_min: lower bound of alpha (optional)
# ep: tolerance 
# M: maximum number of iterations of AMH algorithm
# mc: MCMC length 
# burn: length of burn-in period 
## OUTPUT
# matrix of posterior samples of ('theta', 'tau', 'alpha')

AMH_sampler <- function(x, prior_theta=c(0, sqrt(0.1)), prior_tau=c(0.1, 0.1), 
                       prior_alpha=c(0.1, 0.1), alpha_min=NULL, ep=10^(-8), M=10, mc=5000, burn=1000){
  ## preparation 
  n <- length(x)
  a <- prior_theta[2]
  b <- prior_theta[1]
  c <- prior_tau[1]
  d <- prior_tau[2]
  a0 <- prior_alpha[1]
  b0 <- prior_alpha[2]
  
  ## approximate MH method 
  if(is.null(alpha_min)){
    AMH <- function(u, mu=1){   # sampling function without truncation 
      nn <- length(u)
      R <- sum(log(u))
      S <- sum(u)
      TT <- S / mu - R + nn * log(mu) - nn
      A <- a0 + nn / 2
      B <- b0 + TT
      for(j in 1:M){
        a <- A/B
        A <- a0 - nn * a + nn * a^2 * trigamma(a)
        B <- b0 + (A - a0) / a - nn * log(a) + nn * digamma(a) + TT
        if(abs(a/(A/B)-1) < ep){
          return( c(A, B) )
        }
      }
      return( c(A, B) )
    }
  }else{
    AMH <- function(u, mu=1){   # sampling function under truncation 
      nn <- length(u)
      R <- sum(log(u))
      S <- sum(u)
      TT <- S / mu - R + nn * log(mu) - nn
      A <- a0 + nn / 2
      B <- b0 + TT
      a_new <- (1/B)*Igamma(A+1, B*alpha_min, lower=F)/Igamma(A, B*alpha_min, lower=F)
      for(j in 1:M){
        a <- a_new
        A <- a0 - nn * a + nn * a^2 * trigamma(a)
        B <- b0 + (A - a0) / a - nn * log(a) + nn * digamma(a) + TT
        a_new <- (1/B)*Igamma(A+1, B*alpha_min, lower=F)/Igamma(A, B*alpha_min, lower=F)
        if(abs(a/a_new-1) < ep){
          return( c(A, B) )
        }
      }
      return( c(A, B) )
    }
  }
  
  ## initial values 
  if(is.null(alpha_min)){ 
    al <- 1 
  }else{
    al <- alpha_min + 1 
  }
  w <- rep(1, n)
  th <- ta <- NA
  
  ## matrices to store posterior samples
  Pos_sample <- matrix(NA, mc, 3)
  dimnames(Pos_sample)[[2]] <- c("theta", "tau", "alpha")
  
  ## MCMC iterations
  for(r in 1:mc){
    # tau (scale) 
    tau_shape <- n/2 + c
    tau_rate <- 0.5*(a^2*b^2 + sum(w*x^2) - (a^2*b + sum(w*x))^2/(a^2 + sum(w))) + d
    ta <- 1 / rgamma(1, shape=tau_shape, rate=tau_rate)
    Pos_sample[r, 2] <- ta
    
    # theta (location)
    th_mean <- (a^2*b + sum(w*x)) / (a^2 + sum(w))
    th_sd <- sqrt(ta/(a^2 + sum(w)))
    th <- rnorm(1, mean=th_mean, sd=th_sd)
    Pos_sample[r, 1] <- th
    
    # w (latent variable)
    w <- rgamma(n, shape=al+1/2, rate=al+(x-th)^2/(2*ta))
    
    # alpha (degrees of freedom)
    AB <- AMH(w)
    A <- AB[1]
    B <- AB[2]
    if(is.null(alpha_min)){
      al_proposal <- rgamma(1, shape=A, rate=B)   # proposal 
    }else{
      al_proposal <- al
      try( al_proposal <- rtrunc(1, spec="gamma", a=alpha_min, b=Inf, shape=A, rate=B) )
    }
    log_uniform <- log(runif(1))
    log_ratio <- (a0 - A)*(log(al_proposal) - log(al)) - (al_proposal - al)*(b0 - B) + n*( al_proposal*log(al_proposal) - lgamma(al_proposal) - al*log(al) + lgamma(al) ) + (al_proposal - al)*sum(log(w) - w)
    if(log_uniform <= log_ratio){  # accept & reject
      al <- al_proposal
    }
    Pos_sample[r, 3] <- al
  }
  
  ## Summary
  Pos_sample <- Pos_sample[-(1:burn),]
  return(Pos_sample)
}
