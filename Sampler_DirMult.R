###------------------------------------------------------------------###
###------------------------------------------------------------------###
###                Code implementing MCMC algorithms                 ###
###             for Dirichlet-Multinomial distribution               ###
###------------------------------------------------------------------###
###------------------------------------------------------------------###


## Packages 
library(DirichletReg)
library(GIGrvg)


## Customized function to sample from GIG distribution
myrgig <- function(n, Lambda, Chi, Psi){
  LambdaChiPsi <- cbind(Lambda, Chi, Psi)   
  apply(LambdaChiPsi, 1, function(x) rgig(1, lambda=x[1], chi=x[2], psi=x[3]))
}


## function (PTN distribution)
PTN <- function(p, a, b){
  ta <- sqrt(1/4 + 2*a*p/b^2) + sign(b)/2
  C1 <- ta*abs(b) - b
  C2 <- ta*abs(b)/(2*a)
  TF <- TRUE
  while(TF){
    x <- rgamma(1, shape=p, rate=C1)
    u <- runif(1)
    TF <- (u > exp( -a*(x-C2)^2) )
  }
  return(x)
}



###------------------------------------------------------------------###
###      MCMC algorithm for Dirichlet-Multinomial distribution       ###
###         (data augmentation with normal random variable)          ###
###------------------------------------------------------------------###
## INPUT 
# x: (n,L)-matrix
# prior_alpha: two-dimensional vector of shape and rate of gamma prior for 'alpha' 
# mc: MCMC length 
# burn: length of burn-in period 
## OUTPUT
# matrix of posterior samples of 'alpha'

DAN_sampler <- function(x, prior_alpha=NULL, mc=5000, burn=1000){
  ## preparation 
  L <- dim(x)[2]
  n <- dim(x)[1]
  if(is.null(prior_alpha)){  prior_alpha <- c(1/L, 1) }
  a <- prior_alpha[1]
  b <- prior_alpha[2]
  
  ## function 
  MHratio <- function(new, old){
    (new-1/2)*log(new) - lgamma(new) - new - (old-1/2)*log(old) + lgamma(old) + old
  }
  
  ## initial values 
  al <- rep(1, L)
  
  ## matrices to store posterior samples
  Pos_sample <- matrix(NA, mc, L)
  dimnames(Pos_sample)[[2]] <- paste0("alpha", 1:L)
  
  ## MCMC iterations
  for(r in 1:mc){
    # p (sell probability)
    u <- matrix(rgamma(n*L, shape=t(x)+al, rate=1), L, n)
    p <- u / c(colSums(u)%x%rep(1, L))
    
    # z (latent variable)
    z <- rgamma(n, shape=sum(al), rate=1)
    
    # v (latent variable)
    v <- rgamma(L, shape=n*al, rate=n*al^2)
    
    # rho (latent variable)
    rho <- matrix(rbeta(L*(n-1), shape1=al+(((2:n)-1)/n)%x%rep(1, L), shape2=((n-(2:n)+1)/n)%x%rep(1, L)), L, n-1)
    
    # theta (latent variable)
    Bt <- b - n - sum(log(z)) + rowSums(log(1 / p)) + rowSums(log(1 / rho)) - n - n * log(v)
    B <- pmax(0, Bt) + 1
    A <- B - Bt
    th <- rnorm(L, mean=2*A*al, sd=sqrt(2*A*al))
    
    # eta (latent variable)
    eta <- myrgig(L, Lambda=1/2, Chi=(B*al + th^2/(4*A*al))^2, Psi=1)
    
    # alpha
    alt <- myrgig(L, Lambda=(n+a-1/2)/2, Chi=th^4/(16*eta*A^2), Psi=2*v*n+B^2/eta)
    al_proposal <- sqrt( alt )
    log_uniform <- log( runif(L) )
    I1 <- (log_uniform <= (-2)/(12*n*al_proposal))
    I2 <- (log_uniform[!I1] <= 2*MHratio(new=n*al_proposal[!I1], old=n*al[!I1]))
    al[I1] <- al_proposal[I1]
    al[!I1][I2] <- al_proposal[!I1][I2]
    al[!I1][!I2] <- al[!I1][!I2]
    Pos_sample[r,] <- al
  }
  
  ## Summary
  Pos_sample <- Pos_sample[-(1:burn),]
  return(Pos_sample)
}






###------------------------------------------------------------------###
###      MCMC algorithm for Dirichlet-Multinomial distribution       ###
###        (data augmentation with Poisson random variable)          ###
###------------------------------------------------------------------###
## INPUT 
# x: (n,L)-matrix
# prior_alpha: two-dimensional vector of shape and rate of gamma prior for 'alpha' 
# mc: MCMC length 
# burn: length of burn-in period 
## OUTPUT
# matrix of posterior samples of 'alpha'

DAP_sampler <- function(x, prior_alpha=NULL, mc=5000, burn=1000){
  ## preparation 
  L <- dim(x)[2]
  n <- dim(x)[1]
  if(is.null(prior_alpha)){  prior_alpha <- c(1/L, 1) }
  a <- prior_alpha[1]
  b <- prior_alpha[2]
  
  ## function 
  MHratio <- function(new, old){
    (new-1/2)*log(new) - lgamma(new) - new - (old-1/2)*log(old) + lgamma(old) + old
  }
  
  ## initial values 
  al <- rep(1, L)
  
  ## matrices to store posterior samples
  Pos_sample <- matrix(NA, mc, L)
  dimnames(Pos_sample)[[2]] <- paste0("alpha", 1:L)
  
  ## MCMC iterations
  for(r in 1:mc){
    # p (sell probability)
    u <- matrix(rgamma(n*L, shape=t(x)+al, rate=1), L, n)
    p <- u / c(colSums(u)%x%rep(1, L))
    
    # z (latent variable)
    z <- rgamma(n, shape=sum(al), rate=1)
    
    # v (latent variable)
    v <- rgamma(L, shape=n*al, rate=n*al^2)
    
    # rho (latent variable)
    rho <- matrix(rbeta(L*(n-1), shape1=al+(((2:n)-1)/n)%x%rep(1, L), shape2=((n-(2:n)+1)/n)%x%rep(1, L)), L, n-1)
    
    # ze (latent variable)
    Bt <- b - n - sum(log(z)) + rowSums(log(1/p)) + rowSums(log(1/rho)) - n - n*log(v)
    B <- pmax(0, Bt) + 1
    A <- B - Bt
    ze <- rpois(L, lambda=al*A)
    
    # xi (latent variable)
    xi <- myrgig(L, Lambda=1/2, Chi=B^2*al^2, Psi=1)
    
    # alpha
    alt <- rgamma(L, shape=(ze + n + a)/2, rate=B^2/(2*xi) + n*v)
    al_proposal <- sqrt( alt )
    log_uniform <- log( runif(L) )
    I1 <- (log_uniform <= (-2)/(12*n*al_proposal))
    I2 <- (log_uniform[!I1] <= 2*MHratio(new=n*al_proposal[!I1], old=n*al[!I1]))
    al[I1] <- al_proposal[I1]
    al[!I1][I2] <- al_proposal[!I1][I2]
    al[!I1][!I2] <- al[!I1][!I2]
    Pos_sample[r,] <- al
  }
  
  ## Summary
  Pos_sample <- Pos_sample[-(1:burn),]
  return(Pos_sample)
}




###------------------------------------------------------------------###
###      MCMC algorithm for Dirichlet-Multinomial distribution       ###
###          (data augmentation with PTN random variable)            ###
###------------------------------------------------------------------###
## INPUT 
# x: (n,L)-matrix
# prior_alpha: two-dimensional vector of shape and rate of gamma prior for 'alpha' 
# mc: MCMC length 
# burn: length of burn-in period 
## OUTPUT
# matrix of posterior samples of 'alpha'

DAPT_sampler <- function(x, prior_alpha=NULL, mc=5000, burn=1000){
  ## preparation 
  L <- dim(x)[2]
  n <- dim(x)[1]
  if(is.null(prior_alpha)){  prior_alpha <- c(1/L, 1) }
  a <- prior_alpha[1]
  b <- prior_alpha[2]
  
  ## function 
  MHratio <- function(new, old){
    (new-1/2)*log(new) - lgamma(new) - new - (old-1/2)*log(old) + lgamma(old) + old
  }
  
  ## initial values 
  al <- rep(1, L)
  
  ## matrices to store posterior samples
  Pos_sample <- matrix(NA, mc, L)
  dimnames(Pos_sample)[[2]] <- paste0("alpha", 1:L)
  
  ## MCMC iterations
  for(r in 1:mc){
    # p (sell probability)
    u <- matrix(rgamma(n*L, shape=t(x)+al, rate=1), L, n)
    p <- u / c(colSums(u)%x%rep(1, L))
    
    # z (latent variable)
    z <- rgamma(n, shape=sum(al), rate=1)
    
    # v (latent variable)
    v <- rgamma(L, shape=n*al, rate=n*al^2)
    
    # rho (latent variable)
    rho <- matrix(rbeta(L*(n-1), shape1=al+(((2:n)-1)/n)%x%rep(1, L), shape2=((n-(2:n)+1)/n)%x%rep(1, L)), L, n-1)
    
    # alpha
    p_al <- n + a
    a_al <- n*v
    b_al <- (-1)*rowSums(log(1/p)) + sum(log(z)) + 2*n + n*log(v) - rowSums(log(1/rho)) - b
    pab <- cbind(p_al, a_al, b_al)
    al_proposal <- apply(pab, 1, function(x) PTN(p=x[1], a=x[2], b=x[3]))
    log_uniform <- log( runif(L) )
    I1 <- (log_uniform <= (-2)/(12*n*al_proposal))
    I2 <- (log_uniform[!I1] <= 2*MHratio(new=n*al_proposal[!I1], old=n*al[!I1]))
    al[I1] <- al_proposal[I1]
    al[!I1][I2] <- al_proposal[!I1][I2]
    al[!I1][!I2] <- al[!I1][!I2]
    Pos_sample[r,] <- al
  }
  
  ## Summary
  Pos_sample <- Pos_sample[-(1:burn),]
  return(Pos_sample)
}







###------------------------------------------------------------------###
###      MCMC algorithm for Dirichlet-Multinomial distribution       ###
###              (ERG method by He et al. (2021))                    ###
###------------------------------------------------------------------###
## INPUT 
# x: (n,L)-matrix
# prior_alpha: two-dimensional vector of shape and rate of gamma prior for 'alpha' 
# mc: MCMC length 
# burn: length of burn-in period 
## OUTPUT
# matrix of posterior samples of 'alpha'

ERG_sampler <- function(x, prior_alpha=NULL, M=3, mc=5000, burn=1000){
  ## preparation 
  L <- dim(x)[2]
  n <- dim(x)[1]
  if(is.null(prior_alpha)){  prior_alpha <- c(1/L, 1) }
  a <- prior_alpha[1]
  b <- prior_alpha[2]
  
  ## function 
  ERG <- function(a, b, c, M){
    sum( myrgig(M-1, Lambda=(-3)/2, Chi=1/(2*(1:(M-1))^2), Psi=2*c^2) ) + rgamma(1, shape=a, rate=b)
  }
  
  ## initial values 
  ta <- om <- matrix(NA, L, n)
  eta <- rep(NA, n)
  al <- rep(1, L)
  ga_E <- (-digamma(1))    # constant used in MCMC
  
  ## matrices to store posterior samples
  Pos_sample <- matrix(NA, mc, L)
  dimnames(Pos_sample)[[2]] <- paste0("alpha", 1:L)
  
  ## MCMC iterations
  for(r in 1:mc){
    # ta (unnormalized sell probability)
    ta <- matrix(rgamma(L*n, shape=t(x)+al, rate=1), L, n)
    
    # om (latent variable)
    C1 <- digamma(M + al)
    C2 <- digamma(M)
    b_om <- 2 * al^2 * (C1 - C2) / (C1 - C2 - al * trigamma(M + al))
    a_om <- b_om * (C1 - C2) / (2 * al)
    ABC <- array(NA, dim = c(L, n, 3))
    ABC[,,1] <- outer(a_om, rep(1, n))
    ABC[,,2] <- outer(b_om, rep(1, n))
    ABC[,,3] <- outer(al, rep(1, n))
    om <- apply(ABC, c(1, 2), function(x) ERG(a=x[1], b=x[2], c=x[3], M=M))
    
    # eta (latent variable)
    eta <- rbeta(n, shape1=sum(al), shape2=N)
    
    # alpha
    p_al <- n + a
    a_al <- rowSums(om)
    b_al <- n*ga_E - b + rowSums( log(ta*outer(rep(1, L), eta)) )
    pab <- cbind(p_al, a_al, b_al)
    al <- apply(pab, 1, function(x) PTN(p = x[1], a = x[2], b = x[3]))
    Pos_sample[r, ] <- al
  }
  
  ## Summary
  Pos_sample <- Pos_sample[-(1:burn),]
  return(Pos_sample)
}

