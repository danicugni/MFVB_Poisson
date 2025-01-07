lposterior.si.sigma <- function(param, dati, hyp) {
  #hyp = (alpha, delta)
  #param = (beta, sigma2)
  x <- dati[,-1]
  y <- dati[,1]
  p <- ncol(dati) - 1
  alpha <- hyp[1]
  delta <- hyp[2]
  beta <- param[1:p]
  sigma2 <- param[p+1]
  
  x.beta.y.sum <- as.vector(t(y) %*% (x %*% beta))
  exp.x.beta <- sum(as.vector(exp(x %*% beta)))
  
  if(sigma2 <= 0) {
    return(- Inf)
  }
  else{
    return(-(alpha + p/2 + 1)*log(sigma2) - delta/sigma2 - 
             0.5/sigma2 * as.vector(t(beta) %*% beta) + 
             x.beta.y.sum - exp.x.beta)
  }
}

#MCMC algorithm
MCMC.si.sigma <- function(dati, hyp, R, start, eps) {
  
  start.time <- Sys.time()
  
  n <- nrow(dati)
  p <- ncol(dati) - 1
  accepted <- 0
  
  #============= Matrix of results =============
  out.beta <- matrix(NA, nrow = R, ncol = p)
  out.sigma <- rep(NA, R)
  
  #============= Hyperparameters' setting =============
  alpha <- hyp[1]
  delta <- hyp[2]
  
  #============= Defining parameters =============
  x.beta <- start[1:p]
  x.sigma <- start[p+1]
  
  #============= MCMC algorithm =============
  
  set.seed(1)
  require(mvtnorm)
  for(i in 1:R) {
    
    #drawn sigma2
    x.sigma <- 1/rgamma(1, alpha + p/2, delta + 0.5*
               sum(x.beta*x.beta))
    
    #drawn beta
    xstar.beta <- drop(rmvnorm(1, mean = x.beta, 
                  sigma = eps*solve(post.Jhat.si.sigma[1:p, 1:p])))
    bound <- exp(lposterior.si.sigma(c(xstar.beta, x.sigma), dati,  
             hyp) - lposterior.si.sigma(c(x.beta, x.sigma), dati, 
             hyp))
    if(runif(1) < bound) {
      x.beta <- xstar.beta
      accepted <- accepted + 1
    }
    
    #returning updated values for beta and sigma2
    out.sigma[i] <- x.sigma
    out.beta[i,] <- x.beta
  }
  
  end.time <- Sys.time()
  
  #============= Computing computational effort =============
  
  total.time <- end.time - start.time
  
  return(list(beta.values = out.beta, sigma.values = out.sigma, 
              tasso = accepted/R, 
              tempo = total.time[[1]]))
}