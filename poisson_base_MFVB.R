library(msos)

#Natural fixed-point iteration for the update of optimal
#variational density q.beta
update.local.PMFVB.si.sigma <- function(x, y, x.y.sum, mu.q.beta, 
                                        sigma2.q.beta, 
                                        alpha.q.sigma,
                                        delta.q.sigma, mu.q.1.sigma,
                                        n, p, eps_local, 
                                        maxIter_local, 
                                        Trace_local) {
  #============= Defining history variables =============
  elbo.out.local <- numeric()
  mu.q.beta.out <- matrix(NA, nrow = maxIter_local, ncol = p)
  sigma2.q.beta.out <- array(NA, dim = c(p,p,maxIter_local))
  
  #============= Natural fixed-point iteration =============
  j <- 1
  cond2 <- FALSE
  
  while(cond2 != TRUE & j <= maxIter_local) {
    
    #updating equations of q.mu
    log.exp.mu.q.beta <- as.vector(x %*% as.matrix(mu.q.beta) 
                                   +0.5 * diag(x %*% 
                                    as.matrix(sigma2.q.beta) %*% 
                                    t(x)))
    exp.mu.q.beta <- as.vector(exp(log.exp.mu.q.beta))
    x.exp.mu.q.beta.d1 <- t(x) %*% as.matrix(exp.mu.q.beta)
    nu.q.beta <- -mu.q.beta*(alpha.q.sigma/delta.q.sigma) + 
      x.y.sum - x.exp.mu.q.beta.d1
    sigma2.q.beta <- solve( diag(1,p)*mu.q.1.sigma + t(x) %*%
                              diag(as.vector(exp.mu.q.beta)) %*% x)
    mu.q.beta <- mu.q.beta + sigma2.q.beta %*% nu.q.beta
    
    
    #computing beta-localized lower bound
    entropy.q.beta <- 0.5*p* log(2*pi) + 0.5 * logdet(sigma2.q.beta) 
                      + 0.5 * p
    lprior.beta <- -p/2 * (log(2*pi*delta.q.sigma) - 
                    digamma(alpha.q.sigma))- 0.5*alpha.q.sigma/
                    delta.q.sigma * (sum(diag(sigma2.q.beta)) + 
                    t(mu.q.beta) %*% mu.q.beta)
    logL <- t(y) %*% (x %*% mu.q.beta) - sum(lfactorial(y)) - 
      (t(rep(1,n)) %*% as.matrix(exp.mu.q.beta))
    non.entropy.q.beta <- lprior.beta + logL
    ELBO.local <- entropy.q.beta + non.entropy.q.beta
    
    #updating histories
    mu.q.beta.out[j,] <- mu.q.beta
    sigma2.q.beta.out[,,j] <- sigma2.q.beta
    elbo.out.local <- c(elbo.out.local, ELBO.local)
    
    
    #stopping criterion
    if(j > 1) {
      delta2 <- (abs((elbo.out.local[j] - elbo.out.local[j-1])/
                       elbo.out.local[j-1]))

      if(delta2 < eps_local) cond2 <- 1
      if (Trace_local == 1) {
        print(paste0("iteration:", j, "- local lower bound 
                     increase: ", delta2))
      }
    }
    if (j > maxIter_local) cond2 <- 1
    j <- j + 1
    
  }
  return(list(mu.q.beta = mu.q.beta.out[j-1,], 
              sigma2.q.beta = sigma2.q.beta.out[,,j-1],
              entropy.q.beta = entropy.q.beta, logL = logL))
}


#MFVB algorithm
PMFVB.si.sigma <- function(y, x,  hyp, start, eps_local, eps_global, 
                           maxIter_local, maxIter_global, 
                           Trace_local = 0, Trace_global = 0) {
  
  start.time <- Sys.time()
  
  n <- nrow(x)
  p <- ncol(x)
  x.y.sum <- colSums(y * x)
  
  #============= Hyperparameters' setting =============
  mu.beta <- hyp[[1]]
  alpha.sigma <- hyp[[2]]
  delta.sigma <- hyp[[3]]
  
  #============= Defining history variables =============
  elbo.out.global <- numeric()
  entropy.q.beta.out <- numeric()
  mu.q.beta.out <- matrix(NA, nrow = maxIter_global, ncol = p)
  sigma2.q.beta.out <- array(NA, dim = c(p,p,maxIter_global))
  delta.q.sigma.out <- numeric()
  
  #============= Initializing optimal parameters =============
  mu.q.beta <- start[[1]]
  sigma2.q.beta <- start[[2]]
  alpha.q.sigma <- start[[3]]
  delta.q.sigma <- start[[4]]
  mu.q.1.sigma <- alpha.q.sigma/delta.q.sigma
  
  #============= MFVB algorithm =============
  
  i <- 1
  cond1 <- FALSE
  while(cond1 != TRUE & i <= maxIter_global) { 
    
    #updating q.beta
    q.mu.update <- update.local.PMFVB.si.sigma(x, y, x.y.sum, 
                                               mu.q.beta, 
                                               sigma2.q.beta, 
                                               alpha.q.sigma,
                                               delta.q.sigma, 
                                               mu.q.1.sigma,
                                               n, p, eps_local, 
                                               maxIter_local, 
                                               Trace_local)
    mu.q.beta <- q.mu.update[[1]]
    sigma2.q.beta <- q.mu.update[[2]]
    
    #updating di q.sigma
    alpha.q.sigma <- as.vector(alpha.sigma + p/2)
    delta.q.sigma <- as.vector(delta.sigma + 0.5*
                    (sum(diag(sigma2.q.beta)) + t(mu.q.beta) %*% 
                    mu.q.beta))
    mu.q.1.sigma <- alpha.q.sigma/delta.q.sigma
    
    #computing lower bound
    entropy.q.sigma <- alpha.q.sigma + log(delta.q.sigma * 
                       gamma(alpha.q.sigma)) - (alpha.q.sigma + 1)*
                       digamma(alpha.q.sigma)
    lprior.sigma <-  alpha.sigma*log(delta.sigma) - 
                     lgamma(alpha.sigma) -(alpha.sigma +1)*
                     (log(delta.q.sigma) + (alpha.sigma + 1)*
                     digamma(alpha.q.sigma)) - delta.sigma*
                     mu.q.1.sigma 
    entropy.q.beta <- q.mu.update[[3]]
    logL <- q.mu.update[[4]]
    lprior.beta <- -p/2 * (log(2*pi*delta.q.sigma) - 
                   digamma(alpha.q.sigma))
                   - 0.5*mu.q.1.sigma * 
                   (sum(diag(sigma2.q.beta)) + t(mu.q.beta) %*% 
                   mu.q.beta)
    ELBO.global <- entropy.q.beta + entropy.q.sigma + logL + 
                   lprior.beta  + lprior.sigma
    
    #updating history variables
    mu.q.beta.out[i,] <- mu.q.beta
    sigma2.q.beta.out[,,i] <- sigma2.q.beta
    delta.q.sigma.out <- c(delta.q.sigma.out, delta.q.sigma)
    elbo.out.global <- c(elbo.out.global, ELBO.global)
    
    #stopping criterion
    if(i > 1) {
      delta1 <- (abs((elbo.out.global[i] - elbo.out.global[i-1])/
                       elbo.out.global[i-1]))
      if(delta1 < eps_global) cond1 <- 1
      if (Trace_global == 1) {
        print(paste0("iteration:", i, "- global lower bound 
                     increase: ", delta1))
      }
    }
    if (i > maxIter_global) cond1 <- 1
    i <- i + 1
  }
  
  end.time <- Sys.time()
  
  #============= Computing computational effort =============
  
  total.time <- end.time - start.time
  
  return(list(elbo = elbo.out.global, 
              mu.q.beta = mu.q.beta.out[1:(i-1),],
              sigma2.q.beta =  sigma2.q.beta.out[,,1:(i-1)], 
              alpha.q.sigma = alpha.q.sigma,
              delta.q.sigma = delta.q.sigma.out, iter = i-1, 
              tempo = total.time[[1]]))
} 

