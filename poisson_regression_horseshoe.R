library(msos)

#Natural fixed-point iteration for the update of optimal
#variational density q.beta
update.q.beta.horseshoe <- function(y, x, x.y.sum, mu.q.beta, 
                                    sigma2.q.beta, mu.q.1.lambda2,
                                    mu.q.1.tau2, mu.q.log.tau2, 
                                    mu.q.log.lambda2, n, p, 
                                    eps_beta, maxIter_beta, 
                                    Trace_beta) {
  
  #============= Defining history variables =============
  elbo.out.beta <- numeric()
  mu.q.beta.out <- matrix(NA, nrow = maxIter_beta, ncol = p)
  sigma2.q.beta.out <- array(NA, dim = c(p,p, maxIter_beta))
  
  #============= Natural fixed-point iteration =============
  z <- 1
  cond.beta <- FALSE
  while(cond.beta != TRUE & z <= maxIter_beta) {
      
    #updating equations of q.mu
    exp.q.beta <- exp(x %*% mu.q.beta + 0.5 * diag(x %*% 
        sigma2.q.beta %*% t(x)))
    nu.q.beta <- x.y.sum - t(x) %*% exp.q.beta -
        mu.q.beta*mu.q.1.lambda2*mu.q.1.tau2
    sigma2.q.beta <- - solve(- t(x) %*% 
        diag(as.vector(exp.q.beta)) %*% x - mu.q.1.lambda2*
        mu.q.1.tau2)
    mu.q.beta <- mu.q.beta + sigma2.q.beta %*% nu.q.beta
      
    #computing beta-localized lower bound
    entropy.q.beta <- p/2* log(2*pi) + 0.5 *
      logdet(sigma2.q.beta)  + 0.5 * p
    lprior.beta <- -p/2*log(2*pi) -p/2*mu.q.log.tau2 + 
      0.5*sum(-mu.q.log.lambda2 - (mu.q.beta^2 + 
      diag(sigma2.q.beta))*mu.q.1.lambda2*mu.q.1.tau2)
    logL <- t(y) %*% (x %*% mu.q.beta) - sum(lfactorial(y)) -
      (t(rep(1,n)) %*% exp(x %*% mu.q.beta + 0.5 * 
      diag(x %*% sigma2.q.beta %*% t(x))))
    non.entropy.q.beta <- lprior.beta + logL
    ELBO.beta <- entropy.q.beta + non.entropy.q.beta
    
    #updating histories
    mu.q.beta.out[z,] <- mu.q.beta
    sigma2.q.beta.out[,,z] <- sigma2.q.beta
    elbo.out.beta <- c(elbo.out.beta, ELBO.beta)
    
    #stopping criterion
    if(z > 1) {
      delta.beta <- (abs((elbo.out.beta[z] - elbo.out.beta[z-1])
                         /elbo.out.beta[z-1]))
      
      if(delta.beta < eps_beta) cond.beta <- 1
      if (Trace_beta == 1) {
        print(paste0("iteration:", z, 
              "- beta lower bound increase: ", delta.beta))
      }
    }
    if (z > maxIter_beta) cond.beta <- 1
    z <- z + 1
  }
  
  return(list(mu.q.beta = mu.q.beta.out[z-1,], 
              sigma2.q.beta = sigma2.q.beta.out[,,z-1],
              entropy.q.beta = entropy.q.beta, 
              lprior.beta = lprior.beta, logL = logL))
}

#MFVB algorithm
MFVB.horseshoe <- function(y, x, hyp, start, eps_beta = 1e-4, 
                           eps_ELBO = 1e-4, maxIter_global = 100, 
                           maxIter_beta = 10, Trace_global = 0, 
                           Trace_beta = 0) {
  
  start.time <- Sys.time()
  
  n <- nrow(x)
  p <- ncol(x)
  x.y.sum <- colSums(y * x)
  
  #============= Hyperparameters' setting =============
  mu.beta <- hyp[[1]]
  a.eta <- hyp[[2]]
  b.eta <- hyp[[3]]
  a.tau2 <- hyp[[4]]
  a.nu <- hyp[[5]]
  b.nu <- hyp[[6]]
  a.lambda2 <- hyp[[7]]
  
  
  #============= Defining history variables =============
  elbo.out.global <- numeric()
  mu.q.beta.global.out <- matrix(NA, nrow = maxIter_global, 
                                 ncol = p)
  sigma2.q.beta.global.out <- array(NA, dim = c(p,p,maxIter_global))
  a.q.eta.out <- 1
  b.q.eta.out <- rep(NA, maxIter_global)
  a.q.nu.out <- rep(1, p)
  b.q.nu.out <- matrix(NA, nrow = maxIter_global, ncol = p)
  a.q.tau2.out <- (p+1)/2
  b.q.tau2.out <- matrix(NA, nrow = maxIter_global, ncol = p)
  a.q.lambda2.out <- rep(1, p)
  b.q.lambda2.out <- matrix(NA, nrow = maxIter_global, ncol = p)
  
  #============= Initializing optimal parameters =============
  mu.q.beta <- start[[1]]
  sigma2.q.beta <- start[[2]]
  omega.q.beta <- solve(sigma2.q.beta)
  a.q.eta <- a.q.nu <-  a.q.lambda2 <- 1
  a.q.tau2 <- (p+1)/2
  b.q.eta <- start[[3]]
  b.q.nu <- start[[4]]
  b.q.tau2 <- start[[5]]
  b.q.lambda2 <- start[[6]]
  mu.q.1.tau2 <- a.q.tau2/b.q.tau2
  mu.q.1.lambda2 <- a.q.lambda2/b.q.lambda2
  mu.q.1.nu <- a.q.nu/b.q.nu
  mu.q.log.eta <- log(b.q.eta) - digamma(a.q.eta)
  mu.q.log.nu <- log(b.q.nu) - digamma(a.q.nu)
  mu.q.log.tau2 <- log(b.q.tau2) - digamma(a.q.tau2)
  mu.q.log.lambda2 <- log(b.q.lambda2) - digamma(a.q.lambda2)
  
  
  #============= MFVB algorithm =============
  
  i <- 1
  cond.global <- FALSE
  while(cond.global != TRUE & i <= maxIter_global) { 
    
    
    #updating q.beta
    q.mu.update <- update.q.beta.horseshoe(y, x, x.y.sum, mu.q.beta, 
                                           sigma2.q.beta, 
                                           mu.q.1.lambda2,
                                           mu.q.1.tau2,
                                           mu.q.log.tau2, 
                                           mu.q.log.lambda2,
                                           n, p, eps_beta, 
                                           maxIter_beta, Trace_beta)
    mu.q.beta <- q.mu.update[[1]]
    sigma2.q.beta <- q.mu.update[[2]]
    
    #updating q.eta
    b.q.eta <- 1 + a.q.tau2/b.q.tau2
    mu.q.1.eta <- a.q.eta/b.q.eta
    mu.q.log.eta <- log(b.q.eta) - digamma(a.q.eta)
    
    #updating q.nu
    b.q.nu <- 1 + a.q.lambda2/b.q.lambda2
    mu.q.1.nu <- a.q.nu/b.q.nu
    mu.q.log.nu <- log(b.q.nu) - digamma(a.q.nu)
    
    #updating q.tau2
    b.q.tau2 <- a.q.eta/b.q.eta + 0.5*(sum((diag(sigma2.q.beta) + 
                mu.q.beta^2)*a.q.lambda2/b.q.lambda2))
    mu.q.1.tau2 <- a.q.tau2/b.q.tau2
    mu.q.log.tau2 <- log(b.q.tau2) - digamma(a.q.tau2)
    
    #updating q.lambda2
    if(p == 1)
      b.q.lambda2 <- a.q.nu/b.q.nu + 0.5*(sigma2.q.beta + 
                     mu.q.beta^2)*a.q.tau2/b.q.tau2
    else
      b.q.lambda2 <- a.q.nu/b.q.nu + 0.5*(diag(sigma2.q.beta) + 
                     mu.q.beta^2)*a.q.tau2/b.q.tau2
    mu.q.1.lambda2 <- a.q.lambda2/b.q.lambda2
    mu.q.log.lambda2 <- log(b.q.lambda2) - digamma(a.q.lambda2)
    
    #computing lower bound
    entropy.q.beta <- q.mu.update[[3]]
    
    entropy.q.lambda2 <- sum(a.q.lambda2 + log(b.q.lambda2) +
                         lgamma(a.q.lambda2) - (a.q.lambda2 +1)*
                         digamma(a.q.lambda2))
    
    entropy.q.tau2 <- a.q.tau2 + log(b.q.tau2) +  lgamma(a.q.tau2) -
                      (a.q.tau2 +1)*digamma(a.q.tau2)
    
    entropy.q.eta <- a.q.eta + log(b.q.eta) + lgamma(a.q.eta) -
                     (a.q.eta +1)*digamma(a.q.eta)
    
    entropy.q.nu <- sum(a.q.nu + log(b.q.nu) +  lgamma(a.q.nu) -
                    (a.q.nu +1)*digamma(a.q.nu))
    
    lprior.beta <- -p/2*log(2*pi) -p/2*mu.q.log.tau2 - 
                   0.5*sum(mu.q.log.lambda2)- 
                   0.5*sum((diag(sigma2.q.beta) 
                   + mu.q.beta^2)*mu.q.1.lambda2*mu.q.1.tau2)
    
    logL <- q.mu.update[[5]]
    
    lprior.lambda2 <- -p*lgamma(0.5) - 0.5*sum(0.5*mu.q.log.nu + 
                      3/2*mu.q.log.lambda2 + 
                      mu.q.1.nu*mu.q.1.lambda2)
    
    lprior.tau2 <- -0.5*mu.q.log.eta - lgamma(0.5) 
                   - 3/2*mu.q.log.tau2 - mu.q.1.eta*mu.q.1.tau2
    
    lprior.eta <- -lgamma(0.5) - 3/2*mu.q.log.eta - mu.q.1.eta
    
    lprior.nu <- -p*lgamma(0.5) - sum(3/2*mu.q.log.nu + mu.q.1.nu)
    
    ELBO.global <- entropy.q.beta +  entropy.q.lambda2 + 
                   entropy.q.tau2 + entropy.q.eta + entropy.q.nu + 
                   logL + lprior.beta + lprior.lambda2  + 
                   lprior.tau2 + lprior.eta + lprior.nu
    
    #updating history variables
    elbo.out.global <- c(elbo.out.global, ELBO.global)
    mu.q.beta.global.out[i,] <- mu.q.beta 
    sigma2.q.beta.global.out[,,i] <- sigma2.q.beta
    b.q.eta.out[i] <- b.q.eta
    b.q.nu.out[i,] <- b.q.nu
    b.q.tau2.out[i] <- b.q.tau2
    b.q.lambda2.out[i,] <- b.q.lambda2
    
    #stopping criterion
    if(i > 1) {
      Delta.global <- (abs((elbo.out.global[i] - 
                elbo.out.global[i-1])/elbo.out.global[i-1]))
      if(Delta.global < eps_global) cond.global <- 1
      if (Trace_global == 1) {
        print(paste0("iteration:", i, 
              "- global lower bound increase: ",Delta.global))
      }
    }
    if (i > maxIter_global) cond.global <- 1
    i <- i + 1
  }
  
  end.time <- Sys.time()
  
  #============= Computing computational effort =============
  total.time <- end.time - start.time
  
  return(list(elbo = elbo.out.global, 
              mu.q.beta =as.matrix(mu.q.beta.global.out[1:(i-1),]), 
              sigma2.q.beta=sigma2.q.beta.global.out[,,1:(i-1)], 
              a.q.eta = a.q.eta.out,
              b.q.eta = b.q.eta.out[1:(i-1)], a.q.nu = a.q.nu.out, 
              b.q.nu = as.matrix(b.q.nu.out[1:(i-1),]),
              a.q.tau2 = a.q.tau2.out, b.q.tau2 = b.q.tau2[1:(i-1)],
              a.q.lambda2 = a.q.lambda2.out, 
              b.q.lambda2 = as.matrix(b.q.lambda2.out[1:(i-1),]),
              iter = i-1, tempo = total.time[[1]]))
}
