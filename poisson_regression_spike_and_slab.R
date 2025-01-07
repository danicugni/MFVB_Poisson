library(msos)

#Natural fixed-point iteration for the update of optimal
#variational density q.beta
update.q.beta.spike <- function(y, x, x.y.sum, mu.q.beta, 
                                sigma2.q.beta, mu.q.gamma, 
                                mu.q.1.lambda1, mu.q.log.lambda1, 
                                lambda0,   n, p, eps_beta,
                                maxIter_beta, Trace_beta) {
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
                 (1-mu.q.gamma)*mu.q.beta/lambda0 -
                 mu.q.gamma*mu.q.beta*mu.q.1.lambda1
    sigma2.q.beta <-  solve((1-mu.q.gamma)/lambda0 +
                                 mu.q.gamma*mu.q.1.lambda1+
                                 t(x) %*% diag(as.vector(exp.q.beta)) 
                             %*% x)
    mu.q.beta <- mu.q.beta + sigma2.q.beta %*% nu.q.beta
    
    #computing beta-localized lower bound
    entropy.q.beta <- 0.5*p* log(2*pi) + 0.5 *logdet(sigma2.q.beta)+ 
                      0.5 * p
    lprior.beta <- -p/2*log(2*pi) + sum(-(1-mu.q.gamma)*
                   (log(lambda0)/2 + (mu.q.beta^2 + 
                   diag(sigma2.q.beta))/(2*lambda0)))+ 
                   sum(-mu.q.gamma* 0.5*(mu.q.log.lambda1 + 
                   (mu.q.beta^2 + diag(sigma2.q.beta))/(2)* 
                   mu.q.1.lambda1))
    logL <- t(y) %*% (x %*% mu.q.beta) - sum(lfactorial(y)) -
      (t(rep(1,n)) %*% exp(x %*% mu.q.beta + 0.5 * diag(x %*% 
      sigma2.q.beta %*% t(x))))
    non.entropy.q.beta <- lprior.beta + logL
    ELBO.beta <- entropy.q.beta + non.entropy.q.beta
    
    #updating histories
    mu.q.beta.out[z,] <- mu.q.beta
    sigma2.q.beta.out[,,z] <- sigma2.q.beta
    elbo.out.beta <- c(elbo.out.beta, ELBO.beta)
    
    #stopping criterion
    if(z > 1) {
      delta.beta <- (abs((elbo.out.beta[z] - elbo.out.beta[z-1])/
                    elbo.out.beta[z-1]))
      if(delta.beta < eps_beta) cond.beta <- 1
      if (Trace_beta == 1) {
        print(paste0("iteration:", z, 
              "- beta lower bound increase: ",delta.beta))
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
MFVB.spike.slab <- function(y, x, hyp, start, eps_beta, eps_global, 
                            maxIter_beta, maxIter_gamma, 
                            maxIter_global, Trace_beta = 0, 
                            Trace_global = 0) {
  
  n <- nrow(x)
  p <- ncol(x) 
  x.y.sum <- colSums(y * x)
  
  #============= Hyperparameters' setting =============
  mu.beta <- hyp[[1]][[1]]
  r.lambda1 <- hyp[[2]][[1]]
  delta.lambda1 <- hyp[[2]][[2]]
  a.theta <- hyp[[3]][[1]]
  b.theta <- hyp[[3]][[2]]
  lambda0 <- hyp[[4]]
  
  #============= Defining history variables =============
  elbo.out.global <- numeric()
  mu.q.beta.global.out <- matrix(NA, nrow = maxIter_global, 
                          ncol = p)
  sigma2.q.beta.global.out <- array(NA, dim =c(p,p,maxIter_global))
  mu.q.gamma.global.out <- matrix(NA, nrow = maxIter_global,
                           ncol = p)
  r.q.lambda1.out <- matrix(NA, nrow = maxIter_global, ncol = p)
  delta.q.lambda1.out <- matrix(NA, nrow = maxIter_global, ncol = p)
  a.q.theta.out <- numeric()
  b.q.theta.out <- numeric()
  
  #============= Initializing optimal parameters =============
  mu.q.beta <- start[[1]]
  sigma2.q.beta <- start[[2]]
  omega.q.beta <- solve(sigma2.q.beta)
  mu.q.gamma <- start[[3]]
  w <- log(mu.q.gamma / (1-mu.q.gamma))
  r.q.lambda1 <- start[[4]]
  delta.q.lambda1 <- start[[5]]
  a.q.theta <- start[[6]]
  b.q.theta <- start[[7]]
  mu.q.log.theta <- digamma(a.q.theta) - digamma(a.q.theta + 
                    b.q.theta)
  mu.q.log.1.theta <- digamma(b.q.theta) - digamma(a.q.theta + 
                      b.q.theta)
  mu.q.1.lambda1 <- r.q.lambda1/delta.q.lambda1
  mu.q.log.lambda1 <- log(delta.q.lambda1) - digamma(r.q.lambda1)
  
  
  #============= MFVB algorithm =============
  start.time <- Sys.time()
  
  i <- 1
  cond.global <- FALSE
  while(cond.global != TRUE & i <= maxIter_global) { 
    
    #updating q.beta
    q.mu.update <- update.q.beta.spike(y, x, x.y.sum, mu.q.beta, 
                                       sigma2.q.beta, mu.q.gamma,
                                       mu.q.1.lambda1,
                                       mu.q.log.lambda1,lambda0, 
                                       n, p,eps_beta, 
                                       maxIter_beta, Trace_beta)
    mu.q.beta <- q.mu.update[[1]]
    sigma2.q.beta <- q.mu.update[[2]]
    
    #updating q.gamma
    for(j in 1:p){
      w[j] <- mu.q.log.theta  - mu.q.log.1.theta - 
        (mu.q.beta[j]^2 + diag(sigma2.q.beta)[j])*
        (mu.q.1.lambda1[j]-1/lambda0)/2-
        0.5*(mu.q.log.lambda1[j] - log(lambda0))
      mu.q.gamma[j] <- 1/(1 + exp(-w[j]))
    }
    mu.q.gamma[which(mu.q.gamma == 1.)] <- 0.99
    mu.q.gamma[which(mu.q.gamma == 0.)] <- 0.01
    
    #updating q.lambda
    r.q.lambda1 <- r.lambda1 + mu.q.gamma/2
    if(p == 1)
      delta.q.lambda1 <- delta.lambda1 + 
      mu.q.gamma*(mu.q.beta^2 + sigma2.q.beta)/(2)
    else
      delta.q.lambda1 <- delta.lambda1 + 
      mu.q.gamma*(mu.q.beta^2 + diag(sigma2.q.beta))/(2)
    mu.q.1.lambda1 <- r.q.lambda1/delta.q.lambda1
    mu.q.log.lambda1 <- log(delta.q.lambda1) - digamma(r.q.lambda1)
    
    #updating q.theta
    a.q.theta <- a.theta + sum(mu.q.gamma)
    b.q.theta <- b.theta + p - sum(mu.q.gamma)
    mu.q.log.theta <- digamma(a.q.theta) - digamma(a.q.theta + 
                      b.q.theta)
    mu.q.log.1.theta <- digamma(b.q.theta) - digamma(a.q.theta + 
                        b.q.theta)
    
    #computing lower bound
    entropy.q.beta <- q.mu.update[[3]]
    
    entropy.q.gamma <- - sum(mu.q.gamma * log(mu.q.gamma) + 
                       (1-mu.q.gamma) * log(1-mu.q.gamma))
    
    entropy.q.lambda1 <- sum(-r.q.lambda1*log(delta.q.lambda1) +
                         lgamma(r.q.lambda1)+ delta.q.lambda1*
                         mu.q.1.lambda1+ (r.q.lambda1+1)*
                         mu.q.log.lambda1)
    
    entropy.q.theta <- -(a.q.theta - 1)* mu.q.log.theta -
                       (b.q.theta - 1)*mu.q.log.1.theta + 
                       lbeta(a.q.theta, b.q.theta)
    
    lprior.beta <- -p/2*log(2*pi)+ sum(-(1-mu.q.gamma)*
                   (log(lambda0)/2 + (mu.q.beta^2 + 
                   diag(sigma2.q.beta))/(2*lambda0)))+
                   sum(-(mu.q.gamma)* 0.5*mu.q.log.lambda1 *
                  (mu.q.beta^2 + diag(sigma2.q.beta))/(2)* 
                  mu.q.1.lambda1)
    
    lprior.gamma <- sum((1-mu.q.gamma)*mu.q.log.1.theta+ 
                    mu.q.gamma*mu.q.log.theta) 
    
    lprior.theta <- (a.theta - 1)*mu.q.log.theta + (b.theta - 1)*
                    mu.q.log.1.theta - lbeta(a.theta, b.theta)
    
    lprior.lambda1 <- p*r.lambda1[1]*log(delta.lambda1[1]) - 
                      p*lgamma(r.lambda1[1]) - delta.lambda1[1]*
                      sum(mu.q.1.lambda1)- (r.lambda1[1] + 1)*
                      sum(mu.q.log.lambda1)
    
    logL <- q.mu.update[[5]]
    
    ELBO.global <- entropy.q.beta + entropy.q.gamma + 
                   entropy.q.lambda1 + entropy.q.theta +
                   lprior.beta + lprior.gamma + logL + 
                   lprior.lambda1 + lprior.theta 
    
    #updating history variables
    elbo.out.global <- c(elbo.out.global, ELBO.global)
    mu.q.beta.global.out[i,] <- mu.q.beta 
    sigma2.q.beta.global.out[,,i] <- sigma2.q.beta
    mu.q.gamma.global.out[i,] <- mu.q.gamma
    r.q.lambda1.out[i,] <- r.q.lambda1
    delta.q.lambda1.out[i,] <- delta.q.lambda1
    a.q.theta.out <- c(a.q.theta.out, a.q.theta)
    b.q.theta.out <- c(b.q.theta.out, b.q.theta)
    
    #stopping criterion
    if(i > 1) {
      Delta.global <- (abs((elbo.out.global[i] - 
                      elbo.out.global[i-1])/
                      elbo.out.global[i-1]))
      if(Delta.global < eps_global) cond.global <- 1
      if (Trace_global == 1) {
        print(paste0("iteration:", i, 
              "- global lower bound increase: ", Delta.global))
      }
    }
    if (i > maxIter_global) cond.global <- 1
    i <- i + 1
  }
  
  end.time <- Sys.time()
  
  #============= Computing computational effort =============
  total.time <- end.time - start.time
  
  return(list(elbo = elbo.out.global, 
              mu.q.beta=as.matrix(mu.q.beta.global.out[1:(i-1),]), 
              mu.q.gamma=as.matrix(mu.q.gamma.global.out[1:(i-1),]),
              sigma2.q.beta=sigma2.q.beta.global.out[,,1:(i-1)], 
              a.q.theta = a.q.theta.out,
              b.q.theta = b.q.theta.out, 
              r.q.lambda1=as.matrix(r.q.lambda1.out[1:(i-1),]),
              delta.q.lambda1=as.matrix(delta.q.lambda1.out[1:(i-1),]), 
              iter = i-1, tempo = total.time[[1]]))
}
