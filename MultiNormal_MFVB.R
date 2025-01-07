#Updating equations of q.mu
update.MFVB <- function(y, y.sum, prec0, mu0, n, d, mu.q.prec, 
                        nu.q.prec) {
  
  sigma.q.mu <- solve(n*mu.q.prec + prec0)
  mu.q.mu <- sigma.q.mu %*% (mu.q.prec %*% y.sum + prec0 %*% mu0) 
  y.mean <- t(apply(y, 1, function(x) x - mu.q.mu))
  y.mean.second <- apply(y.mean, 1, function(x) x %*% t(x))
  y.mean.second.sum <- matrix(rowSums(y.mean.second), nrow = d, 
                              ncol = d)
  v.q.prec <- solve(y.mean.second.sum + n*sigma.q.mu + solve(v))
  mu.q.prec <- nu.q.prec * v.q.prec
  return(list(mu.q.mu = mu.q.mu, sigma.q.mu = sigma.q.mu, 
              v.q.prec = v.q.prec, mu.q.prec = mu.q.prec))
}

#MFVB algorithm
multivariate.gaussian.MFVB <- function(y, hyp, start, maxIter = 100, 
                                       eps_ELBO = 1e-4, Trace = 0) {
  
  start.time <- Sys.time()

  n <- nrow(y)
  d <- ncol(y)
  y.sum <- apply(y, 2, sum)
  
  #============= Defining history variables =============
  elbo.out <- numeric()
  mu.q.mu.out <- matrix(NA, nrow = maxIter, ncol = d)
  sigma.q.mu.out <- array(NA, dim = c(d,d,maxIter))
  v.q.prec.out <- array(NA, dim = c(d,d,maxIter))

  #============= Hyperparameters' setting =============
  #hyperparameters (hyp = (mu0, prec0, nu, v))
  mu0 <- hyp[[1]]
  prec0 <- hyp[[2]]
  nu <- hyp[[3]]
  v <- hyp[[4]]
  
  #============= Initializing optimal parameters =============
  mu.q.mu <- start[[1]]
  sigma.q.mu <- solve(start[[2]])
  v.q.prec <- start[[4]]
  mu.q.prec <- start[[3]]*start[[4]]
  
  #constant optimal hyperparameter in the update
  nu.q.prec <- nu + n
  
  #============= MFVB algorithm =============
  require(CholWishart)
  i <- 1
  cond <- FALSE
  while(cond != TRUE & i <= maxIter) {
    
    #updating q.mu
    new.val <- update.MFVB(y, y.sum, prec0, mu0, n, d, mu.q.prec, 
                           nu.q.prec)
    mu.q.mu <-  new.val[[1]]
    sigma.q.mu <- new.val[[2]]
    
    #updating q.prec
    v.q.prec <- new.val[[3]]
    mu.q.prec <- new.val[[4]]
    
    #computing lower bound
    ELBO <- -n*d/2*log(2*pi) - 0.5* as.vector(t(mu.q.mu - mu0) %*% 
      prec0 %*% (mu.q.mu - mu0)) + 0.5*log(det(prec0)) + 
      0.5 * log(det(sigma.q.mu)) - 0.5 * sum(diag(sigma.q.mu * 
      prec0)) + d/2 - nu*d/2*log(2) - nu/2*log(det(v)) - 
      lmvgamma(nu/2, d) + nu.q.prec*d/2*log(2) + nu.q.prec/2*
      log(det(v.q.prec)) + lmvgamma(nu.q.prec/2, 2)
    
    #computing lower bound
    mu.q.mu.out[i,] <- mu.q.mu
    sigma.q.mu.out[,,i] <- sigma.q.mu
    v.q.prec.out[,,i] <- v.q.prec
    elbo.out <- c(elbo.out, ELBO)
    
    
    #stopping criterion
    if (i > 1) {
      #delta1 <- max(abs((parNew-parOld)/parOld))
      delta2 <- (abs((elbo.out[i] - elbo.out[i-1])/elbo.out[i-1]))
      
      #if ((delta1 < eps_Par)&(delta2 < eps_ELBO)) cond <- 1
      if(delta2 < eps_ELBO) cond <- 1
      if (Trace == 1) {
        #print(paste0("iteration: ",i," - parameter variation: ",
        #delta1," - lower bound increase: ",delta2))
        print(paste0("iteration:", i, "- lower bound increase: ", 
                     delta2))
      }
    }
    if (i > maxIter) cond <- 1
    i <- i + 1
  }
  
  end.time <- Sys.time()
  
  #============= Computing computational effort =============
  total.time <- end.time - start.time
  
  return(list(elbo = elbo.out, mu.mu = mu.q.mu.out[1:(i-1),], 
              sigma.mu =  sigma.q.mu.out[,,1:(i-1)], 
              nu.prec = nu.q.prec, v.prec = v.q.prec.out[,,1:(i-1)], 
              iter = i-1, tempo = total.time[[1]]))
}