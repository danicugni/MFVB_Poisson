multivariate.gaussian.MCMC <- function(R, y, hyp, start) {
  
  start.time <- Sys.time()
  
  n <- nrow(y)
  d <- ncol(y)
  y.mean <- apply(y, 2, mean)
  
  #============= Matrix of results =============
  out.mu <- matrix(NA, nrow = R, ncol = d)
  out.prec <- array(NA, dim = c(d,d,R))
  
  #============= Hyperparameters' setting =============
  mu0 <- hyp[[1]]
  prec0 <- hyp[[2]]
  nu <- hyp[[3]]
  v <- hyp[[4]]
  
  #constant hyperparameter in the posterior
  post.nu <- nu + n
  uno.v <- solve(v)
  prec0.nu0 <- prec0 %*% mu0
  
  #============= Defining parameters =============
  xstar.mu <- start[[1]]
  xstar.prec <- start[[2]]
  
  #============= MCMC algorithm =============
  
  require(mvtnorm)
  set.seed(3)
  for(i in 1:R) {
    
    #drawn precision
    #y.mu.mean <- t(apply(y, 1, function(x) x - xstar.mu))
    y.mu.mean <- t(y) - xstar.mu
    #y.mean.second.sum <- t(y.mu.mean) %*% y.mu.mean
    y.mean.second.sum <- y.mu.mean %*% t(y.mu.mean)
    post.v <- solve(uno.v + y.mean.second.sum)
    xstar.prec <- rWishart(n = 1, df = post.nu, Sigma = post.v)[,,1]
    
    #drawn mu
    post.var.mu <- solve(n* xstar.prec + prec0)
    post.mean.mu <- as.vector(post.var.mu %*% (n* xstar.prec 
                                        %*% y.mean + prec0.nu0))
    xstar.mu <- drop(mvtnorm::rmvnorm(n = 1, mean = post.mean.mu, 
                                      sigma = post.var.mu))
    
    #returning updated values for mu and precision
    out.mu[i,] <- xstar.mu
    out.prec[,,i] <- xstar.prec
  }
  
  end.time <- Sys.time()
  
  #============= Computing computational effort =============
  
  total.time <- end.time - start.time
  
  return(list(mu.values = out.mu, prec.values = out.prec, 
              tempo = total.time[[1]]))
}
