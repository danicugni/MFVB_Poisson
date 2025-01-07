#Function for the creation of matrix (n, p) with each row composed 
#by the elements in the production for the update of optimal
#variational densities  q.beta.j and q.gamma.j
prod.elementi.start <- function(x, mu.q.gamma, mu.q.beta, 
                                sigma2.q.beta, n, p) { 
  matrice.prod <- matrix(NA, nrow = n, ncol = p)
  if(p == 1)
    matrice.prod <- as.vector(1 - mu.q.gamma + mu.q.gamma*exp(x *
                    mu.q.beta + sigma2.q.beta * x^2/2))
  else{
    for(j in 1:p)
      matrice.prod[,j] <- 1 - mu.q.gamma[j] + mu.q.gamma[j]*
        exp(x[,j] * mu.q.beta[j] + (sigma2.q.beta[j] * x[,j]^2)/
        2)
  }
  return(as.matrix(matrice.prod))
}

#Function for the change of the column j in the matrix of elements
#in the production for the update of q.beta.j and q.gamma.j
change.prod.elementi.j <- function(x, mu.q.gamma, mu.q.beta, 
                                   sigma2.q.beta, j, p) {
  if(p == 1)
    return(1-mu.q.gamma + mu.q.gamma*exp(x * mu.q.beta + 
           sigma2.q.beta * x^2/2))
  else
    return(1-mu.q.gamma[j] + mu.q.gamma[j]*exp(x[,j] * 
           mu.q.beta[j] + (sigma2.q.beta[j] * x[,j]^2)/2))
}


#Function for the computation of the production in the update of
#optimal variational densities q.beta.j and q.gamma.j
produttoria.q.beta.q.gamma <- function(elementi.prod, j ){
  elementi.prod.partial <- as.matrix(elementi.prod[,-j])
  return(apply(elementi.prod.partial, 1, function(x) prod(x))) 
}


#Natural fixed-point iteration for the update of optimal
#variational density q.beta
update.q.beta.j <- function(y, x, sigma2.beta, mu.q.beta, 
                            sigma2.q.beta, mu.q.gamma,
                          elementi.prod, maxIter_beta, eps_beta,
                          Trace_beta, j, nu.q.beta) {
  
  n <- nrow(x)
  p <- ncol(x)
  
  
  #============= Defining history variables =============
  elbo.out.beta <- c()
  mu.q.beta.out <- rep(NA, maxIter_beta)
  sigma2.q.beta.out <- rep(NA, maxIter_beta)
  nu.q.beta.out <- rep(NA,maxIter_beta)
    
  #============= Natural fixed-point iteration =============
  
    z <- 1
    cond.beta <- FALSE
    if(p == 1)
      produttoria.no.j <- rep(1,n)
    else
      produttoria.no.j<-produttoria.q.beta.q.gamma(elementi.prod,
                        j)
    while(cond.beta != TRUE & z <= maxIter_beta) {
      
      #updating equations of q.mu
      exp.q.beta.j <- mu.q.gamma[j] * exp(x[,j]*mu.q.beta[j] + 
                      sigma2.q.beta[j]*(x[,j]^2)/2)
      nu.q.beta[j] <- -mu.q.beta[j]/sigma2.beta  + t(x[,j] * 
                     mu.q.gamma[j]) %*% y - t(produttoria.no.j) %*% 
                    (x[,j] * exp.q.beta.j)
      sigma2.q.beta[j] <-  as.vector(solve( 1/sigma2.beta + 
                          t(produttoria.no.j) %*% 
                          (x[,j]^2 * exp.q.beta.j)))
      mu.q.beta[j] <- mu.q.beta[j] + sigma2.q.beta[j]*nu.q.beta[j]
      elementi.prod[,j] <- change.prod.elementi.j(x, mu.q.gamma, 
                           mu.q.beta, sigma2.q.beta, j, p) 
        
      #computing beta-localized lower bound
      
      produttoria.total <- apply(elementi.prod, 1, function(x) 
                           prod(x))
        
      entropy.q.beta.j <- 0.5*log(sigma2.q.beta[j]) + 1/2 + 
                          0.5*log(2*pi)
      lprior.beta.j <- - 0.5*log(2*pi) - 1/2*log(sigma2.beta) - 
                       1/(2*sigma2.beta)*(sigma2.q.beta[j] + 
                       mu.q.beta[j]^2)
      if(p == 1){
        logL <- as.vector( t(x * mu.q.gamma * mu.q.beta) %*% y - 
                sum(lfactorial(y)) - sum(elementi.prod))
      } else{
        logL <- as.vector(t(x %*% diag(mu.q.gamma) %*% mu.q.beta) 
                %*% y) - sum(lfactorial(y)) -
                sum(produttoria.total)
      }
      non.entropy.q.beta.j <- lprior.beta.j + logL
        
      ELBO.beta.j <- entropy.q.beta.j + lprior.beta.j + logL 
        
      #updating histories
      mu.q.beta.out[z] <- mu.q.beta[j]
      sigma2.q.beta.out[z] <- sigma2.q.beta[j]
      nu.q.beta.out[z] <- nu.q.beta[j]
      elbo.out.beta <- c(elbo.out.beta, ELBO.beta.j)
        
      #stopping criterion
      if(z > 1) {
        Delta.beta <- (abs((elbo.out.beta[z] - 
                      elbo.out.beta[z-1])/
                      elbo.out.beta[z-1]))
        if(Delta.beta < eps_beta) cond.beta <- 1
        if(Trace_beta == 1) {
          print(paste0("iteration:", z, 
              "- beta lower bound increase: ", Delta.beta))
        }
      }
      if (z > maxIter_beta) cond.beta <- 1
      z <- z + 1
      }
    
  return(list(mu.q.beta.j = mu.q.beta.out[z-1], 
              sigma2.q.beta.j = sigma2.q.beta.out[z-1],
              nu.q.beta.j = nu.q.beta.out[z-1]))
}

#gamma.j-localized component of lower bound
update.gamma.j <- function(mu.q.gamma.j, y, x, alpha.rho, 
                           delta.rho, mu.q.beta, sigma2.q.beta, 
                           elementi.prod, mu.q.log.rho, 
                           mu.q.log.1.rho, mu.q.gamma, j) {
  p <- ncol(x)
  n <- nrow (x)
  
  
  entropy.q.gamma <- - (mu.q.gamma.j * log(mu.q.gamma.j) + 
                     (1-mu.q.gamma.j)*log(1-mu.q.gamma.j))
  
  if(p == 1) { 
    logL <- as.vector(t(x * mu.q.gamma.j * mu.q.beta) %*% y - 
            sum(lfactorial(y)) -sum(1 - mu.q.gamma.j + 
            mu.q.gamma.j*exp(x*mu.q.beta + sigma2.q.beta*x^2/2)))
  } 
  else if(p == 2){
    produttoria.no.j <- produttoria.q.beta.q.gamma(elementi.prod, 
                        j)
    logL <- (t(x[,-j] * mu.q.gamma[-j] * mu.q.beta[-j]) %*% y) +
            (t(x[,j] * mu.q.gamma.j * mu.q.beta[j]) %*% y) - 
            sum(lfactorial(y)) - sum(produttoria.no.j * 
            (1 - mu.q.gamma.j + mu.q.gamma.j*exp(x[,j]*
            mu.q.beta[j] + sigma2.q.beta[j]*(x[,j]^2)/2)))
  } 
  else{
    produttoria.no.j <- produttoria.q.beta.q.gamma(elementi.prod,
                        j)
    logL <- (t(x[,-j] %*% diag(mu.q.gamma[-j]) %*% mu.q.beta[-j]) 
             %*% y) + (t(x[,j] * mu.q.gamma.j * mu.q.beta[j])%*%
             y) - sum(lfactorial(y)) - sum(produttoria.no.j * 
            (1 - mu.q.gamma.j + mu.q.gamma.j*exp(x[,j]*
            mu.q.beta[j] + sigma2.q.beta[j]*(x[,j]^2)/2)))
  }
  
  lprior.gamma <- mu.q.gamma.j*mu.q.log.rho +(1-mu.q.gamma.j)*
                  mu.q.log.1.rho
  
  ELBO.gamma <- entropy.q.gamma + lprior.gamma + logL
  return(ELBO.gamma)
}


#first derivative of gamma.j-localized component of lower bound
der.prime.gamma.j <- function(mu.q.gamma.j, y, x, mu.q.beta, 
                              sigma2.q.beta, mu.q.gamma, 
                              elementi.prod, mu.q.log.rho, 
                              mu.q.log.1.rho, j) {
  p <- ncol(x)
  
  if(p == 1)
    produttoria.no.j <- rep(1, n)
  else
    produttoria.no.j <- produttoria.q.beta.q.gamma(elementi.prod, 
                        j)
  output <- mu.q.log.rho - mu.q.log.1.rho - log(mu.q.gamma.j) - 
            1 + 1/(1-mu.q.gamma.j) + log(1 - mu.q.gamma.j)- 
            mu.q.gamma.j/(1 - mu.q.gamma.j) + t(x[,j] * 
            mu.q.beta[j]) %*% y - sum(produttoria.no.j * 
            (exp(x[,j]*mu.q.beta[j] + sigma2.q.beta[j]*(x[,j]^2)
            /2) -1))
  return(output)
}

#MFVB algorithm
MFVB.norm.bern.beta <- function(y, x, hyp, start, eps_beta, 
                                eps_gamma, eps_global, 
                                maxIter_beta, 
                                maxIter_gamma, maxIter_global, 
                                Trace_beta = 0,Trace_global = 0){
  start.time <- Sys.time()
  
  #============= Hyperparameters' setting =============
  mu.beta <- hyp[[1]][[1]]
  sigma2.beta <- hyp[[1]][[2]]
  alpha.rho <- hyp[[2]][[1]]
  delta.rho <- hyp[[2]][[2]]
  n <- nrow(x)
  p <- ncol(x)
  
  #============= Defining history variables =============
  elbo.out.global <- numeric()
  mu.q.beta.global.out <- matrix(NA, nrow = maxIter_global, 
                          ncol = p)
  sigma2.q.beta.global.out <- matrix(NA, ncol = p,
                              nrow =  maxIter_global)
  mu.q.gamma.out <- matrix(NA, nrow = maxIter_global, ncol = p)
  alpha.q.rho.out <- numeric()
  delta.q.rho.out <- numeric()
  nu.q.beta.global.out <- matrix(NA, nrow = maxIter_global, 
                                 ncol = p)
  
  #============= Initializing optimal parameters =============
  mu.q.beta <- start[[1]]
  sigma2.q.beta <- start[[2]]
  omega.q.beta <- 1/sigma2.q.beta
  mu.q.gamma <- start[[3]]
  alpha.q.rho <- start[[4]]
  delta.q.rho <- start[[5]]
  mu.q.log.rho <- digamma(alpha.q.rho) + digamma(alpha.q.rho + 
                  delta.q.rho)
  mu.q.log.1.rho <- digamma(delta.q.rho) + digamma(alpha.q.rho + 
                    delta.q.rho)
  
  nu.q.beta <- mu.q.beta
  
  elementi.prod <- prod.elementi.start(x, mu.q.gamma, mu.q.beta, 
                                       sigma2.q.beta, n, p)
  
  #============= MFVB algorithm =============
  i <- 1
  cond.global <- FALSE
  
  while(cond.global != TRUE & i <= maxIter_global) { 
    
    #updating q.beta
    for(j in 1:p){
      q.beta.j.update <- update.q.beta.j(y, x, sigma2.beta, 
                         mu.q.beta, sigma2.q.beta, mu.q.gamma,
                         elementi.prod, maxIter_beta[j], 
                         eps_beta[j], Trace_beta[j], j,
                         nu.q.beta)
      mu.q.beta[j] <- q.beta.j.update[[1]]
      sigma2.q.beta[j] <- q.beta.j.update[[2]]
      nu.q.beta[j] <- q.beta.j.update[[3]]
      }
    
    #updating q.gamma
    for(j in 1:p){
      mu.q.gamma[j] <- optim(mu.q.gamma[j],
                             function(z) -update.gamma.j(z, y, x,
                                          alpha.rho, delta.rho,
                                          mu.q.beta,
                                          sigma2.q.beta,
                                          elementi.prod,
                                          mu.q.log.rho,
                                          mu.q.log.1.rho,
                                          mu.q.gamma, j),
                             function(z) -der.prime.gamma.j(z, y,
                                          x, mu.q.beta,
                                          sigma2.q.beta,
                                          mu.q.gamma,
                                          elementi.prod,
                                          mu.q.log.rho,
                                          mu.q.log.1.rho, j),
                             method = "L-BFGS-B",
                             lower = 0.01, upper = 0.99)$par
      elementi.prod[,j] <- change.prod.elementi.j(x, mu.q.gamma,
                           mu.q.beta, sigma2.q.beta, j, p)
      }
    
    #updating q.rho
    alpha.q.rho <- alpha.rho + sum(mu.q.gamma)
    delta.q.rho <- delta.rho + p - sum(mu.q.gamma)
    mu.q.log.rho <- digamma(alpha.q.rho) + digamma(alpha.q.rho + 
                    delta.q.rho)
    mu.q.log.1.rho <- digamma(delta.q.rho) +digamma(alpha.q.rho+ 
                      delta.q.rho)
        
    #computing lower bound
    
    entropy.q.beta <- p/2*log(2*pi)+0.5*sum(log(sigma2.q.beta)) +
                      p/2
    
    entropy.q.gamma <- - sum(mu.q.gamma * log(mu.q.gamma) +
                      (1-mu.q.gamma) * log(1-mu.q.gamma))
    
    entropy.q.rho <- -(alpha.q.rho - 1)* mu.q.log.rho -
                     (delta.q.rho - 1)*mu.q.log.1.rho +
                     lbeta(alpha.q.rho, delta.q.rho)
    
    lprior.beta <- -p/2*log(2*pi) - p/2*log(sigma2.beta) - 
                   1/(2*sigma2.beta)*sum(sigma2.q.beta+ 
                   mu.q.beta^2)
    
    lprior.gamma <- sum(mu.q.gamma*mu.q.log.rho +(1-mu.q.gamma)*
                    mu.q.log.1.rho)
    
    if(p == 1){
      logL <- as.vector(t(x * mu.q.gamma * mu.q.beta) %*% y - 
              sum(lfactorial(y)) - sum(elementi.prod))
    }
    else{
      produttoria.total <- apply(elementi.prod, 1, function(x) 
                           prod(x))
      logL <- as.vector(t(x %*%diag(mu.q.gamma)%*%mu.q.beta) %*% 
              y) - sum(lfactorial(y)) - sum(produttoria.total)
    }

    lprior.rho <- (alpha.rho - 1)*mu.q.log.rho + (delta.rho - 1)*
                  mu.q.log.1.rho - lbeta(alpha.rho, delta.rho)
    
    ELBO.global <- entropy.q.beta + entropy.q.gamma + 
                   entropy.q.rho + lprior.beta + 
                   lprior.gamma + logL + lprior.rho 
    
    #updating history variables
    mu.q.beta.global.out[i,] <- mu.q.beta
    sigma2.q.beta.global.out[i,] <- sigma2.q.beta
    mu.q.gamma.out[i,] <- mu.q.gamma
    alpha.q.rho.out <- c(alpha.q.rho.out, alpha.q.rho)
    delta.q.rho.out <- c(delta.q.rho.out, delta.q.rho)
    elbo.out.global <- c(elbo.out.global, ELBO.global)
    
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
              mu.q.beta = mu.q.beta.global.out[1:(i-1),], 
              mu.q.gamma = mu.q.gamma.out[1:(i-1),],
              sigma2.q.beta=sigma2.q.beta.global.out[1:(i-1),], 
              alpha.q.rho = alpha.q.rho.out, 
              delta.q.rho = delta.q.rho.out, iter = i-1, 
              tempo = total.time[[1]]))
}
