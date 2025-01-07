savs.algorithm <- function(mu.q.beta, X){
  p <- ncol(X)
  mu.q.beta.star <- numeric(p)
  for(j in 1:p){
    mu.j <- 1/mu.q.beta[j]^2
    sum.x.j.2 <- sum(X[,j]^2)
    if(abs(mu.q.beta[j])*sum.x.j.2<= mu.j)
      mu.q.beta.star[j] <- 0
    else
      mu.q.beta.star[j] <- sign(mu.q.beta[j])/sum.x.j.2*
      (abs(mu.q.beta[j])*sum.x.j.2 - mu.j)
  }
  return(mu.q.beta.star)
}
