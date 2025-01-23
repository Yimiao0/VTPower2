calcMargQ <- function(Sigma.marg, mean.t, Sigma.t) {
  J <- length(mean.t)
  inv.Sigma.marg <- solve(Sigma.marg)
  
  Q <- matrix(NA, nrow = 2, ncol = 2)
  Q[1,1] <- sum(inv.Sigma.marg)
  Q[1,2] <- Q[2,1] <- t(rep(1,J))%*%inv.Sigma.marg%*%mean.t
  Q[2,2] <- t(mean.t)%*%inv.Sigma.marg%*%mean.t + sum(diag(inv.Sigma.marg%*%Sigma.t))
  
  return(Q)
}
