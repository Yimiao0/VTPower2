calcCondU <- function(Sigma.y, mean.t, Sigma.t, p) {
  J <- length(mean.t)
  U <- matrix(0, nrow = J+3, ncol = J+3)
  M.mean <- cbind(mean.t, p, p*mean.t)
  
  Sigma.marg <- Sigma.y[2:(J+1),2:(J+1)]
  Sigma.cond <- Sigma.marg - Sigma.y[2:(J+1),1]%*%t(Sigma.y[1,2:(J+1)]) / Sigma.y[1,1]
  inv.Sigma.cond <- solve(Sigma.cond)
  
  U[1:J, 1:J] <- inv.Sigma.cond
  U[1:J, (J+1):(J+3)] <- inv.Sigma.cond%*%M.mean
  U[(J+1):(J+3), 1:J] <- t(U[1:J, (J+1):(J+3)])
  
  middle.mat <- matrix(NA, nrow = 3, ncol = 3)
  middle.mat[1,1] <- t(mean.t)%*%inv.Sigma.cond%*%mean.t + sum(diag(inv.Sigma.cond%*%Sigma.t))
  middle.mat[1,3] <- middle.mat[3,1] <- middle.mat[3,3] <- p*middle.mat[1,1]
  middle.mat[1,2] <- middle.mat[2,1] <- middle.mat[2,3] <- middle.mat[3,2] <- p*t(rep(1,J))%*%inv.Sigma.cond%*%mean.t
  middle.mat[2,2] <- p*sum(inv.Sigma.cond)
  U[(J+1):(J+3), (J+1):(J+3)] <- middle.mat

  return(U)
}
