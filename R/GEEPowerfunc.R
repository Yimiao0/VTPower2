##' Analytical solution of required sample size for longitudinal slope with varying follow-up time points.
##' 
##' @description
##' The GEEPowerAnaly function is to estimate the required sample size for comparison of longitudinal slopes between two groups.
##' Two approaches are implemented in this function: (1) a marginal approach suited for designs without baseline outcome measurement; and (2) a conditional approach adjusting for baseline outcome.
##' 
##' \loadmathjax
##' 
##' Assuming that the treatment assignment \mjeqn{A_i}{} and measurement time points \mjeqn{T_i}{} are independent for all subjects \mjeqn{i}{}, the sample size formula for marginal approach is given by:
##' \mjdeqn{N_M = \frac{(z_{\alpha/2} + z_b)^2}{p(1-p) \delta^2} e_2^\top Q^{-1} e_2}{}
##' where
##' - \mjeqn{\alpha}{} is the type I error.
##' - \mjeqn{b}{} is the type II error, i.e., \mjeqn{b = 1 - power}{}.
##' - \mjeqn{z_{\alpha/2}}{} is the \mjeqn{100(1-\alpha/2)}{}-th percentile of the standard normal distribution, the same applied to \mjeqn{z_b}{}.
##' - \mjeqn{\delta}{} is the group difference in longitudinal slope.
##' - \mjeqn{p}{} is the probability that a subject is assigned to the treatment group.
##' - \mjeqn{e_2}{} is a standard basis vector \mjeqn{(0,1)^\top}{}.
##' - \mjeqn{Q}{} is a \mjeqn{2 \times 2}{} matrix which can be calculate given the covariance matrix of outcomes and the expectation and covariance of the time points. For more information, please refer to the reference.
##' 
##' For conditional approach, the sample size formula is given by:
##' \mjdeqn{N_C = \frac{(z_{\alpha/2} + z_b)^2}{\delta^2} e_{J+3}^\top U^{-1} e_{J+3}}{},
##' where \mjeqn{U^{-1}}{} is a \mjeqn{(J+3) \times (J+3)}{} matrix. For more information of \mjeqn{U}{}, see the reference.
##' 
##' 
##' 
##' @param delta treatment effect on slope
##' @param p probability that a subject is assigned to the treatment group
##' @param alpha type I error
##' @param power the targeted power
##' @param MARGINonly takes logical value `TRUE` or `FALSE`. If `TRUE`, the function will return only the sample size estimated by marginal approach. Otherwise, the function will return both sample sizes estimation by marginal and conditional approach.
##' @param Sigma.y the covariance for post-baseline outcomes excluding baseline outcome if `MARGINonly` is `TRUE`; the covariance for all outcomes including both baseline and post-baseline if `MARGINonly` is `FALSE`
##' @param mean.t the expectation vector of random time points
##' @param Sigma.t the covariance matrix of random time points
##' 
##' 
##' 
##' @return sample size estimated by marginal approach if `MARGINonly` is `TRUE`; sample sizes estimated by marginal and conditional approaches, respectively, if `MARGINonly` is `FALSE`.
##' 
##' 
##' 
##' @references Under review
##' 
##' 
##' 
##' @examples 
##' \dontrun{
##' library(GEEPower)
##' 
##' # For designs with baseline outcome
##' delta <- -0.1
##' 
##' mean.t = c(1,2,3,6)
##' Sigma.t = diag(c(1.5-0.5, 2.5-1.5, 3.5-2.5, 7-5)^2/12)
##' J = length(mean.t)
##' 
##' sigma.y <- 2
##' rho <- 0.5
##' Corr.y <- matrix(NA, nrow = J+1, ncol = J+1)
##' for (i in 1:(J+1)) {
##'   for (j in 1:(J+1)) {
##'     Corr.y[i,j] <- rho^abs(i-j)
##'   }
##' }
##' Sigma.y = sigma.y^2*Corr.y
##' 
##' GEEPowerAnaly(delta, Sigma.y, mean.t, Sigma.t)
##' 
##' 
##' # For designs without baseline outcome
##' delta <- -0.1
##' 
##' mean.t = c(1,2,3,6)
##' Sigma.t = diag(c(1.5-0.5, 2.5-1.5, 3.5-2.5, 7-5)^2/12)
##' J = length(mean.t)
##' 
##' sigma.y <- 2
##' rho <- 0.5
##' Corr.y <- matrix(NA, nrow = J, ncol = J)
##' for (i in 1:J) {
##'   for (j in 1:J) {
##'     Corr.y[i,j] <- rho^abs(i-j)
##'   }
##' }
##' Sigma.y = sigma.y^2*Corr.y
##' 
##' GEEPowerAnaly(delta, Sigma.y, mean.t, Sigma.t, MARGINonly = TRUE)
##' 
##' }
##' 
##' 
##' 
##' 
##' @import mathjaxr
##' @importFrom MASS ginv
##' @export
##' @md





GEEPowerAnaly <- function(delta,
                          Sigma.y,
                          mean.t,
                          Sigma.t,
                          p = 0.5,
                          alpha = 0.05,
                          power = 0.8,
                          MARGINonly = FALSE) {
  
  z.alpha <- qnorm(p = 1 - alpha/2)
  z.b <- qnorm(p = power)
  
  if (MARGINonly) {
    if (nrow(Sigma.y) != length(mean.t)) stop("Sigma.y must be covariance for post-baseline outcomes when MARGINonly is TRUE! Check whether your nrow(Sigma.y) == length(mean.t)!")
    
    Q <- calcMargQ(Sigma.marg = Sigma.y, mean.t, Sigma.t)
    N.marg <- (z.alpha + z.b)^2 * ginv(Q, tol = 1e-30)[2,2] / (p*(1-p)*delta^2)
    names(N.marg) <- "Marginal"
    
    return(N.marg)
    
  } else {
    if (nrow(Sigma.y) != length(mean.t)+1) stop("Sigma.y must be covariance for all outcomes including baseline when MARGINonly is FALSE! Check whether your nrow(Sigma.y) == 1 + length(mean.t)!")
    
    J <- length(mean.t)
    Sigma.marg <- Sigma.y[2:(J+1), 2:(J+1)]
    Q.marg <- calcMargQ(Sigma.marg, mean.t, Sigma.t)
    N.marg <- (z.alpha + z.b)^2 * ginv(Q.marg, tol = 1e-30)[2,2] / (p*(1-p)*delta^2)
    
    U.cond <- calcCondU(Sigma.y, mean.t, Sigma.t, p)
    N.cond <- (z.alpha + z.b)^2 * ginv(U.cond, tol = 1e-30)[J+3,J+3] / delta^2
    
    results <- c(N.marg, N.cond)
    names(results) <- c("Marginal", "Conditional")
    
    return(results)
  }
  
}
