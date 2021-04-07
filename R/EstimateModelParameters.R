library(Matrix)
library(invgamma)

Estimate.Cocluster.Parameters.constraint.trace <- function(x, 
                                          U,
                                          d,
                                          traceDelta = 1000,
                                          mu0 = NULL,
                                          alpha0 = runif(1,1,10),
                                          beta0 = runif(1,1,10),
                                          tau0 = runif(1,1,10),
                                          maxit = 200,
                                          threshold = 1e-4,
                                          print.conv.warning = F){
  n <- nrow(x)
  p <- ncol(x)
  Mu <- Tau <- Xi <- Alpha <- Beta <- numeric(maxit)
  Mu[1] <- cur.mu <- ifelse(is.null(mu0), mean(x), mu0)
  Tau[1] <- cur.tau <- tau0
  Xi[1] <- cur.xi <- traceDelta/p - cur.tau
  Alpha[1] <- cur.alpha <- alpha0
  Beta[1] <- cur.beta <- beta0
  converged <- F
  bl1 <- x %*% U
  bl2 <- matrix(1, n, p) %*% U
  Block1 <- bl1 - cur.mu * bl2
  for(i in 2:maxit){
    #print(i)
    invD <- 1/(cur.tau * d + cur.xi)
    alpha.post <- p/2 + cur.alpha
    beta.post <- diag(Block1 %*% diag(invD) %*% t(Block1))/2 + cur.beta
    expected <- matrix(0, n, 2)
    expected[,1] <- exp(lgamma(alpha.post+1)-log(beta.post)-lgamma(alpha.post))
    expected[,2] <- log(beta.post)-digamma(alpha.post)
    
    # --update mu
    num.mu <- sum(expected[,1] * diag(bl1 %*% diag(invD) %*% t(bl2)))
    den.mu <- sum(expected[,1] * diag(bl2 %*% diag(invD) %*% t(bl2)))
    cur.mu <- num.mu / den.mu
    Block1 <- bl1 - cur.mu * bl2
    
    # --update alpha
    #print(sum(is.nan(expected[,2])))
    routine.alpha <- optim(cur.alpha, function(a){
      if(a <= 0) return(-Inf)
      -(-a*sum(expected[,2]) + (n)*(a*log(cur.beta) - lgamma(a)))})
    if(routine.alpha$convergence != 0){
      stop("Convergence error in alpha!")
    }
    cur.alpha <- routine.alpha$par
    
    # --update beta
    cur.beta <- (n) * cur.alpha / sum(expected[,1])

    # --update tau and xi
    starting.tau <- ifelse(cur.tau < traceDelta/p, cur.tau, runif(1, 1e-7, traceDelta/p))
    routine.tau.xi <- optim(starting.tau, 
                            fn = function(val){
                              if(val <= 0 | val >= traceDelta/p) return(-Inf)
                              xi <- traceDelta/p-val
                              invD.prop <- 1/(val*d+xi)
                              -(
                                n/2*sum(log(invD.prop)) - .5 * sum(expected[,1] * diag(Block1 %*% diag(invD.prop) %*% t(Block1)))
                              )
                            })
    if(routine.tau.xi$convergence != 0){
      stop("Convergence error in tau!")
    }
    cur.tau <- routine.tau.xi$par
    cur.xi <- traceDelta/p-cur.tau
    
    Mu[i] <- cur.mu
    Alpha[i] <- cur.alpha
    Beta[i] <- cur.beta
    Tau[i] <- cur.tau
    Xi[i] <- cur.xi
    if(abs(diff(Mu)[i-1]) < threshold & 
       abs(diff(Alpha)[i-1]) < threshold &
       abs(diff(Beta)[i-1]) < threshold &
       abs(diff(Tau)[i-1]) < threshold &
       abs(diff(Xi)[i-1]) < threshold){
      converged <- T
      break}
  }
  if(converged == F & print.conv.warning) cat("Warning: max-it reached without convergence!")
  return(list(mu = Mu[i],
              alpha = Alpha[i],
              beta = Beta[i],
              tau = Tau[i],
              xi = Xi[i],
              algorithm = list(maxit = i, parameters.at.iterations = list(Mu[1:i], Alpha[1:i], Beta[1:i], Tau[1:i], Xi[1:i])))
  )
}




