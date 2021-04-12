updatePhi_r_marginal <- function(x, Cs, Dist, Mu, Tau, Xi, Alpha, Beta,
                                 phi.old = 1, maxit = 200, eps = 1e-4,
                                 hessian = FALSE) {

  n_k <- as.vector(table(Cs))
  K <- length(n_k)
  routine.phi <- optim(par = phi.old, fn = function(phi){
    if(phi <= 0) return(-Inf)
    eig <- eigen(exp(-Dist/phi))
    val <- 0
    for(k in 1:K){
      block1 <- (x[Cs == k,] - Mu[k]) %*% eig$vec
      alpha.i <- ncol(x)/2 + Alpha[k]
      beta.i <- diag(block1 %*% as(diag(1/(Tau[k]*eig$val + Xi[k])), "sparseMatrix") %*% t(block1))/2 + Beta[k]
      val <- val - n_k[k]/2*sum(log(Tau[k]*eig$val + Xi[k])) - sum(alpha.i * log(beta.i))
    }
    -val
  }, hessian = hessian)
  if(routine.phi$conv != 0) stop("Converge error in Phi!")
  if(!hessian) return(routine.phi$par) else
    return(routine.phi)
}

