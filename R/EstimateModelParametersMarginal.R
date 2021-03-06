Estimate.Cocluster.Parameters.marginal.constraint.trace <- function(x,
                                                                    U,
                                                                    d,
                                                                    mu0,
                                                                    alpha0,
                                                                    beta0,
                                                                    tau0,
                                                                    traceDelta = 5000,
                                                                    maxit = 200,
                                                                    hessian = F,
                                                                    threshold = 1e-4,
                                                                    print.conv.warning = F){
    n <- nrow(x)
    p <- ncol(x)
    Mu <- cur.mu <- mu0
    Tau <- cur.tau <- tau0
    cur.xi <- traceDelta/p - tau0
    Alpha <- cur.alpha <- alpha0
    Beta <- cur.beta <- beta0
    converged <- F
    bl1 <- x %*% U
    bl2 <- matrix(1, n, p) %*% U
    cur.xi <- traceDelta/p - cur.tau
    for(i in 2:maxit){
        # --update mu
        routine.mu <- optim(par = cur.mu, fn = function(mu){
            sum(
                log(
                    diag((bl1 - mu * bl2) %*% diag(1/(cur.tau * d + cur.xi)) %*% t(bl1 - mu * bl2))/2 + cur.beta
                )
            )*(p/2 + cur.alpha)
        }, hessian = hessian)
        if(routine.mu$convergence != 0){
            stop("Convergence error in mu!")
        }
        cur.mu <- routine.mu$par

        # --update alpha
        quadratic <- diag((bl1 - cur.mu * bl2) %*% diag(1/(cur.tau * d + cur.xi)) %*% t(bl1 - cur.mu * bl2))/2
        routine.alpha <- optim(cur.alpha, function(a){
            if(a <= 0) return(-Inf)
            -(
                a*(n*log(cur.beta) - sum(log(quadratic + cur.beta)))-n*(lgamma(a)-lgamma(p/2+a))
            )}, hessian = hessian)
        if(routine.alpha$convergence != 0){
            stop("Convergence error in alpha!")
        }
        cur.alpha <- routine.alpha$par

        # --update beta
        routine.beta <- optim(cur.beta, function(b){
            if(b <= 0) return(-Inf)
            -(n*cur.alpha*log(b)-(p/2+cur.alpha)*sum(log(quadratic+b)))}, hessian = hessian)
        if(routine.beta$convergence != 0){
            stop("Convergence error in beta!")
        }
        cur.beta <- routine.beta$par

        # --update tau and xi
        Block1 <- bl1 - cur.mu * bl2
        starting.tau <- ifelse(cur.tau < traceDelta/p, cur.tau, runif(1, 1e-7, traceDelta/p))
        routine.tau <- optim(starting.tau,
                             fn = function(taup){
                                 if(taup <= 0) return(-Inf)
                                 xip <- traceDelta/p - taup
                                 if(xip <= 0) return(-Inf)
                                 -(
                                     -n/2*sum(log(taup * d + xip)) -
                                         (p/2+cur.alpha) * sum(log(diag(Block1 %*% diag(1/(taup * d + xip)) %*% t(Block1))/2 + cur.beta))
                                 )
                             }, hessian = hessian, control = list(maxit = 1000))
        if(routine.tau$convergence != 0){
            stop("Convergence error in tau!")
        }
        cur.tau <- routine.tau$par
        cur.xi <- traceDelta/p - cur.tau

        Mu[i] <- cur.mu
        Alpha[i] <- cur.alpha
        Beta[i] <- cur.beta
        Tau[i] <- cur.tau
        if(abs(diff(Mu)[i-1]) < threshold &
           abs(diff(Alpha)[i-1]) < threshold &
           abs(diff(Beta)[i-1]) < threshold &
           abs(diff(Tau)[i-1]) < threshold){
            converged <- T
            break}
    }
    if(converged == F & print.conv.warning) cat("Warning: max-it reached without convergence!")
    if(hessian == F) return(list(mu = Mu[i],
                                 alpha = Alpha[i],
                                 beta = Beta[i],
                                 tau = Tau[i],
                                 xi = cur.xi,
                                 algorithm = list(maxit = i, parameters.at.iterations = list(Mu, Alpha, Beta, Tau))))
    else return(list(mu = routine.mu,
                     alpha = routine.alpha,
                     beta = routine.beta,
                     tau = routine.tau))
}

