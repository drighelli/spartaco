main <- function(x, Dist,
                 K,
                 R,
                 traceRatio = 10,
                 max.iter = 10^3,
                 metropolis.iterations = 150,
                 estimate.iterations = 10,
                 input.values = NULL,
                 verbose = F){
    set.seed(seed = NULL)
    if(is.null(input.values)){
        cur.Cs <- best.Ds <- sample(1:K, size = nrow(x), replace = T)
        cur.Ds <- best.Cs <- sample(1:R, size = ncol(x), replace = T)
        cur.phi <- best.phi <- runif(R, 1, 5)
        cur.mu <- best.mu <- matrix(runif(K*R, 1, 10), K, R)
        cur.tau <- best.tau <- matrix(runif(K * R, 1e-7, traceRatio), K, R)
        cur.xi <- best.xi <- traceRatio - best.tau
        cur.alpha <- best.alpha <- matrix(runif(K*R, 1, 3), K, R)
        cur.beta <- best.beta <- matrix(runif(K*R, 1, 3), K, R)} else {
            cur.Cs <- best.Cs <- input.values$best.Cs
            cur.Ds <- best.Ds <- input.values$best.Ds
            cur.phi <- best.phi <- input.values$phi
            cur.mu <- best.mu <- input.values$mu
            cur.tau <- best.tau <- input.values$tau
            cur.xi <- best.xi <- traceRatio - best.tau
            cur.alpha <- best.alpha <- input.values$alpha
            cur.beta <- best.beta <- input.values$beta}
    Cs <- matrix(0, nrow(x), max.iter)
    Cs[,1] <- cur.Cs
    Ds <- matrix(0, ncol(x), max.iter)
    Ds[,1] <- cur.Ds

    Uglob <- list()
    Dglob <- numeric(ncol(x))
    for(r in 1:R){
        eigK <- eigen(exp(-Dist[cur.Ds == r, cur.Ds == r]/cur.phi[r]))
        Uglob[[r]] <- eigK$vec
        Dglob[cur.Ds == r] <- eigK$val
    }

    ll <- -1e+40
    i <- 1
    while(T){
        if(i == max.iter) break
        i <- i+1
        if(verbose) cat(paste("---Iteration",i,"\n"))
        goodK <- sort(unique(cur.Cs))
        goodR <- sort(unique(cur.Ds))
        for(r in goodR){
            traceDelta_r <- traceRatio * sum(cur.Ds == r)
            for(k in goodK){
                if(verbose) cat(paste("r = ",r,", k = ",k,"/",sep=""))
                estimation.parameters <- Estimate.Cocluster.Parameters.marginal.constraint.trace(x = x[cur.Cs == k, cur.Ds == r],
                                                                                                 traceDelta = traceDelta_r,
                                                                                                 U = Uglob[[r]],
                                                                                                 d = Dglob[cur.Ds == r],
                                                                                                 mu0 = cur.mu[k,r],
                                                                                                 alpha0 = cur.alpha[k,r],
                                                                                                 beta0 = cur.beta[k,r],
                                                                                                 tau0 = cur.tau[k,r],
                                                                                                 maxit = estimate.iterations)
                cur.mu[k,r] <- estimation.parameters$mu
                cur.tau[k,r] <- estimation.parameters$tau
                cur.xi[k,r] <- estimation.parameters$xi
                cur.alpha[k,r] <- estimation.parameters$alpha
                cur.beta[k,r] <- estimation.parameters$beta
            }
            cur.phi[r] <- updatePhi_r_marginal(x = x[,cur.Ds == r],
                                               Cs = cur.Cs,
                                               Dist = Dist[cur.Ds == r, cur.Ds == r],
                                               Mu = cur.mu[,r],
                                               Tau = cur.tau[,r],
                                               Xi = cur.xi[,r],
                                               Alpha = cur.alpha[,r],
                                               Beta = cur.beta[,r],
                                               phi.old = cur.phi[r])
            EigenK <- eigen(exp(-Dist[cur.Ds == r, cur.Ds == r]/cur.phi[r]))
            Uglob[[r]] <- EigenK$vec
            Dglob[cur.Ds == r] <- EigenK$val
        }
        if(verbose) cat("\n")
        if(i %% ifelse(K > 1, 2, 1) == 0){
            cur.ds <- MetropolisAllocation(x = x, Uglob = Uglob, Dglob = Dglob,
                                           Cs = cur.Cs, Ds = cur.Ds, Dist = Dist, Mu = cur.mu, Tau = cur.tau, Xi = cur.xi, Alpha = cur.alpha, Beta = cur.beta, Phi = cur.phi,
                                           maxit = metropolis.iterations, min.obs = 5)
            if(verbose) cat(paste("Changed",cur.ds$accepted ,"elements\n"))
            cur.Ds <- cur.ds$Ds
            Uglob <- cur.ds$Uglob
            Dglob <- cur.ds$Dglob
            logL.values <- cur.ds$logL.values
            ll[i] <- cur.ds$logL} else {
                cur.cs <- RowClustering(x = x, Ds = cur.Ds, Mu = cur.mu, Tau = cur.tau, Xi = cur.xi, Alpha = cur.alpha, Beta = cur.beta, Phi = cur.phi, Uglob = Uglob, Dglob = Dglob)
                cur.Cs <- cur.cs$allocation
                goodK <- sort(unique(cur.Cs))
                goodR <- sort(unique(cur.Ds))
                ll[i] <- 0
                for(r in goodR)
                    for(k in goodK){
                        logL.values[k,r] <- logL.Cocluster(x = x[cur.Cs == k, cur.Ds == r],
                                                           Mu = cur.mu[k,r],
                                                           Tau = cur.tau[k,r],
                                                           Xi = cur.xi[k,r],
                                                           Alpha = cur.alpha[k,r],
                                                           Beta = cur.beta[k,r],
                                                           U = Uglob[[r]],
                                                           d = Dglob[cur.Ds == r])
                        ll[i] <- ll[i] + logL.values[k,r]
                    }
            }
        Ds[,i] <- cur.Ds
        Cs[,i] <- cur.Cs

        if(ll[i] == max(ll)){
            best.phi <- cur.phi
            best.mu <- cur.mu
            best.tau <- cur.tau
            best.xi <- cur.xi
            best.alpha <- cur.alpha
            best.beta <- cur.beta
            best.Cs <- cur.Cs
            best.Ds <- cur.Ds
        }

        if(verbose){
            cat(paste("diff(loglikelihood) =",round(diff(ll)[i-1],5),"\n"))
            cat(paste("Row cluster size =", paste(table(cur.Cs), collapse = ", "),"\n"))
            cat(paste("Column cluster size =", paste(table(cur.Ds), collapse = ", "),"\n"))
        }
    }
    return(list(
        phi = best.phi,
        mu = best.mu,
        tau = best.tau,
        xi = best.xi,
        alpha = best.alpha,
        beta = best.beta,
        Cs = best.Cs,
        Ds = best.Ds,
        CS = Cs,
        DS = Ds,
        logL = ll,
        x = x
    ))
}
