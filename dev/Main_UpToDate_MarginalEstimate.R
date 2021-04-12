rm(list=ls())
library(Matrix)
library(invgamma)

functions.directory <- "/Users/andreasottosanti/Dropbox/Post-Doc/Code/SpatialBicluster/NewModel/Functions_ExponentialKernel/"
source(paste(functions.directory,"EstimateModelParametersMarginal.R",sep=""))
source(paste(functions.directory,"EstimatePhi.R",sep=""))
source(paste(functions.directory,"loglikelihood.R",sep=""))
source(paste(functions.directory,"RowAllocation.R",sep=""))
source(paste(functions.directory,"MetropolisAllocation.R",sep=""))


# LOAD THE DATA HERE
#load("/Users/andreasottosanti/Dropbox/Post-Doc/Code/SpatialBicluster/RealData applications/spatialLIBD/Sample 151673/Data/Sample1_V2.Rdata")
x <- as.matrix(sce@assays@data$poisson_deviance_residuals)
coordinates <- cbind(sce$imagerow, sce$imagecol)

Dist <- as.matrix(stats::dist(coordinates))

K <- 3
R <- 3

save.results <- T
save.directory <- "/Users/andreasottosanti/Dropbox/Post-Doc/Code/SpatialBicluster/RealData applications/spatialLIBD/Sample 151673/Results/"
save.to <- paste(save.directory,"Sample1_results_K",
                 K,"_R",R,"_",
                 paste(substr(Sys.time(),1,10),"_",substr(Sys.time(),15,16),"-",substr(Sys.time(),18,19),sep=""),".Rdata",sep="")

cur.Cs <- best.Ds <- sample(1:K, size = nrow(x), replace = T)
cur.Ds <- best.Cs <- sample(1:R, size = ncol(x), replace = T)#kmeans(t(x), centers = R)$clust
traceRatio <- 10
cur.phi <- best.phi <- runif(R, 1, 5)
cur.mu <- best.mu <- matrix(runif(K*R, 1, 10), K, R)
cur.tau <- best.tau <- matrix(runif(K * R, 1e-7, traceRatio), K, R)
cur.xi <- best.xi <- traceRatio - best.tau
cur.alpha <- best.alpha <- matrix(runif(K*R, 1, 3), K, R)
cur.beta <- best.beta <- matrix(runif(K*R, 1, 3), K, R)
max.iter <- 1000
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

metropolis.iterations <- 150
estimate.iterations <- 10
ll <- -1e+40
i <- 1
while(T){
  if(i == max.iter) break
  i <- i+1
  cat(paste("---Iteration",i,"\n"))
  for(r in 1:R){
    traceDelta_r <- traceRatio * sum(cur.Ds == r)
    for(k in 1:K){
      cat(paste("r = ",r,", k = ",k,"/",sep=""))
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
  cat("\n")
  if(i %% ifelse(K > 1, 2, 1) == 0){
    cur.ds <- MetropolisAllocation(x = x, Uglob = Uglob, Dglob = Dglob,
                                   Cs = cur.Cs, Ds = cur.Ds, Dist = Dist, Mu = cur.mu, Tau = cur.tau, Xi = cur.xi, Alpha = cur.alpha, Beta = cur.beta, Phi = cur.phi, 
                                   maxit = metropolis.iterations, min.obs = 5)
    cat(paste("Changed",cur.ds$accepted ,"elements\n"))
    cur.Ds <- cur.ds$Ds
    Uglob <- cur.ds$Uglob
    Dglob <- cur.ds$Dglob
    logL.values <- cur.ds$logL.values
    ll[i] <- cur.ds$logL} else {
      cur.cs <- RowClustering(x = x, Ds = cur.Ds, Mu = cur.mu, Tau = cur.tau, Xi = cur.xi, Alpha = cur.alpha, Beta = cur.beta, Phi = cur.phi, Uglob = Uglob, Dglob = Dglob)
      cur.Cs <- cur.cs$allocation
      for(r in 1:R)
        for(k in 1:K)
          logL.values[k,r] <- logL.Cocluster(x = x[cur.Cs == k, cur.Ds == r], 
                                             Mu = cur.mu[k,r], 
                                             Tau = cur.tau[k,r], 
                                             Xi = cur.xi[k,r], 
                                             Alpha = cur.alpha[k,r], 
                                             Beta = cur.beta[k,r], 
                                             U = Uglob[[r]], 
                                             d = Dglob[cur.Ds == r])
      ll[i] <- sum(logL.values) 
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
  
  cat(paste("diff(loglikelihood) =",round(diff(ll)[i-1],5),"\n"))
  cat(paste("Row cluster size =", paste(table(cur.Cs), collapse = ", "),"\n"))
  cat(paste("Column cluster size =", paste(table(cur.Ds), collapse = ", "),"\n"))
  if(round(diff(ll)[i-1], 5) < 0) print("Decreasing loglikelihood!")
  if(i %% 10 == 0) plot(2:i, ll[2:i])
  if(i %% 20 == 0 & save.results) save.image(save.to)
}

plot(coordinates, col = cur.Ds, pch = 16)

width.sel <- 9
height.sel <- 7
source("/Users/andreasottosanti/Dropbox/Post-Doc/Code/SpatialBicluster/NewModel/Functions/PosteriorSigma2.R")
Sigma2.post <- posteriorSigma2(x = x, Cs = cur.Cs, Ds = cur.Ds, geneNames = rownames(counts(sce)), expDist = exp(-Dist), Mu = cur.mu, Tau = cur.tau, Xi = cur.xi, Alpha = cur.alpha, Beta = cur.beta, Phi = cur.phi)
for(j in 1:ncol(Sigma2.post$Expected.post)){
  df <- data.frame(Genes = 1:nrow(x), y = Sigma2.post$Expected.post[,j], v = Sigma2.post$Inter.low[,j], w = Sigma2.post$Inter.up[,j])
  pdf(file = paste("/Users/andreasottosanti/Dropbox/Post-Doc/Meetings/10 February 2021/Slide/Images/Sigma2_post_",j,".pdf",sep=""),
      width = width.sel, height = height.sel)
  print(ggplot(df, aes(Genes,y))+geom_pointrange(aes(ymin = v, ymax = w))+theme_bw()+ylab(expression(sigma^2))+
          ggtitle(paste("Cluster",j))+theme(plot.title = element_text(hjust = 0.5)))
  dev.off()
}

rownames(counts(sce))[which(Sigma2.post$Expected.post[,2]>.15)]

# ---Bouveiron method
library(blockcluster)

BouvMethod.results <- coclusterContinuous(data = x, nbcocluster = c(K,R))
plot(coordinates, col = BouvMethod.results@colclass+1, main = "Latent Block Model", pch = 16)
table.BouvMethod <- table(BouvMethod.results@colclass, sce$Cluster);table.BouvMethod


# ---kmeans method

kmeans.results <- list()
kmeans.results$rows <- kmeans(x, centers = K)
kmeans.results$columns <- kmeans(t(x), centers = R)
plot(coordinates, col = kmeans.results$columns$cluster, main = "K-means", pch = 16)
table.kmeans <- table(kmeans.results$columns$cluster, sce$Cluster);table.kmeans


# row comparison
table(kmeans.results$rows$cluster, BouvMethod.results@rowclass)
# column comparison
table(kmeans.results$columns$cluster, BouvMethod.results@colclass)


# ---GRAPHS
sce$myCluster <- cur.Ds

width.sel <- 9
height.sel <- 7

Colors <- c("gold", "#e41a1c", "#377eb8", "#4daf4a", "#ff7f00", "#b2df8a", "#a65628",
            "#999999", "black", "grey", "white", "purple")

pdf(file = "/Users/andreasottosanti/Dropbox/Post-Doc/Meetings/10 February 2021/Slide/Images/tissue_layers.pdf", 
    width = width.sel, height = height.sel)
sce_image_clus(sce, sampleid = "151673", clustervar = "layer_guess", 
               spatial = T, title = " - Layer levels", colors = Colors)
dev.off()

pdf(file = "/Users/andreasottosanti/Dropbox/Post-Doc/Meetings/10 February 2021/Slide/Images/tissue_clusters.pdf",
    width = width.sel, height = height.sel)
sce_image_clus(sce, sampleid = "151673", clustervar = "Cluster", 
               spatial = T, title = " - spatialLIBD clustering", colors = Colors)
dev.off()

pdf(file = "/Users/andreasottosanti/Dropbox/Post-Doc/Meetings/10 February 2021/Slide/Images/tissue_myClusters.pdf",
    width = width.sel, height = height.sel)
sce_image_clus(sce, sampleid = "151673", clustervar = "myCluster", 
               spatial = T, title = " - Our method clustering", colors = Colors)
dev.off()

# ---TABLES
tabl <- data.frame(mu = t(cur.mu), rat = t(cur.tau/cur.xi), phi = sqrt(cur.phi))
xtable::xtable(tabl, digits = 3)
