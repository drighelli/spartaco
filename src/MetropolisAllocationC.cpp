#include "RcppArmadillo.h"
//#include "tgmath.h"

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
double logLCoclusterC(arma::Mat<double> x,
    double Mu, double Tau, double Xi, double Alpha,
    double Beta, arma::Mat<double> U, arma::rowvec d)
{
    arma::Mat<double> Block1 = (x - Mu) * U;
    // Block1.brief_print();
    arma::rowvec invD = 1/(Tau*d + Xi);
    double alpha_post_i = x.n_cols/2 + Alpha;
    // perhaps we need a diagonalMatrix instead of a sparseMatrix here
    arma::mat dd = arma::mat(invD.n_elem, invD.n_elem, arma::fill::zeros);
    dd.diag() = invD;
    arma::mat app = ((Block1 * dd) * Block1.t());
    arma::vec beta_post_i = (app.diag()/2) + Beta;
    double pi = 3.141593;
    double ret1 = -((x.n_rows*x.n_cols)/2*log(2*pi));
    double ret2 = (x.n_rows/2)*sum(log(invD));
    double ret3 = x.n_rows*( (Alpha*log(Beta)) -lgamma(Alpha));
    double ret4 = sum(lgamma(alpha_post_i)- (alpha_post_i*log(beta_post_i)));
    double ret = ret1 + ret2 + ret3 + ret4;
    return(ret);
}


// double MetropolisAllocation(x, Cs, Ds,
//     Uglob, // list of matrices of length R, contains the matrices of eigenvec
//     Dglob, // vector of length p, contains the eigenval
//     Dist,
//     Mu, Tau, Xi, Alpha, Beta, Phi, maxit = 10,
//     min.obs = 3, prob.choices = c(1/2,1/2),
//     print.info = FALSE) {
//
//     if(length(prob.choices) != 2) stop("Wrong probabilities input")
//     if(is.null(prob.choices)) prob.choices <- rep(1/3,3)
//     K <- length(table(Cs))
//     R <- length(table(Ds))
//     if(is.vector(Mu)) Mu <- matrix(Mu, K, R)
//     if(is.vector(Tau)) Tau <- matrix(Tau, K, R)
//     if(is.vector(Xi)) Xi <- matrix(Xi, K, R)
//     if(is.vector(Alpha)) Alpha <- matrix(Alpha, K, R)
//     if(is.vector(Beta)) Beta <- matrix(Beta, K, R)
//     D <- as.vector(table(Ds))
//     Ds.init <- Ds
//     accepted <- FALSE
//     accepted <- m.list <- numeric(maxit)
//     logL.values <- matrix(0,K,R)
//     for(r in 1:R){
//     for(k in 1:K){
//     logL.values[k,r] <- logL.Cocluster(x = x[Cs == k, Ds == r], Mu = Mu[k,r], Tau = Tau[k,r], Xi = Xi[k,r], Alpha = Alpha[k,r],
//     Beta = Beta[k,r], U = Uglob[[r]], d = Dglob[Ds == r])
//     }
//     }
//     logL.den <- sum(logL.values)
//     for(iter in seq_len(maxit)) {
//
//     m <- m.list[iter] <- sample(1:3,1, prob = c(.7,.2,.1))#sample(1:4, 1, prob = c(.4, .3, .2, .1))
//     move <- sample(c("sasd","mamd"), size = 1, prob = prob.choices)
//     if(move == "sasd"){
//     if(print.info) cat(paste("Trying sasd with",m,"moves\n"))
//     gr.start <- sample(1:length(D), 1)
//     gr.end <- sample(setdiff(1:length(D), gr.start), 1)
//     j <- sample(1:D[gr.start], m, replace = F)
//     Ds.star <- Ds
//     Ds.star[which(Ds == gr.start)][j] <- gr.end
//     to.be.changed <- unique(c(gr.start, gr.end))
//     logL.values.star <- logL.values
//     Uglob.star <- Uglob
//     Dglob.star <- Dglob
//     for(r in to.be.changed){
//     eigK <- eigen(exp(-Dist[Ds.star == r, Ds.star == r]/Phi[r]))
//     Uglob.star[[r]] <- eigK$vec
//     Dglob.star[Ds.star == r] <- eigK$val
//     for(k in 1:K){
//     logL.values.star[k,r] <- logL.Cocluster(x = x[Cs == k, Ds.star == r], Mu = Mu[k,r], Tau = Tau[k,r], Xi = Xi[k,r], Alpha = Alpha[k,r],
//     Beta = Beta[k,r], U = Uglob.star[[r]], d = Dglob.star[Ds.star == r])
//     }
//     }
//     log.proposal.num <- sum(log(D[gr.start]-0:(m-1)))
//     log.proposal.den <- sum(log(D[gr.end]-1:m))
//     logL.num <- sum(logL.values.star)
//     A <- exp(logL.num + log.proposal.num - logL.den - log.proposal.den)*all(as.vector(table(Ds.star)) >= min.obs)
//     if(runif(1) <= A){
//     Ds <- Ds.star
//     D <- as.vector(table(Ds))
//     logL.values <- logL.values.star
//     logL.den <- logL.num
//     Uglob <- Uglob.star
//     Dglob <- Dglob.star
//     accepted[iter] <- 1}
//     }
//     if(move == "mamd"){
//     if(print.info) cat(paste("Trying mamd with",m,"moves\n"))
//     gr.start <- sample(1:R, m, replace = TRUE)
//     gr.end <- sapply(1:m, function(k) sample(setdiff(1:R, gr.start[k]), 1))
//     q1r <- sapply(1:R, function(r) sum(gr.start == r))
//     q2r <- sapply(1:R, function(r) sum(gr.end == r))
//     j <- sapply(1:R, function(r) ifelse(q1r[r] != 0, return(sample(1:D[r], q1r[r], replace = FALSE)), 0))
//     Ds.star <- Ds
//     for(r in which(q1r != 0)) Ds.star[which(Ds == r)][j[[r]]] <- gr.end[gr.start == r]
//     to.be.changed <- unique(c(gr.start, gr.end))
//     logL.values.star <- logL.values
//     Uglob.star <- Uglob
//     Dglob.star <- Dglob
//     for(r in to.be.changed){
//     eigK <- eigen(exp(-Dist[Ds.star == r, Ds.star == r]/Phi[r]))
//     Uglob.star[[r]] <- eigK$vec
//     Dglob.star[Ds.star == r] <- eigK$val
//     for(k in 1:K){
//     logL.values.star[k,r] <- logL.Cocluster(x = x[Cs == k, Ds.star == r], Mu = Mu[k,r], Tau = Tau[k,r], Xi = Xi[k,r], Alpha = Alpha[k,r],
//     Beta = Beta[k,r], U = Uglob.star[[r]], d = Dglob.star[Ds.star == r])
//     }
//     }
//     log.proposal.num <- sum(log(factorial(q2r)))+sum(sapply(which(q1r != 0), function(r) sum(log(D[r]+0:(q1r[r]-1)))))
//     log.proposal.den <- sum(log(factorial(q1r)))+sum(sapply(which(q2r != 0), function(r) sum(log(D[r]-q1r[r]+1:q2r[r]))))
//     logL.num <- sum(logL.values.star)
//     A <- exp(logL.num + log.proposal.num - logL.den - log.proposal.den)*all(as.vector(table(Ds.star)) >= min.obs)
//     if(runif(1) <= A){
//     Ds <- Ds.star
//     D <- as.vector(table(Ds))
//     logL.values <- logL.values.star
//     logL.den <- logL.num
//     Uglob <- Uglob.star
//     Dglob <- Dglob.star
//     accepted[iter] <- 1}
//     }
//     }
//     results <- list(Ds = Ds, Uglob = Uglob, Dglob = Dglob, logL = logL.den, logL.values = logL.values,
//     accepted = sum(Ds != Ds.init), m.list = m.list)
//     return(results)
// }
