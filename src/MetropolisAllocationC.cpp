//#include "tgmath.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;
#include <stdio.h>
#include <string.h>
// Save on the typing...
typedef std::pair<double, int>  ptype;


// A comparison function to rank values in descending order
bool compare_values(const ptype &p1, const ptype &p2)
{
    return p1.second > p2.second;
}

// Get the top number of observations
// [[Rcpp::export]]
Rcpp::List table_cpp(arma::vec v, bool sort_data = false)
{

    // Create a map
    std::map<double, int> Elt;

    Elt.clear();

    // Fill the map with occurrences per number.
    for (int i = 0; i != v.size(); ++i) {
        Elt[ v[i] ] += 1;
    }

    // Get how many unique elements exist...
    unsigned int n_obs = Elt.size();

    // Switch map to a vector so that we can sort by value
    std::vector< ptype > sorted_Elt(Elt.begin(), Elt.end());

    if(sort_data){
        // Perform the sort with a custom sort function.
        std::sort(sorted_Elt.begin(), sorted_Elt.end(), compare_values);
    }
    // Else, return.

    // Stop here if you do not need to import into R.
    // Why? There is no ability to export a set w/ a pair into R. *cries*

    // Recast for R using Rcpp::*Vectors to avoid another copy)
    Rcpp::NumericVector result_keys(n_obs);
    Rcpp::IntegerVector result_vals(n_obs);

    unsigned int count = 0;

    // Need to use iterators to access objects
    for( std::vector< ptype >::iterator it = sorted_Elt.begin(); it != sorted_Elt.end(); ++it )
    {
        // Move them into split vectors
        result_keys(count) = it->first;
        result_vals(count) = it->second;

        count++;
    }

    return Rcpp::List::create(Rcpp::Named("lengths") = result_vals,
                              Rcpp::Named("values") = result_keys);
}

double logLCoclusterC(arma::mat x,
    double Mu, double Tau, double Xi, double Alpha,
    double Beta, arma::mat U, arma::vec d)
{
    arma::Mat<double> Block1 = (x - Mu) * U;
    arma::vec invD = 1/(Tau*d + Xi);
    // Rcpp::Rcout << "ncols " << x.n_cols << " Alpha: "<< Alpha << "\n";
    double alpha_post_i = (double)x.n_cols/2 + Alpha;
    // perhaps we need a diagonalMatrix instead of a sparseMatrix here
    arma::mat dd = arma::mat(invD.n_elem, invD.n_elem, arma::fill::zeros);
    dd.diag() = invD;
    arma::mat app = ((Block1 * dd) * Block1.t());
    arma::vec beta_post_i = (app.diag()/2) + Beta;
    // beta_post_i.brief_print();
    double pi = 3.141593;
    double ret1 = -(((double)x.n_rows*(double)x.n_cols)/2*log(2*pi));
    double ret2 = ((double)x.n_rows/2)*sum(log(invD));
    double ret3 = (double)x.n_rows*( (Alpha*log(Beta)) -lgamma(Alpha));
    double ret4 = sum(lgamma(alpha_post_i)- (alpha_post_i*log(beta_post_i)));
    double ret = ret1 + ret2 + ret3 + ret4;
    return(ret);
}
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::export]]
double MetropolisAllocationC(NumericMatrix x,
    arma::vec Cs, arma::vec Ds,
    Rcpp::List Uglob, // list of matrices of length R, contains the matrices of eigenvec
    arma::vec Dglob, // vector of length p, contains the eigenval
    arma::mat Dist,
    arma::mat Mu, arma::mat Tau, arma::mat Xi, arma::mat Alpha, arma::mat Beta,
    arma::vec Phi, int maxit = 10, int min_obs = 3,
    Rcpp::NumericVector prob_choices=Rcpp::NumericVector::create(0.5,0.5),
    bool print_info = 0) {

    if ( prob_choices.size() != 2 ) Rcpp::stop("Wrong probabilities input");

    // converting indexes in C++ format, starting from 0
    Cs = Cs -1;
    Ds = Ds -1;
    arma::vec Kvec=unique(Cs); // <----
    arma::vec Rvec=unique(Ds); // <----
    int K=0;
    if (Mu.is_vec())
    {
        K=1;
    } else {
        K=Mu.n_rows;
    }
    int R=Phi.n_elem;
    // int K=Kvec.n_elem;
    // int R=Rvec.n_elem;

    // Check if vectors missing -> how to check if it's a vector?

    arma::vec D = table_cpp(Ds)["lengths"];
    arma::vec Ds_init = Ds;
    bool accepted = 0;
    Rcpp::NumericVector m_list (maxit); //accepted equal to this ones?
    arma::mat logL_values = arma::mat(K, R, arma::fill::zeros);

    arma::vec goodK=arma::sort(Kvec);
    arma::vec goodR=arma::sort(Rvec);


    for (int rit=0; rit<goodR.n_elem; rit++) {
        for (int kit=0; kit<goodK.n_elem; kit++) {
            int r=goodR[rit];
            int k=goodK[kit];
            arma::uvec csk = find(Cs==k);
            arma::uvec dsr = find(Ds==r);
            arma::mat Ur = Uglob[r];

            logL_values(k,r) = logLCoclusterC(x(csk, dsr),
                Mu(k,r), Tau(k,r), Xi(k,r),
                Alpha(k,r), Beta(k,r),
                Ur, Dglob(dsr));

        }
    } // seems to work, check resulting computations
    double logL_den = accu(logL_values);


    for (int iter=0; iter<maxit; iter++) {
        Rcpp::NumericVector m = Rcpp::RcppArmadillo::sample(
            Rcpp::NumericVector::create(1, 2, 3), 1, 0,
            Rcpp::NumericVector::create(0.7, 0.2, 0.1)); // MISSING M.LIST[ITER] <---- is this a matrix?

        Rcpp::CharacterVector move = Rcpp::RcppArmadillo::sample(
            Rcpp::CharacterVector::create("sasd", "mamd"), 1, 0,
            prob_choices);
        // Rcpp::Rcout << move << "\n";
        std::string move1 = Rcpp::as<std::string>(move);
        if ( move1 == "sasd" ) {
            Rcout << move1 << "\n";
            if(print_info) Rcout << "Trying sasd with " << m << " moves\n";

            int gr_start = Rcpp::RcppArmadillo::sample(
                Rcpp::NumericVector::create(arma::linspace(1, D.n_elem)), 1, 0)[1];

            int gr_end = Rcpp::RcppArmadillo::sample(
                Rcpp::NumericVector::create(arma::linspace(1, D.n_elem)), 1, 0)[1];
                setdiff(arma::linspace(1, D.n_elem), gr_start);
            // TO WORK ON
        }
    }


}










