#ifndef DXM_NRMEST_
#define DXM_NRMEST_


#include <RcppArmadillo.h>

long double loglikelihood_nrm(const arma::imat& a, const arma::mat& b, const arma::ivec& ncat, 
				const arma::ivec& pcni, const arma::ivec& pi, const arma::ivec& px, 
				const arma::vec& theta, const arma::vec& mu, const arma::vec& sigma, const arma::ivec& pgroup);

#endif