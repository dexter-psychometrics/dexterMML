#ifndef DXM_PL2EST_
#define DXM_PL2EST_

#include <RcppArmadillo.h>


long double loglikelihood_2pl(const arma::imat& a, const arma::vec& A, const arma::mat& b, const arma::ivec& ncat, 
				const arma::ivec& pni, const arma::ivec& pcni, const arma::ivec& pi, const arma::ivec& px, 
				const arma::vec& theta, const arma::vec& mu, const arma::vec& sigma, const arma::ivec& pgroup);

#endif