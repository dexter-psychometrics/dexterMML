
#ifndef DXM_D2PL_EST_
#define DXM_D2PL_EST_

#include <RcppArmadillo.h>

void estep_2pl_dich(const arma::vec& a, const arma::vec& b, const arma::ivec& pni, const arma::ivec& pcni, const arma::ivec& pi, const arma::ivec& px, 
				const arma::vec& theta, arma::mat& r0, arma::mat& r1, arma::vec& thetabar, arma::vec& sumtheta, arma::vec& sumsig2, const arma::vec& mu, 
				const arma::vec& sigma, const arma::ivec& pgroup, double& ll);
				
#endif