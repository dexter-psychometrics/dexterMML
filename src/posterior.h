#ifndef DXM_POSTERIOR_
#define DXM_POSTERIOR_

#include <RcppArmadillo.h>

void estep(arma::field<arma::mat>& itrace, const arma::ivec& pni, const arma::ivec& pcni, const arma::ivec& pi, const arma::ivec& px, 
				const arma::vec& theta, arma::field<arma::mat>& r, arma::vec& thetabar, arma::vec& sumtheta, arma::vec& sumsig2, const arma::vec& mu, 
				const arma::vec& sigma, const arma::ivec& pgroup, long double& ll);

arma::mat normalized_posterior(arma::field<arma::mat>& itrace, const arma::ivec& pni, const arma::ivec& pcni, const arma::ivec& pi, const arma::ivec& px, 
						const arma::vec& theta, const arma::vec& mu, const arma::vec& sigma, const arma::ivec& pgroup);


long double loglikelihood(arma::field<arma::mat>& itrace, const arma::ivec& pcni, const arma::ivec& pi, const arma::ivec& px, 
							const arma::vec& theta, const arma::vec& mu, const arma::vec& sigma, const arma::ivec& pgroup);

#endif