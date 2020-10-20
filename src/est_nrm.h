
#ifndef DXM_EST_NRM_
#define DXM_EST_NRM_

#include <RcppArmadillo.h>
arma::mat nrm_trace(const arma::vec& theta, const arma::ivec& a, const arma::vec& b, const int ncat, const arma::mat& exp_at);



#endif