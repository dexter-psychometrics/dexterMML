#ifndef DXM_DATA_
#define DXM_DATA_

#include <RcppArmadillo.h>

void persons_ii(const int item1, const int item2, const arma::ivec& ix,
				const arma::ivec& inp, const arma::ivec& icnp, const arma::ivec& ip,
				arma::ivec& persons, arma::ivec& x1, arma::ivec& x2, int& np);

#endif