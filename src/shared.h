#ifndef DXM_SHARED_
#define DXM_SHARED_

#include <RcppArmadillo.h>

inline double SQR(double v){ return v*v; };


arma::vec gaussian_pts(const double mu, const double s, const arma::vec& theta);





#endif