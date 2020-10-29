#ifndef DXM_SHARED_
#define DXM_SHARED_

#include <RcppArmadillo.h>

inline double SQR(const double v){ return v*v; };
inline int SQR(const int v){ return v*v; };
inline double CUB(const double v){ return v*v*v; };
inline int CUB(const int v){ return v*v*v; };

inline int kron(const int a, const int b){return a==b ? 1 : 0; };

arma::vec gaussian_pts(const double mu, const double s, const arma::vec& theta);

arma::field<arma::mat>& field_plus(arma::field<arma::mat>& a, const arma::field<arma::mat>& b); 

arma::field<arma::mat> field_init(const arma::field<arma::mat>& orig);

arma::mat mat_init(const arma::mat& orig);

arma::vec vec_init(const arma::vec& orig);

#pragma omp declare reduction( + : arma::field<arma::mat> : field_plus(omp_out, omp_in)) \
initializer( omp_priv = field_init(omp_orig) )

#pragma omp declare reduction( + : arma::mat : omp_out += omp_in ) \
initializer( omp_priv = mat_init(omp_orig) )

#pragma omp declare reduction( + : arma::vec : omp_out += omp_in ) \
initializer( omp_priv = vec_init(omp_orig) )

#endif