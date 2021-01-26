
#include <RcppArmadillo.h>
#include "minimize.h"
#include "pl2_item.h"

using namespace arma;

// [[Rcpp::export]]
void test_nlm(arma::ivec& a, arma::vec theta, arma::mat& r, const arma::vec& par_in)
{
	ll_pl2 f(a.memptr(), theta.memptr(), r); 
	const double gtol=1e-10;
	int iter=0,err=0;
	double fret;
	vec par = par_in;

	nlm(par, gtol, iter, fret, f, err);
	
	printf("par: %f, %f\nvalue at min: %f\niter: %i\nerror code: %i\n",par[0],par[1],fret, iter,err);
	fflush(stdout);
}
