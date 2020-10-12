
#include <RcppArmadillo.h>

arma::vec gaussian_pts(const double mu, const double s, const arma::vec& theta)
{
	const int nt = theta.n_elem;
	arma::vec out(nt);
	double half = (theta[1] - theta[0])/2;

	for(int i=0; i<nt; i++)
		out[i] = R::pnorm(theta[i]+half,mu,s,true,false) - R::pnorm(theta[i]-half,mu,s,true,false);

	out = out / arma::accu(out);
	
	return out;
};