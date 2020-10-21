
#include <RcppArmadillo.h>

using namespace arma;



arma::vec gaussian_pts(const double mu, const double s, const arma::vec& theta)
{
	const int nt = theta.n_elem;
	arma::vec out(nt);
	double half = (theta[1] - theta[0])/2;

	for(int i=0; i<nt; i++)
		out[i] = R::pnorm(theta[i]+half,mu,s,true,false) - R::pnorm(theta[i]-half,mu,s,true,false);

	out = out / arma::accu(out);
	
	return out;
}


// += for fields with all equal dimensons for use in omp reduction
field<mat>& field_plus(field<mat>& a, const field<mat>& b) 
{
    const int n = b.n_rows;

	for(int i=0; i<n; i++)
		a(i) += b(i);
	
	return a;
}



field<mat> field_init(const field<mat>& orig)
{
	const int n = orig.n_rows;
	field<mat> out(n);

	for(int i=0; i<n; i++)
		out(i) = mat(orig(i).n_rows, orig(i).n_cols, fill::zeros);
	
	return out;
}

mat mat_init(const mat& orig)
{
	return mat(orig.n_rows, orig.n_cols, fill::zeros);
}

vec vec_init(const vec& orig)
{
	return vec(orig.n_elem, fill::zeros);
}

