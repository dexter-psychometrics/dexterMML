
#include <RcppArmadillo.h>
#include "shared.h"

using namespace arma;


arma::vec gaussian_pts(const double mu, const double s, const arma::vec& theta)
{
	const int nt = theta.n_elem;
	arma::vec out(nt);

	for(int i=0; i<nt; i++)
		out[i] = R::dnorm(theta[i],mu,s,false);

	out = out / accu(out);
	
	return out;
}


void scale_theta(const vec& means, const vec& sds, const ivec& gn, const vec& theta_start, vec& theta)
{
  const int N = accu(gn);  
  const double q = accu((gn-1) % square(sds) + gn % square(means));  
  const double mu = accu(means % gn)/N;
  const double sigma = std::sqrt((q-N*SQR(mu))/(N-1));
  
  theta = theta_start * sigma + mu;  
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

imat imat_init(const imat& orig)
{
	return imat(orig.n_rows, orig.n_cols, fill::zeros);
}


void imat_ipl_or(imat& m1, const imat& m2)
{
	const int sz = m1.n_elem;
	const int r4 = sz % 4;
	
	for(int i=0;i<r4;i++) m1[i] = m1[i] || m2[i];
	
	for(int i=r4; i<sz; i+=4)
	{
		m1[i] = m1[i] || m2[i];
		m1[i+1] = m1[i+1] || m2[i+1];
		m1[i+2] = m1[i+2] || m2[i+2];
		m1[i+3] = m1[i+3] || m2[i+3];
	}
}



