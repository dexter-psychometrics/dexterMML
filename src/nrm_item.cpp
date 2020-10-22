#include <RcppArmadillo.h>

using namespace arma;

mat nrm_trace(const vec& theta, const ivec& a, const vec& b, const int ncat, const mat& exp_at)
{
	const int nt = theta.n_elem;
	mat out(nt,ncat);		
	
	vec exp_b = exp(b);
	
	for(int j=0; j< nt; j++)
	{
		double sm=0;
		for(int k=0;k<ncat;k++)
			sm += exp_b[k]*exp_at.at(a[k],j);
		for(int k=0;k<ncat;k++)
			out.at(j,k) = exp_b[k]*exp_at.at(a[k],j)/sm;
	}
		
	return out;
}