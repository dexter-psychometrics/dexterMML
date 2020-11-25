#include <RcppArmadillo.h>


using namespace arma;
using Rcpp::Named;



// to do: min should be >=0
// gives two data sets:
// person, item, score
// item, person, score
// [[Rcpp::export]]
Rcpp::List mat_pre(arma::imat& dat, const int max_score)
{
	const int nit = dat.n_cols, np = dat.n_rows;
	
	imat icat(max_score+1,nit, fill::zeros);
	ivec inp(nit,fill::zeros);
	ivec psum(np,fill::zeros), pni(np,fill::zeros);

// margins
	for(int i=0; i<nit; i++)
	{
		const ivec rsp(dat.colptr(i), np, false, true);

		for(int p=0; p<np; p++)
		{
			if(rsp[p]>=0) // NA test
			{
				psum[p] += rsp[p];				
				inp[i]++;
				pni[p]++;
				icat.at(rsp[p],i)++;			
			}		
		}	
	}

	// cumulative pointers	
	ivec icnp(nit+1), pcni(np+1);
	icnp[0] = 0;
	pcni[0] = 0;
	std::partial_sum(inp.begin(),inp.end(),icnp.begin()+1);
	std::partial_sum(pni.begin(),pni.end(),pcni.begin()+1);

	// response vectors organized two ways
	const int sz = icnp[nit];
	
	ivec ip(sz), ix(sz), pi(sz),px(sz);
	
	ivec pindx(pcni.memptr(), np, true, true);

	for(int i=0; i<nit; i++)
	{
		const ivec rsp(dat.colptr(i), np, false, true);
		int indx = icnp[i];

		for(int p=0; p<np; p++)
		{
			if(rsp[p]>=0)
			{
				ip[indx] = p;
				ix[indx++] = rsp[p];
				pi[pindx[p]] = i;
				px[pindx[p]++] = rsp[p]; 
			}
		}
	}
	
	ivec imax(nit), isum(nit,fill::zeros), ncat(nit, fill::zeros);
	for(int i=0;i<nit;i++)
	{
		for(int k = max_score; k>=0; k--)
			if(icat.at(k,i)>0)
			{
				imax[i] = k;
				break;
			}
		for(int k=0; k<=max_score; k++) if(icat.at(k,i)>0)
		{
			isum[i] += k * icat.at(k,i);
			ncat[i]++;
		}
	}
	
	return Rcpp::List::create(
		Named("pi") = pi, Named("px") = px, Named("ip") = ip, Named("ix") = ix,
		Named("inp") = inp, Named("icnp") = icnp, Named("pni") = pni, Named("pcni") = pcni,
		Named("icat") = icat, Named("ncat") = ncat, Named("imax") = imax, Named("isum") = isum, Named("psum") = psum);
}



// ix and px changed in place if necessary
// returns matrix a
// [[Rcpp::export]]
arma::imat categorize(const arma::ivec& inp, const arma::ivec& pni,
						const arma::ivec& icnp, const arma::ivec& pcni,
						const arma::ivec& ip, const arma::ivec& pi,				
						const arma::imat& icat, const arma::ivec& imax,
						const int max_cat,
						arma::ivec& ix, arma::ivec& px)
{
	
	const int nit = icat.n_cols, np = pni.n_elem;
	imat a(max_cat, nit, fill::zeros);
	imat ai(icat.n_rows, icat.n_cols,fill::zeros);
		
	bool recode = false;
	for(int i=0; i<nit; i++)
	{
		int k=0;
		for(int j=0; j<=imax[i]; j++)
		{
			if(icat.at(j,i) == 0) 
				recode = true;
			else
			{
				a.at(k,i) = j;
				ai.at(j,i) = k;
				k++;
			}
		}			
	}
	
	if(recode)
	{	
#pragma omp parallel for	
		for(int i=0; i<nit; i++)
		{
			for(int j=icnp[i]; j<icnp[i+1]; j++)
				ix[j] = ai.at(ix[j],i);
		}
#pragma omp parallel for	
		for(int p=0; p<np; p++)
		{
			for(int j=pcni[p]; j<pcni[p+1]; j++)
				px[j] = ai.at(px[j],pi[j]);
		}
	}	
	return a;
}


// [[Rcpp::export]]
Rcpp::List design_matrices(const arma::ivec& pni, const arma::ivec& pcni, const arma::ivec& pi, const arma::ivec& pg, const int nit, const int ng)
{
	const int np = pni.n_elem;
	imat item(nit,nit,fill::zeros), group(ng,nit,fill::zeros);
	
	std::vector<bool> bk(nit);
	std::unordered_map<std::vector<bool>, int> booklets;
	
	for(int p=0; p<np; p++)
	{
		std::fill(bk.begin(), bk.end(), false);
		int g = pg[p];
		for(int indx = pcni[p]; indx<pcni[p+1]; indx++)
		{
			bk[pi[indx]] = true;
			group.at(g,pi[indx]) = 1;
		}
		booklets.insert(std::make_pair(bk,p));		
	}
	
	for(auto& iter: booklets )
	{
		int p = iter.second;
		for(int i= pcni[p]; i<pcni[p+1]; i++)
			for(int j=i+1; j<pcni[p+1]; j++)
				item.at(pi[i],pi[j]) = 1;	
	}
	item += item.t();
	item.diag().ones(); 
	return Rcpp::List::create(Named("items")=item, Named("groups")=group);
}
